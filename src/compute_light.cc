#include "compute_light.hh"
#include <random>

std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(0, 1);

inline double clamp(double lo, double hi, double v)
{
    return std::max(lo, std::min(hi, v));
}

Vector reflect(const Vector& incident,
                const Vector& normal)
{
    return incident - 2.0 * normal.dot_product(incident) * normal;
}

Vector get_texture(const Ray &ray, const Texture &texture)
{
    auto pos = ray.tri.uv_pos;

    double u = (1 - ray.u - ray.v) * pos[0][0] + ray.u
                                   * pos[1][0] + ray.v * pos[2][0];
    double v = (1 - ray.u - ray.v) * pos[0][1] + ray.u
                                   * pos[1][1] + ray.v * pos[2][1];
    return texture.get_color(u, v);
}


Vector refract(const Vector& incident,
                const Vector& normal,
                double etai,
                double ior)
{
    double etat = ior;
    double cosi = incident.dot_product(normal);
    cosi = clamp(-1, 1, cosi);

    Vector n = normal;
    if (cosi < 0)
        cosi = -cosi;
    else
    {
        std::swap(etai, etat);
        n = -normal;
    }
    double eta = etai / etat;
    double k = 1 - eta * eta * (1 - cosi * cosi);

    return k < 0 ? Vector(0, 0, 0) : eta * incident + (eta * cosi - sqrt(k)) * n;
}


double fresnel(const Vector &incident,
               const Vector &normal,
               double etat)
{
    double cos_i = incident.dot_product(normal);
    cos_i = clamp(-1, 1, cos_i);

    double etai = 1;

    if (cos_i > 0)
        std::swap(etai, etat);

    double sin_t = etai / etat * sqrt(std::max(0., 1 - cos_i * cos_i));
    // Total internal reflection
    if (sin_t >= 1)
        return 1;
    else
    {
        double cos_t = sqrt(std::max(0., 1 - sin_t * sin_t));
        cos_i = std::abs(cos_i);
        double rs = ((etat * cos_i) - (etai * cos_t)) / ((etat * cos_i) + (etai * cos_t));
        double rp = ((etai * cos_i) - (etat * cos_t)) / ((etai * cos_i) + (etat * cos_t));
        return (rs * rs + rp * rp) / 2;
    }
}

Vector refract_ray(const Scene &scene, const KdTree &tree,
                   const Vector &incident, const Vector &normal,
                   const Vector &inter,
                   double ior, unsigned depth) // get refraction + reflection
{

    double kr = fresnel(incident, normal, ior);
    bool outside = incident.dot_product(normal) < 0;
    Vector bias = BIAS * normal;

    Vector refraction_color;
    if (kr < 1)
    {
        Vector refract_dir = refract(incident, normal, 1, ior).norm_inplace();
        Vector refract_ori = outside ? inter - bias : inter + bias;
        Ray ray_refr(refract_ori, refract_dir);
        ray_refr.ni = ior;

        refraction_color = cast_ray(scene, ray_refr, tree, depth + 1);
    }

    Vector reflect_dir = reflect(incident, normal).norm_inplace();
    Vector reflect_ori = outside ? inter + bias : inter - bias;
    Ray ray_refl(reflect_ori, reflect_dir);

    Vector reflection_color = cast_ray(scene, ray_refl, tree, depth + 1);
    return reflection_color * kr + refraction_color * (1 - kr);
}

Vector ray_to_light(const std::vector<Light*> lights,
                    const Ray &ray, const KdTree &tree)
{
    for (const auto *light : lights) // send ray in every light
    {
        Ray r(ray.o, -light->dir);
        if (!tree.search_inter(r))
            return light->color;
    }
    return Vector(0.f, 0.f, 0.f);
}

Vector cast_ray(const Scene &scene,
                Ray &ray, const KdTree &tree,
                unsigned char depth)
{
    if (depth >= scene.depth) // max depth bi-path tracing
        return ray_to_light(scene.lights, ray, tree);

    double dist = -1;
    if (tree.search(ray, dist))
    {
        const auto material = scene.map.at(scene.mat_names[ray.tri.id]);
        const Vector inter = ray.o + ray.dir * dist;
        Vector normal = ray.tri.normal[0] * (1 - ray.u - ray.v)
          + ray.tri.normal[1] * ray.u +  ray.tri.normal[2] * ray.v;
        normal.norm_inplace();

        Vector direct_color = direct_light(scene, material, ray,
                                           tree, inter, normal, depth);
        Vector indirect_color;

        if (material.illum == 4) //transparence
        {
            Vector refr = reflect(ray.dir, normal);

            Vector origin = inter + refr * BIAS;
            Ray r(origin, refr);

            indirect_color = cast_ray(scene, r, tree, depth + 1) * 0.8; //Fixme
        }
        else if (material.illum == 5)
        {
            indirect_color = refract_ray(scene, tree, ray.dir,
                                         normal, inter,
                                         material.ni, depth) * material.kd;
        }
        else
            indirect_color = indirect_light(scene, tree,
                                            inter, normal, material, ray,
                                            depth);

        return direct_color / M_PI + indirect_color;
    }
    return ray_to_light(scene.lights, ray, tree);
}

void create_coordinate_system(const Vector &N, Vector &Nt, Vector &Nb)
{
    if (std::fabs(N[0]) > std::fabs(N[1]))
        Nt = Vector(N[2], 0, -N[0]) * (1 / sqrtf(N[0] * N[0] + N[2] * N[2]));
    else
        Nt = Vector(0, -N[2], N[1]) * (1 / sqrtf(N[1] * N[1] + N[2] * N[2]));
    Nb = N.cross_product(Nt);
}

Vector uniform_sample_hemisphere(double r1, double r2)
{
    double sinTheta = sqrtf(1 - r1 * r1);
    double phi = 2 * M_PI * r2;
    double x = sinTheta * cosf(phi);
    double z = sinTheta * sinf(phi);
    return Vector(x, r1, z);
  //  return Vector(x, z, r1);
}

float chiGGX(float v)
{
    return v > 0 ? 1 : 0;
}

double GGX_Distribution(const Vector &n, const Vector &h, double alpha)
{
    double NoH = n.dot_product(h);
    double alpha2 = alpha * alpha;
    double NoH2 = NoH * NoH;
    //double den = NoH2 * alpha2 + (1 - NoH2);
    double den = NoH2 * (alpha2 - 1) + 1; // http://graphicrants.blogspot.0om/2013/08/specular-brdf-reference.html
    return (chiGGX(NoH) * alpha2) / ( M_PI * den * den );
}

// https://computergraphics.stackexchange.com/questions/4394/path-tracing-the-cook-torrance-brdf
double g_cook_torrance(const Vector &normal, const Vector &view, const Vector &h,
                       const Vector &sample)
{
    double nh_2 = 2 * normal.dot_product(h);
    double nw = normal.dot_product(view);
    double oh = view.dot_product(h);

    double mid = nh_2 * nw / oh;
    double right = nh_2 * sample.dot_product(normal) / oh;

    return std::min(1.0, std::min(mid, right));
}

double schlick(double n1, double n2, double cos_i)
{
    double ro = (n1 - n2) / (n1 + n2);
    ro *= ro;

    return ro + (1 - ro) * pow((1 - cos_i), 5);
}

double schlick_fresnel(double n1, double n2, const Vector &normal,
        const Vector &incident)
{
    double r0 = (n1-n2) / (n1+n2);
    r0 *= r0;
    double cosX = -normal.dot_product(incident);
    if (n1 > n2)
    {
        double n = n1/n2;
        double sinT2 = n*n*(1.0-cosX*cosX);
        // Total internal reflection
        if (sinT2 > 1.0)
            return 1.0;
        cosX = sqrt(1.0-sinT2);
    }
    double x = 1.0-cosX;
    double ret = r0+(1.0-r0)*x*x*x*x*x;

    double kr = 0.8;
    // adjust reflect multiplier for object reflectivity
    ret = (kr + (1.0- kr) * ret);
    return ret;
}

double beckman(const Vector &normal, const Vector &h, double m)
{
    double nh = normal.dot_product(h);
    double left = 1.0 / (M_PI * (m * m) * (pow(nh, 4)));

    double e = exp((pow(nh, 2) - 1) / (m * m * nh * nh)); // tan = (nh ** 2 - 1) / nh ** 2

    return left * e;
}

double fresnel_transmision(const Vector &incident, const Vector &m)
{
    double nt = 1.5 * 1.5;
    double ni = 1.0 * 1.0;

    double c = incident.dot_product(m);
    double g = sqrt(nt / ni - 1 + c * c);

    double left = ((g - c) * (g - c)) / ((g + c) * (g + c));
    left /= 2;

    double right = pow((c * (g + c) - 1), 2);
    right /= pow(c * (g - c) + 1, 2);
    right += 1;

    return left * right;
}

Vector indirect_light(const Scene &scene,
                      const KdTree &tree, const Vector &inter,
                      const Vector &normal,
                      const Material &material,
                      const Ray &ray,
                      unsigned char depth)
{
    const unsigned nb_ray = scene.nb_ray / ((depth + 1) * 2);
    Vector nt;
    Vector nb;
    Vector indirect_color;
    create_coordinate_system(normal, nt, nb);

   // const Vector vo = inter - scene.cam_pos;
    const Vector vo = inter - scene.cam_pos;
    //const Vector vo = -ray.dir;

    Vector diffuse;
    Vector spec;

    for (unsigned i = 0; i  < nb_ray; ++i)
    {
        double r1 = distribution(generator); // r1 = incident.dot_product(normal)
        double r2 = distribution(generator);

        Vector sample = uniform_sample_hemisphere(r1, r2);
        Vector incident(sample[0] * nb[0] + sample[1] * normal[0] + sample[2] * nt[0],
                        sample[0] * nb[1] + sample[1] * normal[1] + sample[2] * nt[1],
                        sample[0] * nb[2] + sample[1] * normal[2] + sample[2] * nt[2]);

        Vector origin = inter + incident * BIAS;
        Ray ray(origin, incident);

        Vector li = cast_ray(scene, ray, tree, depth + 1);
        diffuse += li * r1;

        Vector hr = vo + incident;
        hr /= hr.get_dist();

        // double d = GGX_Distribution(normal, hr, 0.2); // roughness 0.2
        double d = beckman(normal, hr, 0.2); // roughness 0.2
        double g = g_cook_torrance(normal, vo, hr, incident);
        double f = schlick(1, 1.5, incident.dot_product(hr));
        //     double f = fresnel_transmision(incident, hr);  // broken

        double fr = d * f * g / (4 * incident.dot_product(normal) * vo.dot_product(normal));

        /*
        if (material.ni > 1.2) // not ideal refraction
        {
            Vector ht = -1 * (1 * incident + 1.5 * vo); // ni fix to 1 no 1.5
            ht /= ht.get_dist();

            double ft = incident.dot_product(ht) * vo.dot_product(ht);
            ft /= (incident.dot_product(normal) * (vo.dot_product(normal)));

            double dt = beckman(normal, ht, 0.2);
            double gt = g_cook_torrance(normal, vo, ht, incident);
            double fresnel_t = schlick(1, 1.5, incident.dot_product(ht));

            double right = (1.5 * 1.5) * (1 - fresnel_t) * gt * dt;
            right /= pow((1.0 * (incident.dot_product(ht)) + 1.5 * (vo.dot_product(ht))), 2);

            ft *= right;
            fr += ft;
        }*/
        spec += M_PI * 2 * li * fr * r1; // p = 1 / 2PI
 // p = 1 / 2PI//        spec += ((M_PI / 2 * li * d * f * g) / (normal.dot_product(vo)));
    }

    auto kd_map = scene.map_text.find(material.kd_name); // apply diffuse light
    if (kd_map == scene.map_text.end()) // case not texture
        indirect_color = 2 * diffuse * material.kd; // lambert
    else
    {
        const Vector &text = get_texture(ray, kd_map->second); // case texture
        indirect_color = 2 * diffuse * text;
    }

    spec *= material.ks;
    indirect_color += spec; // apply specular
    indirect_color /= (double)(nb_ray);

    return indirect_color;
}


Vector direct_light(const Scene &scene, const Material &material,
                    const Ray &ray, const KdTree &tree,
                    const Vector &inter, const Vector &normal,
                    int depth)
{
    Vector color;

    if (material.illum == 4) //transparence
    {
        //doule kr = fresnel(light.dir, normal, material.ni);

        Vector reflect_dir = reflect(ray.dir, normal);
        reflect_dir.norm_inplace();

        Vector origin = inter + reflect_dir * BIAS;
        Ray r(origin, reflect_dir);

        color += cast_ray(scene, r, tree, depth + 1) * 0.8;
        return color;
    }
    else if (material.illum == 5) // refraction
    {
        color = refract_ray(scene, tree, ray.dir,
                                     normal, inter, material.ni, depth);

    }

    double rat;
    for (const auto *light : scene.lights) // diffuse light
    {
        Vector L = light->compute_light(inter, tree, rat);

        double diff = 0.f;
        if (rat > 0)
        {
            diff = L.dot_product(normal);
            if (diff < 0)
                diff = 0;
        }

        double spec = 0;
        if (!(L[0] == 0 && L[1] == 0 && L[2] == 0))
        {
            Vector R = reflect(L, normal);
            R.norm_inplace();

            double spec_coef = ray.dir.dot_product(R);
            if (spec_coef < 0)
                spec_coef = 0;
            spec = pow(spec_coef, material.ns);
            if (spec < 0)
                spec = 0;
        }
        if (material.illum < 4 && diff)
        {
            auto kd_map = scene.map_text.find(material.kd_name);
            if (kd_map != scene.map_text.end())
            {
                const Vector &text = get_texture(ray, kd_map->second);
                color += light->color *  text * diff * rat;
            }
            else
                color += light->color * material.kd * diff * rat;
        }
        if (material.illum != 1 && spec)
            color += (light->color * spec * material.ks);
    }

    auto ka_map = scene.map_text.find(material.ka_name);
    if (ka_map != scene.map_text.end())
    {
        const Vector &text = get_texture(ray, ka_map->second);
        color += text * scene.a_light;
    }
    else
        color += material.ka * scene.a_light;

    return color;
}
