#include "compute_light.hh"


inline
double clamp(double lo, double hi, double v)
{
    return std::max(lo, std::min(hi, v));
}

Vector reflect(const Vector& incident,
                const Vector& normal)
{
    return incident - (normal * (2.0 * normal.dot_product(incident)));
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
    double cosi = incident.dot_product(normal);
    cosi = clamp(-1, 1, cosi);

    Vector n = normal;
    if (cosi < 0)
        cosi = -cosi;
    else
    {
        std::swap(etai, ior);
        n = -1 * normal;
    }
    double eta = etai / ior;
    double k = 1 - eta * eta * (1 - cosi * cosi);

    return k < 0 ? Vector(0, 0, 0) : eta * incident + (eta * cosi -sqrt(k)) * n;
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
    // As a consequence of the conservation of energy, transmittance is given by:
    // kt = 1 - kr;
}

Vector cast_ray(const Scene &scene,
                Ray &ray, const KdTree &tree,
                unsigned char depth)
{
    if (depth > 1) // max depth
    {
        for (const auto *light : scene.lights) // send ray in every light
        {
            Ray r(ray.o, -light->dir);
            if (!tree.search_inter(r))
                return light->color;
        }
        return Vector(0.f, 0.f, 0.f);
    }
    double dist = -1;
    if (tree.search(ray, dist))
    {
        const auto material = scene.map.at(scene.mat_names[ray.tri.id]);
        const Vector inter = ray.o + ray.dir * dist;
        Vector normal = ray.tri.normal[0] * (1 - ray.u - ray.v)
          + ray.tri.normal[1] * ray.u +  ray.tri.normal[2] * ray.v;
        normal.norm_inplace();

        /*
        Vector direct_color = direct_light(scene, material, ray,
                                           tree, inter, normal, depth);
        */
        Vector indirect_color = indirect_light(scene, tree,
                                               inter, normal, depth);

      //  Vector res = (direct_color  + 2 * indirect_color );
      //  Vector res = direct_color + indirect_color * M_PI / 2;
        
        auto kd_map = scene.map_text.find(material.kd_name);
        if (kd_map == scene.map_text.end()) // case not texture 
            return indirect_color * 2 * material.kd;
        else
        {
            const Vector &text = get_texture(ray, kd_map->second); // case texture 
            return indirect_color * 2 * text;
        }
    }

    for (const auto *light : scene.lights)
    {
        if (light->dir.dot_product(ray.dir) < 0) // FIXME assume hit directional
            return light->color;
    }

    return Vector(0.f, 0.f, 0.f);
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

#include <random>
std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(0, 1);

Vector indirect_light(const Scene &scene,
                      const KdTree &tree, const Vector &inter,
                      const Vector &normal, unsigned char depth)
{
    const unsigned nb_ray = 24;
    Vector nt;
    Vector nb;
    Vector indirect_color;
    create_coordinate_system(normal, nt, nb);
    //constexpr double inv_pdf = 2 * M_PI;
    //const double pdf = 1.0 / (2 * M_PI);

    for (unsigned i = 0; i  < nb_ray; ++i)
    {
        double r1 = distribution(generator);
        double r2 = distribution(generator);
        Vector sample = uniform_sample_hemisphere(r1, r2);
        Vector sample_world(sample[0] * nb[0] + sample[1] * normal[0] + sample[2] * nt[0],
                            sample[0] * nb[1] + sample[1] * normal[1] + sample[2] * nt[1],
                            sample[0] * nb[2] + sample[1] * normal[2] + sample[2] * nt[2]);
        //sample_world.norm_inplace();

        Vector origin = inter + sample_world * 0.001; // bias
        Ray ray(origin, sample_world);
        indirect_color += cast_ray(scene, ray, tree, depth + 1) * r1; // * inv_pdf;
    }
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
        //double kr = fresnel(light.dir, normal, material.ni);

        Vector refr = reflect(ray.dir, normal);
        double bias = 0.001;

        Vector origin = inter + normal * bias;
        Ray r(origin, refr);

        color += cast_ray(scene, r, tree, depth + 1) * 0.8;
        return color;
    }
    else if (false && material.illum == 5)
    {
        double kr = fresnel(ray.dir, normal, material.ni);
        bool outside = ray.dir.dot_product(normal) < 0;
        Vector bias = 0.001 * normal;

        Vector refraction_color;
        if (kr < 1)
        {
            Vector refraction_direction = refract(ray.dir, normal, material.ni, material.ni).norm_inplace();
            Vector refract_ray_orig = outside ? inter - bias : inter + bias;
            Ray ray_refr(refract_ray_orig, refraction_direction);
            ray_refr.ni = material.ni;

            return refraction_color = cast_ray(scene, ray_refr, tree, depth + 1);
        }

        Vector reflect_direction = reflect(ray.dir, normal).norm_inplace();
        Vector reflection_ray_orig = outside ? inter + bias : inter - bias;
        Ray ray_refr(reflection_ray_orig, reflect_direction);

        Vector reflection_color = cast_ray(scene, ray_refr, tree, depth + 1);

        return reflection_color * kr + refraction_color * (1 - kr) * 0.9;
    }

    double rat;
    for (const auto *light : scene.lights)
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
