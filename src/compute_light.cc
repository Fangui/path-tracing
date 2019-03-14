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

Vector refract(const Vector& incident,
                const Vector& normal,
                double ior)
{
    double cosi = incident.dot_product(normal);
    cosi = clamp(-1, 1, cosi);
    double etai = 1;

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
               double ior)
{
    double cosi = incident.dot_product(normal);
    cosi = clamp(-1, 1, cosi);

    double etai = 1;

    if (cosi > 0)
        std::swap(etai, ior);

    double sint = etai / ior * sqrt(std::max(0., 1 - cosi * cosi));
    // Total internal reflection
    if (sint >= 1)
        return 1;
    else {
        double cost = sqrt(std::max(0., 1 - sint * sint));
        cosi = std::abs(cosi);
        double Rs = ((ior * cosi) - (etai * cost)) / ((ior * cosi) + (etai * cost));
        double Rp = ((etai * cosi) - (ior * cost)) / ((etai * cosi) + (ior * cost));
        return (Rs * Rs + Rp * Rp) / 2;
    }
    // As a consequence of the conservation of energy, transmittance is given by:
    // kt = 1 - kr;
}

Vector cast_ray(const Scene &scene, 
                Ray &ray, const KdTree &tree,
                unsigned char depth)
{
    if (depth > 4) // max depth
        return Vector(0.f, 0.f, 0.f);

    double dist = -1;
    if (tree.search(ray, dist))
    {
        const auto material = scene.map.at(scene.mat_names[ray.tri.id]);
        const Vector inter = ray.o + ray.dir * dist;
        Vector normal = ray.tri.normal[0] * (1 - ray.u - ray.v) 
          + ray.tri.normal[1] * ray.u +  ray.tri.normal[2] * ray.v; //FIXME
        normal.norm_inplace();

        Vector direct_color = direct_light(scene, material, ray, 
                                           tree, inter, normal, depth);

        return direct_color;
        Vector indirect_color = indirect_light(scene, tree,
                                               inter, normal, depth);

        Vector res = (direct_color / M_PI  + 2 * indirect_color)  * 0.2; //FIXME albedo
        return res;

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
    const unsigned nb_ray = 4;
    Vector nt;
    Vector nb;
    Vector indirect_color;
    create_coordinate_system(normal, nt, nb);
    constexpr double inv_pdf = 2 * M_PI;

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
        indirect_color += cast_ray(scene, ray, tree, depth + 1) * r1 * inv_pdf;
    }
    indirect_color /= (double)nb_ray;

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
        r.ni = material.ni;

        color += cast_ray(scene, r, tree, depth + 1) * 0.8;
        return color;
    }
    else if (material.illum == 5)
    {
        Vector refl = reflect(ray.dir, normal).norm_inplace();
        Vector refr = refract(ray.dir, normal, material.ni).norm_inplace();

        double bias = 0.001;
        if (refr.dot_product(normal) < 0)
            bias = -bias;

        Vector refr_ray_o = inter + normal * bias;
        Ray r(refr_ray_o, refl);
        Ray r_refr(refr_ray_o, refr);

        r.ni = material.ni;
        r_refr.ni = material.ni;

        Vector refl_color = cast_ray(scene, r, tree, depth + 1);
        Vector refr_color = cast_ray(scene, r_refr, tree, depth + 1);
        double kr = fresnel(ray.dir, normal, material.ni);
        return refl_color * kr + refr_color * (1 - kr); /*


         double m = ray.ni / material.ni;
         double c = (-1 * normal).dot_product(ray.dir);

         Vector r = m * ray.dir + (m * c - sqrt(1 - m * m * (1 - c * c))) * normal;
         Vector o = ray.o  + 0.001 * r;

         Ray refr(o, r);
         refr.ni = material.ni;

         return cast_ray(scene, refr, tree, depth + 1) * 0.8;*/

    }

    for (const auto &light : scene.lights)
    {
        Vector L = light.dir * -1;

        Vector o_shadow = inter + L * 0.001; // biais
        Ray shadow_ray(o_shadow, L);

        double diff = 0.f;
        if (!tree.search_inter(shadow_ray))
        {
            diff = L.dot_product(normal);
            if (diff < 0)
                diff = 0;
        }
        Vector R = reflect(L, normal);
        R.norm_inplace();

        double spec_coef = ray.dir.dot_product(R);
        if (spec_coef < 0)
            spec_coef = 0;
        double spec = pow(spec_coef, material.ns);

        if (diff)
        {
            auto kd_map = scene.map_kd.find(material.kd_name);
            if (kd_map != scene.map_kd.end())
            {
                auto pos = ray.tri.uv_pos;
                const auto &map = kd_map->second;

                const Vector &text = map.get_color(pos[0][0], pos[0][1]) * (1 - ray.u - ray.v) 
                                   + map.get_color(pos[1][0], pos[1][1]) * ray.u 
                                   + map.get_color(pos[2][0], pos[2][1]) * ray.v;
                //const Vector &text = map.get_color(pos[0][0], pos[0][1]);
                color += light.color *  text * diff;
            }
            else
                color += light.color * material.kd * diff;
        }
        if (material.illum != 1)
            color += (light.color * spec * material.ks);
    }

    color += material.ka * scene.a_light;

    if (color[0] > 1)
        color[0] = 1;

    if (color[1] > 1)
        color[1] = 1;

    if (color[2] > 1)
        color[2] = 1;

    return color;
}
