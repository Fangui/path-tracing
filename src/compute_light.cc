#include "compute_light.hh"

Vector reflect(const Vector& incident,
               const Vector& normal)
{
    return incident - (normal * (2.0 * normal.dot_product(incident)));
}

Vector cast_ray(const Scene &scene, 
                Ray &ray, const KdTree &tree,
                unsigned char depth)
{
    if (depth > 2) // max depth
        return Vector(0.f, 0.f, 0.f);

    float dist = -1;
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

        Vector res = direct_color * (1. / M_PI) + indirect_color * 2; //FIXME albedo
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

Vector uniform_sample_hemisphere(float r1, float r2)
{
    float sinTheta = sqrtf(1 - r1 * r1);
    float phi = 2 * M_PI * r2;
    float x = sinTheta * cosf(phi);
    float z = sinTheta * sinf(phi);
    return Vector(x, r1, z);
  //  return Vector(x, z, r1); 
}

#include <random> 
std::default_random_engine generator;
std::uniform_real_distribution<float> distribution(0, 1); 

Vector indirect_light(const Scene &scene,
                      const KdTree &tree, const Vector &inter,
                      const Vector &normal, unsigned char depth)
{
    const unsigned nb_ray = 128;
    Vector nt;
    Vector nb;
    Vector indirect_color;
    create_coordinate_system(normal, nt, nb);
    constexpr float inv_pdf = 2 * M_PI;

    for (unsigned i = 0; i  < nb_ray; ++i)
    {
        float r1 = distribution(generator);
        float r2 = distribution(generator);
        Vector sample = uniform_sample_hemisphere(r1, r2);
        Vector sample_world(sample[0] * nb[0] + sample[1] * normal[0] + sample[2] * nt[0],
                            sample[0] * nb[1] + sample[1] * normal[1] + sample[2] * nt[1],
                            sample[0] * nb[2] + sample[1] * normal[2] + sample[2] * nt[2]);
        //sample_world.norm_inplace();

        Vector origin = inter + sample_world * 0.001; // bias
        Ray ray(origin, sample_world);
        indirect_color += cast_ray(scene, ray, tree, depth + 1) * r1 * inv_pdf;
    }
    indirect_color *= (1. / (float)nb_ray);

    return indirect_color;
}

Vector direct_light(const Scene &scene, const Material &material,
                    const Ray &ray, const KdTree &tree, 
                    const Vector &inter, const Vector &normal,
                    int depth)
{
    Vector color;
    for (const auto &light : scene.lights)
    {
        Vector L = light.dir * -1;

        Vector o_shadow = inter + L * 0.001; // biais
        Ray shadow_ray(o_shadow, L);

        float diff = 0.f;
        bool b = tree.search_inter(shadow_ray);
        if (!b)
        {
            diff = L.dot_product(normal);
            if (diff < 0)
                diff = 0;
        }
        Vector R = reflect(L, normal);
        R.norm_inplace();

        float spec_coef = ray.dir.dot_product(R);
        if (spec_coef < 0)
            spec_coef = 0;
        float spec = pow(spec_coef, material.ns);

        if (material.illum == 4) //transparence
        {
            if (b)
                continue;
            Vector origin = inter + normal * 0.001;
            Vector ref = reflect(light.dir, normal);
            Ray r(origin, ref);

            color +=  cast_ray(scene, r, tree, depth + 1);
        }
        else if (diff)
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
