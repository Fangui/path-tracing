#include "compute_light.hh"

Vector reflect(const Vector& incident,
               const Vector& normal)
{
    return incident - (normal * (2.0 * normal.dot_product(incident)));
}

Vector direct_light(const Material &material, const Scene &scene,
                    const Ray &ray, const KdTree &tree, 
                    const Vector& inter)
{
    Vector color;
    for (const auto &light : scene.lights)
    {
        Vector L = light.dir * -1;
        L.norm_inplace();

        Vector o_shadow = inter + L * 0.001; // biais
        Ray shadow_ray(o_shadow, L);

        float diff = 0.f;
        Vector normal_n = ray.tri.normal[0] * (1 - ray.u - ray.v) 
                        + ray.tri.normal[1] * ray.u +  ray.tri.normal[2] * ray.v; //FIXME
        normal_n.norm_inplace();

        if (!tree.search_inter(shadow_ray))
        {
            diff = L.dot_product(normal_n);
            if (diff < 0)
                diff = 0;
        }
        Vector R = reflect(L, normal_n);
        R.norm_inplace();

        float spec_coef = ray.dir.norm().dot_product(R);
        if (spec_coef < 0)
            spec_coef = 0;
        float spec = pow(spec_coef, material.ns);

        if (diff)
            color += light.color * material.kd * diff;
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
