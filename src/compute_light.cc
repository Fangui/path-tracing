#include "compute_light.hh"

Vector reflect(const Vector& incident,
               const Vector& normal)
{
    return incident - (normal * 2.0 * normal.dot_product(incident));
}

Vector direct_light(const Material &material, const Scene &scene,
                    const Ray &ray, const KdTree &tree)
{
    Vector color;
    for (const auto &light : scene.lights)
    {
        Vector L = light.dir * -1;
        L.norm_inplace();

        Vector o_shadow = ray.o + L * 0.001;
        Ray shadow_ray(o_shadow, L);

        if (tree.search_inter(shadow_ray))
            return Vector(0, 0, 0);

        Vector normal_n = ray.tri.normal[0].norm(); //FIXME
        float diff = L.dot_product(normal_n);
        if (diff < 0)
            diff = 0;

        Vector R = L - normal_n * (2 * normal_n.dot_product(L));
     //   Vector R = reflect(L, normal_n);
        R.norm_inplace();

        float spec_coef = ray.dir.norm().dot_product(R);
        if (spec_coef < 0)
            spec_coef = 0;
        float spec = pow(spec_coef, material.ns);

        color += light.color * material.kd * diff;
        if (material.illum != 1)
            color += (light.color * spec * material.ks);
    }

    color += material.ka * scene.a_light;

    if (color[0] > 1)
        color.set_x(1);

    if (color[1] > 1)
        color.set_y(1);

    if (color[2] > 1)
        color.set_z(1);

    return color;
}
