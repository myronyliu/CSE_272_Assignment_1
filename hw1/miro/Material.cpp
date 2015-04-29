#define _USE_MATH_DEFINES
#include "Material.h"

Material::Material()
{
}

Material::~Material()
{
}

Vector3
Material::shade(const Ray&, const HitInfo&, const Scene&) const
{
    return Vector3(1.0f, 1.0f, 1.0f);
}

Vector3
Material::randReflect(const Ray&, const HitInfo&) const{
    double phi = 2.0 * M_PI*(double)rand() / (double)RAND_MAX;
    double theta = acos(2.0*(double)rand() / RAND_MAX - 1.0);
    return Vector3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
}

float
Material::BRDF(const Ray&, const HitInfo&) const{
    return 1.0 / M_PI;
}
