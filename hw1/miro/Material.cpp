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
Material::powerPerPatchPerSolidAngle(const Vector3& normal, const Vector3& direction) const { return m_powerPerPatch / (4.0*M_PI); }

vec3pdf
Material::randReflect(const Ray&, const HitInfo&) const{
    double phi = 2.0 * M_PI*((double)rand() / RAND_MAX);
    double theta = acos(2.0*(double)rand() / RAND_MAX - 1.0);
    return vec3pdf(Vector3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)), 1.0 / (4.0*M_PI));
}

vec3pdf
Material::randEmit(const Vector3& n) const{
    double phi = 2.0 * M_PI*(double)rand() / (double)RAND_MAX;
    double theta = acos((double)rand() / RAND_MAX);
    Vector3 v = Vector3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
    Vector3 z = n.normalized();
    float a = dot(Vector3(1, 0, 0), z);
    float b = dot(Vector3(0, 1, 0), z);
    Vector3 y;
    if (fabs(a) < fabs(b)) y = Vector3(1, 0, 0).orthogonal(z).normalize();
    else y = Vector3(0, 1, 0).orthogonal(z).normalize();
    Vector3 x = cross(y, z).normalize();
    return vec3pdf(v[0] * x + v[1] * y + v[2] * z, 1.0 / (2.0*M_PI));
}