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
Material::radiance(const Vector3& normal, const Vector3& direction) const { return m_powerPerArea / (4.0*M_PI); }


vec3pdf Material::randReflect(const Vector3& ray, const Vector3& hit) const {
    double phi = 2.0 * M_PI*(double)rand() / (double)RAND_MAX;
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

Vector3 Material::sum_L_cosTheta_dOmega() const {
    Vector3 integral = 0;
    float divisions_theta = 100;
    float divisions_phi = 100;
    float dPhi = 2.0*M_PI / divisions_phi;
    float dTheta = M_PI / divisions_theta;
    for (float i = 0; i < divisions_theta; i++){
        float Theta0 = i*dTheta;
        float Theta1 = Theta0 + dTheta;
        float Theta = (Theta0 + Theta1) / 2.0;
        for (float j = 0; j < divisions_phi; j++){
            double Phi = 2.0*M_PI*dPhi;
            double d_CosThetaOmega = dPhi*(cos(2 * Theta0) - cos(2 * Theta1)) / 4.0;
            integral += radiance(Vector3(0, 0, 1), Vector3(sin(Theta)*cos(Phi), sin(Theta)*sin(Phi), cos(Theta)))*d_CosThetaOmega;
        }
    }
    return integral;
}