#define _USE_MATH_DEFINES
#include "Material.h"
#include "Ray.h"
#include "Scene.h"
#include <algorithm>
#include <random>

Material::Material() : m_interacting(true) {}
Material::Material(const bool& interacting) : m_interacting(interacting) {}
Material::~Material() {}


Vector3
Material::radiance(const Vector3& normal, const Vector3& direction) const { return m_powerPerArea / (4.0*M_PI); }


vec3pdf Material::randReflect(const Vector3& ray, const Vector3& hit, const bool& isFront) const {
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

Vector3
Material::shade(const Ray& ray, const HitInfo& hit, const Scene& scene, const bool& isFront) const
{
    Ray rayLight;
    HitInfo hitLight;
    Vector3 L = Vector3(0.0f, 0.0f, 0.0f);

    const Vector3 viewDir = -ray.d; // d is a unit vector

    const PointLights *plightlist = scene.pointLights();
    // loop over all of the POINT lights
    PointLights::const_iterator plightIter;
    for (plightIter = plightlist->begin(); plightIter != plightlist->end(); plightIter++)
    {
        PointLight* pLight = *plightIter;
        Vector3 l = pLight->position() - hit.P;
        rayLight.o = hit.P;
        rayLight.d = l.normalized();
        Vector3 brdf = BRDF(rayLight.d, hit.N, -ray.d);
        if (brdf == 0) continue;
        if (scene.trace(hitLight, rayLight, 0.0001, l.length())) continue;

        // the inverse-squared falloff
        float falloff = l.length2();

        // get the diffuse component
        float nDotL = std::max(0.0f, dot(hit.N, l));
        Vector3 result = pLight->color();

        L += nDotL / falloff * pLight->wattage() *brdf * result;
    }

    const AreaLights *alightlist = scene.areaLights();
    // loop over all of the lights
    AreaLights::const_iterator alightIter;
    for (alightIter = alightlist->begin(); alightIter != alightlist->end(); alightIter++)
    {
        AreaLight* aLight = *alightIter;
        vec3pdf vp = aLight->randPt();
        Vector3 l = vp.v - hit.P; // shoot a shadow ray to a random point on the area light
        rayLight.o = hit.P;
        rayLight.d = l.normalized();

        Vector3 brdf = BRDF(rayLight.d, hit.N, -ray.d);
        if (brdf == 0) continue;
        // if the shadow ray hits the "backside of the light" continue to the next area light
        if (!aLight->intersect(hitLight, rayLight)){
            //printf("front-side of light not visible\n");
            continue;
        }
        // if the shadow ray is occluded by another (hence the "skip") object continue the next light
        if (scene.trace(hitLight, rayLight, aLight, 0.0001, l.length())){
            //printf("random point on light occluded\n");
            continue;
        }

        // the inverse-squared falloff
        float falloff = l.length2();

        float nDotL = dot(hit.N, l);
        Vector3 result = aLight->color();

        L += std::max(0.0f, dot(hitLight.N, -l))* nDotL / falloff* aLight->wattage() / aLight->area() *brdf * result / (vp.p);
    }

    // add the ambient component
    L += m_ka;

    return L;
}