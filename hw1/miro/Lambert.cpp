#define _USE_MATH_DEFINES
#include "Lambert.h"
#include "Ray.h"
#include "Scene.h"
#include <algorithm>

Lambert::Lambert(const Vector3 & kd, const Vector3 & ka) :
    m_kd(kd), m_ka(ka)
{

}

Lambert::~Lambert()
{
}

Vector3
Lambert::shade(const Ray& ray, const HitInfo& hit, const Scene& scene) const
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
        //printf("POINT LIGHT\n");
        PointLight* pLight = *plightIter;
        Vector3 l = pLight->position() - hit.P;
        rayLight.o = hit.P;
        rayLight.d = l.normalized();
        if (scene.trace(hitLight, rayLight, 0.0001, l.length())) continue;

        // the inverse-squared falloff
        float falloff = l.length2();

        // normalize the light direction
        l /= sqrt(falloff);

        // get the diffuse component
        float nDotL = dot(hit.N, l);
        Vector3 result = pLight->color();
        result *= m_kd;

        L += std::max(0.0f, nDotL / falloff * pLight->wattage() / PI) * result;
    }

    const AreaLights *alightlist = scene.areaLights();
    // loop over all of the lights
    AreaLights::const_iterator alightIter;
    for (alightIter = alightlist->begin(); alightIter != alightlist->end(); alightIter++)
    {
        AreaLight* aLight = *alightIter;
        Vector3 l = aLight->randPt() - hit.P;
        rayLight.o = hit.P;
        rayLight.d = l.normalized();
        if (scene.trace(hitLight, rayLight, aLight, 0.0001, l.length())) {
            //printf("good, skipped\n");
            continue;
        }

        // the inverse-squared falloff
        float falloff = l.length2();

        // normalize the light direction
        l /= sqrt(falloff);

        // get the diffuse component
        float nDotL = dot(hit.N, l);
        Vector3 result = aLight->color(); // add checks for back front on/off later
        result *= m_kd;

        L += std::max(0.0f, nDotL / falloff * aLight->wattage() / PI) * result;
    }
    
    // add the ambient component
    L += m_ka;
    
    return L;
}

Vector3 Lambert::randReflect(const Ray& ray, const HitInfo& hit) const{
    //return Vector3(1, 2, 3);
    double phi = 2.0 * M_PI*(double)rand() / (double)RAND_MAX;
    double theta = acos((double)rand() / RAND_MAX);
    Vector3 v = Vector3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
    //printf("%f %f %f\n", v[0], v[1], v[2]);
    // generate a basis with surface normal hit.N as the z-axis
    Vector3 z = hit.N.normalized();
    float a = dot(Vector3(1, 0, 0), z);
    float b = dot(Vector3(0, 1, 0), z);
    Vector3 y;
    if (fabs(a) < fabs(b)) y = Vector3(1, 0, 0).orthogonal(z).normalize();
    else y = Vector3(0, 1, 0).orthogonal(z).normalize();
    Vector3 x = cross(y, z).normalize();
    //printf("  %f %f %f\n", x[0], x[1], x[2]);
    //printf("  %f %f %f\n", y[0], y[1], y[2]);
    //printf("  %f %f %f\n", z[0], z[1], z[2]);
    return v[0] * x + v[1] * y + v[2] * z;
}

float Lambert::BRDF(const Ray& ray, const HitInfo& hit) const{
    return 1.0 / M_PI;
}