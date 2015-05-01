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

float Lambert::BRDF(const Ray& in, const HitInfo& hit, const Ray& out) const { return 1.0 / M_PI; }

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
        vec3pdf vp = aLight->randPt();
        Vector3 l = vp.v - hit.P; // shoot a shadow ray to a random point on the area light
        rayLight.o = hit.P;
        rayLight.d = l.normalized();
        
        // if the shadow ray hits the "backside of the light" continue to the next area light
        if (!aLight->intersect(hitLight, rayLight)) continue;
        // if the shadow ray is occluded by another (hence the "skip") object continue the next light
        if (scene.trace(hitLight, rayLight, aLight, 0.0001, l.length())) continue;

        // the inverse-squared falloff
        float falloff = l.length2();

        // normalize the light direction
        l /= sqrt(falloff);

        // get the diffuse component
        float nDotL = dot(hit.N, l);
        Vector3 result = aLight->color(); // add checks for back front on/off later
        result *= m_kd;

        L += std::max(0.0f, nDotL / falloff * aLight->wattage() / aLight->area() / PI) * result / (vp.p);
    }
    
    // add the ambient component
    L += m_ka;
    
    return L;
}

vec3pdf Lambert::randReflect(const Ray& ray, const HitInfo& hit) const{
    double phi = 2.0 * M_PI*(double)rand() / (double)RAND_MAX;
    double theta = acos((double)rand() / RAND_MAX);
    Vector3 v = Vector3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
    // generate a basis with surface normal hit.N as the z-axis
    Vector3 z = hit.N.normalized();
    float a = dot(Vector3(1, 0, 0), z);
    float b = dot(Vector3(0, 1, 0), z);
    Vector3 y;
    if (fabs(a) < fabs(b)) y = Vector3(1, 0, 0).orthogonal(z).normalize();
    else y = Vector3(0, 1, 0).orthogonal(z).normalize();
    Vector3 x = cross(y, z).normalize();
    return vec3pdf(v[0] * x + v[1] * y + v[2] * z, 1 / (2.0*M_PI));
}

vec3pdf Lambert::randEmit(const Vector3& n) const {
    double u = (double)rand() / (double)RAND_MAX;
    double v = 2.0 * M_PI*(double)rand() / (double)RAND_MAX;
    Vector3 d = Vector3(cos(v)*sqrt(u), sin(v)*sqrt(u), sqrt(1 - u));
    Vector3 z = n.normalized();
    float a = dot(Vector3(1, 0, 0), z);
    float b = dot(Vector3(0, 1, 0), z);
    Vector3 y;
    if (fabs(a) < fabs(b)) y = Vector3(1, 0, 0).orthogonal(z).normalize();
    else y = Vector3(0, 1, 0).orthogonal(z).normalize();
    Vector3 x = cross(y, z).normalize();
    return vec3pdf(d[0] * x + d[1] * y + d[2] * z, 1.0 / sqrt(1 - u));
}