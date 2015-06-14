#define _USE_MATH_DEFINES
#include "RefractiveInterface.h"
#include "Ray.h"
#include "Scene.h"
#include <algorithm>
#include <random>

RefractiveInterface::RefractiveInterface(const Vector3 & kr, const Vector3 & kt, const Vector3 & ka, const float& n0, const float& n1) : Material(), m_kr(kr), m_kt(kt), m_ka(ka), m_n0(n0), m_n1(n1) {}
RefractiveInterface::~RefractiveInterface() {}

Vector3 RefractiveInterface::BRDF(const Vector3& in, const Vector3& normal0, const Vector3& out, const bool& isFront) const {
    // normal is assumed to be facing the FRONT of the interface.
    // First check to see that all three vectors are in the same plane
    Vector3 normal = normal0;
    if (isFront == false) normal = -normal0;
    bool coPlanarCheck = false;
    if (fabs(dot(in, cross(normal, out).normalize())) < 0.000001 ||
        fabs(dot(normal, cross(out, in).normalize())) < 0.000001 ||
        fabs(dot(out, cross(in, normal).normalize())) < 0.000001) coPlanarCheck = true;
    if (coPlanarCheck == false) {
        std::cout << "not coplanar" << std::endl;
        return 0;
    }
    ////////////////////////////////////////////////////////////////////////////////////
    float nIn = m_n0;
    float nOut = m_n1;
    float cosIn = dot(normal, in);
    float cosOut = dot(normal, out);
    float sinIn = sqrt(1.0f - cosIn*cosIn); // sin's are all positive by definition of sqrt(...)
    float sinOut = sqrt(1.0f - cosOut*cosOut);
    if (cosIn < 0) nIn = m_n1;
    if (cosOut > 0) nOut = m_n0;
    float sinCrit = fmin(m_n0, m_n1) / fmax(m_n0, m_n1); // critical sin(theta) for total internal reflection
    if (cosIn*cosOut < 0) { // rays are on opposite sides of interface (TRANSMITTED)
        float nSinIn = nIn*sinIn;
        float nSinOut = nOut*sinOut;
        if (fmin(nSinIn, nSinOut) / fmax(nSinIn, nSinOut) > 0.99999) return (m_kt[0] + m_kt[1] + m_kt[2]) / 3;
        else return 0;
    }
    else if (cosIn*cosOut > 0) { // rays are on the same side of interface (REFLECTED)
        if (fmin(sinIn, sinOut) / fmax(sinIn, sinOut) > 0.99999) return (m_kr[0] + m_kr[1] + m_kr[2]) / 3;
        else return 0;
    }
    else { // at least one of the rays is parallel to the interface
        if (cosIn == 0) {
            if (cosOut == 0) return 1; // special case of ray skimming the surface (no energy lost)
            if (nIn > nOut) return 0;
            if (fmin(sinCrit, sinOut) / fmax(sinCrit, sinOut) > 0.99999) return (m_kt[0] + m_kt[1] + m_kt[2]) / 3;
            return 0;
        }
        else if (cosOut == 0) {
            if (nIn < nOut) return 0;
            if (fmin(sinCrit, sinIn) / fmax(sinCrit, sinIn) > 0.99999) return (m_kt[0] + m_kt[1] + m_kt[2]) / 3;
            return 0;
        }
    }
    return 0;
}

vec3pdf RefractiveInterface::randReflect(const Vector3& in, const Vector3& normal0, const bool& isFront) const{
    //return vec3pdf(-in, 1);
    if (m_n0 == m_n1) return vec3pdf(-in, 1);
    Vector3 normal = normal0;
    if (isFront == false) normal = -normal0;
    float cosIn = dot(in, normal);
    float sinIn = sqrt(1 - cosIn*cosIn);
    float nIn = m_n0;
    float nOpp = m_n1;
    if (cosIn < 0) {
        nIn = m_n1;
        nOpp = m_n0;
    }
    float sinCrit = fmin(m_n0, m_n1) / fmax(m_n0, m_n1);
    if (nIn > nOpp && sinIn > sinCrit) return vec3pdf(2 * dot(normal, in)*normal - in, 1);
    float sinOut = sinIn*nIn / m_n0;
    float thetaOut = asin(sinOut);
    Vector3 rotAxis = cross(normal, in).normalize();
    if (nIn == m_n1) rotAxis = -rotAxis;
    return vec3pdf(normal.rotated(thetaOut, rotAxis), 1);//*/
}

Vector3
RefractiveInterface::shade(const Ray& ray, const HitInfo& hit, const Scene& scene, const bool& isFront) const
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
        Vector3 brdf = BRDF(rayLight.d, hit.N, -ray.d, isFront);
        if (brdf == 0) continue;
        if (scene.trace(hitLight, rayLight, 0.0001, l.length())) continue;

        // the inverse-squared falloff
        float falloff = l.length2();

        float nDotL = fabs(dot(hit.N, l));
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

        Vector3 brdf = BRDF(rayLight.d, hit.N, -ray.d, isFront);
        if (brdf == 0) continue;
        // if the shadow ray hits the "backside of the light" continue to the next area light
        if (!aLight->intersect(hitLight, rayLight)){
            continue;
        }
        // if the shadow ray is occluded by another (hence the "skip") object continue the next light
        if (scene.trace(hitLight, rayLight, aLight, 0.0001, l.length())){
            continue;
        }

        // the inverse-squared falloff
        float falloff = l.length2();

        float nDotL = fabs(dot(hit.N, l));
        Vector3 result = aLight->color();

        L += std::max(0.0f, dot(hitLight.N, -l))* 0.0f, nDotL / falloff*
            aLight->wattage() / aLight->area() *brdf * result / (vp.p);
    }

    // add the ambient component
    L += m_ka;

    return L;
}