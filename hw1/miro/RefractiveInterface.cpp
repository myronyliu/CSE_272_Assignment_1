#define _USE_MATH_DEFINES
#include "RefractiveInterface.h"
#include "Ray.h"
#include "Scene.h"
#include <algorithm>
#include <random>

RefractiveInterface::RefractiveInterface(const Vector3 & kr, const Vector3 & kt, const Vector3 & ka, const float& n) :
m_kr(kr), m_kt(kt), m_ka(ka), m_n(n) {}
RefractiveInterface::~RefractiveInterface() {}

// in and normal should be on the same side when this is called
Vector3 RefractiveInterface::BRDF(const Vector3& in, const Vector3& normal, const Vector3& out, const bool& isFront) const {
    // First check to see that all three vectors are in the same plane
    if (!coplanar(in, normal, out)){
        return Vector3(0, 0, 0);
    }

    Vector3 sinVecIn = cross(in, normal);
    Vector3 sinVecOut = cross(out, normal);

    if (dot(sinVecIn, sinVecOut) > 0) {
        return Vector3(0, 0, 0);
    }
    
    float cosIn = dot(in, normal);
    float cosOut = dot(out, normal);

    //return Vector3(1, 1, 1) / fabs(cosOut);

    if (cosOut == 0) return Vector3(0, 0, 0); // this case cannot be handled due to division by zero

    if (fabs((cosOut - cosIn) / cosIn) < 0.000001) return m_kr / fabs(cosOut); // reflected

    float sinIn = sinVecIn.length();
    float sinOut = sinVecOut.length();
    float sinIn_sinOut = sinIn / sinOut;

    float nOut_nIn; // ratio of indices of refraction n_out/n_in 
    if (isFront) nOut_nIn = m_n;
    else nOut_nIn = 1 / m_n; // this is also sin(theta_critical)

    if (fabs((sinIn_sinOut - nOut_nIn) / nOut_nIn) < 0.000001) return m_kt / cosOut;
    return Vector3(0, 0, 0);
}

vec3pdf RefractiveInterface::randReflect(const Vector3& in, const Vector3& normal, const bool& isFront) const{
    //return vec3pdf(-in, fabs(dot(in,normal)));
    float cosIn = dot(in, normal);
    float sinIn = sqrt(1 - cosIn*cosIn);

    float nOut_nIn;
    if (isFront) nOut_nIn = m_n;
    else nOut_nIn = 1 / m_n;

    Vector3 outR = 2 * dot(normal, in)*normal - in;

    //return vec3pdf(outR, 1);

    if (nOut_nIn < 1 && sinIn>nOut_nIn) {
        return vec3pdf(outR, 1);
    }

    // otherwise, both transmission and reflection are possible

    float sinOut = sinIn / nOut_nIn;
    float cosOut = sqrt(1 - sinOut*sinOut);
    Vector3 rotAxis = cross(normal, in);
    Vector3 outT = -normal.rotated(asin(sinOut), rotAxis);

    float rhoPara = (nOut_nIn*cosIn - cosOut) / (nOut_nIn*cosIn + cosOut);
    float rhoPerp = (cosIn - nOut_nIn*cosOut) / (cosIn + nOut_nIn*cosOut);

    float probR = (rhoPara*rhoPara + rhoPerp*rhoPerp) / 2;

    float rn = (float)rand() / RAND_MAX;
    if (rn < probR) {
        return vec3pdf(outR, probR);
    }
    else {
        return vec3pdf(outT, 1 - probR);
    }
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