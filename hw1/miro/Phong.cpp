#define _USE_MATH_DEFINES
#include "Phong.h"
#include "Ray.h"
#include "Scene.h"
#include <algorithm>
#include <random>

Phong::Phong(const Vector3 & kd, const Vector3 & ks, const float n, const Vector3 & v) :
m_kd(kd), m_ks(ks), m_n(n), m_v(v)
{
    m_v.normalize();
}

Phong::~Phong()
{
}

float Phong::BRDF(const Vector3& in, const Vector3& normal, const Vector3& out, const bool& isFront) const {
    if (dot(normal, out) < 0) return 0;
    else {
        double diffuse = (float)rand() / RAND_MAX;
        if (diffuse < (m_kd[0] + m_kd[1] + m_kd[2]) / 3) {
            return 1.0f / M_PI;
        }
        else {
            Vector3 R = 2.0 * dot(normal, in) * normal - in;
            return  (m_n+1) / (2.0 * M_PI) * pow(dot(R, out), m_n);
        }
    }
}

Vector3
Phong::shade(const Ray& ray, const HitInfo& hit, const Scene& scene, const bool& isFront) const
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
        float brdf = BRDF(rayLight.d, hit.N, -ray.d);
        if (brdf == 0) continue;
        if (scene.trace(hitLight, rayLight, 0.0001, l.length())) continue;

        // the inverse-squared falloff
        float falloff = l.length2();

        // normalize the light direction
        l /= sqrt(falloff);

        // get the diffuse component
        float nDotL = dot(hit.N, l);
        Vector3 result = pLight->color();

        L += std::max(0.0f, nDotL / falloff * pLight->wattage()*brdf) * result;
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

        float brdf = BRDF(rayLight.d, hit.N, -ray.d);
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

        // normalize the light direction
        l /= sqrt(falloff);

        // get the diffuse component
        float nDotL = dot(hit.N, l);
        Vector3 result = aLight->color();

        L += std::max(0.0f, dot(hitLight.N, -l))*
            std::max(0.0f, nDotL / falloff*
            aLight->wattage() / aLight->area()*brdf) * result / (vp.p);
    }
    
    // add the ambient component
    L += m_ka;
    
    return L;
}//*/

vec3pdf Phong::randReflect(const Vector3& in, const Vector3& normal, const bool& isFront) const{
    Vector3 out = 2 * dot(normal, in)*normal - in;
    return vec3pdf(out, 1);
}

vec3pdf Phong::randEmit(const Vector3& n) const {
    return vec3pdf(n, 1); // like a laserpointer?
}

Vector3 Phong::radiance(const Vector3& normal, const Vector3& direction) const {
    if (dot(normal, direction) < 0) return Vector3(0, 0, 0);
    return m_powerPerArea / (2.0*M_PI);
}

Vector3 Phong::sum_L_cosTheta_dOmega() const {
    return 2 * M_PI*radiance(Vector3(0, 0, 1), Vector3(0, 0, 1));
}
