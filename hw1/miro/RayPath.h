#ifndef RAYPATH_H_INCLUDED
#define RAYPATH_H_INCLUDED

#include <vector>
#include "Ray.h"
#include "Light.h"

class RayPath {
public:
    RayPath() {}
    RayPath(Ray rayInit) { m_ray.push_back(rayInit); }

    std::vector<HitInfo> m_hit;
    std::vector<float> m_cosF;
    std::vector<float> m_cosB;
    std::vector<Vector3> m_estimator; // CUMULATIVE estimator (can also be amplified, since we do Russian Roulette)

    std::vector<Ray> m_ray;
    std::vector<float> m_length2;
    std::vector<float> m_prob;  // CUMULATIVE probability DENSITY in area space

};

class EyePath : public RayPath
{
public:
    EyePath() : RayPath() {}
    EyePath(Ray rayInit) : RayPath(rayInit) {}
};

class LightPath : public RayPath {
public:
    LightPath() : RayPath(), m_light(nullptr) {}
    LightPath(Ray rayInit) : RayPath(rayInit), m_light(nullptr) {}

    HitInfo hit(const int& i) {
        if (i < 0) return m_lightHit;
        else return m_hit[i];
    }

    float m_originProb; // the probability density of sampling the first point on the light (not including the emitted direction PDF, which is in m_prob[0])
    HitInfo m_lightHit; // the hitinfo on the light
    Light* m_light;
};

#endif // RAYPATH_H_INCLUDED
