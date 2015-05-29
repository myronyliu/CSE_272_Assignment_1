#ifndef CSE168_RAYPATH_H_INCLUDED
#define CSE168_RAYPATH_H_INCLUDED

#include <vector>
#include "Ray.h"
#include "Light.h"

class RayPath
{
public:
    RayPath(Ray rayInit) : m_light(nullptr) {
        m_ray.push_back(rayInit);
    };

    // brdf is aligned with hit
    std::vector<HitInfo> m_hit;
    std::vector<float> m_brdf;

    // the following are all aligned. ray_i corresponds to the ray emitted off of hit_i
    std::vector<Ray> m_ray;
    std::vector<float> m_decay; // CUMULATIVE decay in the flux
    std::vector<float> m_prob;  // CUMULATIVE probability DENSITY (can be greater than 1)
    std::vector<float> m_cosB;
    std::vector<float> m_cosF;
    std::vector<float> m_length2;

    Light * m_light;
};

#endif // CSE168_RAYPATH_H_INCLUDED