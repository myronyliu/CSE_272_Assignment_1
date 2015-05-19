#ifndef CSE168_RAYPATH_H_INCLUDED
#define CSE168_RAYPATH_H_INCLUDED

#include <vector>
#include "Ray.h"
#include "Light.h"

class RayPath
{
public:
    RayPath(Ray rayInit) : m_rayInit(rayInit), m_light(nullptr) {
        m_rays.push_back(rayInit);
    };

    Ray m_rayInit;
    std::vector<Ray> m_rays;
    std::vector<HitInfo> m_hits;
    std::vector<float> m_fluxDecay;
    std::vector<float> m_probs;
    std::vector<float> m_brdfs;

    Light * m_light;
};

#endif // CSE168_RAYPATH_H_INCLUDED