#ifndef CSE168_RAYPATH_H_INCLUDED
#define CSE168_RAYPATH_H_INCLUDED

#include "Ray.h"
#include <vector>

class RayPath
{
public:
    RayPath(Ray rayInit) : m_rayInit(rayInit) { m_rays.push_back(rayInit); };

    Ray m_rayInit;
    Vector3 m_normalInit;
    std::vector<Ray> m_rays;
    std::vector<HitInfo> m_hits;
    std::vector<float> m_probs;
};

#endif // CSE168_RAYPATH_H_INCLUDED