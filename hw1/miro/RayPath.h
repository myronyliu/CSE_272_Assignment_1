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

    // the following are all aligned.
    // for LIGHTpath: ray_i corresponds to the ray emitted off of hit_i
    // for   EYEpath: ray_i corresponds to the ray incident on hit_i
    // Hence for LIGHTpath, there will be one more hit than rays
    // This is following LAFORTUNE's paper
    std::vector<Ray> m_ray;
    std::vector<float> m_cosB;
    std::vector<float> m_cosF;
    std::vector<float> m_length2;
    std::vector<float> m_decay; // CUMULATIVE decay in the flux
    std::vector<float> m_prob;  // CUMULATIVE probability DENSITY (can be greater than 1)
                                // eyePath.m_prob[k] is dimensionless (or 1/Steradians if you like)
                                // lightPath.m_prob[k] has units of 1/A
                                // cosB[k]/length2[k] is the conversion factor
    // TODO: maybe we should just remove m_decay. Since the first term is 1, we can just use m_decay[k]=m_prob[k]/m_prob[0]

    Light* m_light;
};

#endif // CSE168_RAYPATH_H_INCLUDED