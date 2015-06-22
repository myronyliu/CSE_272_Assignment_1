#ifndef BVH_H_INCLUDED
#define BVH_H_INCLUDED

#include "Miro.h"
#include "Object.h"
#include "AreaLight.h"

class BVH
{
public:
    void build(Objects* objs);
    void build(AreaLights* alights);

    bool intersect(HitInfo& result, const Ray& ray,
        float tMin = 0.0f, float tMax = MIRO_TMAX) const;

    // Same as above but skips an Object
    // Skipping is used in shade(...) for instance.
    // See Lambert.cpp where we skip the current areaLight of interest
    bool intersect(HitInfo& result, const Ray& ray, const Object* skip,
        float tMin = 0.0f, float tMax = MIRO_TMAX) const;

protected:
    Objects* m_objects;
    AreaLights* m_areaLights;
};

#endif // BVH_H_INCLUDED
