#ifndef CSE168_BVH_H_INCLUDED
#define CSE168_BVH_H_INCLUDED

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

    bool intersect(HitInfo& result, const Ray& ray, const Object* skip,
        float tMin = 0.0f, float tMax = MIRO_TMAX) const;

    bool intersectLights(HitInfo& results, const Ray& ray,
        float tMin = 0.0f, float tMax = MIRO_TMAX) const;

protected:
    Objects* m_objects;
    AreaLights* m_areaLights;
};

#endif // CSE168_BVH_H_INCLUDED
