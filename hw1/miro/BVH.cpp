#include "BVH.h"
#include "Ray.h"
#include "Console.h"

void
BVH::build(Objects * objs)
{
    // construct the bounding volume hierarchy
    m_objects = objs;
}

void
BVH::build(AreaLights* aLights) {
    m_areaLights = aLights;
}


bool
BVH::intersect(HitInfo& minHit, const Ray& ray, float tMin, float tMax) const
{
    // Here you would need to traverse the BVH to perform ray-intersection
    // acceleration. For now we just intersect every object.

    bool hit = false;
    HitInfo tempMinHit;
    minHit.t = MIRO_TMAX;

    for (size_t i = 0; i < m_objects->size(); ++i)
    {
        if ((*m_objects)[i]->intersect(tempMinHit, ray, tMin, tMax))
        {
            hit = true;
            if (tempMinHit.t < minHit.t)
                minHit = tempMinHit;
        }
    }
    /*if (hit == false) return hit;
    if (minHit.object->material()->isInteracting() == false) {
        Ray newRay(minHit.P,minHit.object->randReflect(-ray.d, minHit.N, minHit.P).v); // ugh, what should i do with the pdf? right now, the only interface we have is deterministic so it's fine
        intersect(minHit, newRay, tMin, tMax);
    }//*/
    return hit;
}

bool
BVH::intersect(HitInfo& minHit, const Ray& ray, const Object* skip, float tMin, float tMax) const
{
    // Here you would need to traverse the BVH to perform ray-intersection
    // acceleration. For now we just intersect every object.

    bool hit = false;
    HitInfo tempMinHit;
    minHit.t = MIRO_TMAX;
    
    for (size_t i = 0; i < m_objects->size(); ++i)
    {
        if ((*m_objects)[i] == skip) continue;
        if ((*m_objects)[i]->intersect(tempMinHit, ray, tMin, tMax))
        {
            hit = true;
            if (tempMinHit.t < minHit.t)
                minHit = tempMinHit;
        }
    }
    /*if (hit == false) return hit;
    if (minHit.object->material()->isInteracting() == false) {
        Ray newRay(minHit.P, minHit.object->randReflect(-ray.d, minHit.N, minHit.P).v);
        intersect(minHit, newRay, tMin, tMax);
    }//*/
    return hit;
}