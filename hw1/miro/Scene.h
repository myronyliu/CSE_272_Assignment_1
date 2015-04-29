#ifndef CSE168_SCENE_H_INCLUDED
#define CSE168_SCENE_H_INCLUDED

#include "Miro.h"
#include "Object.h"
#include "PointLight.h"
#include "AreaLight.h"
#include "Light.h"
#include "BVH.h"

class Camera;
class Image;

class Scene
{
public:
    void addObject(Object* pObj)        {m_objects.push_back(pObj);}
    const Objects* objects() const      {return &m_objects;}

    void addPointLight(PointLight* pObj)     {m_pointLights.push_back(pObj);}
    const PointLights* pointLights() const        {return &m_pointLights;}
    void addAreaLight(AreaLight* pObj) {
        m_objects.push_back(pObj);
        m_areaLights.push_back(pObj);
    }
    const AreaLights* areaLights() const        { return &m_areaLights; }

    void preCalc();
    void openGL(Camera *cam);

    void raytraceImage(Camera *cam, Image *img);
    void pathtraceImage(Camera *cam, Image *img);
    void biDitraceImage(Camera *cam, Image *img);
    bool trace(HitInfo& minHit, const Ray& ray,
        float tMin = 0.0f, float tMax = MIRO_TMAX) const;
    bool trace(HitInfo& minHit, const Ray& ray, const Object* skip,
        float tMin = 0.0f, float tMax = MIRO_TMAX) const;
    Vector3 recursiveTrace(HitInfo& hit, const Ray& ray, int bounces, int maxbounces);
    Vector3 recursiveTrace(const Ray& ray, int bounces, int maxbounces);

protected:
    Objects m_objects;
    BVH m_bvh;
    PointLights m_pointLights;
    AreaLights m_areaLights;
};

extern Scene * g_scene;

#endif // CSE168_SCENE_H_INCLUDED
