#ifndef CSE168_SCENE_H_INCLUDED
#define CSE168_SCENE_H_INCLUDED

#include "Miro.h"
#include "Object.h"
#include "PointLight.h"
#include "AreaLight.h"
#include "Light.h"
#include "BVH.h"

struct Splat{
    float x = 0 ;
    float y = 0 ;
    Vector3 rad = Vector3(0,0,0);
};

class Camera;
class Image;

class Scene
{
public:
    void setSamplesPerPix(int i) { m_samplesPerPix = i; }
    int samplesPerPix() { return m_samplesPerPix; }
    void setPhotonSamples(int i) { m_photonSamples = i; }
    int photonSamples() { return m_photonSamples; }
    void setMaxBounces(int i) { m_maxBounces = i; }
    int maxBounces() { return m_maxBounces; }
    void addObject(Object* pObj)        {m_objects.push_back(pObj);}
    const Objects* objects() const      {return &m_objects;}

    void addPointLight(PointLight* pObj) {m_pointLights.push_back(pObj);}
    const PointLights* pointLights() const {return &m_pointLights;}
    void addAreaLight(AreaLight* pObj) {
        m_objects.push_back(pObj);
        m_areaLights.push_back(pObj);
    }
    const AreaLights* areaLights() const        { return &m_areaLights; }

    void preCalc();
    void openGL(Camera *cam);

    void raytraceImage(Camera *cam, Image *img);
    void pathtraceImage(Camera *cam, Image *img);
    void photontraceImage(Camera *cam, Image *img);
    void biditraceImage(Camera *cam, Image *img);
    bool trace(HitInfo& minHit, const Ray& ray,
        float tMin = 0.0f, float tMax = MIRO_TMAX) const;
    bool trace(HitInfo& minHit, const Ray& ray, const Object* skip,
        float tMin = 0.0f, float tMax = MIRO_TMAX) const;
    Vector3 recursiveTrace_fromEye(const Ray& ray, int bounces, int maxbounces);
    Splat recursiveTrace_fromLight(Camera* cam, Image *img, const Ray& ray, int bounces, int maxbounces);

protected:
    Objects m_objects;
    BVH m_bvh;
    int m_samplesPerPix = 1000;
    int m_photonSamples = 100000;
    int m_maxBounces = 20;
    PointLights m_pointLights;
    AreaLights m_areaLights;
};

extern Scene * g_scene;

#endif // CSE168_SCENE_H_INCLUDED
