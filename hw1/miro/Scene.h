#ifndef CSE168_SCENE_H_INCLUDED
#define CSE168_SCENE_H_INCLUDED

#include "Miro.h"
#include "Object.h"
#include "PointLight.h"
#include "AreaLight.h"
#include "Light.h"
#include "BVH.h"
#include "RayPath.h"
#include "PhotonMap.h"

class Camera;
class Image;

struct LightPDF {
    Light* l;
    float p;
};

class Scene
{
public:
    void setEmittedPhotonsPerLight(std::vector<int> s) { m_emittedPhotonsPerLight = s; }
    void setSamplesPerPix(int i) { m_samplesPerPix = fmax(0, i); }
    void setBidiSamplesPerPix(int i) { m_bidiSamplesPerPix = fmax(0, i); }
    int samplesPerPix() { return m_samplesPerPix; }
    int bidiSamplesPerPix() { return m_bidiSamplesPerPix; }
    void setPhotonSamples(int i) { m_photonSamples = fmax(0, i); }
    int photonSamples() { return m_photonSamples; }
    void setMaxBounces(int i) { m_maxBounces = fmax(0, i); }
    int maxBounces() { return m_maxBounces; }
    void setMaxEyePaths(int i) { m_maxEyePaths = fmax(0, i); }
    void setMaxLightPaths(int i) { m_maxLightPaths = fmax(0, i); }
    int maxEyePaths() { return m_maxEyePaths; }
    int maxLightPaths() { return m_maxLightPaths; }
    void addObject(Object* pObj) { m_objects.push_back(pObj); }
    const Objects* objects() const { return &m_objects; }
    void setSamplingHeuristic(float p) { m_samplingHeuristic = fmin(fmax(0.0f, p), 1.0f); }
    float samplingHeuristic() { return m_samplingHeuristic; }

    void setPreview(bool preview) { m_preview = preview;  }
    bool preview(){ return m_preview; }

    void addPointLight(PointLight* pObj, int s = 1024) {
        m_lights.push_back(pObj);
        m_pointLights.push_back(pObj);
        m_emittedPhotonsPerLight.push_back(s);
    }
    const PointLights* pointLights() const {return &m_pointLights;}
    void addAreaLight(AreaLight* pObj, int s = 1024) {
        m_lights.push_back(pObj);
        m_objects.push_back(pObj);
        m_areaLights.push_back(pObj);
        m_emittedPhotonsPerLight.push_back(s);
    }
    const AreaLights* areaLights() const        { return &m_areaLights; }

    void preCalc();
    void openGL(Camera *cam);

    void raytraceImage(Camera *cam, Image *img);
    void pathtraceImage(Camera *cam, Image *img);
    void photontraceImage(Camera *cam, Image *img);
    void biditraceImage(Camera *cam, Image *img);
    void unifiedpathtraceImage(Camera *cam, Image *img);

    bool trace(HitInfo& minHit, const Ray& ray,
        float tMin = 0.0f, float tMax = MIRO_TMAX) const;
    bool trace(HitInfo& minHit, const Ray& ray, const Object* skip,
        float tMin = 0.0f, float tMax = MIRO_TMAX) const;
    Vector3 recursiveTrace_fromEye(const Ray& ray, int bounces, int maxbounces);
    // trace a ray through the scene and return an image with accumlated pixel values from that single photon
    void tracePhoton(Camera *cam, Image *img, const LightPDF& lightAndProb, const RayPDF& rayAndProb);

    LightPDF randLightByWattage();


    EyePath randEyePath(float i, float j, Camera* cam, Image* img, const int& bounces = -1);
    LightPath randLightPath(Light* light = NULL, const int& bounces = -1);
    void bounceRayPath(RayPath &, const int& paths);
    Vector3 bidiFlux(int i, int j, LightPath lightPath, EyePath eyePath);
    Vector3 uniFlux(const int& i, const int& j, const LightPath& lightPath, const EyePath& eyePath, PhotonMap* photonMap, const bool& explicitConnection = true, const int& nLightPaths = 1, const std::vector<PhotonDeposit>& photons = std::vector<PhotonDeposit>(0));
    Vector3 uniFluxDE(const int& j, const EyePath& eyePath, PhotonMap* photonMap, const int& nLightPaths = 1);

    std::pair<PhotonMap*, std::vector<LightPath*>> generatePhotonMap();
    std::pair<Vector3, Vector3> axisAlignedBounds();

protected:
    Objects m_objects;
    BVH m_bvh;
    int m_samplesPerPix = 100;
    int m_bidiSamplesPerPix = 50;
    int m_photonSamples = 100000000;
    int m_maxBounces = 20;
    int m_maxEyePaths = 1;
    int m_maxLightPaths = 0;
    int m_nGatheredPhotons = 32;
    float m_photonGatheringRadius = 0.2f; // radius for gathering photons in the vicinity

    std::vector<int> m_emittedPhotonsPerLight;
    Lights m_lights;
    PointLights m_pointLights;
    AreaLights m_areaLights;
    float m_samplingHeuristic = 0.5f; // probability of sampling BRDF. Complement is for sampling light
    bool m_preview = false;
};

extern Scene * g_scene;

#endif // CSE168_SCENE_H_INCLUDED
