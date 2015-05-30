#include <math.h>
#include "Miro.h"
#include "Scene.h"
#include "Camera.h"
#include "Image.h"
#include "Console.h"

#include "ParallelogramLight.h"
#include "PointLight.h"
#include "Sphere.h"
#include "Parallelogram.h"
#include "TriangleMesh.h"
#include "Triangle.h"
#include "Lambert.h"
#include "Mirror.h"
#include "RefractiveInterface.h"
#include "MiroWindow.h"

void
makeRoomScene(){
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;
    g_image->resize(218,218);

    // set up the camera
    g_camera->setEye(Vector3(0, -4, 1));
    g_camera->setLookAt(Vector3(0, -3, 1));
    g_camera->setUp(Vector3(0, 0, 1));
    g_camera->setFOV(40);

    g_scene->setPreview(true);
    g_scene->setSamplesPerPix(256);
    g_scene->setBidiSamplesPerPix(4);
    g_scene->setMaxBounces(100);
    g_scene->setMaxEyePaths(1);
    g_scene->setMaxLightPaths(0);
    g_scene->setPhotonSamples(10000000);

    // create room geometry

    Lambert* mat = new Lambert(Vector3(1.0f, 1.0f, 1.0f));
    mat->setKd(0.8f);
    Mirror* mir = new Mirror(Vector3(1.0f, 1.0f, 1.0f));
    mir->setKs(0.8f);
    Lambert* coverMat = new Lambert(Vector3(0.0f, 0.0f, 0.0f));
    coverMat->setKd(0.0f);
    RefractiveInterface* waterMat = new RefractiveInterface(Vector3(1.0f, 1.0f, 1.0f), Vector3(1.0f, 1.0f, 1.0f));
    waterMat->setRefractiveIndexFront(1.0f);
    waterMat->setRefractiveIndexBack(1.5f);

    Parallelogram * wall_F = new Parallelogram(Vector3(0, 1, 1), Vector3(1, 0, 0), Vector3(0, 0, 1), 1, 1); // far
    Parallelogram * wall_L = new Parallelogram(Vector3(-1, 0, 1), Vector3(0, 1, 0), Vector3(0, 0, 1), 1, 1); // left
    Parallelogram * wall_R = new Parallelogram(Vector3(1, 0, 1), Vector3(0, 0, 1), Vector3(0, 1, 0), 1, 1); // right
    Parallelogram * wall_T = new Parallelogram(Vector3(0, 0, 2), Vector3(0, 1, 0), Vector3(1, 0, 0), 1, 1); // top
    Parallelogram * wall_B = new Parallelogram(Vector3(0, 0, 0), Vector3(1, 0, 0), Vector3(0, 1, 0), 1, 1); // bottom
    Parallelogram * water_T = new Parallelogram(Vector3(0, 0, 1), Vector3(1, 0, 0), Vector3(0, 1, 0), 1, 1);
    Parallelogram * water_N = new Parallelogram(Vector3(0, -1, 0.5), Vector3(1, 0, 0), Vector3(0, 0, 1), 1, 0.5);

    ParallelogramLight * light = new ParallelogramLight(Vector3(0, 0, 1.98), Vector3(1, 0, 0), Vector3(0, 1, 0), 0.1, 0.1);
    Parallelogram * cover = new Parallelogram(Vector3(0, 0, 1.98), Vector3(0, 1, 0), Vector3(1, 0, 0), 0.1, 0.1);
    cover->disableBack();
    
    wall_T->setMaterial(mat);
    wall_B->setMaterial(mat);
    wall_L->setMaterial(mat);
    wall_R->setMaterial(mat);
    wall_F->setMaterial(mir);
    cover->setMaterial(coverMat);
    water_T->setMaterial(waterMat);
    water_N->setMaterial(waterMat);
    

    //light->flip(); cover->flip(); light->setWattage(2);

    // add objects to scene
    g_scene->addObject(wall_B);
    g_scene->addObject(wall_F);
    g_scene->addObject(wall_L);
    g_scene->addObject(wall_R);
    g_scene->addObject(wall_T);
    g_scene->addObject(cover);
    g_scene->addAreaLight(light, 1000000);
    //g_scene->addObject(water_T);
    //g_scene->addObject(water_N);

    // let objects do pre-calculations if needed
    g_scene->preCalc();

    /*PhotonMap* photonMap = g_scene->generatePhotonMapTest();
    std::vector<PhotonDeposit> photons = photonMap->getPhotons();
    for (int i = 0; i < photons.size(); i++) {
        Sphere* sphere = new Sphere();
        sphere->setCenter(photons[i].m_location);
        sphere->setRadius(0.1);
        g_scene->addObject(sphere);
    }
    float eps = 1;
    int n = 3;
    std::vector<PhotonDeposit> nearbyPhotons;
    nearbyPhotons = photonMap->getNearestPhotons(Vector3(-eps, -eps, 1 - eps), n);
    nearbyPhotons = photonMap->getNearestPhotons(Vector3(-eps, -eps, 1 + eps), n);
    nearbyPhotons = photonMap->getNearestPhotons(Vector3(-eps, eps, 1 - eps), n);
    nearbyPhotons = photonMap->getNearestPhotons(Vector3(-eps, eps, 1 + eps), n);
    nearbyPhotons = photonMap->getNearestPhotons(Vector3(eps, -eps, 1 - eps), n);
    nearbyPhotons = photonMap->getNearestPhotons(Vector3(eps, -eps, 1 + eps), n);
    nearbyPhotons = photonMap->getNearestPhotons(Vector3(eps, eps, 1 - eps), n);
    nearbyPhotons = photonMap->getNearestPhotons(Vector3(eps, eps, 1 + eps), n);
    //*/
}

int
main(int argc, char*argv[])
{
    // create a scene
    makeRoomScene();
    //unsigned int cw;
    //_controlfp_s(&cw, 0, 0);
    //cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_ZERODIVIDE|
    //        EM_DENORMAL|EM_INVALID);
    //unsigned int cwOriginal;
    //_controlfp_s(&cwOriginal,cw, _MCW_EM);
    MiroWindow miro(&argc, argv);
    miro.mainLoop();

    return 0; // never executed
}
