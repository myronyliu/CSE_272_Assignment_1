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
#include "MiroWindow.h"

void
makeRoomScene(){
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;
    g_image->resize(256, 256);

    // set up the camera
    g_camera->setBGColor(Vector3(1.0f, 1.0f, 1.0f));
    g_camera->setEye(Vector3(0, -4, 1));
    g_camera->setLookAt(Vector3(0, -3, 1));
    g_camera->setUp(Vector3(0, 0, 1));
    g_camera->setFOV(40);

    g_scene->setSamplesPerPix(100);
    g_scene->setMaxBounces(20);
    g_scene->setPhotonSamples(100000);

    // create room geometry

    Material* mat = new Lambert(Vector3(0.8f, 0.8f, 0.8f));
    Material* coverMat = new Lambert(Vector3(0.0f, 0.0f, 0.0f));
    
    Parallelogram * wall_F = new Parallelogram(Vector3(0, 1, 1), Vector3(1, 0, 0), Vector3(0, 0, 1), 1, 1); // far
    Parallelogram * wall_L = new Parallelogram(Vector3(-1, 0, 1), Vector3(0, 1, 0), Vector3(0, 0, 1), 1, 1); // left
    Parallelogram * wall_R = new Parallelogram(Vector3(1, 0, 1), Vector3(0, 0, 1), Vector3(0, 1, 0), 1, 1); // right
    Parallelogram * wall_T = new Parallelogram(Vector3(0, 0, 2), Vector3(0, 1, 0), Vector3(1, 0, 0), 1, 1); // top
    Parallelogram * wall_B = new Parallelogram(Vector3(0, 0, 0), Vector3(1, 0, 0), Vector3(0, 1, 0), 1, 1); // bottom

    //Parallelogram * wall_N = new Parallelogram(Vector3(0, -1, 1), Vector3(0, 0, 1), Vector3(1, 0, 0), 1, 1); // near
    //wall_N->setMaterial(mat);
    //g_scene->addObject(wall_N);

    Parallelogram * cover = new Parallelogram(Vector3(0, 0, 1.98), Vector3(0, 1, 0), Vector3(1, 0, 0), 0.1, 0.1);
    ParallelogramLight * light = new ParallelogramLight(Vector3(0, 0, 1.98), Vector3(1, 0, 0), Vector3(0, 1, 0), 0.1, 0.1);

    cover->disableBack();

    wall_T->setMaterial(mat);
    wall_B->setMaterial(mat);
    wall_L->setMaterial(mat);
    wall_R->setMaterial(mat);
    wall_F->setMaterial(mat);
    cover->setMaterial(coverMat);

    light->flip(); cover->flip();
    light->setWattage(2);

    // add objects to scene
    g_scene->addObject(wall_B);
    g_scene->addObject(wall_F);
    g_scene->addObject(wall_L);
    g_scene->addObject(wall_R);
    g_scene->addObject(wall_T);
    g_scene->addObject(cover);
    g_scene->addAreaLight(light);

    /*int w = g_image->width();
    int h = g_image->height();
    for (float i = 0; i < w; i+=16){
        for (float j = 0; j < h; j+=16){
            float sa = g_camera->pixelSolidAngle(i,j,w,h);
            Parallelogram* para = new Parallelogram(Vector3(2*i/w, 2*j/w, sa/2), Vector3(0, 0, 1), Vector3(0, 0, 1), sa, sa);
            g_scene->addObject(para);
        }
    }*/

    // let objects do pre-calculations if needed
    g_scene->preCalc();
}

int
main(int argc, char*argv[])
{
    // create a scene
    makeRoomScene();

    MiroWindow miro(&argc, argv);
    miro.mainLoop();

    return 0; // never executed
}
