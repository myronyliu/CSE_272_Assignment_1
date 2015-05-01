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
makeSpiralScene()
{
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;

    g_image->resize(512, 512);
    
    // set up the camera
    g_camera->setBGColor(Vector3(1.0f, 1.0f, 1.0f));
    g_camera->setEye(Vector3(-5, 2, 3));
    g_camera->setLookAt(Vector3(0, 0, 0));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(45);

    // create and place a point light source
    PointLight * light = new PointLight;
    light->setPosition(Vector3(-3, 15, 3));
    //light->setPosition(-Vector3(-5, 2, 3));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(1000);
    g_scene->addPointLight(light);

    // create a spiral of spheres
    Material* mat = new Lambert(Vector3(1.0f, 0.0f, 0.0f));
    const int maxI = 200;
    const float a = 0.15f;
    for (int i = 1; i < maxI; ++i)
    {
        float t = i/float(maxI);
        float theta = 4*PI*t;
        float r = a*theta;
        float x = r*cos(theta);
        float y = r*sin(theta);
        float z = 2*(2*PI*a - r);
        Sphere * sphere = new Sphere;
        sphere->setCenter(Vector3(x,y,z));
        sphere->setRadius(r/10);
        sphere->setMaterial(mat);
        g_scene->addObject(sphere);
    }
    
    // let objects do pre-calculations if needed
    g_scene->preCalc();
}

void
makeBunnyScene()
{
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;

    g_image->resize(128, 128);
    
    // set up the camera
    g_camera->setBGColor(Vector3(0.0f, 0.0f, 0.2f));
    g_camera->setEye(Vector3(-2, 3, 5));
    g_camera->setLookAt(Vector3(-.5, 1, 0));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(45);

    // create and place a point light source
    PointLight * light = new PointLight;
    light->setPosition(Vector3(-3, 15, 3));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(500);
    g_scene->addPointLight(light);

    Material* mat = new Lambert(Vector3(1.0f));

    TriangleMesh * bunny = new TriangleMesh;
    bunny->load("bunny.obj");
    
    // create all the triangles in the bunny mesh and add to the scene
    for (int i = 0; i < bunny->numTris(); ++i)
    {
        Triangle* t = new Triangle;
        t->setIndex(i);
        t->setMesh(bunny);
        t->setMaterial(mat); 
        g_scene->addObject(t);
    }
    
    // create the floor triangle
    TriangleMesh * floor = new TriangleMesh;
    floor->createSingleTriangle();
    floor->setV1(Vector3(  0, 0,  10));
    floor->setV2(Vector3( 10, 0, -10));
    floor->setV3(Vector3(-10, 0, -10));
    floor->setN1(Vector3(0, 1, 0));
    floor->setN2(Vector3(0, 1, 0));
    floor->setN3(Vector3(0, 1, 0));
    
    Triangle* t = new Triangle;
    t->setIndex(0);
    t->setMesh(floor);
    t->setMaterial(mat); 
    g_scene->addObject(t);
    
    // let objects do pre-calculations if needed
    g_scene->preCalc();
}

void
makeTestScene() {
    Material* mat = new Lambert(Vector3(1.0f, 0.0f, 0.0f));
    Material* matLight = new Lambert(Vector3(1.0, 1.0, 1.0));
    mat->setEmittance(0);
    mat->setEmitted(Vector3(0, 0, 0));
    matLight->setEmittance(1.0);
    matLight->setEmitted(Vector3(1, 1, 1));

    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;

    g_image->resize(128, 512);

    // set up the camera
    g_camera->setBGColor(Vector3(1.0f, 1.0f, 1.0f));
    g_camera->setEye(Vector3(7, 0, 1));
    g_camera->setLookAt(Vector3(0, 0, 1));
    g_camera->setUp(Vector3(0, 0, 1));
    g_camera->setFOV(70);

    // create and place a point light source
    PointLight * plight = new PointLight;
    plight->setPosition(Vector3(0, 0, 10));
    plight->setColor(Vector3(1, 1, 1));
    plight->setWattage(1000);
    //g_scene->addPointLight(plight);

    ParallelogramLight * light = new ParallelogramLight(Vector3(0, 0, 5), Vector3(0, 1, 0), Vector3(1, 0, 0), 0.5, 0.5);
    light->setWattage(100);
    light->setColor(Vector3(1, 0, 0));
    light->setMaterial(matLight);
    g_scene->addAreaLight(light);

    // create occluding spheres
    Sphere* top = new Sphere;
    Sphere* bot = new Sphere;
    top->setCenter(Vector3(0, 0, 1));
    bot->setCenter(Vector3(0, 0, -1));
    top->setRadius(1);
    bot->setRadius(1);
    top->setMaterial(mat);
    bot->setMaterial(mat);
    g_scene->addObject(top);
    g_scene->addObject(bot);

    ParallelogramLight* p = new ParallelogramLight(Vector3(0, 0, 0.5), Vector3(1, 0, 0), Vector3(0, 1, 0.2), 2.5, 2.5);
    //g_scene->addAreaLight(p);
    for (int i = 0; i < 1000; i++){
        Sphere * s = new Sphere;
        s->setCenter(p->randPt().v);
        printf("( %f, %f, %f )\n", s->center()[0], s->center()[1], s->center()[2]);
        s->setRadius(0.02);
        //g_scene->addObject(s);
    }
    
    /*HitInfo hit;
    Ray ray;
    hit.N = Vector3(0, 0, -1);
    for (int i = 0; i < 1000; i++) {
        Vector3 v = mat->randReflect(ray, hit);
        Sphere * sphere = new Sphere;
        sphere->setCenter(v);
        sphere->setRadius(0.02);
        sphere->setMaterial(mat);
        g_scene->addObject(sphere);
    }
    hit.N = Vector3(0, 0, 1);
    for (int i = 0; i < 1000; i++) {
        Vector3 v = mat->randReflect(ray, hit);
        Sphere * sphere = new Sphere;
        sphere->setCenter(v);
        sphere->setRadius(0.02);
        sphere->setMaterial(mat);
        g_scene->addObject(sphere);
    }
    */
    // let objects do pre-calculations if needed
    g_scene->preCalc();
}

void
makeWalls() {
}

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

    // create and place a point light source
    PointLight * plight = new PointLight;
    plight->setPosition(Vector3(0, 0, 1.99));
    plight->setColor(Vector3(1, 1, 1));
    plight->setWattage(1);
    //g_scene->addPointLight(plight);

    // create room geometry

    Material* mat = new Lambert(Vector3(1.0f, 1.0f, 1.0f));
    mat->setEmittance(0.0);
    mat->setEmitted(Vector3(0.0, 0.0, 0.0));

    Parallelogram * wall_F = new Parallelogram(Vector3(0, 1, 1), Vector3(1, 0, 0), Vector3(0, 0, 1), 1, 1); // far
    Parallelogram * wall_L = new Parallelogram(Vector3(-1, 0, 1), Vector3(0, 1, 0), Vector3(0, 0, 1), 1, 1); // left
    Parallelogram * wall_R = new Parallelogram(Vector3(1, 0, 1), Vector3(0, 0, 1), Vector3(0, 1, 0), 1, 1); // right
    Parallelogram * wall_T = new Parallelogram(Vector3(0, 0, 2), Vector3(0, 1, 0), Vector3(1, 0, 0), 1, 1); // top
    Parallelogram * wall_B = new Parallelogram(Vector3(0, 0, 0), Vector3(1, 0, 0), Vector3(0, 1, 0), 1, 1); // bottom

    Parallelogram * cover = new Parallelogram(Vector3(0, 0, 1.98), Vector3(0, 1, 0), Vector3(1, 0, 0), 0.1, 0.1);
    cover->disableBack();

    wall_T->setMaterial(mat);
    wall_B->setMaterial(mat);
    wall_L->setMaterial(mat);
    wall_R->setMaterial(mat);
    wall_F->setMaterial(mat);
    cover->setMaterial(mat);
    /*Parallelogram * sideN = new Parallelogram(Vector3(0, -0.5, 2.4), Vector3(1, 0, 0), Vector3(0, 0, 1), 0.5, 0.5); // near
    Parallelogram * sideF = new Parallelogram(Vector3(0, 0.5, 2.4), Vector3(1, 0, 0), Vector3(0, 0, 1), 0.5, 0.5); // far
    Parallelogram * sideL = new Parallelogram(Vector3(-0.5, 0, 2.4), Vector3(0, 1, 0), Vector3(0, 0, 1), 0.5, 0.5); // left
    Parallelogram * sideR = new Parallelogram(Vector3(0.5, 0, 2.4), Vector3(0, 0, 1), Vector3(0, 1, 0), 0.5, 0.5); // right
    sideN->setMaterial(mat);
    sideF->setMaterial(mat);
    sideL->setMaterial(mat);
    sideR->setMaterial(mat);
    g_scene->addObject(sideN);
    g_scene->addObject(sideF);
    g_scene->addObject(sideL);
    g_scene->addObject(sideR);*/

    ParallelogramLight * light = new ParallelogramLight(Vector3(0, 0, 1.98), Vector3(1, 0, 0), Vector3(0, 1, 0), 0.1, 0.1);
    light->flip(); cover->flip();
    light->setWattage(10000);
    light->setColor(Vector3(1, 1, 1));

    // add objects to scene
    g_scene->addObject(wall_B);
    g_scene->addObject(wall_F);
    g_scene->addObject(wall_L);
    g_scene->addObject(wall_R);
    g_scene->addObject(wall_T);
    g_scene->addObject(cover);
    g_scene->addAreaLight(light);


    Parallelogram* imgPlane = new Parallelogram;
    *imgPlane = g_camera->imagePlane(256, 256);
    //g_scene->addObject(imgPlane);

    for (int i = 0; i < 1000; i++){
        raypdf rr = light->randRay();
        Sphere* sphere = new Sphere;
        sphere->setRadius(0.001);
        sphere->setCenter(rr.r.o+rr.r.d);
        //g_scene->addObject(sphere);
    }

    // let objects do pre-calculations if needed
    g_scene->preCalc();
}


int
main(int argc, char*argv[])
{
    // create a scene
    //makeBunnyScene();
    //makeSpiralScene();
    //makeTestScene();
    makeRoomScene();

    MiroWindow miro(&argc, argv);
    miro.mainLoop();

    return 0; // never executed
}

