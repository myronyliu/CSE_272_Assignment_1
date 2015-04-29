#define _USE_MATH_DEFINES
#include "Miro.h"
#include "Scene.h"
#include "Camera.h"
#include "Image.h"
#include "Console.h"

Scene * g_scene = 0;

void
Scene::openGL(Camera *cam)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    cam->drawGL();

    // draw objects
    for (size_t i = 0; i < m_objects.size(); ++i)
        m_objects[i]->renderGL();

    glutSwapBuffers();
}

void
Scene::preCalc()
{
    Objects::iterator it;
    for (it = m_objects.begin(); it != m_objects.end(); it++)
    {
        Object* pObject = *it;
        pObject->preCalc();
    }
    PointLights::iterator lit;
    for (lit = m_pointLights.begin(); lit != m_pointLights.end(); lit++)
    {
        PointLight* pLight = *lit;
        pLight->preCalc();
    }

    m_bvh.build(&m_objects);
}

void
Scene::raytraceImage(Camera *cam, Image *img)
{
    Ray ray;
    HitInfo hitInfo;
    Vector3 shadeResult;

    // loop over all pixels in the image
    for (int j = 0; j < img->height(); ++j)
    {
        for (int i = 0; i < img->width(); ++i)
        {
            ray = cam->eyeRay(i, j, img->width(), img->height());
            if (!trace(hitInfo, ray)) continue;
            shadeResult = hitInfo.material->shade(ray, hitInfo, *this);
            img->setPixel(i, j, shadeResult);
        }
        img->drawScanline(j);
        glFinish();
        printf("Rendering Progress: %.3f%%\r", j / float(img->height())*100.0f);
        fflush(stdout);
    }
    printf("Rendering Progress: 100.000%\n");
    debug("done Raytracing!\n");
}

Vector3 Scene::recursiveTrace(const Ray& ray, int bounces, int maxbounces) {
    if (bounces >= maxbounces) return Vector3(0, 0, 0);
    //printf("BOUNCE %i\n", bounces);
    HitInfo hit;
    if (!trace(hit, ray)) {
        //printf("ray left scene\n");
        return Vector3(0, 0, 0);
    }
    /*
    printf("  new hit point:   ( %f, %f, %f )\n", hit.P[0], hit.P[1], hit.P[2]);
    printf("  current ray pos: ( %f, %f, %f )\n", ray.o[0], ray.o[1], ray.o[2]);
    printf("  current ray dir: ( %f, %f, %f )\n", ray.d[0], ray.d[1], ray.d[2]);
    */

    double rn = (double)rand() / (double)RAND_MAX;
    double em = hit.material->getEmittance();
    if (rn <= em) {
        //printf("terminated in %i bounces\n", bounces);
        return (1.0 / em)*hit.material->getEmitted();
    }
    else {
        Vector3 newDir = hit.material->randReflect(ray, hit); // pick BRDF weighted random direction
        //printf("%f, %f, %f\n", hit.N[0], hit.N[1], hit.N[2]);
        Ray newRay = ray;
        newRay.o = hit.P;
        newRay.d = newDir;
        float brdf = hit.material->BRDF(ray, hit);
        float cos = dot(hit.N, newDir);
        //float cos = hitInfo.N[0];
        Vector3 dir = hit.material->shade(ray, hit, *this);
        //printf("  brdf: %f\n", brdf);
        //printf("  cos : %f\n", cos);
        //printf("gathered irradiance: ( %f, %f, %f )\n", dir[0], dir[1], dir[2]);
        //HitInfo newhit;
        return
            dir + (1.0 / (1.0 - em))*(2.0*M_PI)*brdf*cos*
            recursiveTrace(newRay, bounces + 1, maxbounces);
    }
}

void
Scene::pathtraceImage(Camera *cam, Image *img)
{
    Ray ray;
    HitInfo hitInfo;
    Vector3 shadeResult;
    Vector3 newDir;
    int samples = 10;
    
    // loop over all pixels in the image
    for (int j = 0; j < img->height(); ++j){
    //for (int j = img->height() / 4; j < img->height(); j += 2 * img->height()){
        for (int i = 0; i < img->width(); ++i){
        //for (int i = img->width() / 2; i < img->width(); i += 2 * img->width()){
            Vector3 pixSum = Vector3(0.0, 0.0, 0.0);
            ray = cam->eyeRay(i, j, img->width(), img->height());
            if (!trace(hitInfo, ray)) continue; // don't bother if the first shot misses
            for (int k = 0; k < samples; k++){
                pixSum += recursiveTrace(ray, 0, 10);
                Vector3 p = pixSum;
                //printf("( %f, %f, %f )\n", p[0], p[1], p[2]);
            }
            img->setPixel(i, j, pixSum / (double)samples);
        }
        img->drawScanline(j);
        glFinish();
        printf("Rendering Progress: %.3f%%\r", j / float(img->height())*100.0f);
        fflush(stdout);
    }

    printf("Rendering Progress: 100.000%\n");
    debug("done Raytracing!\n");
}


bool
Scene::trace(HitInfo& minHit, const Ray& ray, float tMin, float tMax) const
{
    return m_bvh.intersect(minHit, ray, tMin, tMax);
}

bool
Scene::trace(HitInfo& minHit, const Ray& ray, const Object* skip, float tMin, float tMax) const
{
    return m_bvh.intersect(minHit, ray, skip, tMin, tMax);
}
