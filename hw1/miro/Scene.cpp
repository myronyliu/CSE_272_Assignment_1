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

Vector3 Scene::recursiveTrace_fromEye(const Ray& ray, int bounces, int maxbounces) {
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

    double rn = (double)(1 + rand()) / (double)(1 + RAND_MAX);
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
        float brdf = hit.material->BRDF(ray, hit, newRay);
        float cos = dot(hit.N, newDir);
        //float cos = hitInfo.N[0];
        Vector3 dir = hit.material->shade(ray, hit, *this); // gathered direct lighting
        //printf("  brdf: %f\n", brdf);
        //printf("  cos : %f\n", cos);
        //printf("gathered irradiance: ( %f, %f, %f )\n", dir[0], dir[1], dir[2]);
        //HitInfo newhit;
        return
            dir + (1.0 / (1.0 - em))*(2.0*M_PI)*brdf*cos*
            recursiveTrace_fromEye(newRay, bounces + 1, maxbounces);
    }
}

void
Scene::pathtraceImage(Camera *cam, Image *img)
{
    HitInfo hitInfo;
    
    // loop over all pixels in the image
    for (int j = 0; j < img->height(); ++j){
    //for (int j = img->height() / 4; j < img->height(); j += 2 * img->height()){
        for (int i = 0; i < img->width(); ++i){
        //for (int i = img->width() / 2; i < img->width(); i += 2 * img->width()){
            Ray ray00 = cam->eyeRay((float)i - 0.5, (float)j - 0.5, img->width(), img->height());
            Ray ray01 = cam->eyeRay((float)i - 0.5, (float)j + 0.5, img->width(), img->height());
            Ray ray10 = cam->eyeRay((float)i + 0.5, (float)j - 0.5, img->width(), img->height());
            Ray ray11 = cam->eyeRay((float)i + 0.5, (float)j + 0.5, img->width(), img->height());
            if (!trace(hitInfo, ray00) &&
                !trace(hitInfo, ray01) &&
                !trace(hitInfo, ray10) &&
                !trace(hitInfo, ray11)) continue;
            Vector3 pixSum = Vector3(0.0, 0.0, 0.0);
            for (int k = 0; k < m_samplesPerPix; k++){
                Ray ray = cam->eyeRayJittered(i, j, img->width(), img->height());
                if (!trace(hitInfo, ray)) continue;
                pixSum += recursiveTrace_fromEye(ray, 0, m_maxBounces);
                //Vector3 p = pixSum;
                //printf("( %f, %f, %f )\n", p[0], p[1], p[2]);
            }
            img->setPixel(i, j, pixSum / (double)m_samplesPerPix);
        }
        img->drawScanline(j);
        glFinish();
        printf("Rendering Progress: %.3f%%\r", j / float(img->height())*100.0f);
        fflush(stdout);
    }

    printf("Rendering Progress: 100.000%\n");
    debug("done Raytracing!\n");
}

Splat Scene::recursiveTrace_fromLight(Camera *cam, Image *img, const Ray& ray, int bounces, int maxbounces) {
    Splat splat;
    HitInfo hit, hitImg;
    Parallelogram imgPlane = cam->imagePlane(img->width(),img->height());
    bool objIntersected = trace(hit, ray);
    if (imgPlane.intersect(hitImg,ray)) {
        if (!objIntersected || hit.t > hitImg.t) {
            std::vector<float> pix = cam->imgProject(hitImg.P, img->width(), img->height());
            splat.x = pix[0];
            splat.y = pix[1];
            splat.rad = hitImg.material->shade(ray, hitImg, *this);
        }
    }
    // otherwise ray is not going towards image plane
    if (bounces >= maxbounces) return Vector3(0, 0, 0);
    if (!objIntersected) return Vector3(0, 0, 0);
    double rn = (double)(1 + rand()) / (double)(1 + RAND_MAX);
    double em = hit.material->getEmittance();
    if (rn <= em) return (1.0 / em)*hit.material->getEmitted();
    Vector3 newDir = hit.material->randReflect(ray, hit); // pick BRDF weighted random direction
    Ray newRay = ray;
    newRay.o = hit.P;
    newRay.d = newDir;
    float brdf = hit.material->BRDF(ray, hit, newRay);
    float cos = dot(hit.N, newDir);
    Vector3 dir = hit.material->shade(ray, hit, *this);
    return
        dir + (1.0 / (1.0 - em))*(2.0*M_PI)*brdf*cos*
        recursiveTrace_fromLight(newRay, bounces + 1, maxbounces);
    }
}


void
Scene::photontraceImage(Camera *cam, Image *img)
{
    HitInfo hit;

    // loop over all point lights
    // loop over all area lights
    for (int aL = 0; aL < m_areaLights.size(); aL++) {
        for (int p = 0; p < m_photonSamples; p++) {
            Ray ray = m_areaLights[aL]->randRay();
            //recursiveTrace_fromLight(ray, 0, m_maxBounces);
        }
    }
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
