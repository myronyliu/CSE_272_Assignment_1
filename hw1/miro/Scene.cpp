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
    HitInfo hit;
    if (!trace(hit, ray)) {
        return Vector3(0, 0, 0);
    }
    //Vector3 down(0, 0, -1);
    //if (2.0 - hit.P[2] < 0.00001 && hit.N == down && bounces>0) printf("ceiling hit on bounce\n");
    double rn = (double)(rand()) / RAND_MAX;
    double em = hit.material->emittance();
    //printf("%f%\r", em);
    if (em == 1.0 || rn < em) {
        Vector3 rad = hit.material->powerPerPatchPerSolidAngle(hit.N, ray.o - hit.P);
        if (bounces == 0) return (1.0 / em)*rad;
        else return Vector3(0, 0, 0);
    }
    rn = (double)rand() / RAND_MAX;
    vec3pdf vp;
    if (m_samplingHeuristic == 1 || rn < m_samplingHeuristic) { // sample BRDF
        vp = hit.material->randReflect(ray, hit);// pick BRDF weighted random direction
    }
    else { // sample AreaLight
        int numAL = m_areaLights.size();
        int randALind = std::fmin(floor((double)numAL* rand() / RAND_MAX), numAL - 1);
        AreaLight* randAL = m_areaLights[randALind];
        vp = randAL->randPt();
        Vector3 lightPt = vp.v;
        vp.v -= hit.P;
        vp.p *= vp.v.length2() / fabs(dot(vp.v.normalized(),randAL->normal(lightPt))); // convert PDF from 1/A to 1/SA
    }
    Vector3 newDir = vp.v;
    Ray newRay;
    newRay.o = hit.P;
    newRay.d = newDir;
    float brdf = hit.material->BRDF(ray, hit, newRay);
    float cos = dot(hit.N, newDir);
    Vector3 gather = hit.material->shade(ray, hit, *this); // gathered direct lighting
    return
        gather + (1.0 / (1.0 - em))/vp.p*brdf*cos*
        recursiveTrace_fromEye(newRay, bounces + 1, maxbounces);
}

void
Scene::pathtraceImage(Camera *cam, Image *img)
{
    HitInfo hitInfo;
    
    // loop over all pixels in the image
    for (int j = 0; j < img->height(); ++j){
        for (int i = 0; i < img->width(); ++i){
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

void Scene::tracePhoton(Camera *cam, vector<vector<Vector3>>& img, const Light& light, const raypdf& rp) {
    float w = img.size();
    float h = img[0].size();
    Vector3 pix;
    HitInfo hit, hitTry;
    Ray ray, newRay, rayToEye;
    ray = rp.r;
    float power = light.wattage() / rp.p; // Monte Carlo sampling so divide by PDF of the random ray
    for (int bounces = 0; bounces < m_maxBounces; bounces++) {
        if (!trace(hit, ray)) {
            //printf("bounce %i left scene\n",bounces);
            //printf(" %f %f %f\n", ray.d[0], ray.d[1], ray.d[2]);
            //printf(" %f %f %f\n", ray.o[0], ray.o[1], ray.o[2]);
            //std::cout << ray.d << std::endl;
            return; // ray left scene
        }double rn = (double)rand() / RAND_MAX;
        double em = hit.material->emittance();
        if (em == 1.0 || rn < em) {
            //printf("bounce %i absorbed\n", bounces);
            return; // photon was absorbed
        }
        vec3pdf vp = hit.material->randReflect(ray, hit); // pick BRDF weighted random direction
        Vector3 newDir = vp.v;
        newRay.o = hit.P;
        newRay.d = newDir;
        float brdf_toNew = hit.material->BRDF(ray, hit, newRay); // BRDF for bouncing to new random direction
        float cos_toNew = dot(hit.N, newDir); // cosine between new and current ray directions
        rayToEye.o = hit.P;
        rayToEye.d = (cam->eye() - hit.P).normalize();
        if (!trace(hitTry, rayToEye)) { // check if anything is occluding the eye from current hitpoint
            pix = cam->imgProject(hit.P, w, h); // find the pixel the onto which the current hitpoint projects
            int x = round(pix[0]);
            int y = round(pix[1]);
            if (pix[2]>0 && x > -1 && x<w && y>-1 && y < h) { // check that the pixel is within the viewing window
                //std::cout << "updating pixel" << std::endl;
                float brdf_toEye = hit.material->BRDF(ray, hit, rayToEye); // BRDF between current and shadow ray
                float cos_toEye = dot(hit.N, rayToEye.d); // cosine between shadow ray and current ray
                Vector3 gather = hit.material->shade(ray, hit, *this);
                img[x][y] += power*gather*brdf_toEye*cos_toEye;
            }
            else printf("impossible: %f %f %f \n,",pix[0],pix[1],pix[2]);
        }
        // update power and ray in anticipation for next bounce
        power *= brdf_toNew*cos_toNew;// / vp.p;
        ray = newRay;
    }
}


void
Scene::photontraceImage(Camera *cam, Image *img)
{
    int w = img->width();
    int h = img->height();
    vector<vector<Vector3>> im;
    for (int i = 0; i < w; i++) im.push_back(vector<Vector3>(h, Vector3(0, 0, 0)));
    // loop over all point lights
    // loop over all area lights
    for (unsigned int aL = 0; aL < m_areaLights.size(); aL++) {
        for (int p = 0; p < m_photonSamples; p++) {
            Light* light = m_areaLights[aL];
            raypdf rp = light->randRay();
            //printf("%f %f %f\n", rp.r.d[0], rp.r.d[1], rp.r.d[2]);
            //printf("%f %f %f\n", rp.r.o[0], rp.r.o[1], rp.r.o[2]);
            tracePhoton(cam, im, *light, rp);
            if (p % 100 == 0) {
                printf("Rendering Progress for AreaLight %i: %.3f%%\r", aL, p / float(m_photonSamples)*100.0f);
                fflush(stdout);
            }
        }
    }
    printf("Rendering Progress: 100.000%\n");
    debug("done Raytracing!\n");
    for (int i = 0; i < w; i++){
        for (int j = 0; j < h; j++){
            im[i][j] /= (float)m_photonSamples;
            im[i][j] /= cam->pixelCosine(i, j, w, h)*cam->pixelSolidAngle(i, j, w, h);
            //printf("%f %f %f\n", im[i][j][0], im[i][j][1], im[i][j][2]);
            img->setPixel(i, j, im[i][j]);
        }
    }
    img->draw();
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
