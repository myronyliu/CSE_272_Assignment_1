#define _USE_MATH_DEFINES
#include "Miro.h"
#include "Scene.h"
#include "Camera.h"
#include "Image.h"
#include "Console.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <ppl.h>

#include <algorithm>

using namespace std;
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
            shadeResult = hitInfo.object->material()->shade(ray, hitInfo, *this);
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
    double rn = (double)(rand() / RAND_MAX);
    double em = hit.object->material()->emittance();
    if (em == 1.0 || rn < em) {
        Vector3 rad = hit.object->material()->radiance(hit.N, ray.o - hit.P);
        if (bounces == 0) return (1.0 / em)*rad;
        else return Vector3(0, 0, 0);
    }
    rn = (double)rand() / RAND_MAX;
    vec3pdf vp;
    if (m_samplingHeuristic == 1 || rn < m_samplingHeuristic) { // sample BRDF
        vp = hit.object->randReflect(-ray.d, hit.N, hit.P); // pick BRDF weighted random direction
    }
    else { // sample AreaLight
        int numAL = m_areaLights.size();
        int randALind = std::fmin(floor((double)numAL* rand() / RAND_MAX), numAL - 1);
        AreaLight* randAL = m_areaLights[randALind];
        vp = randAL->randPt();
        Vector3 lightPt = vp.v;
        vp.v -= hit.P;
        vp.p *= vp.v.length2() / fabs(dot(vp.v.normalized(), randAL->normal(lightPt))); // convert PDF from 1/A to 1/SA
    }
    Vector3 newDir = vp.v.normalized();
    Ray newRay;
    newRay.o = hit.P;
    newRay.d = newDir;
    float brdf = hit.object->BRDF(-ray.d, hit.N, newRay.d, hit.P);
    float cos = fabs(dot(hit.N, newDir)); // changed this to fabs for transmissible materials such as RefractiveInterface
    Vector3 gather = hit.object->shade(ray, hit, *this, hit.P); // gathered direct lighting
    if (hit.object->material()->isInteracting() == false) return gather + (1.0 / (1.0 - em)) / vp.p*brdf*cos*recursiveTrace_fromEye(newRay, bounces, maxbounces);
    else return gather + (1.0 / (1.0 - em)) / vp.p*brdf*cos*recursiveTrace_fromEye(newRay, bounces + 1, maxbounces);
}

void
Scene::pathtraceImage(Camera *cam, Image *img)
{
    HitInfo hitInfo;
    
    std::ofstream plotfile;
    plotfile.open("pathtraceplot.txt");
    
    int integrationStart = glutGet(GLUT_ELAPSED_TIME);
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
                //if (dot(hitInfo.N, Vector3(0, 0, -1)) > 0) { pixSum += Vector3(1.0, 0.0, 0.0);  }
                pixSum += recursiveTrace_fromEye(ray, 0, m_maxBounces) / (double)m_samplesPerPix;
                if (j == img->width() / 2 && i == img->height() / 2){
                    plotfile << pixSum[0] / (k + 1) * m_samplesPerPix << std::endl;
                }
            }
            img->setPixel(i, j, pixSum);
        }
        if (preview())
        {
            img->drawScanline(j);
        }
        glFinish();
        printf("Rendering Progress: %.3f%%\r", j / float(img->height())*100.0f);
        fflush(stdout);
    }
    if (!preview())
    {
        img->draw();
    }
    glFinish();

    printf("Rendering Progress: 100.000%\n");
    debug("done Raytracing!\n");
    int integrationEnd = glutGet(GLUT_ELAPSED_TIME);
    std::cout << "Rendering took " << ((integrationEnd - integrationStart) / 1000.0f) << "s" << std::endl;

    debug("done Raytracing!\n");
    plotfile.close();
}

void Scene::tracePhoton(Camera *cam, Image *img, const LightPDF& lp, const RayPDF& rp) {
    float w = img->width();
    float h = img->height();
    Light* light = lp.l;
    Vector3 pix;
    HitInfo hit, tryHitEye;
    Ray rayIn, rayOut, rayToEye;
    rayOut = rp.m_ray;
    Vector3 power = light->wattage()*light->color() / lp.p;
    //Vector3 power = light.wattage()*light.color()*dot(rayOut.d,light.normal(rayOut.o)) / rp.p / light.area(); // Monte Carlo sampling so divide by PDF of the random ray
    // The following is for the initial emission (to make the light visible)
    rayToEye.o = rayOut.o;
    rayToEye.d = (cam->eye() - rayToEye.o).normalize();
    // Add direct light to pixel
    if (!trace(tryHitEye, rayToEye) && dot(rayToEye.d, light->normal(rayToEye.o)) > 0) { // check if anything is occluding the eye from current hitpoint
        pix = cam->imgProject(rayToEye.o, w, h); // find the pixel the onto which the current hitpoint projects
        int x = round(pix[0]);
        int y = round(pix[1]);
        if (pix[2]>0 && x > -1 && x<w && y>-1 && y < h) { // check that the pixel is within the viewing window
            Vector3 initValue = light->material()->radiance(light->normal(rayToEye.o), rayToEye.d)*light->area();
            initValue *= dot(rayToEye.d, light->normal(rayToEye.o));
            initValue *= cam->pixelCosine(pix[0], pix[1], w, h);
            initValue /= light->material()->sum_L_cosTheta_dOmega();
            initValue /= (rayToEye.o - cam->eye()).length2();
            img->setPixel(x, y, img->getPixel(x, y) + power*initValue);
        }
    }
    // now for the bounces
    rayIn = rayOut;
    for (int bounces = 0; bounces < m_maxBounces; bounces++) {
        if (!trace(hit, rayIn)) {
            return; // ray left scene
        }
        double rn = (double)rand() / RAND_MAX;
        double em = hit.object->material()->emittance();
        rayToEye.o = hit.P;
        rayToEye.d = (cam->eye() - hit.P).normalize();
        float brdf = hit.object->material()->BRDF(-rayIn.d, hit.N, rayToEye.d); // BRDF between eye-ray and incoming-ray
        if (!trace(tryHitEye, rayToEye)) { // check if anything is occluding the eye from current hitpoint
            pix = cam->imgProject(hit.P, w, h); // find the pixel the onto which the current hitpoint projects
            int x = round(pix[0]);
            int y = round(pix[1]);
            //img->setPixel(x, y, Vector3(1.0, 0.0, 0.0)); return;
            if (dot(Vector3(0, 0, -1), hit.N) != 0) {}
            if (pix[2]>0 && x >= 0 && x<w && y >= 0 && y < h) { // check that the pixel is within the viewing window
                float cos0 = dot(hit.N, rayToEye.d);
                float lengthSqr = (cam->eye() - hit.P).length2();
                float cosAlpha = cam->pixelCosine(pix[0], pix[1], w, h);
                if (em == 1.0 || rn < 0.2) {
                    img->setPixel(x, y, img->getPixel(x, y) + (1 / 0.2) * power*brdf*cos0*cosAlpha / lengthSqr);
                    return; // photon was absorbed
                } // otherwise photon will be reflected
                else
                {
                    img->setPixel(x, y, img->getPixel(x, y) + (1 / 0.8) * power*brdf*cos0*cosAlpha / lengthSqr);
                }
            }
        }
        vec3pdf vp = hit.object->material()->randReflect(-rayIn.d, hit.N); // pick BRDF.cos weighted random direction
        Vector3 dirOut = vp.v;
        rayOut.o = hit.P;
        rayOut.d = dirOut;
        rayIn = rayOut;
        power = power * brdf;
    }
}


void
Scene::photontraceImage(Camera *cam, Image *img)
{
    int w = img->width();
    int h = img->height();
    std::ofstream plotfile;
    plotfile.open("photontraceplot.txt");
    
    int integrationStart = glutGet(GLUT_ELAPSED_TIME);
    for (int p = 0; p < m_photonSamples; p++) { // shoot a photon...
        LightPDF lp = randLightByWattage(); // ... off of a random light (I don't think we need the PDF here)
        Light* light = lp.l;
        RayPDF rp = light->randRay();
        tracePhoton(cam, img, lp, rp);
        if (p % 100 == 0) {
            printf("Rendering Progress: %.3f%%\r", p / float(m_photonSamples)*100.0f);
            fflush(stdout);
        }
        if (((float)p / m_photonSamples)<0.8 && p>0 && p % (m_photonSamples / 3) == 0)
        {
            img->draw();
        }
        if (p % 1000000 == 0) {
            plotfile << img->getPixel(w/2, h/2)[0] *m_photonSamples << std::endl;
    }
    }
    for (int i = 0; i < w; i++){
        for (int j = 0; j < h; j++){
            Vector3 pix = img->getPixel(i, j);
            img->setPixel(i, j,
                pix / (cam->pixelCosine(i, j, w, h) * cam->pixelSolidAngle(i, j, w, h) * m_photonSamples));
        }
    }
    img->draw();
    glFinish();
    printf("Rendering Progress: 100.000%\n");

    int integrationEnd = glutGet(GLUT_ELAPSED_TIME);
    std::cout << "Rendering took " << ((integrationEnd - integrationStart) / 1000.0f) << "s" << std::endl;

    debug("done Photontracing!\n");
    plotfile.close();
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

LightPDF Scene::randLightByWattage() {
    LightPDF rlbw;
    int n = m_lights.size();
    vector<float> wattageConcattage(n + 1, 0);
    for (int i = 0; i < n; i++){
        float w = m_lights[i]->wattage();
        wattageConcattage[i + 1] = wattageConcattage[i] + w;
    }
    float r = wattageConcattage[n] * (float)rand() / RAND_MAX;
    for (int i = 0; i < n; i++){ // find the interval in which r lies and return the light along with PDF;
        if (r < wattageConcattage[i] || r > wattageConcattage[i + 1]) continue;
        rlbw.l = m_lights[i];
        rlbw.p = m_lights[i]->wattage() / wattageConcattage[n];
        return rlbw;
    }
}


void
Scene::biditraceImage(Camera *cam, Image *img)
{
    int w = img->width();
    int h = img->height();
    HitInfo hitInfo;

    std::ofstream plotfile;
    plotfile.open("biditraceplot.txt");

    int integrationStart = glutGet(GLUT_ELAPSED_TIME);

    //for (int y = h-1; y >-1; y--)
    for (int y = 0; y < h; y++)
    {
        for (int x = 0; x < w; x++)
        {
            Ray ray00 = cam->eyeRay((float)x - 0.5, (float)y - 0.5, img->width(), img->height());
            Ray ray01 = cam->eyeRay((float)x - 0.5, (float)y + 0.5, img->width(), img->height());
            Ray ray10 = cam->eyeRay((float)x + 0.5, (float)y - 0.5, img->width(), img->height());
            Ray ray11 = cam->eyeRay((float)x + 0.5, (float)y + 0.5, img->width(), img->height());
            if (!trace(hitInfo, ray00) &&
                !trace(hitInfo, ray01) &&
                !trace(hitInfo, ray10) &&
                !trace(hitInfo, ray11)) continue;
            Vector3 fluxSumOverSamples(0, 0, 0);

            Concurrency::parallel_for(0, m_bidiSamplesPerPix, [&](int k) {
                EyePath eyePath = randEyePath(x, y, cam, img);
                if (eyePath.m_hit.size() == 0) return;
                LightPath lightPath = randLightPath();
                Vector3 fluxSum = estimateFlux(0, 0, lightPath, eyePath);
                for (int i = 0; i <= lightPath.m_hit.size(); i++) {
                    for (int j = 1; j <= eyePath.m_hit.size(); j++) {
                        fluxSum += estimateFlux(i, j, lightPath, eyePath);
                    }
                }
                fluxSumOverSamples += fluxSum;
                if (y == h / 2 && x == w / 2){
                    plotfile << fluxSum[0] / (k + 1) / M_PI / 0.04 << std::endl;
                }
            });
            img->setPixel(x, y, fluxSumOverSamples / m_bidiSamplesPerPix);
        }
        if (preview())
        {
            img->drawScanline(y);
        }
        glFinish();
        printf("Rendering Progress: %.3f%%\r", y / float(img->height())*100.0f);
        fflush(stdout);
    }
    if (!preview())
    {
        img->draw();
    }
    glFinish();

    printf("Rendering Progress: 100.000%\n");

    int integrationEnd = glutGet(GLUT_ELAPSED_TIME);
    std::cout << "Rendering took " << ((integrationEnd - integrationStart) / 1000.0f) << "s" << std::endl;

    debug("Done Bidi Pathtracing!\n");
    plotfile.close();

}

EyePath Scene::randEyePath(float x, float y, Camera* cam, Image* img) {
    EyePath eyePath(cam->eyeRay(x, y, img->width(), img->height()));
    eyePath.m_prob.push_back(1);
    bounceRayPath(eyePath, m_maxEyePaths);
    return eyePath;
}
LightPath Scene::randLightPath() {
    HitInfo hit;
    LightPDF lp = randLightByWattage();
    Light* light = lp.l;

    RayPDF rp = light->randRay();
    Ray rayInit = rp.m_ray;
    HitInfo lightHit(0.0f, rayInit.o, light->normal(rayInit.o), light);
    LightPath lightPath(rayInit);
    lightPath.m_prob.push_back(rp.m_dProb);

    // lightPath specific values
    lightPath.m_light = light;
    lightPath.m_lightHit = lightHit;
    lightPath.m_originProb = rp.m_oProb;

    bounceRayPath(lightPath, m_maxLightPaths);
    return lightPath;
}

void Scene::bounceRayPath(RayPath & raypath, const int& paths) {
    if (paths < 1) return;
    HitInfo hit;
    Ray rayInit = raypath.m_ray[0];
    if (!trace(hit, rayInit)) return;
    raypath.m_hit.push_back(hit);
    raypath.m_cosB.push_back(std::max(0.0f, dot(hit.N, -rayInit.d)));
    raypath.m_length2.push_back((hit.P - rayInit.o).length2());

    int bounce = 1;
    while (bounce < paths)
    {
        Ray lastRay = raypath.m_ray.back();
        HitInfo lastHit = raypath.m_hit.back();

        vec3pdf vp = lastHit.object->material()->randReflect(-lastRay.d, lastHit.N);
        Ray newRay(lastHit.P, vp.v);
        if (!trace(hit, newRay)) return;
        float brdf = lastHit.object->material()->BRDF(-lastRay.d, lastHit.N, vp.v);
        float cos = std::max(0.0f, dot(lastHit.N, newRay.d));
        float cosPrime = std::max(0.0f, dot(hit.N, -newRay.d));
        float distSqr = (hit.P - lastHit.P).length2();

        raypath.m_hit.push_back(hit);
        raypath.m_brdf.push_back(brdf);
        raypath.m_cosF.push_back(cos);
        raypath.m_cosB.push_back(cosPrime);

        raypath.m_ray.push_back(newRay);
        raypath.m_length2.push_back(distSqr);
        raypath.m_prob.push_back(brdf*cos*raypath.m_prob.back());
        float decay = brdf*cos;
        if (raypath.m_decay.size() != 0) decay*=raypath.m_decay.back();
        raypath.m_decay.push_back(decay);

        bounce++;
    }
}

Vector3 Scene::estimateFlux(int i, int j, LightPath lightPath, EyePath eyePath) {
    if (i == 0 && j == 0) return eyePath.m_hit[0].object->material()->radiance(eyePath.m_hit[0].N, -eyePath.m_ray[0].d);
    HitInfo hit;
    float intersectEpsilon = 0.00001;
    // The following is for the explicit connection
    HitInfo hit_E = eyePath.m_hit[j - 1];
    HitInfo hit_L;
    if (i == 0) hit_L = lightPath.m_lightHit;
    else hit_L = lightPath.m_hit[i - 1];
    Material* mat_E = hit_E.object->material();
    Material* mat_L = hit_L.object->material();
    Ray shadow_LtoE(hit_L.P, (hit_E.P - hit_L.P).normalize()); // shadow ray from Light path to Eye path
    float shadowLength2 = (hit_E.P - hit_L.P).length2();
    if (trace(hit, shadow_LtoE, intersectEpsilon, sqrt(shadowLength2) - intersectEpsilon)) return Vector3(0, 0, 0);
    if (i == 0) {
        if (!lightPath.m_light->intersect(hit, Ray(hit_E.P, (lightPath.m_lightHit.P - hit_E.P).normalize()))) return Vector3(0, 0, 0);
    }
    float cosF_L = std::max(0.0f, dot(shadow_LtoE.d, hit_L.N));
    if (cosF_L == 0) return Vector3(0, 0, 0);
    float cosB_E = std::max(0.0f, dot(-shadow_LtoE.d, hit_E.N));
    if (cosB_E == 0) return Vector3(0, 0, 0);
    float brdf_E = mat_E->BRDF(-shadow_LtoE.d, hit_E.N, -eyePath.m_ray[j - 1].d);
    if (brdf_E == 0) return Vector3(0, 0, 0);
    float brdf_L = 1;
    if (i > 0) brdf_L = mat_L->BRDF(shadow_LtoE.d, hit_L.N, -lightPath.m_ray[i - 1].d);
    if (brdf_L == 0) return Vector3(0, 0, 0);
    float dProb_EtoL = cosB_E*brdf_E;
    float dProb_LtoE = cosF_L*brdf_L;
    ///////////////////////////////////////////////////////////////////////////////////////////////
    vector<float> probF(i + j + 1); // forward from light to eye (excludes the const emission probability, since it's just a constant)
    vector<float> probB(i + j + 1); // backward from eye to light
    vector<float> length2(i + j + 1);
    vector<float> cosF(i + j);
    vector<float> cosB(i + j);
    vector<float> brdf(i + j);
    for (int k = 0; k < i; k++) {
        probF[k] = lightPath.m_prob[k];
        length2[k] = lightPath.m_length2[k];
        cosB[k] = lightPath.m_cosB[k];
        if (k == i - 1) break;
        cosF[k] = lightPath.m_cosF[k];
        brdf[k] = lightPath.m_brdf[k];
    }
    if (i > 0) {
        cosF[i - 1] = cosF_L;
        brdf[i - 1] = brdf_L;
    }
    length2[i] = shadowLength2;
    probB[i] = dProb_EtoL*eyePath.m_prob[j - 1];
    probF[i] = dProb_LtoE;
    if (i > 0) probF[i] *= probF[i - 1];
    cosF[i] = eyePath.m_cosB[j - 1];
    cosB[i] = cosB_E;
    brdf[i] = brdf_E;
    for (int k = i + 1; k < i + j; k++) {
        int u = i + j - k;
        cosF[k] = eyePath.m_cosB[u-1]; // note the swap
        cosB[k] = eyePath.m_cosF[u-1];
        brdf[k] = eyePath.m_brdf[u-1];
        length2[k] = eyePath.m_length2[u];
        probB[k] = eyePath.m_prob[u];
        probF[k] = probF[k - 1] * brdf[k] * cosF[k];
    }
    length2[i + j] = eyePath.m_length2[0];
    probB[i + j] = eyePath.m_prob[0];
    probF[i + j] = probF[i + j - 1] * brdf[i + j - 1] * cosF[i + j - 1];
    for (int k = i - 1; k > -1; k--) probB[k] = probB[k + 1] * brdf[k + 1] * cosB[k + 1];
    //////////////////////////////////////////////////////////////////////////
    float formFactor = cosF_L * cosB_E / shadowLength2; // form factor
    Vector3 flux = lightPath.m_light->wattage()*brdf_L* brdf_E * formFactor;
    float prob = probB[i + 1] * (cosF[i] / length2[i + 1]);
    if (i > 0) prob *= probF[i - 1] * (cosB[i - 1] / length2[i - 1]);
    if (i > 1) flux *= lightPath.m_decay[i - 2];
    if (j > 1) flux *= eyePath.m_decay[j - 2];
    float probSum = 0;
    for (int k = 0; k < i + j; k++) {
        float p = probB[k + 1] * (cosF[k] / length2[k + 1]);
        if (k > 0) p *= probF[k - 1] * (cosB[k - 1] / length2[k - 1]);
        probSum += p;
    }
    flux *= prob / probSum;
    return flux;
}


PhotonMap* Scene::generatePhotonMapTest() {
    SequentialPhotonMap spm;
    Vector3 power(100, 100, 100);
    spm.addPhoton(PhotonDeposit(power, Vector3(-1, -1, 0), 1));
    spm.addPhoton(PhotonDeposit(power, Vector3(-1, -1, 2), 1));
    spm.addPhoton(PhotonDeposit(power, Vector3(-1, 1, 0), 1));
    spm.addPhoton(PhotonDeposit(power, Vector3(-1, 1, 2), 1));
    spm.addPhoton(PhotonDeposit(power, Vector3(1, -1, 0), 1));
    spm.addPhoton(PhotonDeposit(power, Vector3(1, -1, 2), 1));//*/
    spm.addPhoton(PhotonDeposit(power, Vector3(1, 1, 0), 1));
    spm.addPhoton(PhotonDeposit(power, Vector3(1, 1, 2), 1));
    PhotonMap* pm = new PhotonMap;
    pm->buildBalancedTree(spm);
    return pm;
}

PhotonMap* Scene::generatePhotonMap() {
    HitInfo hit;
    SequentialPhotonMap spm;
    for (int i = 0; i < m_lights.size(); i++) {
        Light* light = m_lights[i];
        int nPhotons = m_emittedPhotonsPerLight[i];
        Vector3 photonPower = light->wattage() / nPhotons;
        for (int j = 0; j < nPhotons; j++) {
            RayPDF rp = light->randRay();
            Ray ray = rp.m_ray;
            PhotonDeposit photon(photonPower, ray.o, 1);
            photon.m_prob = rp.m_oProb*rp.m_dProb;
            spm.addPhoton(photon);
            while (true) {
                if (!trace(hit, ray)) break;
                ray.o = hit.P;
                ray.d = hit.object->material()->randReflect(-ray.d, hit.N).v;
                photon.m_location = hit.P;
                //photon.m_prob // TODO: put in the probability for photon propagation
                Vector3 reflectance = hit.object->material()->reflectance();
                // Surface is same color as incident photon
                if (reflectance[0] / reflectance[1] == photon.m_power[0] / photon.m_power[1] &&
                    reflectance[1] / reflectance[2] == photon.m_power[1] / photon.m_power[2]) {
                    float rn = (float)rand() / RAND_MAX;
                    if (reflectance[0] == 1 || rn < reflectance[0]) spm.addPhoton(photon);
                    else break;
                }
                else {
                    float rn0 = (float)rand() / RAND_MAX;
                    float rn1 = (float)rand() / RAND_MAX;
                    float rn2 = (float)rand() / RAND_MAX;
                    if (!(reflectance[0] == 1 || rn0 < reflectance[0])) photon.m_power[0] = 0;
                    if (!(reflectance[1] == 1 || rn1 < reflectance[1])) photon.m_power[1] = 0;
                    if (!(reflectance[2] == 1 || rn2 < reflectance[2])) photon.m_power[2] = 0;
                    if (photon.m_power[0] > 0 || photon.m_power[1] > 0 || photon.m_power[2] > 0) spm.addPhoton(photon);
                    else break;
                }
                // TODO: handle cases where two of the channels match, so we don't get rainbow colors all over the place
            }
        }
    }
    PhotonMap* pm = new PhotonMap;
    pm->buildBalancedTree(spm);
    //std::vector<PhotonDeposit> asdf = pm->getPhotons();
    return pm;
}

Vector3 Scene::estimateFlux(int i, int j, LightPath lightPath, EyePath eyePath, PhotonMap* photonMap) {
    if (i == 0 && j == 0) return eyePath.m_hit[0].object->material()->radiance(eyePath.m_hit[0].N, -eyePath.m_ray[0].d);
    HitInfo hit;
    float intersectEpsilon = 0.00001;
    // The following is for the explicit connection
    HitInfo hit_E = eyePath.m_hit[j - 1];
    HitInfo hit_L;
    if (i == 0) hit_L = lightPath.m_lightHit;
    else hit_L = lightPath.m_hit[i - 1];
    Material* mat_E = hit_E.object->material();
    Material* mat_L = hit_L.object->material();
    Ray shadow_LtoE(hit_L.P, (hit_E.P - hit_L.P).normalize()); // shadow ray from Light path to Eye path
    float shadowLength2 = (hit_E.P - hit_L.P).length2();
    if (trace(hit, shadow_LtoE, intersectEpsilon, sqrt(shadowLength2) - intersectEpsilon)) return Vector3(0, 0, 0);
    if (i == 0) {
        if (!lightPath.m_light->intersect(hit, Ray(hit_E.P, (lightPath.m_lightHit.P - hit_E.P).normalize()))) return Vector3(0, 0, 0);
    }
    float cosF_L = std::max(0.0f, dot(shadow_LtoE.d, hit_L.N));
    if (cosF_L == 0) return Vector3(0, 0, 0);
    float cosB_E = std::max(0.0f, dot(-shadow_LtoE.d, hit_E.N));
    if (cosB_E == 0) return Vector3(0, 0, 0);
    float brdf_E = mat_E->BRDF(-shadow_LtoE.d, hit_E.N, -eyePath.m_ray[j - 1].d);
    if (brdf_E == 0) return Vector3(0, 0, 0);
    float brdf_L = 1;
    if (i > 0) brdf_L = mat_L->BRDF(shadow_LtoE.d, hit_L.N, -lightPath.m_ray[i - 1].d);
    if (brdf_L == 0) return Vector3(0, 0, 0);
    float dProb_EtoL = cosB_E*brdf_E;
    float dProb_LtoE = cosF_L*brdf_L;
    ///////////////////////////////////////////////////////////////////////////////////////////////
    vector<float> probF(i + j + 1); // forward from light to eye (excludes the const emission probability, since it's just a constant)
    vector<float> probB(i + j + 1); // backward from eye to light
    vector<float> length2(i + j + 1);
    vector<float> cosF(i + j);
    vector<float> cosB(i + j);
    vector<float> brdf(i + j);
    for (int k = 0; k < i; k++) {
        probF[k] = lightPath.m_prob[k];
        length2[k] = lightPath.m_length2[k];
        cosB[k] = lightPath.m_cosB[k];
        if (k == i - 1) break;
        cosF[k] = lightPath.m_cosF[k];
        brdf[k] = lightPath.m_brdf[k];
    }
    if (i > 0) {
        cosF[i - 1] = cosF_L;
        brdf[i - 1] = brdf_L;
    }
    length2[i] = shadowLength2;
    probB[i] = dProb_EtoL*eyePath.m_prob[j - 1];
    probF[i] = dProb_LtoE;
    if (i > 0) probF[i] *= probF[i - 1];
    cosF[i] = eyePath.m_cosB[j - 1];
    cosB[i] = cosB_E;
    brdf[i] = brdf_E;
    for (int k = i + 1; k < i + j; k++) {
        int u = i + j - k;
        cosF[k] = eyePath.m_cosB[u - 1]; // note the swap
        cosB[k] = eyePath.m_cosF[u - 1];
        brdf[k] = eyePath.m_brdf[u - 1];
        length2[k] = eyePath.m_length2[u];
        probB[k] = eyePath.m_prob[u];
        probF[k] = probF[k - 1] * brdf[k] * cosF[k];
    }
    length2[i + j] = eyePath.m_length2[0];
    probB[i + j] = eyePath.m_prob[0];
    probF[i + j] = probF[i + j - 1] * brdf[i + j - 1] * cosF[i + j - 1];
    for (int k = i - 1; k > -1; k--) probB[k] = probB[k + 1] * brdf[k + 1] * cosB[k + 1];
    //////////////////////////////////////////////////////////////////////////
    RadiusDensityPhotons rdp = photonMap->radiusDensityPhotons(hit_E.P, 1);
    Vector3 flux = brdf_L*cosF_L*brdf_E*rdp.m_density;
    if (j > 1) flux *= eyePath.m_decay[j - 2];
    if (i > 1) flux *= lightPath.m_decay[i - 2];
    float prob = probB[i + 1] * (cosF[i] / length2[i + 1])*probF[i] * (cosB[i] / length2[i]);
    float probSum = 0;
    for (int k = 0; k < i + j; k++) {
        float p = probB[k + 1] * (cosF[k] / length2[k + 1])*probF[k] * (cosB[k] / length2[k]);
        probSum += p;
    }
    flux *= prob / probSum;
    return flux;
}

void
Scene::unifiedpathtraceImage(Camera *cam, Image *img) {
    int w = img->width();
    int h = img->height();
    HitInfo hitInfo;

    std::ofstream plotfile;
    plotfile.open("unifiedpathtraceplot.txt");

    int integrationStart = glutGet(GLUT_ELAPSED_TIME);

    //PhotonMap* photonMap = generatePhotonMap();
    PhotonMap* photonMap = generatePhotonMapTest();

    for (int y = h - 1; y > -1; y--)
    //for (int y = 0; y < h; y++)
    {
        for (int x = 0; x < w; x++)
        {
            Ray ray00 = cam->eyeRay((float)x - 0.5, (float)y - 0.5, img->width(), img->height());
            Ray ray01 = cam->eyeRay((float)x - 0.5, (float)y + 0.5, img->width(), img->height());
            Ray ray10 = cam->eyeRay((float)x + 0.5, (float)y - 0.5, img->width(), img->height());
            Ray ray11 = cam->eyeRay((float)x + 0.5, (float)y + 0.5, img->width(), img->height());
            if (!trace(hitInfo, ray00) &&
                !trace(hitInfo, ray01) &&
                !trace(hitInfo, ray10) &&
                !trace(hitInfo, ray11)) continue;
            Vector3 fluxSumOverSamples(0, 0, 0);

            Concurrency::parallel_for(0, m_bidiSamplesPerPix, [&](int k) {
                EyePath eyePath = randEyePath(x, y, cam, img);
                if (eyePath.m_hit.size() == 0) return;
                LightPath lightPath = randLightPath();
                Vector3 fluxSum = estimateFlux(0, 0, lightPath, eyePath);
                for (int i = 0; i <= lightPath.m_hit.size(); i++) {
                    for (int j = 1; j <= eyePath.m_hit.size(); j++) {
                        fluxSum += estimateFlux(i, j, lightPath, eyePath, photonMap);
                    }
                }
                fluxSumOverSamples += fluxSum;
                if (y == h / 2 && x == w / 2){
                    plotfile << fluxSum[0] / (k + 1) / M_PI / 0.04 << std::endl;
                }
            });
            img->setPixel(x, y, fluxSumOverSamples / m_bidiSamplesPerPix);
        }
        if (preview())
        {
            img->drawScanline(y);
        }
        glFinish();
        printf("Rendering Progress: %.3f%%\r", y / float(img->height())*100.0f);
        fflush(stdout);
    }
    if (!preview())
    {
        img->draw();
    }
    glFinish();

    printf("Rendering Progress: 100.000%\n");

    int integrationEnd = glutGet(GLUT_ELAPSED_TIME);
    std::cout << "Rendering took " << ((integrationEnd - integrationStart) / 1000.0f) << "s" << std::endl;

    debug("Done Unified Pathtracing!\n");
    plotfile.close();
}
