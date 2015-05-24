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
    double rn = (double)(rand() / RAND_MAX);
    double em = hit.material->emittance();
    if (em == 1.0 || rn < em) {
        Vector3 rad = hit.material->radiance(hit.N, ray.o - hit.P);
        if (bounces == 0) return (1.0 / em)*rad;
        else return Vector3(0, 0, 0);
    }
    rn = (double)rand() / RAND_MAX;
    vec3pdf vp;
    if (m_samplingHeuristic == 1 || rn < m_samplingHeuristic) { // sample BRDF
        vp = hit.material->randReflect(-ray.d, hit.N); // pick BRDF weighted random direction
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
    Vector3 newDir = vp.v;
    Ray newRay;
    newRay.o = hit.P;
    newRay.d = newDir;
    float brdf = hit.material->BRDF(-ray.d, hit.N, newRay.d);
    float cos = std::max(0.0f, dot(hit.N, newDir));
    Vector3 gather = hit.material->shade(ray, hit, *this); // gathered direct lighting
    return
        gather + (1.0 / (1.0 - em)) / vp.p*brdf*cos*
        recursiveTrace_fromEye(newRay, bounces + 1, maxbounces);
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
    rayOut = rp.r;
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
        if (pix[2] > 0 && x > -1 && x<w && y>-1 && y < h) { // check that the pixel is within the viewing window
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
        double em = hit.material->emittance();
        rayToEye.o = hit.P;
        rayToEye.d = (cam->eye() - hit.P).normalize();
        float brdf = hit.material->BRDF(-rayIn.d, hit.N, rayToEye.d); // BRDF between eye-ray and incoming-ray
        if (!trace(tryHitEye, rayToEye)) { // check if anything is occluding the eye from current hitpoint
            pix = cam->imgProject(hit.P, w, h); // find the pixel the onto which the current hitpoint projects
            int x = round(pix[0]);
            int y = round(pix[1]);
            //img->setPixel(x, y, Vector3(1.0, 0.0, 0.0)); return;
            if (dot(Vector3(0, 0, -1), hit.N) != 0) {}
            if (pix[2]>0 && x >= 0 && x < w && y >= 0 && y < h) { // check that the pixel is within the viewing window
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
        vec3pdf vp = hit.material->randReflect(-rayIn.d, hit.N); // pick BRDF.cos weighted random direction
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
        if (((float)p / m_photonSamples) < 0.8 && p > 0 && p % (m_photonSamples / 3) == 0)
        {
            img->draw();
        }
        if (p % 1000000 == 0) {
            plotfile << img->getPixel(w / 2, h / 2)[0] * m_photonSamples / (p + 1) << std::endl;
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

    float W = 0.5;

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
            vector<Vector3> fluxSums;
            Vector3 fluxSum(0, 0, 0);

            //Concurrency::parallel_for(0, bidiSamplesPerPix(), [&](int k){
            for (int k = 0; k < m_bidiSamplesPerPix; k++) {
                RayPath eyePath = randEyePath(x, y, cam, img);
                if (eyePath.m_hits.size() == 0) {
                    return;
                }
                RayPath lightPath = randLightPath();

                int dj = eyePath.m_hits.size();
                int di = lightPath.m_hits.size();
                fluxSum += eyePath.m_hits[0].material->radiance(eyePath.m_hits[0].N, -eyePath.m_rays[0].d);
                bool lightFlag = false;
                if (fluxSum[0] > 0 || fluxSum[1] > 0 || fluxSum[2] > 0) lightFlag = true;
                for (int i = 0; i < di; i++){
                    for (int j = 1; j < dj; j++) {
                        fluxSum += estimateFlux(i, j, eyePath, lightPath);
                        //printf("%f %f %f\n", fluxSum[0], fluxSum[1], fluxSum[2]);
                        //system("PAUSE");
                    }
                }
                if (y == h / 2 && x == w / 2){
                    plotfile << fluxSum[0] / (k + 1) / M_PI / 0.04 << std::endl;
                }
            }
            //});
            img->setPixel(x, y, fluxSum / bidiSamplesPerPix());
            //img->setPixel(x, y, fluxSum);
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

RayPath Scene::randEyePath(float x, float y, Camera* cam, Image* img) {
    Ray rInit = cam->eyeRay(x, y, img->width(), img->height());
    RayPath raypath(rInit);
    raypath.m_probs.push_back(1.0f);
    raypath.m_fluxDecay.push_back(1.0f);

    RayPath asdf = generateRayPath(raypath);
    for (int i = 0; i < asdf.m_rays.size() - 1; i++) {
        if (asdf.m_rays[i].o == asdf.m_rays[i + 1].o) {
            printf("eye ray ORIGIN %i %i same out of %i\n", i, i + 1, asdf.m_rays.size());
            if (asdf.m_rays[i].d == asdf.m_rays[i + 1].d) {
                printf("eye ray %i %i same out of %i\n", i, i + 1, asdf.m_rays.size());
            }
        }
    }
    return asdf;
}

RayPath Scene::randLightPath() {
    HitInfo hit;
    LightPDF lp = randLightByWattage();
    Light* light = lp.l;

    RayPDF rp = light->randRay();
    RayPath raypath(rp.r);

    raypath.m_light = light;
    HitInfo first_hit(0.0f, rp.r.o, light->normal(rp.r.o));
    first_hit.material = light->material();
    raypath.m_hits.push_back(first_hit);
    raypath.m_fluxDecay.push_back(1.0f);
    if (trace(hit, rp.r)) {
        float cos = std::max(0.0f, dot(first_hit.N, rp.r.d));
        float cosPrime = std::max(0.0f, dot(hit.N, -rp.r.d));
        float distSqr = (hit.P - rp.r.o).length2();
        raypath.m_probs.push_back(rp.p*cosPrime / distSqr);
    }
    else raypath.m_probs.push_back(1.0f);
    return generateRayPath(raypath);
}

RayPath Scene::generateRayPath(RayPath & raypath) {
    HitInfo hit;
    if (!trace(hit, raypath.m_rayInit)) return raypath;
    raypath.m_hits.push_back(hit);

    int bounce = 0;
    while (bounce < m_maxPaths)
    {
        Ray lastRay = raypath.m_rays.back();
        HitInfo lastHit = raypath.m_hits.back();
        vec3pdf vp = lastHit.material->randReflect(-lastRay.d, lastHit.N);
        Ray newRay(lastHit.P, vp.v);
        if (!trace(hit, newRay)) return raypath;
        float brdf = lastHit.material->BRDF(-lastRay.d, lastHit.N, vp.v);
        float cos = std::max(0.0f, dot(lastHit.N, newRay.d));
        float cosPrime = std::max(0.0f, dot(hit.N, -newRay.d));
        float distSqr = (hit.P - lastHit.P).length2();
        raypath.m_rays.push_back(newRay);
        raypath.m_hits.push_back(hit);
        raypath.m_probs.push_back(brdf*cos*cosPrime * raypath.m_probs.back() / distSqr);
        raypath.m_fluxDecay.push_back(brdf*cos*raypath.m_fluxDecay.back());

        bounce++;
    }
    return raypath;
}

Vector3 Scene::estimateFlux(int i, int j, RayPath eyePath, RayPath lightPath) {
    HitInfo hit;
    float intersectEpsilon = 0.00001;
    // The following is for the explicit connection
    HitInfo hitE = eyePath.m_hits[j - 1];
    HitInfo hitL = lightPath.m_hits[i];
    const Material* matE = hitE.material;
    const Material* matL = hitL.material;
    Ray LshadowE(hitL.P, (hitE.P - hitL.P).normalize()); // shadow ray from Light path to Eye path
    float LshadowE_length2 = (hitE.P - hitL.P).length2();
    if (trace(hit, LshadowE, intersectEpsilon, sqrt(LshadowE_length2) - intersectEpsilon)) return Vector3(0, 0, 0);
    if (i == 0) {
        Vector3 lightPoint = lightPath.m_hits[0].P;
        if (!lightPath.m_light->intersect(hit, Ray(hitE.P, (lightPoint - hitE.P).normalize()))) return Vector3(0, 0, 0);
    }
    float LcosE = std::max(0.0f, dot(LshadowE.d, hitL.N));
    float EcosL = std::max(0.0f, dot(-LshadowE.d, hitE.N));
    float brdfE = matE->BRDF(-LshadowE.d, hitE.N, -eyePath.m_rays[j - 1].d);
    float brdfL = 1;
    if (i>0) brdfL = matL->BRDF(LshadowE.d, hitL.N, -lightPath.m_rays[i - 1].d);
    float EdecayL = EcosL*brdfE; // NOT cumulative
    float LdecayE = LcosE*brdfL;
    float EprobL = EdecayL*LcosE / LshadowE_length2;
    float LprobE = LdecayE*EcosL / LshadowE_length2;
    ///////////////////////////////////////////////////////////////////////////////////////////////
    vector<float> probF(i + j + 1); // forward from light to eye
    vector<float> probB(i + j + 1); // backward from eye to light
    vector<float> decayF(i + j + 1);
    vector<float> decayB(i + j + 1);
    vector<float> brdf(i + j + 1);
    vector<float> cosF(i + j + 1);
    vector<float> cosB(i + j + 1);
    vector<float> length2(i + j + 1);
    for (int k = 0; k < i; k++) {
        probF[k] = lightPath.m_probs[k];
        decayF[k] = lightPath.m_fluxDecay[k];
        const Material* m = lightPath.m_hits[k].material;
        HitInfo h = lightPath.m_hits[k];
        HitInfo hNext = lightPath.m_hits[k + 1];
        Vector3 n = h.N;
        Vector3 nNext = hNext.N;
        Vector3 out = lightPath.m_rays[k].d;
        if (k == 0) brdf[k] = 1;
        else {
            Vector3 in = -lightPath.m_rays[k - 1].d;
            brdf[k] = m->BRDF(in, n, out);
        }
        cosF[k] = std::max(0.0f, dot(n, out));
        cosB[k] = std::max(0.0f, dot(nNext, -out));
        length2[k] = (hNext.P - h.P).length2();
        if (length2[k] == 0) {
            printf("location 1\n");
            printf("%i\n", k);
            printf("%f %f %f\n", h.P[0], h.P[1], h.P[2]);
            printf("%f %f %f\n", hNext.P[0], hNext.P[1], hNext.P[2]);
            system("PAUSE");
        }
    }
    probF[i] = LprobE;
    decayF[i] = LdecayE;
    if (i > 0){
        decayF[i] *= decayF[i - 1];
        probF[i] *= probF[i - 1];
    }
    cosF[i] = LcosE;
    cosB[i] = EcosL;
    brdf[i] = brdfL;
    brdf[i + 1] = brdfE;
    length2[i] = LshadowE_length2;
    if (length2[i] == 0) {
        cout << "location 2" << endl;
        system("PAUSE");
    }
    for (int k = i + 1; k < i + j + 1; k++) {
        int u = i + j - k;
        probB[k] = eyePath.m_probs[u];
        decayB[k] = eyePath.m_fluxDecay[u];
        const Material* m = eyePath.m_hits[u].material;
        Vector3 n = eyePath.m_hits[u].N;
        Vector3 out = -eyePath.m_rays[u].d;
        cosF[k] = std::max(0.0f, dot(n, out));
        if (u < j - 1) {
            Vector3 in = eyePath.m_rays[u + 1].d;
            brdf[k] = m->BRDF(in, n, out);
        }
        if (u == 0) cosB[k] = 1.0f;
        else {
            Vector3 nNext = eyePath.m_hits[u - 1].N;
            cosB[k] = std::max(0.0f, dot(nNext, -out));
        }
        if (u == j - 1) length2[k] = (eyePath.m_hits[u].P - eyePath.m_rays[u].o).length2();
        else length2[k] = (eyePath.m_rays[u + 1].o - eyePath.m_rays[u].o).length2();
        /*if (length2[k] == 0) {
        printf("location 3\n");
        printf("%i %i\n", u, j-1);
        printf("%f %f %f\n", eyePath.m_hits[u].P[0], eyePath.m_hits[u].P[1], eyePath.m_hits[u].P[2]);
        printf("%f %f %f\n", eyePath.m_rays[u].o[0], eyePath.m_rays[u].o[1], eyePath.m_rays[u].o[2]);
        system("PAUSE");
        }*/
        decayF[k] = decayF[k - 1] * brdf[k] * cosF[k];
        probF[k] = probF[k - 1] * brdf[k] * cosF[k] * cosB[k] / length2[k];
    }
    for (int k = i; k > -1; k--) {
        probB[k] = probB[k + 1] * brdf[k + 1] * cosB[k];
        decayB[k] = decayB[k + 1] * brdf[k + 1] * cosB[k] * cosF[k] / length2[k];
    }
    //////////////////////////////////////////////////////////////////////////
    /*printf("%i %i\n", i, j);
    for (int k = 0; k < i + j + 1; k++) {
    printf("%f %f %f %f %f %f %f %f\n", cosF[k], cosB[k], brdf[k], length2[k], decayF[k], decayB[k], probF[k], probB[k]);
    }//*/
    Vector3 fluxSum(0, 0, 0);
    float probSum = 0;
    for (int k = 0; k < i + j; k++) {
        if (fabs(length2[k]) < 0.00001) continue;
        float form = cosF[k] * cosB[k] / length2[k]; // form factor
        if (form > 100000) {
            printf("\n%f\n", form);
            printf("%f\n", cosF[k]);
            printf("%f\n", cosB[k]);
            printf("%f\n", length2[k]);
            system("PAUSE");
        }//*/
        Vector3 flux = lightPath.m_light->wattage()*brdf[k] * brdf[k + 1] * decayB[k + 1] * form;
        float prob = probF[k + 1];
        if (k > 0) {
            prob *= probB[k - 1];
            flux *= decayF[k - 1];
        }
        probSum += prob;
        fluxSum += flux*prob;
    }
    fluxSum /= probSum;
    //printf("%f %f %f\n",fluxSum[0],fluxSum[1] < 0,fluxSum[2]);
    if (fluxSum[0] < 0 || fluxSum[1] < 0 || fluxSum[2] < 0) printf("negative flux\n");
    if (fluxSum[0] == -INFINITY || fluxSum[1] == -INFINITY || fluxSum[2] == -INFINITY) printf("negative infinite flux\n");
    if (fluxSum[0] == INFINITY || fluxSum[1] == INFINITY || fluxSum[2] == INFINITY) printf("positive infinite flux\n");
    return fluxSum;
}

PhotonMap Scene::generatePhotonMap() {
    HitInfo hit;
    PhotonMap map;
    for (int i = 0; i < m_lights.size(); i++) {
        Light* light = m_lights[i];
        int nPhotons = m_emittedPhotonsPerLight[i];
        Vector3 photonPower = light->wattage() / nPhotons;
        for (int j = 0; j < nPhotons; j++) {
            Ray ray = light->randRay().r;
            PhotonDeposit photon(photonPower, ray.o);
            map.push_back(photon);
            while (true) {
                if (!trace(hit, ray)) break;
                ray.o = hit.P;
                ray.d = hit.material->randReflect(-ray.d, hit.N).v;
                photon.m_location = hit.P;
                Vector3 reflectance = hit.material->reflectance();
                // Surface is same color as incident photon
                if (reflectance[0] / reflectance[1] == photon.m_power[0] / photon.m_power[1] &&
                    reflectance[1] / reflectance[2] == photon.m_power[1] / photon.m_power[2]) {
                    float rn = (float)rand() / RAND_MAX;
                    if (reflectance[0] == 1 || rn < reflectance[0]) map.push_back(photon);
                    else break;
                }
                else {
                    float rn0 = (float)rand() / RAND_MAX;
                    float rn1 = (float)rand() / RAND_MAX;
                    float rn2 = (float)rand() / RAND_MAX;
                    if (!(reflectance[0] == 1 || rn0 < reflectance[0])) photon.m_power[0] = 0;
                    if (!(reflectance[1] == 1 || rn1 < reflectance[1])) photon.m_power[1] = 0;
                    if (!(reflectance[2] == 1 || rn2 < reflectance[2])) photon.m_power[2] = 0;
                    if (photon.m_power[0] > 0 || photon.m_power[1] > 0 || photon.m_power[2] > 0) map.push_back(photon);
                    else break;
                }
                // TODO: handle cases where two of the channels match, so we don't get rainbow colors all over the place
            }
        }
    }
}

Vector3 Scene::estimateFlux(int i, int j, RayPath eyePath, RayPath lightPath, PhotonMap photonMap) {
    Vector3 flux = Vector3(0, 0, 0);
    HitInfo h;
    float intersectEpsilon = 0.00001;
    Ray rayE = eyePath.m_rays[j - 1];
    HitInfo hitE = eyePath.m_hits[j - 1];
    HitInfo hitL = lightPath.m_hits[i];
    const Material* matE = hitE.material;
    const Material* matL = hitL.material;
    float xDot = dot(Vector3(1, 0, 0), hitE.N);
    float yDot = dot(Vector3(0, 0, 1), hitE.N);
    Vector3 yAxis;
    if (fabs(yDot) > fabs(xDot)) yAxis = Vector3(0, 1, 0).orthogonal(hitE.N).normalize();
    else yAxis = Vector3(1, 0, 0).orthogonal(hitE.N).normalize();
    Vector3 xAxis = cross(yAxis, hitE.N);
    float r = sqrt(m_photonMapRadius*rand() / RAND_MAX);
    float phi = 2.0f * PI*rand() / RAND_MAX;
    Vector3 perturbation = r*cos(phi)*xAxis + r*sin(phi)*yAxis;
    Vector3 hitPointE_perturbed = hitE.P + perturbation;
    Ray LshadowE(hitL.P, (hitPointE_perturbed - hitL.P).normalize());
    float LshadowE_length2 = (hitPointE_perturbed - hitL.P).length2();

    if (trace(h, LshadowE, intersectEpsilon, sqrt(LshadowE_length2) - intersectEpsilon)) return flux;
    if (i == 0) {
        if (!lightPath.m_light->intersect(h, Ray(hitE.P + perturbation,-LshadowE.d))) return flux;
        float f = matE->BRDF(-LshadowE.d, hitE.N, -perturbation.normalized()) * matE->BRDF(perturbation.normalized(), hitE.N, -rayE.d);
        flux = lightPath.m_light->wattage()*lightPath.m_light->rayPDF(LshadowE)*f*eyePath.m_fluxDecay[j - 1];
    }
    else {
        float f = matE->BRDF(-LshadowE.d, hitE.N, -perturbation.normalized()) * matE->BRDF(perturbation.normalized(), hitE.N, -rayE.d);
        float brdfL = matL->BRDF(LshadowE.d, hitL.N, perturbation.normalized());
        flux = lightPath.m_light->wattage()*lightPath.m_fluxDecay[i] * brdfL*eyePath.m_fluxDecay[j - 1];
    }
    return flux;
}