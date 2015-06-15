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
    //for (int j = img->height() - 1; j > -1; j--)
    {
        for (int i = 0; i < img->width(); ++i)
        {
            ray = cam->eyeRay(i, j, img->width(), img->height());
            if (!trace(hitInfo, ray)) {
                continue;
            }
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
    Vector3 brdf = hit.object->BRDF(-ray.d, hit.N, newRay.d, hit.P);
    float cos = fabs(dot(hit.N, newDir)); // changed this to fabs for transmissible materials such as RefractiveInterface
    Vector3 gather = hit.object->shade(ray, hit, *this, hit.P); // gathered direct lighting
    return gather + (1.0 / (1.0 - em)) / vp.p*brdf*cos*recursiveTrace_fromEye(newRay, bounces + 1, maxbounces);
}

void
Scene::pathtraceImage(Camera *cam, Image *img)
{
    HitInfo hitInfo;
    int integrationStart = glutGet(GLUT_ELAPSED_TIME);
    
       // loop over all pixels in the image
    //for (int j = 0; j < img->height(); ++j)
    for (int j = img->height() - 1; j > -1; j--)
    {
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
                pixSum += recursiveTrace_fromEye(ray, 0, m_maxBounces) / (double)m_samplesPerPix;
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
        double reflectance = hit.object->material()->reflectance()[0];
        rayToEye.o = hit.P;
        rayToEye.d = (cam->eye() - hit.P).normalize();
        Vector3 brdf = hit.object->material()->BRDF(-rayIn.d, hit.N, rayToEye.d); // BRDF between eye-ray and incoming-ray
        if (!trace(tryHitEye, rayToEye)) { // check if anything is occluding the eye from current hitpoint
            pix = cam->imgProject(hit.P, w, h); // find the pixel the onto which the current hitpoint projects
            int x = round(pix[0]);
            int y = round(pix[1]);
            if (dot(Vector3(0, 0, -1), hit.N) != 0) {}
            if (pix[2]>0 && x >= 0 && x<w && y >= 0 && y < h) { // check that the pixel is within the viewing window
                float cos0 = dot(hit.N, rayToEye.d);
                float lengthSqr = (cam->eye() - hit.P).length2();
                float cosAlpha = cam->pixelCosine(pix[0], pix[1], w, h);
                if (reflectance == 1 || rn < reflectance) {
                    img->setPixel(x, y, img->getPixel(x, y) + (1.0f / reflectance) * power*brdf*cos0*cosAlpha / lengthSqr);
                }
                else {
                    img->setPixel(x, y, img->getPixel(x, y) + (1.0f / (1.0f - reflectance)) * power*brdf*cos0*cosAlpha / lengthSqr);
                    return; // photon was absorbed
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

    Vector3 floorPoint = cam->imgProject(Vector3(0, 0, 0), w, h);
    int floorPointPixel[] = { round(floorPoint[0]), round(floorPoint[1]) };
    float pixelFactor = cam->pixelCosine(floorPoint[0], floorPoint[1], w, h) * cam->pixelSolidAngle(floorPoint[0], floorPoint[1], w, h);
    printf("center-point of floor is at pixel ( %i , %i )\n", floorPointPixel[0], floorPointPixel[1]);

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
        if (((float)p / m_photonSamples)<0.8 && p>0 && p % (m_photonSamples / 3) == 0) img->draw();
        if (p % 1000000 == 0) {
            Vector3 dataPoint = img->getPixel(floorPointPixel[0], floorPointPixel[1]) / (p*pixelFactor);
            plotfile << dataPoint[0] << std::endl;
        }
    }
    for (int i = 0; i < w; i++) {
        for (int j = 0; j < h; j++) {
            Vector3 pix = img->getPixel(i, j);
            img->setPixel(i, j, pix / (cam->pixelCosine(i, j, w, h) * cam->pixelSolidAngle(i, j, w, h) * m_photonSamples));
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
    float random = (float)rand() / RAND_MAX;
    float r = wattageConcattage[n] * random;
    for (int i = 0; i < n; i++){ // find the interval in which r lies and return the light along with PDF;
        if (r < wattageConcattage[i] || r > wattageConcattage[i + 1]) continue;
        rlbw.l = m_lights[i];
        rlbw.p = m_lights[i]->wattage() / wattageConcattage[n];
        return rlbw;
    }
}

pair<Vector3, Vector3> Scene::axisAlignedBounds() {
    if (m_objects.size() == 0) return pair<Vector3, Vector3>(Vector3(0, 0, 0), Vector3(0, 0, 0));
    pair<Vector3, Vector3> objBounds = m_objects[0]->axisAlignedBounds();
    Vector3 minBounds = objBounds.first;
    Vector3 maxBounds = objBounds.second;
    for (int i = 1; i < m_objects.size(); i++) {
        objBounds = m_objects[i]->axisAlignedBounds();
        Vector3 xyz = objBounds.first;
        Vector3 XYZ = objBounds.second;
        if (xyz.x < minBounds.x) minBounds.x = xyz.x;
        if (xyz.y < minBounds.y) minBounds.y = xyz.y;
        if (xyz.z < minBounds.z) minBounds.z = xyz.z;
        if (XYZ.x > maxBounds.x) maxBounds.x = XYZ.x;
        if (XYZ.y > maxBounds.y) maxBounds.y = XYZ.y;
        if (XYZ.z > maxBounds.z) maxBounds.z = XYZ.z;
    }
    return pair<Vector3, Vector3>(minBounds, maxBounds);
}

Vector3 Scene::bidiRadiance(int i, int j, LightPath lightPath, EyePath eyePath) {
    if (i == 0 && j == 0) return eyePath.m_hit[0].object->material()->radiance(eyePath.m_hit[0].N, -eyePath.m_ray[0].d);
    if (i < 0 || j < 1) return Vector3(0, 0, 0);

    vector<float> probF;
    vector<float> probB;
    Vector3 estimatorLink;
    if (!forwardBackwardProbs(i, j, lightPath, eyePath, true, probF, probB, estimatorLink)) return Vector3(0, 0, 0);

    Vector3 flux = lightPath.m_light->wattage()*estimatorLink;
    if (i > 1) flux *= lightPath.m_estimator[i - 2];
    if (j > 1) flux *= eyePath.m_estimator[j - 2];

    float probSum = 0;
    for (int k = 0; k < i + j; k++) {
        float probPI = 1;
        probPI *= probB[k + 1];
        if (k > 0) probPI *= probF[k - 1];

        if (k == i) flux *= probPI;

        probSum += probPI;
    }
    if (isfinite(1 / probSum)) {
        return flux / probSum;
    }
    else return Vector3(0, 0, 0);
}

void Scene::biditraceImage(Camera *cam, Image *img) {
    int w = img->width();
    int h = img->height();
    HitInfo hitInfo;
    int integrationStart = glutGet(GLUT_ELAPSED_TIME);

    for (int y = h-1; y >-1; y--)
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
                Vector3 fluxSum = bidiRadiance(0, 0, lightPath, eyePath);
                for (int i = 0; i <= lightPath.m_hit.size(); i++) {
                    for (int j = 1; j <= eyePath.m_hit.size(); j++) {
                        fluxSum += bidiRadiance(i, j, lightPath, eyePath);
                    }
                }
                fluxSumOverSamples += fluxSum;
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

}



EyePath Scene::randEyePath(float x, float y, Camera* cam, Image* img, const int& bounces) {
    HitInfo hit;
    Ray rayInit = cam->eyeRay(x, y, img->width(), img->height());
    EyePath eyePath(rayInit);
    if (!trace(hit, rayInit)) return eyePath;
    else {
        eyePath.m_prob.push_back(1);
        eyePath.m_hit.push_back(hit);
        eyePath.m_cosB.push_back(std::max(0.0f, dot(hit.N, -rayInit.d)));
        eyePath.m_length2.push_back((hit.P - rayInit.o).length2());
        if (hit.object->material()->reflectance()[0] == 0) return eyePath;
    }

    if (bounces < 0) bounceRayPath(eyePath, m_maxEyePaths);
    else bounceRayPath(eyePath, bounces);
    return eyePath;
}

LightPath Scene::randLightPath(Light* lightInput, const int& bounces) {
    HitInfo hit;
    Light* light;
    if (lightInput == NULL) { light = randLightByWattage().l; }
    else light = lightInput;

    RayPDF rp = light->randRay();
    Ray rayInit = rp.m_ray;
    HitInfo lightHit(0.0f, rayInit.o, light->normal(rayInit.o), light);
    LightPath lightPath(rayInit);

    // lightPath specific values
    lightPath.m_light = light;
    lightPath.m_lightHit = lightHit;
    lightPath.m_originProb = rp.m_oProb;

    if (!trace(hit, rayInit)) return lightPath;
    else {
        float cosPrime = std::max(0.0f, dot(hit.N, -rayInit.d));
        float distSqr = (hit.P - lightHit.P).length2();
        lightPath.m_hit.push_back(hit);
        lightPath.m_cosB.push_back(cosPrime);
        lightPath.m_length2.push_back(distSqr);
        lightPath.m_prob.push_back(rp.m_dProb*cosPrime / distSqr);
        if (hit.object->material()->reflectance()[0] == 0) return lightPath;
    }

    if (bounces < 0) bounceRayPath(lightPath, m_maxLightPaths);
    else bounceRayPath(lightPath, bounces);
    return lightPath;
}

void Scene::bounceRayPath(RayPath & rayPath, const int& maxBounces) {
    if (maxBounces < 1) return;
    
    HitInfo newHit;
    int bounce = 1;
    bool terminate = false;
    while (bounce < maxBounces)
    {
        Ray lastRay = rayPath.m_ray.back();
        HitInfo lastHit = rayPath.m_hit.back();

        // Russian Roulette
        Vector3 lastReflectanceRGB = lastHit.object->material()->reflectance();
        float lastReflectance = (lastReflectanceRGB[0] + lastReflectanceRGB[1] + lastReflectanceRGB[2]) / 3;
        float rn = (float)rand() / RAND_MAX;
        if (rn > lastReflectance) terminate = true; // terminate upon this next hit

        vec3pdf vp = lastHit.object->material()->randReflect(-lastRay.d, lastHit.N);
        Ray newRay(lastHit.P, vp.v);
        if (!trace(newHit, newRay)) return;
        Vector3 newReflectanceRGB = newHit.object->material()->reflectance();
        float newReflectance = (newReflectanceRGB[0] + newReflectanceRGB[1] + newReflectanceRGB[2]) / 3;
        if (newReflectance == 0) return;

        Vector3 brdf = lastHit.object->material()->BRDF(-lastRay.d, lastHit.N, newRay.d);
        float cos = std::max(0.0f, dot(lastHit.N, newRay.d));
        float cosPrime = std::max(0.0f, dot(newHit.N, -newRay.d));
        float distSqr = (newHit.P - lastHit.P).length2();

        Vector3 estimator;
        if (!terminate)
            estimator = brdf*cos / vp.p / lastReflectance;
        else
            estimator = brdf*cos / vp.p / (1.0f - lastReflectance);
        if (rayPath.m_estimator.size() != 0) estimator *= rayPath.m_estimator.back();

        rayPath.m_estimator.push_back(estimator);
        rayPath.m_hit.push_back(newHit);
        rayPath.m_cosF.push_back(cos);
        rayPath.m_cosB.push_back(cosPrime);

        rayPath.m_ray.push_back(newRay);
        rayPath.m_length2.push_back(distSqr);

        float prob;
        if (!terminate) {
            prob = rayPath.m_prob.back()*vp.p*(cosPrime / distSqr)*lastReflectance;
        }
        else {
            prob = rayPath.m_prob.back()*vp.p*(cosPrime / distSqr)*(1.0f - lastReflectance);
        }
        rayPath.m_prob.push_back(prob);

        if (terminate) return;

        bounce++;
    }
}

pair<PhotonMap*, vector<LightPath*>> Scene::generatePhotonMap() {
    int nPaths = 0;
    for (int i = 0; i < m_lights.size(); i++) nPaths += m_emittedPhotonsPerLight[i];
    vector<LightPath*> paths(nPaths);

    SequentialPhotonMap spm;
    int pathCount = 0;
    for (int i = 0; i < m_lights.size(); i++) {
        Light* light = m_lights[i];
        int nPhotons = m_emittedPhotonsPerLight[i];
        Vector3 photonPower = light->wattage() / nPhotons;
        int printStep = fmax(1, nPhotons / 100);
        for (int j = 0; j < nPhotons; j++) {
            if (j % printStep == 0 || j == nPhotons - 1) printf("Bouncing photon %i/%i light %i __________\r", j, nPhotons, i);
            LightPath* path = new LightPath;
            *path = randLightPath(light, m_maxLightPaths);
            //if (path->m_hit.size() > 1) {
            //    cout << "oh nos\n" << endl;
            //}
            paths[pathCount] = path;
            pathCount++;
            //spm.addPhoton(PhotonDeposit(photonPower, path, -1));
            for (int k = 0; k < path->m_hit.size(); k++) spm.addPhoton(PhotonDeposit(photonPower, path, k));
        }
    }
    printf("\nPopulating photon Map with %i photons...\n", spm.nPhotons());
    return pair<PhotonMap*, vector<LightPath*>>(spm.buildBalancedTree(0, true), paths);
}

Vector3 Scene::uniRadiance(const int& i, const int& j, const LightPath& lightPath, const EyePath& eyePath, PhotonMap* photonMap,
    const bool& explicitConnection, const int& nLightPaths, const Vector3& density, const float& radiusInput)
{
    if (i == 0 && j == 0) return eyePath.m_hit[0].object->material()->radiance(eyePath.m_hit[0].N, -eyePath.m_ray[0].d);
    if (i < 0 || j < 1) return Vector3(0, 0, 0);

    vector<float> probF;
    vector<float> probB;
    Vector3 estimatorLink;
    if (!forwardBackwardProbs(i, j, lightPath, eyePath, explicitConnection, probF, probB, estimatorLink)) return Vector3(0, 0, 0);

    Vector3 flux;
    float diskArea = M_PI*m_photonGatheringRadius*m_photonGatheringRadius;
    if (explicitConnection == true) {
        flux = lightPath.m_light->wattage()*estimatorLink;
        if (i > 1) flux *= lightPath.m_estimator[i - 2];
        if (j > 1) flux *= eyePath.m_estimator[j - 2];
    }
    else {
        flux = density*estimatorLink;
        if (i > 0) flux *= lightPath.m_estimator[i - 1]; // note the extra segment i-1 as opposed to i-2
        if (j > 1) flux *= eyePath.m_estimator[j - 2];
    }

    //return flux;

    float probSum = 0;
    for (int k = 0; k < i + j; k++) {
        float probPI = 1;
        probPI *= probB[k + 1];
        if (k > 0) probPI *= probF[k - 1];
        float probDE = probF[k] * probB[k + 1] * diskArea;

        probPI = 0; // just for testing photonmapping only

        if (k == i) {
            if (explicitConnection == true) flux *= probPI;
            else flux *= probDE;
        }

        probSum += probPI + probDE*nLightPaths;
    }
    if (isfinite(1 / probSum)) {
        return flux / probSum;
    }
    else {
        return Vector3(0, 0, 0);
    }
}


Vector3 Scene::uniRadianceDE(const int& j, const EyePath& eyePath, PhotonMap* photonMap, const int& nLightPaths) {
    HitInfo hit_E = eyePath.m_hit[j - 1];
    vector<PhotonDeposit> photons = photonMap->getPhotons(hit_E.P, m_photonGatheringRadius);
    Vector3 density = 0;
    float diskArea = M_PI*m_photonGatheringRadius*m_photonGatheringRadius;
    for (int k = 0; k < photons.size(); k++) density += photons[k].m_power;
    density /= diskArea;

    //cout << density[0] << endl;

    Vector3 flux(0,0,0);

    int smallPhotons = 0;

    for (auto & photon : photons) {
        int i = photon.m_hitIndex;
        // usually with explicit connection we would make i = photon.m_hitIndex +1
        // we do one less than usual because we make the explicit connection (for weighting) with the second to last hit
        // the extra vertex (which is treated as not actually existent and merged with it's neighbor) is handled in uniRadiance(...)
        LightPath lightPath = *photon.m_lightPath;
        Vector3 fluxAccum = uniRadiance(i, j, lightPath, eyePath, photonMap, false, nLightPaths, density);
        flux += fluxAccum;
        if (fluxAccum[0] < 0.000001) {
            smallPhotons++;
        }
        //return flux;
    }
    //cout << smallPhotons << "/" << photons.size() << endl;
    return flux;
}

void
Scene::unifiedpathtraceImage(Camera *cam, Image *img) {
    int w = img->width();
    int h = img->height();
    HitInfo hitInfo;
    
    int integrationStart = glutGet(GLUT_ELAPSED_TIME);

    pair<PhotonMap*,vector<LightPath*>> mapAndPaths = generatePhotonMap();
    int nLightPaths = mapAndPaths.second.size();

    for (int y = h - 1; y > -1; y--)
    //for (int y = 0; y < h; y++)
    {
        //for (int x = w / 2 - 1; x < w / 2 + 1; x++)
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

            /*Ray ray = cam->eyeRay(x, y, img->width(), img->height());
            trace(hitInfo, ray00);
            vector<PhotonDeposit> photons = mapAndPaths.first->getPhotons(hitInfo.P, m_photonGatheringRadius);
            Vector3 density(0, 0, 0);
            for (int k = 0; k<photons.size(); k++) {
                density += photons[k].m_power;
            }
            density /= M_PI*m_photonGatheringRadius*m_photonGatheringRadius;
            img->setPixel(x, y, density / (2 * M_PI));//*/

            Concurrency::parallel_for(0, m_bidiSamplesPerPix, [&](int k) {
                EyePath eyePath = randEyePath(x, y, cam, img);
                if (eyePath.m_hit.size() == 0) return;
                int randIndex = (int)nLightPaths*((float)rand() / RAND_MAX);
                int lightPathIndex = min(nLightPaths - 1, randIndex );
                LightPath lightPath = *mapAndPaths.second[lightPathIndex];
                Vector3 fluxSum = uniRadiance(0, 0, lightPath, eyePath, mapAndPaths.first, true, nLightPaths);
                for (int j = 1; j <= eyePath.m_hit.size(); j++) {
                    fluxSum += uniRadianceDE(j, eyePath, mapAndPaths.first, nLightPaths);
                    for (int i = 0; i <= lightPath.m_hit.size(); i++) {
                        //fluxSum += uniRadiance(i, j, lightPath, eyePath, mapAndPaths.first, true, nLightPaths);
                    }
                }
                fluxSumOverSamples += fluxSum;
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
}




bool Scene::forwardBackwardProbs(const int& i, const int& j, const LightPath& lightPath, const EyePath& eyePath, const bool& explicitConnection,
    vector<float>& probF, vector<float>& probB, Vector3& estimatorLink)
{
    HitInfo hit;
    float intersectEpsilon = 0.00001;

    // First, compute the shadow connection FROM lightPath TO eyePath
    HitInfo hit_E = eyePath.m_hit[j - 1];
    HitInfo hit_L;
    if (i == 0) hit_L = lightPath.m_lightHit;
    else hit_L = lightPath.m_hit[i - 1];
    Material* mat_E = hit_E.object->material();
    Material* mat_L = hit_L.object->material();
    Ray shadow_LtoE(hit_L.P, (hit_E.P - hit_L.P).normalize());
    float shadowLength2 = (hit_E.P - hit_L.P).length2();

    // The following are checks for whether or not the shadow connection is valid
    if (trace(hit, shadow_LtoE, intersectEpsilon, sqrt(shadowLength2) - intersectEpsilon)) return false;
    if (i == 0 && !lightPath.m_light->intersect(hit, Ray(hit_E.P, (lightPath.m_lightHit.P - hit_E.P).normalize()))) return false;
    float cosF_L = std::max(0.0f, dot(shadow_LtoE.d, hit_L.N));
    if (cosF_L == 0) return false;
    float cosB_E = std::max(0.0f, dot(-shadow_LtoE.d, hit_E.N));
    if (cosB_E == 0) return false;
    Vector3 brdf_E = mat_E->BRDF(-shadow_LtoE.d, hit_E.N, -eyePath.m_ray[j - 1].d);
    if (brdf_E == 0) return false;
    Vector3 brdf_L = 1;
    if (i > 0) brdf_L = mat_L->BRDF(shadow_LtoE.d, hit_L.N, -lightPath.m_ray[i - 1].d);
    if (brdf_L == 0) return false;

    // The following are the probabilities (not cumulative) of making the transition from eyepath to lightpath, or visa versa
    float dProb_EtoL = mat_E->reflectPDF(-eyePath.m_ray[j - 1].d, hit_E.N, -shadow_LtoE.d)*(cosF_L / shadowLength2);
    float dProb_LtoE = 0;
    if (i == 0) dProb_LtoE = mat_L->emitPDF(hit_L.N, lightPath.m_ray[0].d)*(cosB_E / shadowLength2);
    else dProb_LtoE = mat_L->reflectPDF(-lightPath.m_ray[i - 1].d, hit_L.N, shadow_LtoE.d)*(cosB_E / shadowLength2);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // At this point, we have computed all the relevant quantities for the shadow connection. Proceed to build the probability sequences //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    probF.resize(i + j + 1);        // forward from light to eye (excludes the emission origin probability, since it's just a constant 1/lightArea)
    probB.resize(i + j + 1);        // backward from eye to light
    vector<float> length2(i + j);   // The following three are for converting to AREA space PDFs
    vector<float> cosF(i + j);      // cosF[k] and cosB[k] are the forward/backward angles at hit[k]
    vector<float> cosB(i + j);

    // For the first half most quantities can just be copied from lightPath
    for (int k = 0; k < i; k++) {
        probF[k] = lightPath.m_prob[k];
        length2[k] = lightPath.m_length2[k];
        cosB[k] = lightPath.m_cosB[k];
        if (k == i - 1) break;
        cosF[k] = lightPath.m_cosF[k];  // The last cosF on the lightPath is undefined
    }

    // Next, fill in the quantities at the shadow connection
    if (i > 0) { 
        cosF[i - 1] = cosF_L;
    }
    length2[i] = shadowLength2;
    probB[i] = dProb_EtoL*eyePath.m_prob[j - 1];
    probF[i] = dProb_LtoE;
    if (i > 0) probF[i] *= probF[i - 1];
    cosF[i] = eyePath.m_cosB[j - 1];
    cosB[i] = cosB_E;

    // For the latter half, most quantities can be copied from eyePath
    for (int k = i + 1; k < i + j; k++) {
        int u = i + j - k;                  // hit i+1 on the cocatenated path corresponds to the u=j-1 hit on the eyePath
        cosF[k] = eyePath.m_cosB[u - 1];    // note the swap (B on the lightpath becomes F on the concatenated path)
        cosB[k] = eyePath.m_cosF[u - 1];
        length2[k] = eyePath.m_length2[u];
        probB[k] = eyePath.m_prob[u];
        float reverseReflectProb = 0;
        if (u == j - 1) { // This needs to handled specially since the incoming direction is the shadow ray
            reverseReflectProb = eyePath.m_hit[u].object->material()->reflectPDF(-shadow_LtoE.d, eyePath.m_hit[u].N, -eyePath.m_ray[u].d);
        } else {
            reverseReflectProb = eyePath.m_hit[u].object->material()->reflectPDF(eyePath.m_ray[u+1].d, eyePath.m_hit[u].N, -eyePath.m_ray[u].d);
        }
        probF[k] = probF[k - 1] * reverseReflectProb * (cosB[k] / length2[k]);
    }
    probF[i + j] = eyePath.m_prob[0];   // For now, let's just assume we don't do photon splatting ...
    probB[i + j] = probB[i + j - 1];    // ... so the probability for the camRay (both F and B) is just 1

    // Finally, fill in the missing values for probB in the front half of the concatenated path
    for (int k = i - 1; k > -1; k--) {
        float reverseReflectProb = 0;
        if (k == i - 1) { // As before, this one is special due to the outgoing direction being shadow ray
            reverseReflectProb = lightPath.m_hit[k].object->material()->reflectPDF(shadow_LtoE.d, lightPath.m_hit[k].N, -lightPath.m_ray[k].d);
        }
        else {
            reverseReflectProb = lightPath.m_hit[k].object->material()->reflectPDF(lightPath.m_ray[k + 1].d, lightPath.m_hit[k].N, -lightPath.m_ray[k].d);
        }
        if (k == 0) { // The first probB[0] has pdf conversion cos(theta) defined by the hitpoint on the light
            probB[k] = probB[k + 1] * reverseReflectProb * (std::max(0.0f, dot(lightPath.m_lightHit.N, lightPath.m_ray[0].d) / length2[0]));
        }
        else {
            probB[k] = probB[k + 1] * reverseReflectProb * (cosF[k - 1] / length2[k]);
        }
    }
    
    // Fill in the contribution to the estimator from the excluded shadow-connection
    // Note: the units for these two are NOT the same (these are computed only for convenience)
    // For the first, one still needs to multiply by light-wattage
    // For the latter, one still needs to multiply by the photon-power-density
    // This is, of course, in addition to the accumlated estimator values
    if (explicitConnection == true) {
        estimatorLink = brdf_L*brdf_E*(cosF_L * cosB_E / shadowLength2); // in paranthesis is just the form factor
    }
    else {
        estimatorLink = mat_E->BRDF(-lightPath.m_ray[i].d, hit_E.N, -eyePath.m_ray[j - 1].d);
    }

    // The following is a paranoid check. TODO: remove this when we are sure our indexing is correct
    for (int index = 0; index < i + j; index++) {
        if (probF[index] == 0.0f || probB[index] == 0.0f || cosF[index] == 0.0f || cosB[index] == 0.0f || length2[index] == 0.0f) {
            cout << "Missed a value" << endl;
        }
    }
    return true;
}
