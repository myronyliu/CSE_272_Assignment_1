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
        vp.p *= vp.v.length2() / fabs(dot(vp.v.normalized(),randAL->normal(lightPt))); // convert PDF from 1/A to 1/SA
    }
    Vector3 newDir = vp.v;
    Ray newRay;
    newRay.o = hit.P;
    newRay.d = newDir;
    float brdf = hit.material->BRDF(-ray.d, hit.N, newRay.d);
    float cos = dot(hit.N, newDir);
    Vector3 gather = hit.material->shade(ray, hit, *this); // gathered direct lighting
    return
            gather + (1.0 / (1.0 - em))/vp.p*brdf*cos*0.8*
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
                //if (dot(hitInfo.N, Vector3(0, 0, -1)) > 0) { pixSum += Vector3(1.0, 0.0, 0.0);  }
                pixSum += recursiveTrace_fromEye(ray, 0, m_maxBounces) / (double)m_samplesPerPix;
            }
            img->setPixel(i, j, pixSum);
        }
        img->drawScanline(j);
        glFinish();
        printf("Rendering Progress: %.3f%%\r", j / float(img->height())*100.0f);
        fflush(stdout);
    }

    printf("Rendering Progress: 100.000%\n");
    debug("done Raytracing!\n");
}

void Scene::tracePhoton(Camera *cam, Image *img, const Light& light, const RayPDF& rp) {
    float w = img->width();
    float h = img->height();
    Vector3 pix;
    HitInfo hit, tryHitEye;
    Ray rayIn, rayOut, rayToEye;
    rayOut = rp.r;
    Vector3 power = light.wattage()*light.color() / rp.p; // Monte Carlo sampling so divide by PDF of the random ray
    // The following is for the initial emission (to make the light visible)
/*    rayToEye.o = rayOut.o;
    rayToEye.d = (cam->eye() - rayToEye.o).normalize();
    if (!trace(tryHitEye, rayToEye)) { // check if anything is occluding the eye from current hitpoint
        pix = cam->imgProject(rayToEye.o, w, h); // find the pixel the onto which the current hitpoint projects
        int x = round(pix[0]);
        int y = round(pix[1]);
        if (pix[2]>0 && x > -1 && x<w && y>-1 && y < h) { // check that the pixel is within the viewing window
            Vector3 initValue = light.material()->radiance(light.normal(rayOut.o), rayToEye.d)*light.area();
            initValue *= dot(light.normal(rayOut.o), rayToEye.d);
            initValue *= cam->pixelCosine(pix[0], pix[1], w, h);
            initValue /= light.material()->sum_L_cosTheta_dOmega();
            initValue /= (rayOut.o - cam->eye()).length2();
            img[x][y] += power*initValue;
        }
    }*/
    // now for the bounces
    rayIn = rayOut;
    for (int bounces = 0; bounces < m_maxBounces; bounces++) {
        if (!trace(hit, rayIn)) {
            return; // ray left scene
        }
        double rn = (double)rand() / RAND_MAX;
        double em = hit.material->emittance();
        if (em == 1.0 || rn < 0.2) {
            return; // photon was absorbed
        } // otherwise photon will be reflected
        rayToEye.o = hit.P;
        rayToEye.d = (cam->eye() - hit.P).normalize();
        if (!trace(tryHitEye, rayToEye)) { // check if anything is occluding the eye from current hitpoint
            pix = cam->imgProject(hit.P, w, h); // find the pixel the onto which the current hitpoint projects
            int x = round(pix[0]);
            int y = round(pix[1]);
            if (pix[2]>0 && x > -1 && x<w && y>-1 && y < h) { // check that the pixel is within the viewing window
                float brdf = hit.material->BRDF(rayToEye.d, hit.N, -rayIn.d); // BRDF between eye-ray and incoming-ray
                float cos = dot(rayToEye.d, hit.N); // cosine between eye-ray and surface-normal
                float lengthSqr = (cam->eye() - hit.P).length2();
                float cosAlpha = cam->pixelCosine(pix[0], pix[1], w, h);
                img->setPixel(x, y, img->getPixel(x, y) + power/m_photonSamples*brdf*cos*cosAlpha / lengthSqr);
            }
        }
        vec3pdf vp = hit.material->randReflect(-rayIn.d, hit.N); // pick BRDF.cos weighted random direction
        Vector3 dirOut = vp.v;
        rayOut.o = hit.P;
        rayOut.d = dirOut;
        rayIn = rayOut;
        power /= vp.p; // Monte Carlo
    }
}


void
Scene::photontraceImage(Camera *cam, Image *img)
{
    int w = img->width();
    int h = img->height();
    for (int p = 0; p < m_photonSamples; p++) { // shoot a photon...
        LightPDF lp = randLightByWattage(); // ... off of a random light (I don't think we need the PDF here)
        Light* light = lp.l;
        RayPDF rp = light->randRay();
        tracePhoton(cam, img, *light, rp);
        if (p % 100 == 0) {
            printf("Rendering Progress: %.3f%%\r", p / float(m_photonSamples)*100.0f);
            fflush(stdout);
        }
        if (p>0 && p % (m_photonSamples / 5) == 0)
        {
            img->draw();
        }
    }
    printf("Rendering Progress: 100.000%\n");
    debug("done Photontracing!\n");
    for (int i = 0; i < w; i++){
        for (int j = 0; j < h; j++){
            Vector3 pix = img->getPixel(i, j);
            img->setPixel(i, j,
                pix / (cam->pixelCosine(i, j, w, h)*cam->pixelSolidAngle(i, j, w, h)));
        }
    }
    img->draw();
    printf("Rendering Progress: 100.000%\n");
    debug("done Raytracing!\n");
    glFinish();
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