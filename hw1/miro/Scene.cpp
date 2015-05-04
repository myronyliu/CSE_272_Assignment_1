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
        double em = hit.material->emittance();
        rayToEye.o = hit.P;
        rayToEye.d = (cam->eye() - hit.P).normalize();
        float brdf = hit.material->BRDF(-rayIn.d, hit.N, rayToEye.d); // BRDF between eye-ray and incoming-ray
        if (!trace(tryHitEye, rayToEye)) { // check if anything is occluding the eye from current hitpoint
            pix = cam->imgProject(hit.P, w, h); // find the pixel the onto which the current hitpoint projects
            int x = round(pix[0]);
            int y = round(pix[1]);
            //img->setPixel(x, y, Vector3(1.0, 0.0, 0.0)); return;
            if (dot(Vector3(0, 0, -1), hit.N) != 0) {  }
            if (pix[2]>0 && x >= 0 && x<w && y>=0 && y < h) { // check that the pixel is within the viewing window
                float cos0 = dot(hit.N, rayToEye.d);
                float lengthSqr = (cam->eye() - hit.P).length2();
                float cosAlpha = cam->pixelCosine(pix[0], pix[1], w, h);
                if (em == 1.0 || rn < 0.2) {
                    img->setPixel(x, y, img->getPixel(x, y) + (1 /0.2) * power*brdf*cos0*cosAlpha / lengthSqr);
                    return; // photon was absorbed
                } // otherwise photon will be reflected
                else
                {
                    img->setPixel(x, y, img->getPixel(x, y) + (1/ 0.8) * power*brdf*cos0*cosAlpha / lengthSqr);
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
    for (int p = 0; p < m_photonSamples; p++) { // shoot a photon...
        LightPDF lp = randLightByWattage(); // ... off of a random light (I don't think we need the PDF here)
        Light* light = lp.l;
        RayPDF rp = light->randRay();
        tracePhoton(cam, img, lp, rp);
        if (p % 100 == 0) {
            printf("Rendering Progress: %.3f%%\r", p / float(m_photonSamples)*100.0f);
            fflush(stdout);
        }
        if (((float) p/m_photonSamples)<0.8 && p>0 && p % (m_photonSamples / 3) == 0)
        {
            img->draw();
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
    printf("Rendering Progress: 100.000%\n");
    debug("done Photontracing!\n");
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


void
Scene::biditraceImage(Camera *cam, Image *img)
{
    int w = img->width();
    int h = img->height();
    HitInfo hitInfo;

    float W = 0.5;

    for (int y = 0; y < h; y++)
    //for (int y = h-1; y >= 0; y--)
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
            Vector3 fluxSum(0, 0, 0);
            for (int k = 0; k < bidiSamplesPerPix(); k++){
                RayPath eyePath = randEyePath(x, y, cam, img);
                if (eyePath.m_hits.size() == 0) {
                    continue;
                }
                RayPath lightPath = randLightPath();

                int dj = eyePath.m_hits.size();
                int di = lightPath.m_hits.size();
                vector<Vector3> fixedLengthFlux(di + dj - 1, Vector3(0, 0, 0));
                vector<float> count(di + dj - 1, 0);

                fixedLengthFlux[0] = eyePath.m_hits[0].material->radiance(eyePath.m_hits[0].N, -eyePath.m_rays[0].d);
                count[0] = 1;
                for (unsigned int i = 0; i < di; i++){
                    for (unsigned int j = 1; j < dj; j++) {
                        fixedLengthFlux[i + j] += estimateFlux(i, j, eyePath, lightPath);
                        count[i + j] += 1;
                    }
                }
                for (int i = 0; i < di + dj - 1; i++) fluxSum += fixedLengthFlux[i] / count[i];
            }
            float cosSA = cam->pixelCosine(x, y, img->width(), img->height())*cam->pixelSolidAngle(x, y, img->width(), img->height());
            img->setPixel(x, y, fluxSum / bidiSamplesPerPix()/cosSA);
        }
        img->drawScanline(y);
        glFinish();
        printf("Rendering Progress: %.3f%%\r", y / float(img->height())*100.0f);
        fflush(stdout);
    }

    printf("Rendering Progress: 100.000%\n");
    debug("Done Bidi Pathtracing!\n");

}

RayPath Scene::randEyePath(float x, float y, Camera* cam, Image* img) {
    RayPath raypath(cam->eyeRay(x, y, img->width(), img->height()));
    return generateRayPath(raypath);
}

RayPath Scene::randLightPath() {
    LightPDF lp = randLightByWattage();
    Light* light = lp.l;

    RayPDF rp = light->randRay();
    RayPath raypath(rp.r);
    raypath.m_fluxDecay.push_back(1.0);

    raypath.m_probs.push_back(1.0); /////////

    HitInfo first_hit(0.0f, rp.r.o, light->normal(rp.r.o));
    raypath.m_hits.push_back(first_hit);
    first_hit.material = light->material();
    raypath.m_light = light;

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
        float cos = dot(lastHit.N, newRay.d);
        raypath.m_rays.push_back(newRay);
        raypath.m_hits.push_back(hit);
        raypath.m_probs.push_back(raypath.m_probs.back() * vp.p);
        if (raypath.m_fluxDecay.size() == 0) raypath.m_fluxDecay.push_back(brdf*cos); // for eyePath
        else raypath.m_fluxDecay.push_back(brdf*cos*raypath.m_fluxDecay.back()); // for lightPath

        bounce++;
    }
    return raypath;
}

Vector3 Scene::estimateFlux(int i, int j, RayPath eyePath, RayPath lightPath) {
    Vector3 flux = Vector3(0, 0, 0);
    HitInfo h;
    
    if (i == 0){
        Vector3 lightPoint = lightPath.m_hits[0].P;
        Ray rayEye = eyePath.m_rays[j-1];
        HitInfo hit = eyePath.m_hits[j-1];
        const Material* mat = hit.material;
        Ray rayShadow(lightPoint, (hit.P - lightPoint).normalize()); // visibility ray FROM eye point TO light
        float shadowLength2 = (hit.P - lightPoint).length2();
        if (trace(h, rayShadow, 0.00001, sqrt(shadowLength2) - 0.00001)) { // plus epsilon because we want to intersect cover
            return flux;
        }
        if (dot(lightPath.m_hits[0].N, rayShadow.d) <= 0) {
            return flux; // hits backside of light
        }
        //float brdfChain = pow(mat->BRDF(rayShadow.d, hit.N, -rayEye.d), j-1);
        float brdf = mat->BRDF(rayShadow.d, hit.N, -rayEye.d);
        float form = dot(lightPath.m_hits[0].N, -rayShadow.d)*dot(hit.N, rayShadow.d) / shadowLength2;
        flux = lightPath.m_light->wattage() * brdf *  form;
        if (j > 1){
            flux *= eyePath.m_fluxDecay[j - 2];
        }
    }
    else
    {
        HitInfo hitj = eyePath.m_hits[j - 1];
        HitInfo hiti = lightPath.m_hits[i];
        const Material* mati = hiti.material;
        const Material* matj = hitj.material;
        Ray rayShadow(hitj.P, (hiti.P - hitj.P).normalize());
        float shadowLength2 = (hiti.P - hitj.P).length2();
        if (!trace(h, rayShadow, 0.00001, sqrt(shadowLength2) - 0.00001)) {
            //float brdfChain = pow(mati->BRDF(-rayShadow.d, hiti.N, -lightPath.m_rays[i - 1].d), i + j - 2);
            float brdfi = mati->BRDF(-rayShadow.d, hiti.N, -lightPath.m_rays[i - 1].d);
            float brdfj = matj->BRDF(rayShadow.d, hitj.N, -eyePath.m_rays[j - 1].d);
            float form = dot(-rayShadow.d, hiti.N)*dot(rayShadow.d, hitj.N) / shadowLength2;
            flux = lightPath.m_light->wattage() * brdfi * brdfj * form;
            flux *= lightPath.m_fluxDecay[i - 1] * eyePath.m_fluxDecay[j - 1];
        }
    }
    return flux;
}
