#include <stdio.h>
#include <stdlib.h>
#include "Miro.h"
#include "Camera.h"
#include "Image.h"
#include "Scene.h"
#include "Console.h" 
#include "OpenGL.h"

Camera * g_camera = 0;

static bool firstRayTrace = true; 

const float HalfDegToRad = DegToRad/2.0f;

Camera::Camera() :
    m_bgColor(0,0,0),
    m_renderer(RENDER_OPENGL),
    m_eye(0,0,0),
    m_viewDir(0,0,-1),
    m_up(0,1,0),
    m_lookAt(FLT_MAX, FLT_MAX, FLT_MAX),
    m_fov((45.)*(PI/180.))
{
    calcLookAt();
}


Camera::~Camera()
{

}


void
Camera::click(Scene* pScene, Image* pImage)
{
    calcLookAt();
    static bool firstRayTrace = false;

    if (m_renderer == RENDER_OPENGL)
    {
        glDrawBuffer(GL_BACK);
        pScene->openGL(this);
        firstRayTrace = true;
    }
    else if (m_renderer == RENDER_RAYTRACE)
    {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        glDrawBuffer(GL_FRONT);
        if (firstRayTrace)
        {
            pImage->clear(bgColor());
            pScene->raytraceImage(this, g_image);
            firstRayTrace = false;
        }
        
        g_image->draw();
    }
    else if (m_renderer == RENDER_PHOTONTRACE)
    {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        glDrawBuffer(GL_FRONT);
        if (firstRayTrace)
        {
            pImage->clear(bgColor());
            pScene->photontraceImage(this, g_image);
            firstRayTrace = false;
        }

        g_image->draw();
    }
    else if (m_renderer == RENDER_PATHTRACE)
    {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        glDrawBuffer(GL_FRONT);
        if (firstRayTrace)
        {
            pImage->clear(bgColor());
            pScene->pathtraceImage(this, g_image);
            firstRayTrace = false;
        }

        g_image->draw();
    }
}


void
Camera::calcLookAt()
{
    // this is true when a "lookat" is not used in the config file
    if (m_lookAt.x != FLT_MAX)
    {
        setLookAt(m_lookAt);
        m_lookAt.set(FLT_MAX, FLT_MAX, FLT_MAX);
    }
}


void
Camera::drawGL()
{
    // set up the screen with our camera parameters
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(fov(), g_image->width()/(float)g_image->height(),
                   0.01, 10000);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    Vector3 vCenter = eye() + viewDir();
    gluLookAt(eye().x, eye().y, eye().z,
              vCenter.x, vCenter.y, vCenter.z,
              up().x, up().y, up().z);
}

Ray
Camera::eyeRay(float x, float y, int imageWidth, int imageHeight)
{
    // first compute the camera coordinate system 
    // ------------------------------------------

    // wDir = e - (e+m_viewDir) = -m_vView
    const Vector3 wDir = Vector3(-m_viewDir).normalize(); 
    const Vector3 uDir = cross(m_up, wDir).normalize(); 
    const Vector3 vDir = cross(wDir, uDir);    



    // next find the corners of the image plane in camera space
    // --------------------------------------------------------

    const float aspectRatio = (float)imageWidth/(float)imageHeight; 


    const float top     = tan(m_fov*HalfDegToRad); 
    const float right   = aspectRatio*top; 

    const float bottom  = -top; 
    const float left    = -right; 



    // transform x and y into camera space 
    // -----------------------------------

    const float imPlaneUPos = left   + (right - left)*((x+0.5f)/(float)imageWidth); 
    const float imPlaneVPos = bottom + (top - bottom)*((y+0.5f)/(float)imageHeight); 

    return Ray(m_eye, (imPlaneUPos*uDir + imPlaneVPos*vDir - wDir).normalize());
}

Ray Camera::eyeRayJittered(float x, float y, int imageWidth, int imageHeight){
    float dx = 0.5 - (double)rand() / (double)RAND_MAX;
    float dy = 0.5 - (double)rand() / (double)RAND_MAX;
    return eyeRay(x + dx, y + dy, imageWidth, imageHeight);
}

Parallelogram Camera::imagePlane(int imageWidth, int imageHeight){
    float W = imageWidth;
    float H = imageHeight;
    const Vector3 wDir = Vector3(-m_viewDir).normalize();
    const Vector3 uDir = cross(m_up, wDir).normalize();
    const Vector3 vDir = cross(wDir, uDir);
    const float aspectRatio = W / H;
    const float top = tan(m_fov*HalfDegToRad);
    const float right = aspectRatio*top;
    const float bottom = -top;
    const float left = -right;
    const float imPlaneUPos00 = left + (right - left)*(0.5f / W);
    const float imPlaneVPos00 = bottom + (top - bottom)*(0.5f / H);
    const float imPlaneUPos01 = left + (right - left)*(0.5f / W);
    const float imPlaneVPos01 = bottom + (top - bottom)*(((H - 1) + 0.5f) / H);
    const float imPlaneUPos10 = left + (right - left)*(((W - 1) + 0.5f) / W);
    const float imPlaneVPos10 = bottom + (top - bottom)*(0.5f / H);
    const float imPlaneUPos11 = left + (right - left)*(((W - 1) + 0.5f) / W);
    const float imPlaneVPos11 = bottom + (top - bottom)*(((H - 1) + 0.5f) / H);
    Vector3 c00 = imPlaneUPos00*uDir + imPlaneVPos00*vDir - wDir;
    Vector3 c01 = imPlaneUPos01*uDir + imPlaneVPos01*vDir - wDir;
    Vector3 c10 = imPlaneUPos10*uDir + imPlaneVPos10*vDir - wDir;
    Vector3 c11 = imPlaneUPos11*uDir + imPlaneVPos11*vDir - wDir;
    Parallelogram p;
    p.setCenter(m_eye + (c00 + c11) / 2.0);
    p.setVecX((c10 - c00).normalize());
    p.setVecX((c01 - c00).normalize());
    p.disableFront();
    p.enableBack();
    Material* mat = new Lambert(Vector3(1, 1, 1));
    p.setMaterial(mat);
    return p;
}

std::vector<float>
Camera::imgProject(const Vector3& pt,int imageWidth, int imageHeight){
    float W = imageWidth;
    float H = imageHeight;
    const Vector3 wDir = Vector3(-m_viewDir).normalize();
    const Vector3 uDir = cross(m_up, wDir).normalize();
    const Vector3 vDir = cross(wDir, uDir);
    const float aspectRatio = W / H;
    const float T = tan(m_fov*HalfDegToRad);
    const float R = aspectRatio*T;
    const float B = -T;
    const float L = -R;
    float u = dot(pt, uDir);
    float v = dot(pt, vDir);
    float w = dot(pt, -wDir); // this should be positive
    u /= w;
    v /= w;
    float px = W*(u - L) / (R - L) - 0.5;
    float py = H*(v - B) / (T - B) - 0.5;
    std::vector<float> pt2D = { px, py };
    return pt2D;
}