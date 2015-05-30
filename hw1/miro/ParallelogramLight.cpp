#define _USE_MATH_DEFINES
#include "ParallelogramLight.h"

ParallelogramLight::ParallelogramLight(const Vector3& center,
    const Vector3& vecX, const Vector3& vecY,
    const float& spanX, const float& spanY)
{
    m_center = center;
    m_vecX = vecX;
    m_vecY = vecY;
    m_spanX = spanX;
    m_spanY = spanY;
    m_front = true; // one-sided light
    m_back = false;
    m_color = Vector3(1, 1, 1);
    m_wattage = 100;
    Material* mat = new Lambert(Vector3(1.0, 1.0, 1.0));
    mat->setEmittance(1.0);
    mat->setPowerPerArea(m_wattage*m_color / area());
    setMaterial(mat);
}

vec3pdf ParallelogramLight::randPt() const {
    double rx = 1.0 - 2.0*(double)rand() / (double)RAND_MAX;
    double ry = 1.0 - 2.0*(double)rand() / (double)RAND_MAX;
    Vector3 vx = rx*spanX()*vecX();
    Vector3 vy = ry*spanY()*vecY();
    double area = 4.0*m_spanX*m_spanY*cross(m_vecX, m_vecY).length();
    return vec3pdf(m_center + vx + vy, 1.0 / area);
}

RayPDF ParallelogramLight::randRay() const {
    vec3pdf o = randPt();
    vec3pdf d = m_material->randEmit(cross(m_vecX,m_vecY).normalize());
    return RayPDF(Ray(o.v, d.v),o.p,d.p);
}

float ParallelogramLight::rayPDF(const Ray& ray) const {
    Vector3 cr = ray.o - m_center; // vector from center to ray origin
    Vector3 vecY = m_vecY;
    Vector3 YperpX = vecY.orthogonal(m_vecX).normalize();
    Vector3 cr_x = cr - (dot(cr, YperpX) / dot(m_vecY, YperpX))*m_vecY; // find componenets along skew-basis
    Vector3 cr_y = cr - cr_x;
    if (cr_x.length() > m_spanX || cr_y.length() > m_spanY) return 0; // ray originates from outside of parallelogram

    return m_material->emitPDF(cross(m_vecX, m_vecY).normalize(), ray.d) / area();
}


void
ParallelogramLight::renderGL() {
    glColor3f(1, 1, 1);
    glPushMatrix();
    glTranslatef(m_center.x, m_center.y, m_center.z);
    Vector3 pt00 = -m_spanX*m_vecX - m_spanY*m_vecY;
    Vector3 pt01 = -m_spanX*m_vecX + m_spanY*m_vecY;
    Vector3 pt10 = m_spanX*m_vecX - m_spanY*m_vecY;
    Vector3 pt11 = m_spanX*m_vecX + m_spanY*m_vecY;
    float s = 1.0 / 4.0;
    Vector3 n = cross(m_vecX,m_vecY).normalize()*sqrt(m_spanX*m_spanX + m_spanY*m_spanY)*s;
    int a = 4;
    for (int i = 0; i < a; i++){
        for (int j = 0; j < a; j++){
            Vector3 p = (2.0*i / (a - 1.0) - 1.0)*m_spanX*m_vecX + (2.0*j / (a - 1.0) - 1.0)*m_spanY*m_vecY;
            Vector3 q = n + 1.5*p;
            glBegin(GL_LINES);
            glVertex3f(p[0], p[1], p[2]);
            glVertex3f(q[0], q[1], q[2]);
            glEnd();
        }
    }
    glBegin(GL_QUADS);
    glVertex3f(pt00[0], pt00[1], pt00[2]);
    glVertex3f(pt01[0], pt01[1], pt01[2]);
    glVertex3f(pt11[0], pt11[1], pt11[2]);
    glVertex3f(pt10[0], pt10[1], pt10[2]);
    glEnd();
    glPopMatrix();
}