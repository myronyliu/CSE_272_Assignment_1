#define _USE_MATH_DEFINES
#include "QuadLight.h"

QuadLight::QuadLight(PolygonMesh * m, unsigned int i)
{
    m_mesh = m;
    m_index = i;
    m_front = true; // one-sided light
    m_back = false;
    m_color = Vector3(1, 1, 1);
    m_wattage = 100;
    Lambert * mat = new Lambert(Vector3(1.0, 1.0, 1.0));
    mat->setKd(0);
    mat->setEmittance(1.0);
    mat->setPowerPerArea(m_wattage*m_color / area());
    setMaterial(mat);
}


RayPDF QuadLight::randRay() const {
    vec3pdf o = randPt();
    vec3pdf d = m_material->randEmit(normal(o.v));
    return RayPDF(Ray(o.v, d.v),o.p,d.p);
}

float QuadLight::rayPDF(const Ray& ray) const {
    return m_material->emitPDF(normal(ray.o), ray.d) / area();
}

void
QuadLight::renderGL() {
    Vector3 center = (A() + C()) / 2;
    glColor3f(1, 1, 1);
    glPushMatrix();
    //glTranslatef(center[0], center[1], center[2]);
    float s = 1.0 / 4.0;
    Vector3 n = m_mesh->normal(m_mesh->quadNormalIndex(m_index).m_a)*sqrt(area())*s;

    glBegin(GL_LINES);
    glVertex3f(center[0], center[1], center[2]);
    glVertex3f(center[0] + n[0], center[1] + n[1], center[2] + n[2]);
    glEnd();

    glBegin(GL_QUADS);
    glVertex3f(A()[0], A()[1], A()[2]);
    glVertex3f(B()[0], B()[1], B()[2]);
    glVertex3f(C()[0], C()[1], C()[2]);
    glVertex3f(D()[0], D()[1], D()[2]);
    glEnd();
    glPopMatrix();
}