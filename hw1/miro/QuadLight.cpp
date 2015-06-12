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