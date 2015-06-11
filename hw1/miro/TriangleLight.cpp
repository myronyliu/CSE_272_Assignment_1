#define _USE_MATH_DEFINES
#include "TriangleLight.h"

TriangleLight::TriangleLight(TriangleMesh * m, unsigned int i)
{
    m_mesh = m;
    m_index = i;
    m_front = true; // one-sided light
    m_back = false;
    m_color = Vector3(1, 1, 1);
    m_wattage = 100;
    Lambert * mat = new Lambert(Vector3(1.0, 1.0, 1.0));
    mat->setKd(1);
    mat->setEmittance(1.0);
    mat->setPowerPerArea(m_wattage*m_color / area());
    setMaterial(mat);
}


RayPDF TriangleLight::randRay() const {
    vec3pdf o = Triangle::randPt();
    vec3pdf d = m_material->randEmit(Triangle::normal(o.v));
    return RayPDF(Ray(o.v, d.v),o.p,d.p);
}

float TriangleLight::rayPDF(const Ray& ray) const {
    Vector3 b = barycentric(ray.o);
    if (b[0] < 0 || b[1] < 0 || b[2] < 0 || b[0]>1 || b[1]>1 || b[2]>1) return 0;
    else return m_material->emitPDF(Triangle::normal(ray.o), ray.d) / Triangle::area();
}