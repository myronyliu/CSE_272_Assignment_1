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
    m_front = true; // if you want back-side of light, create another light with opposing normal
    m_back = false;
    m_color = Vector3(1, 1, 1);
    m_wattage = 100;
    Material* mat = new Lambert(Vector3(1.0, 1.0, 1.0));
    mat->setEmittance(1.0);
    mat->setEmitted(Vector3(1.0, 0.0, 0.0));
    setMaterial(mat);
}


Vector3 ParallelogramLight::randPt() const {
    double rx = 1.0 - 2.0*(double)rand() / (double)RAND_MAX;
    double ry = 1.0 - 2.0*(double)rand() / (double)RAND_MAX;
    Vector3 vx = rx*spanX()*vecX();
    Vector3 vy = ry*spanY()*vecY();
    return m_center + vx + vy;
}