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
    mat->setEmitted(Vector3(1.0, 1.0, 1.0));
    setMaterial(mat);
}

Vector3 ParallelogramLight::randPt() const {
    double rx = 1.0 - 2.0*(double)rand() / (double)RAND_MAX;
    double ry = 1.0 - 2.0*(double)rand() / (double)RAND_MAX;
    Vector3 vx = rx*spanX()*vecX();
    Vector3 vy = ry*spanY()*vecY();
    return m_center + vx + vy;
}

/*
void
Parallelogram::renderGL() {
    glColor3f(1, 1, 1);
    glPushMatrix();
    glTranslatef(m_center.x, m_center.y, m_center.z);
    Vector3 pt00 = -m_spanX*m_vecX - m_spanY*m_vecY;
    Vector3 pt01 = -m_spanX*m_vecX + m_spanY*m_vecY;
    Vector3 pt10 = m_spanX*m_vecX - m_spanY*m_vecY;
    Vector3 pt11 = m_spanX*m_vecX + m_spanY*m_vecY;
    float s = 1.0 / 4.0;
    Vector3 n = normal()*sqrt(m_spanX*m_spanX + m_spanY*m_spanY)*s;
    for (int i = 0; i < 5; i++){
        for (int j = 0; j < 5; j++){
            Vector3 p = (2.0*i / 5.0 - 1.0)*m_spanX*m_vecX + (2.0*i / 5.0 - 1.0)*m_spanY*m_vecY;
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
}*/