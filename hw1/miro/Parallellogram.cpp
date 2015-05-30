#include "Parallelogram.h"
#include "Ray.h"
#include "Console.h"

Parallelogram::Parallelogram() {}
Parallelogram::Parallelogram(
    const Vector3& center,
    const Vector3& vecX, const Vector3& vecY,
    const float& spanX, const float& spanY)
{
    m_center = center;
    m_vecX = vecX.normalized();
    m_vecY = vecY.normalized();
    m_spanX = spanX;
    m_spanY = spanY;
    m_front = true;
    m_back = true;
}
Parallelogram::Parallelogram(
    const Vector3& center,
    const Vector3& vecX, const Vector3& vecY,
    const float& spanX, const float& spanY,
    const bool& front, const bool& back)
{
    m_center = center;
    m_vecX = vecX.normalized();
    m_vecY = vecY.normalized();
    m_spanX = spanX;
    m_spanY = spanY;
    m_front = front;
    m_back = back;
}

Parallelogram::~Parallelogram()
{
}

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
    float a = 0.9;
    Vector3 n00 = a*n - (1 - a)*pt00*s;
    Vector3 n01 = a*n - (1 - a)*pt01*s;
    Vector3 n10 = a*n - (1 - a)*pt10*s;
    Vector3 n11 = a*n - (1 - a)*pt11*s;
    glBegin(GL_QUADS);
    glVertex3f(pt00[0], pt00[1], pt00[2]);
    glVertex3f(pt01[0], pt01[1], pt01[2]);
    glVertex3f(pt11[0], pt11[1], pt11[2]);
    glVertex3f(pt10[0], pt10[1], pt10[2]);
    glEnd();
    glBegin(GL_LINES);
    glVertex3f(0, 0, 0);
    glVertex3f(n[0], n[1], n[2]);
    glEnd();
    glBegin(GL_TRIANGLES);
    glVertex3f(n[0], n[1], n[2]);
    glVertex3f(n00[0], n00[1], n00[2]);
    glVertex3f(n11[0], n11[1], n11[2]);
    glVertex3f(n[0], n[1], n[2]);
    glVertex3f(n01[0], n01[1], n01[2]);
    glVertex3f(n10[0], n10[1], n10[2]);
    glEnd();
    glPopMatrix();
}

bool
Parallelogram::intersect(HitInfo& result, const Ray& ray,
                  float tMin, float tMax)
{
    Vector3 RtoP_perp = (m_center - ray.o).orthogonal(m_vecX, m_vecY); // perpendicular displacement vector from ray to plane
    if (dot(RtoP_perp, normal()) >= 0 && m_back == false) return false;
    if (dot(RtoP_perp, normal()) <= 0 && m_front == false) return false;
    Vector3 dir = ray.d;
    Vector3 step_perp = dir.orthogonal(m_vecX, m_vecY); // perpendicular (relative to plane) step for ray
    if (dot(RtoP_perp, step_perp) <= 0 || step_perp.length() == 0.0) return false;
    float nSteps = RtoP_perp.length() / step_perp.length();
    Vector3 r = ray.o + nSteps*ray.d; // location where the ray hits in the (extended) plane

    Vector3 cr = r - m_center; // vector from center to r
    Vector3 YperpX = m_vecY.orthogonal(m_vecX).normalized();
    Vector3 cr_x = cr - (dot(cr, YperpX) / dot(m_vecY, YperpX))*m_vecY; // find componenets along skew-basis
    Vector3 cr_y = cr - cr_x;
    if (cr_x.length() > m_spanX || cr_y.length() > m_spanY) return false; // ray hits outside of plane
    if (nSteps<tMin || nSteps>tMax) return false;

    Vector3 eps = 0.000001*(ray.o - r).normalize();
    eps = Vector3(0, 0, 0);
    result.t = nSteps;
    result.N = -step_perp.normalized();
    result.P = r + eps;
    result.object = this;
    return true;
}

std::pair<Vector3, Vector3> Parallelogram::axisAlignedBounds() {
    Vector3 diagP = m_spanX*m_vecX + m_spanY*m_vecY;
    Vector3 diagM = m_spanX*m_vecX - m_spanY*m_vecY;
    Vector3 dCenter;
    dCenter.x = fmax(fabs(diagP.x), fabs(diagM.x));
    dCenter.y = fmax(fabs(diagP.y), fabs(diagM.y));
    dCenter.z = fmax(fabs(diagP.z), fabs(diagM.z));
    return std::pair<Vector3, Vector3>(m_center - dCenter, m_center + dCenter);
}