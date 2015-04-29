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
    m_vecX = vecX;
    m_vecY = vecY;
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
    m_vecX = vecX;
    m_vecY = vecY;
    m_spanX = spanX;
    m_spanY = spanY;
    m_front = front;
    m_back = back;
}

Parallelogram::~Parallelogram()
{
}

void
Parallelogram::renderGL()
{
    glColor3f(1, 1, 1);
    glPushMatrix();
    glTranslatef(m_center.x, m_center.y, m_center.z);
    Vector3 pt00 = -m_spanX*m_vecX - m_spanY*m_vecY;
    Vector3 pt01 = -m_spanX*m_vecX + m_spanY*m_vecY;
    Vector3 pt10 = m_spanX*m_vecX - m_spanY*m_vecY;
    Vector3 pt11 = m_spanX*m_vecX + m_spanY*m_vecY;
    glBegin(GL_QUADS);
    glVertex3f(pt00[0], pt00[1], pt00[2]);
    glVertex3f(pt01[0], pt01[1], pt01[2]);
    glVertex3f(pt11[0], pt11[1], pt11[2]);
    glVertex3f(pt10[0], pt10[1], pt10[2]);
    glEnd();
    //glutWireSphere(sqrt(m_spanX*m_spanX+m_spanY*m_spanY)/100, 20, 20);
    glPopMatrix();
}

bool
Parallelogram::intersect(HitInfo& result, const Ray& ray,
                  float tMin, float tMax)
{
    if (m_back == false && m_front == false) return false;
    Vector3 RtoP_perp = (m_center - ray.o).orthogonal(m_vecX, m_vecY); // perpendicular displacement vector from ray to plane
    if (dot(RtoP_perp, normal()) >= 0 && m_back == false) return false;
    if (dot(RtoP_perp, normal()) <= 0 && m_front == false) return false;
    Vector3 dir = ray.d;
    Vector3 step_perp = dir.orthogonal(m_vecX, m_vecY); // perpendicular (relative to plane) step for ray
    if (dot(RtoP_perp, step_perp) < 0) return false;
    float nSteps = RtoP_perp.length() / step_perp.length();
    Vector3 r = ray.o + nSteps*ray.d; // location where the ray hits in the (extended) plane

    float r_x = dot(r-m_center, m_vecX.normalized()); // find components along skew-basis
    float r_y = dot(r-m_center, m_vecY.normalized());

    if (fabs(r_x) > m_spanX || fabs(r_y) > m_spanY) return false; // ray hits outside of plane
    if (nSteps<tMin || nSteps>tMax) return false;

    result.t = nSteps; 
    result.N = -step_perp.normalized();
    //result.P = r;
    result.P = r + (ray.o-result.P)*0.000001;
    result.material = this->m_material; 
    
    return true;
}