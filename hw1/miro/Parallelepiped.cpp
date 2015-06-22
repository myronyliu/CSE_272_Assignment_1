#include "Parallelepiped.h"
#include "Ray.h"
#include "Matrix3x3.h"
#include "Parallelogram.h"
#include "Console.h"

Parallelepiped::Parallelepiped(
    const Vector3& center,
    const Vector3& vecX, const Vector3& vecY, const Vector3& vecZ,
    const float& spanX, const float& spanY, const float& spanZ) :
    m_center(center),
    m_vecX(vecX.normalized()),
    m_vecY(vecY.normalized()),
    m_vecZ(vecZ.normalized()),
    m_spanX(spanX),
    m_spanY(spanY),
    m_spanZ(spanZ)
{
    if (dot(cross(vecX, vecY), vecZ) < 0) m_vecZ = -vecZ;
}


vec3pdf Parallelepiped::randPt() const {
    double rx = 1.0 - 2.0*(double)rand() / (double)RAND_MAX;
    double ry = 1.0 - 2.0*(double)rand() / (double)RAND_MAX;
    Vector3 vx = rx*spanX()*vecX();
    Vector3 vy = ry*spanY()*vecY();
    double area = 4.0*m_spanX*m_spanY*cross(m_vecX, m_vecY).length();
    return vec3pdf(m_center + vx + vy, 1.0 / area);
}

Vector3 Parallelepiped::normal(const Vector3& v) const {
    Vector3 r = v - m_center;
    Vector3 rAligned = Matrix3x3(m_vecX, m_vecY, m_vecZ).invert()*r;
    float dx = fabs(rAligned.x) - 1.0f;
    float dy = fabs(rAligned.x) - 1.0f;
    float dz = fabs(rAligned.x) - 1.0f;
    bool interiorFlag = true;
    if (dx >= 0 || dy >= 0 || dz >= 0) interiorFlag = false;
    dx = fabs(dx);
    dy = fabs(dy);
    dz = fabs(dz);
    Vector3 n;
    if (dx < dy && dx < dz) {
        if (rAligned.x < 0) n = -cross(m_vecY, m_vecZ).normalize();
        else n = cross(m_vecY, m_vecZ).normalize();
    }
    else if (dy < dz && dy < dx) {
        if (rAligned.y < 0) n = -cross(m_vecZ, m_vecX).normalize();
        else n = cross(m_vecZ, m_vecX).normalize();
    }
    else {
        if (rAligned.z < 0) n = -cross(m_vecX, m_vecY).normalize();
        else n = cross(m_vecX, m_vecY).normalize();
    }
    if (interiorFlag == true) n = -n;
    return n;
}

void
Parallelepiped::renderGL() {
    glColor3f(1, 1, 1);
    glPushMatrix();
    glTranslatef(m_center.x, m_center.y, m_center.z);
    Vector3 pt000 = -m_spanX*m_vecX - m_spanY*m_vecY - m_spanZ*m_vecZ;
    Vector3 pt001 = -m_spanX*m_vecX - m_spanY*m_vecY + m_spanZ*m_vecZ;
    Vector3 pt010 = -m_spanX*m_vecX + m_spanY*m_vecY - m_spanZ*m_vecZ;
    Vector3 pt011 = -m_spanX*m_vecX + m_spanY*m_vecY + m_spanZ*m_vecZ;
    Vector3 pt100 = m_spanX*m_vecX - m_spanY*m_vecY - m_spanZ*m_vecZ;
    Vector3 pt101 = m_spanX*m_vecX - m_spanY*m_vecY + m_spanZ*m_vecZ;
    Vector3 pt110 = m_spanX*m_vecX + m_spanY*m_vecY - m_spanZ*m_vecZ;
    Vector3 pt111 = m_spanX*m_vecX + m_spanY*m_vecY + m_spanZ*m_vecZ;

    glBegin(GL_QUADS);

    glVertex3f(pt000[0], pt000[1], pt000[2]);
    glVertex3f(pt001[0], pt001[1], pt001[2]);
    glVertex3f(pt011[0], pt011[1], pt011[2]);
    glVertex3f(pt010[0], pt010[1], pt010[2]);

    glVertex3f(pt100[0], pt100[1], pt100[2]);
    glVertex3f(pt101[0], pt101[1], pt101[2]);
    glVertex3f(pt111[0], pt111[1], pt111[2]);
    glVertex3f(pt110[0], pt110[1], pt110[2]);

    glVertex3f(pt000[0], pt000[1], pt000[2]);
    glVertex3f(pt001[0], pt001[1], pt001[2]);
    glVertex3f(pt101[0], pt101[1], pt101[2]);
    glVertex3f(pt100[0], pt100[1], pt100[2]);

    glVertex3f(pt010[0], pt010[1], pt010[2]);
    glVertex3f(pt011[0], pt011[1], pt011[2]);
    glVertex3f(pt111[0], pt111[1], pt111[2]);
    glVertex3f(pt110[0], pt110[1], pt110[2]);

    glVertex3f(pt000[0], pt000[1], pt000[2]);
    glVertex3f(pt010[0], pt010[1], pt010[2]);
    glVertex3f(pt110[0], pt110[1], pt110[2]);
    glVertex3f(pt100[0], pt100[1], pt100[2]);

    glVertex3f(pt001[0], pt001[1], pt001[2]);
    glVertex3f(pt011[0], pt011[1], pt011[2]);
    glVertex3f(pt111[0], pt111[1], pt111[2]);
    glVertex3f(pt101[0], pt101[1], pt101[2]);

    glEnd();

    glPopMatrix();
}

bool
Parallelepiped::intersect(HitInfo& result, const Ray& ray, float tMin, float tMax)
{
    Parallelogram faces[6];

    faces[0] = Parallelogram(m_center - m_spanX*m_vecX, m_vecZ, m_vecY, m_spanZ, m_spanY);
    faces[1] = Parallelogram(m_center + m_spanX*m_vecX, m_vecY, m_vecZ, m_spanY, m_spanZ);

    faces[2] = Parallelogram(m_center - m_spanY*m_vecY, m_vecX, m_vecZ, m_spanX, m_spanZ);
    faces[3] = Parallelogram(m_center + m_spanY*m_vecY, m_vecZ, m_vecX, m_spanZ, m_spanX);

    faces[4] = Parallelogram(m_center - m_spanZ*m_vecZ, m_vecY, m_vecX, m_spanY, m_spanX);
    faces[5] = Parallelogram(m_center + m_spanZ*m_vecZ, m_vecX, m_vecY, m_spanX, m_spanY);

    HitInfo hit;
    for (int i = 0; i < 6; i++) {
        if (faces[i].intersect(hit, ray, tMin, tMax)) {
            result.t = hit.t;
            result.P = hit.P;
            result.N = hit.N;
            result.object = this;
            for (int j = i; j < 6; j++) {
                if (!faces[j].intersect(hit, ray, tMin, tMax)) continue;
                if (hit.t < result.t) {
                    result.t = hit.t;
                    result.P = hit.P;
                    result.N = hit.N;
                }
            }
            return true;
        }
    }
    return false;
}

std::pair<Vector3, Vector3> Parallelepiped::axisAlignedBounds() {
    return std::pair<Vector3, Vector3>(Vector3(0, 0, 0), Vector3(0, 0, 0));
}