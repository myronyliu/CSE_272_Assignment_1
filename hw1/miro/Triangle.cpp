#include "Triangle.h"
#include "Ray.h"


Triangle::Triangle(TriangleMesh * m, unsigned int i) : m_mesh(m), m_index(i) {}


Triangle::~Triangle() {}


void
Triangle::renderGL()
{
    TriangleMesh::TupleI3 ti3 = m_mesh->vIndex(m_index);
    const Vector3 & v0 = m_mesh->vertex(ti3.m_x); //vertex a of triangle
    const Vector3 & v1 = m_mesh->vertex(ti3.m_y); //vertex b of triangle
    const Vector3 & v2 = m_mesh->vertex(ti3.m_z); //vertex c of triangle
    
    glColor3f(1, 1, 1);
    glPushMatrix();
    glBegin(GL_TRIANGLES);
    glVertex3f(v0.x, v0.y, v0.z);
    glVertex3f(v1.x, v1.y, v1.z);
    glVertex3f(v2.x, v2.y, v2.z);
    glEnd();
    glPopMatrix();
}

bool
Triangle::intersect(HitInfo& result, const Ray& r,float tMin, float tMax)
{
    TriangleMesh::TupleI3 ti3 = m_mesh->vIndex(m_index);
    Vector3 A = m_mesh->vertex(ti3.m_x); //vertex a of triangle
    Vector3 B = m_mesh->vertex(ti3.m_y); //vertex b of triangle
    Vector3 C = m_mesh->vertex(ti3.m_z); //vertex c of triangle
    Vector3 AB = B - A;
    Vector3 AC = C - A;
    Vector3 ABxAC = cross(AB, AC); // un-normalized surface normal

    float det = -dot(r.d, ABxAC);
    if (det == 0) return false;

    Vector3 b = r.o - A;

    float alpha = -dot(r.d, cross(b, AC)) / det;
    float beta = -dot(r.d, cross(AB, b)) / det;

    if (alpha < 0 || beta < 0 || alpha + beta>1) return false;

    float nSteps = dot(b, ABxAC) / det;
    if (nSteps<tMin || nSteps>tMax) return false;

    result.t = nSteps;
    result.N = ABxAC.normalized();
    if (dot(r.d, result.N) > 0) result.N = -result.N;
    result.P = r.o + result.t*r.d;
    result.object = this;
    return true;
}
