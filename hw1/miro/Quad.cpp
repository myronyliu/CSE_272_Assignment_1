#include "Quad.h"
#include "Ray.h"

#include "Parallelogram.h"


Quad::Quad(PolygonMesh * m, unsigned int i, const bool& front, const bool& back) : m_mesh(m), m_index(i), m_front(front), m_back(back) {}

Quad::~Quad() {}

vec3pdf Quad::randPt() const {
    float r1 = (float)rand() / RAND_MAX;
    float r2 = (float)rand() / RAND_MAX;
    return vec3pdf(r1*(B() - A()) + r2*(D() - A()), 1.0f / area());
}

Vector3 Quad::normal(const Vector3& ptInput) const {
    return cross(B() - A(), D() - A()).normalize(); // blargh, just do this for now

    /*Vector3 v0 = m_mesh->vertex(m_mesh->quadVertexIndex(m_index).m_w);
    Vector3 v1 = m_mesh->vertex(m_mesh->quadVertexIndex(m_index).m_x);
    Vector3 v2 = m_mesh->vertex(m_mesh->quadVertexIndex(m_index).m_y);
    Vector3 v3 = m_mesh->vertex(m_mesh->quadVertexIndex(m_index).m_z);

    Vector3 n0 = m_mesh->normal(m_mesh->quadNormalIndex(m_index).m_w);
    Vector3 n1 = m_mesh->normal(m_mesh->quadNormalIndex(m_index).m_x);
    Vector3 n2 = m_mesh->normal(m_mesh->quadNormalIndex(m_index).m_y);
    Vector3 n3 = m_mesh->normal(m_mesh->quadNormalIndex(m_index).m_z);

    Vector3 crossNormal = cross(v1 - v0, v3 - v0);
    Vector3 pt = pt - dot(ptInput - v0, crossNormal)*crossNormal;

    Vector3 A00*/
}

void
Quad::renderGL()
{
    PolygonMesh::TupleI4 ti4 = m_mesh->quadVertexIndex(m_index);
    Vector3 v0 = m_mesh->vertex(ti4.m_a); //vertex a of Quad
    Vector3 v1 = m_mesh->vertex(ti4.m_b); //vertex b of Quad
    Vector3 v2 = m_mesh->vertex(ti4.m_c); //vertex c of Quad
    Vector3 v3 = m_mesh->vertex(ti4.m_d); //vertex c of Quad
    
    glColor3f(1, 1, 1);
    glPushMatrix();
    glBegin(GL_QUADS);
    glVertex3f(v0.x, v0.y, v0.z);
    glVertex3f(v1.x, v1.y, v1.z);
    glVertex3f(v2.x, v2.y, v2.z);
    glVertex3f(v3.x, v3.y, v3.z);
    glEnd();
    glPopMatrix();
}

bool
Quad::intersect(HitInfo& result, const Ray& r,float tMin, float tMax)
{
    PolygonMesh::TupleI4 ti4 = m_mesh->quadVertexIndex(m_index);
    Vector3 A = m_mesh->vertex(ti4.m_a); //vertex a of Quad
    Vector3 B = m_mesh->vertex(ti4.m_b); //vertex b of Quad
    Vector3 C = m_mesh->vertex(ti4.m_c); //vertex c of Quad
    Vector3 D = m_mesh->vertex(ti4.m_d); //vertex c of Quad
    
    Vector3 center = A + C / 2;
    Vector3 vecX = (B - A).normalize();
    Vector3 vecY = (D - A).normalize();
    float spanX = (B - A).length() / 2;
    float spanY = (D - A).length() / 2;

    Parallelogram parallelogram(center, vecX, vecY, spanX, spanY, m_front, m_back);
    return parallelogram.intersect(result, r, tMin, tMax);
}
