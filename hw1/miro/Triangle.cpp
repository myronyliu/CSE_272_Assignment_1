#include "Triangle.h"
#include "Ray.h"


Triangle::Triangle(PolygonMesh * m, unsigned int i, const bool& front, const bool& back) : m_mesh(m), m_index(i), m_front(front), m_back(back) {}

Triangle::~Triangle() {}

vec3pdf Triangle::randPt() const {
    float r1 = (float)rand() / RAND_MAX;
    float r2 = (float)rand() / RAND_MAX;
    Vector3 A = corner0();
    Vector3 B = corner1();
    Vector3 C = corner2();
    Vector3 pt = (1 - sqrt(r1))*A + (sqrt(r1)*(1 - r2))*B + r2*sqrt(r1)*C;
    return vec3pdf(pt, 1.0 / Triangle::area());
}

Vector3 Triangle::barycentric(const Vector3& ptInput) const {
    // get 3d coordinates of the vertices along with "cross-product normal"
    PolygonMesh::TupleI3 ti3 = m_mesh->triangleVertexIndex(m_index);
    Vector3 v1 = m_mesh->vertex(ti3.m_x);
    Vector3 v2 = m_mesh->vertex(ti3.m_y);
    Vector3 v3 = m_mesh->vertex(ti3.m_z);
    Vector3 crossNormal = cross(v2 - v1, v3 - v1).normalize();
    // project ptInput onto the plane of the triangle
    Vector3 pt = ptInput - dot(ptInput - v1, crossNormal)*crossNormal;
    // project everything onto one of the axis-aligned planes by dropping a coordinate
    int projectionAxis = std::max({ crossNormal[0], crossNormal[1], crossNormal[2] }); // the axis to be dropped
    int xAxis = (projectionAxis + 1) % 3;
    int yAxis = (projectionAxis + 2) % 3;
    float x = pt[xAxis];
    float x1 = v1[xAxis];
    float x2 = v2[xAxis];
    float x3 = v3[xAxis];
    float y = pt[yAxis];
    float y1 = v1[yAxis];
    float y2 = v2[yAxis];
    float y3 = v3[yAxis];
    // compute the barycentric coordinates (b1,b2,b3)
    float det = (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3);
    float b1 = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)) / det;
    float b2 = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)) / det;
    float b3 = 1 - b1 - b2;
    return Vector3(b1, b2, b3);
}

Vector3 Triangle::normal(const Vector3& pt) const {
    Vector3 n0 = m_mesh->normal(m_mesh->triangleNormalIndex(m_index).m_x);
    Vector3 n1 = m_mesh->normal(m_mesh->triangleNormalIndex(m_index).m_y);
    Vector3 n2 = m_mesh->normal(m_mesh->triangleNormalIndex(m_index).m_z);
    Vector3 b = barycentric(pt);
    return b[0] * n0 + b[1] * n1 + b[2] * n2;
}

void
Triangle::renderGL()
{
    PolygonMesh::TupleI3 ti3 = m_mesh->triangleVertexIndex(m_index);
    Vector3 v0 = m_mesh->vertex(ti3.m_x); //vertex a of triangle
    Vector3 v1 = m_mesh->vertex(ti3.m_y); //vertex b of triangle
    Vector3 v2 = m_mesh->vertex(ti3.m_z); //vertex c of triangle
    
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
    PolygonMesh::TupleI3 ti3 = m_mesh->triangleVertexIndex(m_index);
    Vector3 A = m_mesh->vertex(ti3.m_x); //vertex a of triangle
    Vector3 B = m_mesh->vertex(ti3.m_y); //vertex b of triangle
    Vector3 C = m_mesh->vertex(ti3.m_z); //vertex c of triangle
    Vector3 AB = B - A;
    Vector3 AC = C - A;
    Vector3 ABxAC = cross(AB, AC); // un-normalized surface normal

    float det = -dot(r.d, ABxAC);
    if (det == 0) return false;
    if (m_back == false && det < 0) return false;
    if (m_front == false && det > 0) return false;

    Vector3 b = r.o - A;

    float alpha = -dot(r.d, cross(b, AC)) / det;
    float beta = -dot(r.d, cross(AB, b)) / det;

    if (alpha < 0 || beta < 0 || alpha + beta>1) return false;

    float nSteps = dot(b, ABxAC) / det;
    if (nSteps<tMin || nSteps>tMax) return false;

    result.t = nSteps;
    if (det > 0) result.N = ABxAC.normalized();
    else result.N = -ABxAC.normalized();
    result.P = r.o + result.t*r.d;
    result.object = this;
    return true;
}
