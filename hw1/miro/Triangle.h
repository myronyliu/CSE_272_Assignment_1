#ifndef CSE168_TRIANGLE_H_INCLUDED
#define CSE168_TRIANGLE_H_INCLUDED

#include "PlanarObject.h"
#include "TriangleMesh.h"
#include "Vector3.h"

/*
    The Triangle class stores a pointer to a mesh and an index into its
    triangle array. The mesh stores all data needed by this Triangle.
*/
class Triangle : public virtual Object
{
public:
    Triangle(TriangleMesh * m = 0, unsigned int i = 0);
    virtual ~Triangle();

    void setIndex(unsigned int i) {m_index = i;}
    void setMesh(TriangleMesh* m) { m_mesh = m; }
    inline Vector3 corner0() const { return m_mesh->vertex(m_mesh->vIndex(m_index).m_x); }
    inline Vector3 corner1() const { return m_mesh->vertex(m_mesh->vIndex(m_index).m_y); }
    inline Vector3 corner2() const { return m_mesh->vertex(m_mesh->vIndex(m_index).m_z); }
    virtual Vector3 normal() const { return cross(corner1() - corner0(), corner2() - corner0()); }

    virtual void renderGL();
    virtual bool intersect(HitInfo& result, const Ray& ray,
                           float tMin = 0.0f, float tMax = MIRO_TMAX);

    virtual Vector3 shade(const Ray& ray, const HitInfo& hit, const Scene& scene, const Vector3& point = Vector3(0, 0, 0)) const {
        bool isFront = dot(hit.N, normal()) > 0;
        return m_material->shade(ray, hit, scene, isFront);
    }
    virtual float BRDF(const Vector3& in, const Vector3& n, const Vector3& out, Vector3& point = Vector3(0, 0, 0)) const {
        bool isFront = dot(n, normal()) > 0;
        return m_material->BRDF(in, n, out, isFront);
    }
    virtual vec3pdf randReflect(const Vector3& in, const Vector3& n, const Vector3& point = Vector3(0, 0, 0)) const  {
        bool isFront = dot(n, normal()) > 0;
        return m_material->randReflect(in, n, isFront);
    }
    
protected:
    TriangleMesh* m_mesh;
    unsigned int m_index;
};

#endif // CSE168_TRIANGLE_H_INCLUDED
