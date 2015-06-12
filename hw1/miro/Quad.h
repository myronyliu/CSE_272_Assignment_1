#ifndef CSE168_Quad_H_INCLUDED
#define CSE168_Quad_H_INCLUDED

#include "PlanarObject.h"
#include "PolygonMesh.h"
#include "Vector3.h"

/*
    The Quad class stores a pointer to a mesh and an index into its
    Quad array. The mesh stores all data needed by this Quad.
*/
class Quad : public virtual Object
{
public:
    Quad(PolygonMesh * m = 0, unsigned int i = 0, const bool& front = true, const bool& back = true);
    virtual ~Quad();

    void setIndex(unsigned int i) {m_index = i;}
    void setMesh(PolygonMesh* m) { m_mesh = m; }
    inline Vector3 corner0() const { return m_mesh->vertex(m_mesh->quadVertexIndex(m_index).m_w); }
    inline Vector3 corner1() const { return m_mesh->vertex(m_mesh->quadVertexIndex(m_index).m_x); }
    inline Vector3 corner2() const { return m_mesh->vertex(m_mesh->quadVertexIndex(m_index).m_y); }
    inline Vector3 corner3() const { return m_mesh->vertex(m_mesh->quadVertexIndex(m_index).m_z); }
    virtual Vector3 normal(const Vector3& pt) const;
    virtual float area() const { return cross(corner1() - corner0(), corner2() - corner0()).length(); }
    virtual vec3pdf randPt() const;

    virtual void disableBack() { m_back = false; }
    virtual void enableBack() { m_back = true; }
    virtual void disableFront() { m_front = false; }
    virtual void enableFront() { m_front = true; }
    virtual void enable() { m_front = true; m_back = true; }
    virtual void disable() { m_front = false; m_back = false; }

    virtual void renderGL();
    virtual bool intersect(HitInfo& result, const Ray& ray,
                           float tMin = 0.0f, float tMax = MIRO_TMAX);

    virtual Vector3 shade(const Ray& ray, const HitInfo& hit, const Scene& scene, const Vector3& point = Vector3(0, 0, 0)) const {
        bool isFront = dot(hit.N, normal(point)) > 0;
        return m_material->shade(ray, hit, scene, isFront);
    }
    virtual float BRDF(const Vector3& in, const Vector3& n, const Vector3& out, Vector3& point = Vector3(0, 0, 0)) const {
        bool isFront = dot(n, normal(point)) > 0;
        return m_material->BRDF(in, n, out, isFront);
    }
    virtual vec3pdf randReflect(const Vector3& in, const Vector3& n, const Vector3& point = Vector3(0, 0, 0)) const  {
        bool isFront = dot(n, normal(point)) > 0;
        return m_material->randReflect(in, n, isFront);
    }
    
protected:
    PolygonMesh* m_mesh;
    unsigned int m_index;
    bool m_back;
    bool m_front;
};

#endif // CSE168_Quad_H_INCLUDED
