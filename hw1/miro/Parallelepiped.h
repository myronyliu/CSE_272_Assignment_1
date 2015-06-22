#ifndef PARALLELEPIPED_H_INCLUDED
#define PARALLELEPIPED_H_INCLUDED

#include "Vector3.h"
#include "Object.h"

class Parallelepiped : public virtual Object
{
public:
    Parallelepiped() {}
    Parallelepiped(const Vector3& center,
        const Vector3& vecX, const Vector3& vecY, const Vector3& vecZ,
        const float& spanX, const float& spanY, const float& spanZ);
    virtual ~Parallelepiped() {}

    void setCenter(const float x, const float y, const float z) { m_center = Vector3(x, y, z); }
    void setCenter(const Vector3& v) { m_center = v; }
    void setSpanX(const float f) { m_spanX = f; }
    void setSpanY(const float f) { m_spanY = f; }
    void setVecX(const float x, const float y, const float z) { m_vecX = Vector3(x, y, z).normalize(); }
    void setVecY(const float x, const float y, const float z) { m_vecY = Vector3(x, y, z).normalize(); }
    void setVecX(const Vector3 v) { m_vecX = v.normalized(); }
    void setVecY(const Vector3 v) { m_vecY = v.normalized(); }

    const Vector3 center() const { return m_center; }
    const float spanX() const { return m_spanX; }
    const float spanY() const { return m_spanY; }
    const float spanZ() const { return m_spanZ; }
    const Vector3 vecX() const { return m_vecX; }
    const Vector3 vecY() const { return m_vecY; }
    const Vector3 vecZ() const { return m_vecZ; }
    virtual Vector3 normal(const Vector3& v) const;
    virtual float area() const { return 8.0f*(
        m_spanX*m_spanY*cross(m_vecX, m_vecY).length()+
        m_spanY*m_spanZ*cross(m_vecY, m_vecZ).length()+
        m_spanZ*m_spanX*cross(m_vecZ, m_vecX).length());
    }

    virtual void renderGL();
    virtual bool intersect(HitInfo& result, const Ray& ray,
                           float tMin = 0.000001f, float tMax = MIRO_TMAX);

    virtual vec3pdf randPt() const;

    virtual Vector3 shade(const Ray& ray, const HitInfo& hit, const Scene& scene, const Vector3& point) const {
        bool isFront = dot(hit.N, normal(point)) > 0;
        return m_material->shade(ray, hit, scene, isFront);
    }
    virtual Vector3 BRDF(const Vector3& in, const Vector3& n, const Vector3& out, Vector3& point) const {
        bool isFront = dot(n, normal(point)) > 0;
        return m_material->BRDF(in, n, out, isFront);
    }
    virtual vec3pdf randReflect(const Vector3& in, const Vector3& n, const Vector3& point) const  {
        bool isFront = dot(n, normal(point)) > 0;
        return m_material->randReflect(in, n, isFront);
    }

    virtual std::pair<Vector3, Vector3> axisAlignedBounds();

protected:
    Vector3 m_center;
    float m_spanX;
    float m_spanY;
    float m_spanZ;
    Vector3 m_vecX; // note these not be perpendicular unless one desires a rectangle
    Vector3 m_vecY;
    Vector3 m_vecZ;

};

#endif // PARALLELEPIPED_H_INCLUDED
