#ifndef SPHERE_H_INCLUDED
#define SPHERE_H_INCLUDED

#include "Vector3.h"
#include "Object.h"

class Sphere : public Object
{
public:
    Sphere();
    Sphere(const Vector3& center, const float& radius) : m_center(center), m_radius(radius) {}
    virtual ~Sphere();

    void setCenter(const Vector3& v)    {m_center = v;}
    void setRadius(const float f)       {m_radius = f;}

    const Vector3& center() const       {return m_center;}
    float radius() const                {return m_radius;}

    virtual void renderGL();
    virtual bool intersect(HitInfo& result, const Ray& ray,
                           float tMin = 0.0f, float tMax = MIRO_TMAX);

    virtual Vector3 normal(const Vector3& v) const { return (v - m_center).normalize(); }

    virtual Vector3 shade(const Ray& ray, const HitInfo& hit, const Scene& scene, const Vector3& point = Vector3(0, 0, 0)) const {
        bool isFront = dot(hit.N, normal(point)) > 0;
        return m_material->shade(ray, hit, scene, isFront);
    }
    virtual Vector3 BRDF(const Vector3& in, const Vector3& n, const Vector3& out, Vector3& point = Vector3(0, 0, 0)) const {
        bool isFront = dot(n, normal(point)) > 0;
        return m_material->BRDF(in, n, out, isFront);
    }
    virtual vec3pdf randReflect(const Vector3& in, const Vector3& n, const Vector3& point = Vector3(0, 0, 0)) const  {
        bool isFront = dot(n, normal(point)) > 0;
        return m_material->randReflect(in, n, isFront);
    }


protected:
    Vector3 m_center;
    float m_radius;
};

#endif // SPHERE_H_INCLUDED
