#ifndef CSE168_OBJECT_H_INCLUDED
#define CSE168_OBJECT_H_INCLUDED

#include <vector>
#include "Miro.h"
#include "Ray.h"
#include "Material.h"

class Object
{
public:
    Object() {}
    virtual ~Object() {}

    void setMaterial(Material* m) {m_material = m;}
    Material* material() const { return m_material; }

    virtual void renderGL() {}
    virtual void preCalc() {}
    virtual Vector3 normal(const Vector3& v) const { return Vector3(0, 0, 1); }

    virtual bool intersect(HitInfo& result, const Ray& ray,
                           float tMin = 0.0f, float tMax = MIRO_TMAX) = 0;

    virtual Vector3 shade(const Ray& ray, const HitInfo& hit, const Scene& scene, const Vector3& point = Vector3(0, 0, 0)) const { return Vector3(0, 0, 0); }
    virtual float BRDF(const Vector3& in, const Vector3& normal, const Vector3& out, Vector3& point = Vector3(0, 0, 0)) const { return 0; }
    virtual vec3pdf randReflect(const Vector3& in, const Vector3& normal, const Vector3& point = Vector3(0,0,0)) const { return vec3pdf(); }

    virtual std::pair<Vector3, Vector3> axisAlignedBounds() { return std::pair<Vector3, Vector3>(Vector3(0, 0, 0), Vector3(0, 0, 0)); }
protected:
    Material* m_material;
};

typedef std::vector<Object*> Objects;

#endif // CSE168_OBJECT_H_INCLUDED
