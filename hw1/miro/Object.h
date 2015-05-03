#ifndef CSE168_OBJECT_H_INCLUDED
#define CSE168_OBJECT_H_INCLUDED

#include <vector>
#include "Miro.h"
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
protected:
    Material* m_material;
};

typedef std::vector<Object*> Objects;

#endif // CSE168_OBJECT_H_INCLUDED
