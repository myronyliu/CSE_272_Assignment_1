#ifndef CSE168_MATERIAL_H_INCLUDED
#define CSE168_MATERIAL_H_INCLUDED

#include <cmath>
#include "Miro.h"
#include "Vector3.h"

struct vec3pdf {
    Vector3 v; // random vector
    double p; // PDF for that random vector
    vec3pdf(Vector3 V, double P) { v = V; p = P; }
};

class Material
{
public:
    Material();
    virtual ~Material();

    virtual void preCalc() {}
    
    virtual Vector3 shade(const Ray& ray, const HitInfo& hit,
                          const Scene& scene) const;

    virtual vec3pdf randEmit(const Vector3& n) const;
    virtual vec3pdf randReflect(const Ray& ray, const HitInfo& hit) const;
    virtual float BRDF(const Ray& in, const HitInfo& hit, const Ray& out) const { return 0; }

    float emittance() const { return m_emittance; }
    Vector3 powerPerPatch() const { return m_powerPerPatch; }
    void setEmittance(const float& e) { m_emittance = e; }
    void setPowerPerPatch(const Vector3& c) { m_powerPerPatch = c; }
    virtual Vector3 powerPerPatchPerSolidAngle(const Vector3&normal, const Vector3& direction) const;

protected:
    Vector3 m_powerPerPatch = Vector3(0, 0, 0); // emitted power per dx.dy patch area
    float m_emittance = 0; // probability of reflecting light
};

#endif // CSE168_MATERIAL_H_INCLUDED
