#ifndef CSE168_MATERIAL_H_INCLUDED
#define CSE168_MATERIAL_H_INCLUDED

#include <cmath>
#include "Miro.h"
#include "Vector3.h"

struct vec3pdf {
    Vector3 v; // random vector
    double p; // PDF for that random vector
    vec3pdf() { v = Vector3(0, 0, 0), p = 1; }
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
    virtual vec3pdf randReflect(const Vector3& ray, const Vector3& hit) const;
    virtual float BRDF(const Vector3& in, const Vector3& normal, const Vector3& out) const { return 0; }

    virtual Vector3 reflectance() const { return Vector3(0,0,0); }
    float emittance() const { return m_emittance; }
    Vector3 powerPerArea() const { return m_powerPerArea; }
    void setEmittance(const float& e) { m_emittance = e; }
    void setPowerPerArea(const Vector3& c) { m_powerPerArea = c; }
    // I think the following function is just radiance (Maybe we'll rename it earlier to make things less confusing);
    virtual Vector3 radiance(const Vector3&normal, const Vector3& direction) const;
    virtual Vector3 sum_L_cosTheta_dOmega() const;

protected:
    Vector3 m_powerPerArea = Vector3(0, 0, 0); // emitted power per unit dx.dy patch area
    float m_emittance = 0; // probability of emitting light
};

#endif // CSE168_MATERIAL_H_INCLUDED
