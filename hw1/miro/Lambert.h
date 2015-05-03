#ifndef CSE168_LAMBERT_H_INCLUDED
#define CSE168_LAMBERT_H_INCLUDED

#include "Material.h"

class Lambert : public Material
{
public:
    Lambert(const Vector3 & kd = Vector3(1),
            const Vector3 & ka = Vector3(0));
    virtual ~Lambert();

    const Vector3 & kd() const {return m_kd;}
    const Vector3 & ka() const {return m_ka;}

    void setKd(const Vector3 & kd) {m_kd = kd;}
    void setKa(const Vector3 & ka) {m_ka = ka;}

    virtual void preCalc() {}
    
    virtual Vector3 shade(const Ray& ray, const HitInfo& hit,
                          const Scene& scene) const;

    // Returns a random direction for an emitted photon given a surface normal
    virtual vec3pdf randEmit(const Vector3& n) const;
    // Generates a random ray in the upper hemisphere according the BRDF
    virtual vec3pdf randReflect(const Vector3& ray, const Vector3& hit) const;
    // BRDF
    virtual float BRDF(const Vector3& in, const Vector3& normal, const Vector3& out) const;
    virtual Vector3 radiance(const Vector3& normal, const Vector3& direction) const;
    virtual Vector3 sum_L_cosTheta_dOmega() const;

protected:
    Vector3 m_kd;
    Vector3 m_ka;
};

#endif // CSE168_LAMBERT_H_INCLUDED
