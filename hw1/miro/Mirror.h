#ifndef CSE168_MIRROR_H_INCLUDED
#define CSE168_MIRROR_H_INCLUDED

#include "Material.h"

class Mirror : public Material
{
public:
    Mirror(const Vector3 & ks = Vector3(1),
        const Vector3 & ka = Vector3(0));
    virtual ~Mirror();

    const Vector3 & ks() const { return m_ks; }
    const Vector3 & ka() const { return m_ka; }

    void setKs(const Vector3 & ks) { m_ks = ks; }
    void setKa(const Vector3 & ka) { m_ka = ka; }

    virtual void preCalc() {}

    virtual Vector3 shade(const Ray& ray, const HitInfo& hit,
        const Scene& scene) const;

    // Returns a random direction for an emitted photon given a surface normal
    virtual vec3pdf randEmit(const Vector3& n) const;
    // Generates a random ray in the upper hemisphere according the BRDF*cos
    virtual vec3pdf randReflect(const Vector3& in, const Vector3& normal) const;
    virtual float BRDF(const Vector3& in, const Vector3& normal, const Vector3& out) const;
    virtual Vector3 radiance(const Vector3& normal, const Vector3& direction) const;
    virtual Vector3 sum_L_cosTheta_dOmega() const;

protected:
    Vector3 m_ks;
    Vector3 m_ka;
};

#endif // CSE168_MIRROR_H_INCLUDED
