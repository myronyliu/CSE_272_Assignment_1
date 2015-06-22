#ifndef PHONG_H_INCLUDED
#define PHONG_H_INCLUDED

#include "Material.h"

class Phong : public Material
{
public:
    Phong(const Vector3 & kd = Vector3(1),
        const Vector3 & ks = Vector3(0),
        const float n = 1.0f);
    virtual ~Phong();

    const Vector3 & kd() const { return m_kd; }
    const Vector3 & ks() const { return m_ks; }
    const float n() const { return m_n; }

    void setKd(const Vector3 & kd) { m_kd = kd; }
    void setKs(const Vector3 & ks) { m_ks = ks; }
    void setN(const float n) { m_n = n;}

    virtual void preCalc() {}

    virtual Vector3 shade(const Ray& ray, const HitInfo& hit, const Scene& scene, const bool& isFront = true) const;

    // Returns a random direction for an emitted photon given a surface normal
    virtual vec3pdf randEmit(const Vector3& n) const;
    virtual float emitPDF(const Vector3& n, const Vector3& v) const;
    // Generates a random ray in the upper hemisphere according the BRDF*cos
    virtual vec3pdf randReflect(const Vector3& in, const Vector3& normal,const bool& isFront = true) const;
    virtual Vector3 BRDF(const Vector3& in, const Vector3& normal, const Vector3& out, const bool& isFront = true) const;
    virtual Vector3 radiance(const Vector3& normal, const Vector3& direction) const;
    virtual Vector3 sum_L_cosTheta_dOmega() const;

    virtual Vector3 reflectance() const { return m_kd + m_ks; }
    virtual Vector3 transmittance() const { return Vector3(0, 0, 0); }
protected:
    Vector3 m_kd;
    Vector3 m_ks;
    float m_n;
};

#endif // PHONG_H_INCLUDED
