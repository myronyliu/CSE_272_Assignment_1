#ifndef CSE168_REFRACTIVEINTERFACE_H_INCLUDED
#define CSE168_REFRACTIVEINTERFACE_H_INCLUDED

#include "Material.h"

class RefractiveInterface : public Material
{
public:
    RefractiveInterface(const Vector3& kr = Vector3(1), const Vector3& kt = Vector3(1), const Vector3& ka = Vector3(0), const float& n0 = 1, const float& n1 = 1);
    virtual ~RefractiveInterface();

    const Vector3 & kr() const { return m_kr; }
    const Vector3 & kt() const { return m_kt; }
    const Vector3 & ka() const { return m_ka; }

    void setKr(const Vector3 & kr) { m_kr = kr; }
    void setKt(const Vector3 & kt) { m_kt = kt; }
    void setKa(const Vector3 & ka) { m_ka = ka; }
    void setRefractiveIndexFront(const float& n0) { m_n0 = n0; }
    void setRefractiveIndexBack(const float& n1) { m_n1 = n1; }

    virtual void preCalc() {}

    virtual Vector3 shade(const Ray& ray, const HitInfo& hit, const Scene& scene, const bool& isFront = true) const;

    // Returns a random direction for an emitted photon given a surface normal
    virtual vec3pdf randEmit(const Vector3& n) const { return vec3pdf(Vector3(0, 0, 0), 1); }
    virtual float emitPDF(const Vector3& n, const Vector3& v) const { return 0; }
    // Generates a random ray in the upper hemisphere according the BRDF*cos
    virtual vec3pdf randReflect(const Vector3& in, const Vector3& normal, const bool& isFront) const;
    virtual float BRDF(const Vector3& in, const Vector3& normal, const Vector3& out, const bool& isFront) const;
    virtual Vector3 radiance(const Vector3& normal, const Vector3& direction) const { return Vector3(0, 0, 0); }
    virtual Vector3 sum_L_cosTheta_dOmega() const { return Vector3(0, 0, 0); }

    virtual Vector3 reflectance() const { return m_kr; }
    virtual Vector3 transmittance() const { return m_kt; }
protected:
    float m_n0; // index of refraction of the medium on the FRONT side
    float m_n1;
    Vector3 m_kt;
    Vector3 m_kr;
    Vector3 m_ka;
};

#endif // CSE168_REFRACTIVEINTERFACE_H_INCLUDED
