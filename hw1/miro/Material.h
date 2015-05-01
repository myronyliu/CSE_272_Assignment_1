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

    float getEmittance() const { return emittance; }
    Vector3 getEmitted() const { return emitted; }
    void setEmittance(const float& e) { emittance = e; }
    void setEmitted(const Vector3& c) { emitted = c; }

protected:
    Vector3 emitted = Vector3(0,0,0); // emitted RGB color
    float emittance = 0; // probability of reflecting light
};

#endif // CSE168_MATERIAL_H_INCLUDED