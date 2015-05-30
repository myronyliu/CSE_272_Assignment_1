#ifndef CSE168_LIGHT_H_INCLUDED
#define CSE168_LIGHT_H_INCLUDED

#include <vector>
#include "Vector3.h"
#include "Object.h"
#include "Ray.h"

struct RayPDF {
    Ray m_ray;
    float m_oProb;
    float m_dProb;
    RayPDF() : m_ray(Ray(Vector3(0, 0, 0), Vector3(0, 0, 0))), m_oProb(0), m_dProb(0) {}
    RayPDF(Ray ray, float oProb, float dProb) : m_ray(ray), m_oProb(oProb), m_dProb(dProb) {}
};

class Light: public virtual Object
{
public:
    void setWattage(float f);
    float wattage() const { return m_wattage; }
    void setColor(Vector3 v);
    Vector3 color() const { return m_color; }
    void preCalc() {} // use this if you need to

    //virtual bool intersect(HitInfo& result, const Ray& ray, float tMin = 0.0f, float tMax = MIRO_TMAX) { return false; }
    virtual float area() const { return 1; }
    virtual RayPDF randRay() const { return RayPDF(); }
    virtual float rayPDF(const Ray& ray) const { return 0; }
    virtual vec3pdf randPt() const { return vec3pdf(Vector3(0, 0, 0), 1); }

protected:
    Vector3 m_color;
    float m_wattage;
};

typedef std::vector<Light*> Lights;

#endif // CSE168_LIGHT_H_INCLUDED
