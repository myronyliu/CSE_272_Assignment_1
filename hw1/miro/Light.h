#ifndef CSE168_LIGHT_H_INCLUDED
#define CSE168_LIGHT_H_INCLUDED

#include <vector>
#include "Vector3.h"
#include "Object.h"
#include "Ray.h"

struct RayPDF {
    Ray r;
    float p;
    RayPDF() { r = Ray(Vector3(0, 0, 0), Vector3(0, 0, 0)); p = 0; }
    RayPDF(Ray R, float P) { r = R; p = P; }
};

class Light: public virtual Object
{
public:
    void setWattage(float f) { m_wattage = f; }
    float wattage() const { return m_wattage; }
    void setColor(Vector3 v);
    Vector3 color() const { return m_color; }
    void preCalc() {} // use this if you need to

    //virtual bool intersect(HitInfo& result, const Ray& ray, float tMin = 0.0f, float tMax = MIRO_TMAX) { return false; }
    virtual float area() const { return 1; }
    virtual RayPDF randRay() const { return RayPDF(); }
    virtual vec3pdf randPt() const { return vec3pdf(Vector3(0, 0, 0), 1); }
    

protected:
    Vector3 m_color;
    float m_wattage;
};

typedef std::vector<Light*> Lights;

#endif // CSE168_LIGHT_H_INCLUDED
