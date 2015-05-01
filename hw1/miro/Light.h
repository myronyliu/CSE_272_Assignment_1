#ifndef CSE168_LIGHT_H_INCLUDED
#define CSE168_LIGHT_H_INCLUDED

#include <vector>
#include "Vector3.h"
#include "Object.h"
#include "Ray.h"

struct raypdf {
    Ray r;
    double p;
    raypdf(Ray R, double P) { r = R; p = P; }
};

class Light: public virtual Object
{
public:
    void setWattage(float f) { m_wattage = f; }
    float wattage() const { return m_wattage; }
    void setColor(Vector3 v);
    Vector3 color() const { return m_color; }
    void preCalc() {} // use this if you need to

    virtual raypdf randRay() const { return raypdf(Ray(), 1); }
    virtual vec3pdf randPt() const { return vec3pdf(Vector3(0, 0, 0), 1); }

protected:
    Vector3 m_color;
    float m_wattage;
};

typedef std::vector<Light*> Lights;

#endif // CSE168_LIGHT_H_INCLUDED
