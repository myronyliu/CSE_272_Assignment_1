#ifndef CSE168_LIGHT_H_INCLUDED
#define CSE168_LIGHT_H_INCLUDED

#include <vector>
#include "Vector3.h"
#include "Object.h"

class Light: public virtual Object
{
public:
    void setWattage(float f) { m_wattage = f; }
    float wattage() const { return m_wattage; }
    void setColor(Vector3 v) { m_color = v; }
    Vector3 color() const { return m_color; }
    void preCalc() {} // use this if you need to

    virtual Vector3 randPt() const { return Vector3(0, 0, 0); }

protected:
    Vector3 m_color;
    float m_wattage;
};

typedef std::vector<Light*> Lights;

#endif // CSE168_LIGHT_H_INCLUDED
