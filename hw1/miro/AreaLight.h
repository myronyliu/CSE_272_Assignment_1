#ifndef CSE168_AREALIGHT_H_INCLUDED
#define CSE168_AREALIGHT_H_INCLUDED

#include <vector>
#include "Vector3.h"
#include "Light.h"

class AreaLight : public Light
{
public:
    //void setColorFront(const Vector3& v) { m_colorFront = v; }
    //void setColorBack(const Vector3& v) { m_colorBack = v; }
    //void setWattageFront(float f) { m_wattageFront = f; }
    //void setWattageBack(float f) { m_wattageBack= f; }

    //float wattageFront() const { return m_wattageFront; }
    //float wattageBack() const { return m_wattageBack; }
    //const Vector3 & colorFront() const { return m_colorFront; }
    //const Vector3 & colorBack() const { return m_colorBack; }
    
    //virtual bool backEnabled() const { return true; }
    //virtual bool frontEnabled() const { return true; }
    virtual Vector3 randPt() const { return Vector3(0, 0, 0); }

    void preCalc() {} // use this if you need to

protected:
    //Vector3 m_colorBack;
    //Vector3 m_colorFront;
    //float m_wattageBack;
    //float m_wattageFront;
};

typedef std::vector<AreaLight*> AreaLights;

#endif // CSE168_AREALIGHT_H_INCLUDED
