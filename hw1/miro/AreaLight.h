#ifndef CSE168_AREALIGHT_H_INCLUDED
#define CSE168_AREALIGHT_H_INCLUDED

#include <vector>
#include "Vector3.h"
#include "Light.h"

// AreaLights are one-sided
// They cannot be intersected from behind, nor does any light come off the back-side
// For two-sided area lights, create a new arealight with the same geometry
class AreaLight : public Light
{
public:
    virtual float area() const { return 0; }

    virtual Vector3 radiance() { return Vector3(0, 0, 0); }
    void preCalc() {} // use this if you need to

protected:
};

typedef std::vector<AreaLight*> AreaLights;

#endif // CSE168_AREALIGHT_H_INCLUDED
