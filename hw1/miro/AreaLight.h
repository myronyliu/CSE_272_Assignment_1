#ifndef AREALIGHT_H_INCLUDED
#define AREALIGHT_H_INCLUDED

#include <vector>
#include "Vector3.h"
#include "Light.h"

// AreaLights are one-sided
// They cannot be intersected from behind, nor does any light come off the back-side
// For two-sided area lights, create a new arealight with the same geometry
class AreaLight : public Light
{
public:

    void preCalc() {} // use this if you need to

protected:
};

typedef std::vector<AreaLight*> AreaLights;

#endif // AREALIGHT_H_INCLUDED
