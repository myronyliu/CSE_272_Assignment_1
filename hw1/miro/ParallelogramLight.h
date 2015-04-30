#ifndef CSE168_PARALLELOGRAMLIGHT_H_INCLUDED
#define CSE168_PARALLELOGRAMLIGHT_H_INCLUDED

#include <vector>
#include "Vector3.h"
#include "Parallelogram.h"
#include "AreaLight.h"
#include "Lambert.h" // this is here because we default the light to "diffuse"

class ParallelogramLight : public Parallelogram, public AreaLight
{
public:
    ParallelogramLight();
    ParallelogramLight(
        const Vector3& center,
        const Vector3& vecX, const Vector3& vecY,
        const float& spanX, const float& spanY);
    //~ParallelogramLight();

    virtual void disableFront() { printf("disableFront() not available for lights\n"); }
    virtual void disableBack() { printf("disableBack() not available for lights\n"); }
    virtual void enableFront() { printf("enableFront() not available for lights\n"); }
    virtual void enableBack() { printf("enableBack() not available for lights\n"); }
    //virtual void renderGL();

    virtual Vector3 randPt() const;
    void preCalc() {} // use this if you need to

protected:

};

typedef std::vector<ParallelogramLight*> ParallelogramLights;

#endif // CSE168_PARALLELOGRAMLIGHT_H_INCLUDED
