#ifndef PARALLELOGRAMLIGHT_H_INCLUDED
#define PARALLELOGRAMLIGHT_H_INCLUDED

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
    
    virtual void renderGL();

    virtual RayPDF randRay() const;
    virtual float rayPDF(const Ray& ray) const;
    void preCalc() {} // use this if you need to

protected:

};

typedef std::vector<ParallelogramLight*> ParallelogramLights;

#endif // PARALLELOGRAMLIGHT_H_INCLUDED
