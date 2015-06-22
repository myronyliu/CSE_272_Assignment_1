#ifndef QUADLIGHT_H_INCLUDED
#define QUADLIGHT_H_INCLUDED

#include <vector>
#include "Vector3.h"
#include "Quad.h"
#include "AreaLight.h"
#include "Lambert.h" // this is here because we default the light to "diffuse"

class QuadLight : public Quad, public AreaLight
{
public:
    QuadLight();
    QuadLight(PolygonMesh * m = 0, unsigned int i = 0);

    virtual RayPDF randRay() const;
    virtual float rayPDF(const Ray& ray) const;
    void preCalc() {} // use this if you need to
    void renderGL();

protected:

};

typedef std::vector<QuadLight*> QuadLights;

#endif // QUADLIGHT_H_INCLUDED
