#ifndef TRIANGLELIGHT_H_INCLUDED
#define TRIANGLELIGHT_H_INCLUDED

#include <vector>
#include "Vector3.h"
#include "Triangle.h"
#include "AreaLight.h"
#include "Lambert.h" // this is here because we default the light to "diffuse"

class TriangleLight : public Triangle, public AreaLight
{
public:
    TriangleLight();
    TriangleLight(PolygonMesh * m = 0, unsigned int i = 0);

    virtual RayPDF randRay() const;
    virtual float rayPDF(const Ray& ray) const;
    void preCalc() {} // use this if you need to

protected:

};

typedef std::vector<TriangleLight*> TriangleLights;

#endif // TRIANGLELIGHT_H_INCLUDED
