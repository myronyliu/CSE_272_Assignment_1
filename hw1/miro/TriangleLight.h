#ifndef CSE168_TRIANGLELIGHT_H_INCLUDED
#define CSE168_TRIANGLELIGHT_H_INCLUDED

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

    virtual void disableFront() { printf("disableFront() not available for lights\n"); }
    virtual void disableBack() { printf("disableBack() not available for lights\n"); }
    virtual void enableFront() { printf("enableFront() not available for lights\n"); }
    virtual void enableBack() { printf("enableBack() not available for lights\n"); }

    //virtual Vector3 normal(const Vector3& pt) const;
    //virtual float area() const { return cross(corner1() - corner0(), corner2() - corner0()).length() / 2; }
    virtual RayPDF randRay() const;
    virtual float rayPDF(const Ray& ray) const;
    void preCalc() {} // use this if you need to

protected:

};

typedef std::vector<TriangleLight*> TriangelLights;

#endif // CSE168_TRIANGLELIGHT_H_INCLUDED
