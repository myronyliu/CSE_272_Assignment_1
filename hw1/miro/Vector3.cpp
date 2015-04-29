#include "Vector3.h"

Vector3 Vector3::project(const Vector3 & v) {
    Vector3 u = v.normalized();
    return dot(*this, u)*u;
}

Vector3 Vector3::project(const Vector3 & u, const Vector3 & v) {
    return *this - orthogonal(u, v);
}

Vector3 Vector3::orthogonal(const Vector3 & v) { // return the component that is perpendicular to v
    return *this - project(v);
}

Vector3 Vector3::orthogonal(const Vector3 & u, const Vector3 &v) { // return the component that is perpendicular to the plane spanned by u,v
    return project(cross(u,v));
}

Vector3 Vector3::inBasis(const Vector3 & u, const Vector3 & v, const Vector3& w){
    return Vector3(dot(*this, u), dot(*this, v), dot(*this, w));
}