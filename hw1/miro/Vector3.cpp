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

bool coplanar(const Vector3 & uIn, const Vector3& vIn, const Vector3& wIn) {
    Vector3 u = uIn.normalized();
    Vector3 v = vIn.normalized();
    Vector3 w = vIn.normalized();
    if (fabs(dot(u, v)) > 0.999999) return true;
    if (fabs(dot(v, w)) > 0.999999) return true;
    if (fabs(dot(w, u)) > 0.999999) return true;
    if (dot(u, cross(v, w).normalize()) < 0.000001) return true;
    if (dot(v, cross(w, u).normalize()) < 0.000001) return true;
    if (dot(w, cross(u, v).normalize()) < 0.000001) return true;
    return false;
}