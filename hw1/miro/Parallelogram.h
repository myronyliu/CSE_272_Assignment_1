#ifndef CSE168_PARALLELOGRAM_H_INCLUDED
#define CSE168_PARALLELOGRAM_H_INCLUDED

#include "Vector3.h"
#include "Object.h"

class Parallelogram : public virtual Object
{
public:
    Parallelogram();
    Parallelogram(
        const Vector3& center,
        const Vector3& vecX, const Vector3& vecY,
        const float& spanX, const float& spanY);
    Parallelogram(
        const Vector3& center,
        const Vector3& vecX, const Vector3& vecY,
        const float& spanX, const float& spanY,
        const bool& front, const bool& back);
    virtual ~Parallelogram();

    void setCenter(const float x, const float y, const float z) { m_center = Vector3(x, y, z); }
    void setCenter(const Vector3& v) { m_center = v; }
    void setSpanX(const float f) { m_spanX = f; }
    void setSpanY(const float f) { m_spanY = f; }
    void setVecX(const float x, const float y, const float z) { m_vecX = Vector3(x, y, z).normalize(); }
    void setVecY(const float x, const float y, const float z) { m_vecY = Vector3(x, y, z).normalize(); }
    void setVecX(const Vector3 v) { m_vecX = v.normalized(); }
    void setVecY(const Vector3 v) { m_vecY = v.normalized(); }
    void flip() {
        Vector3 v = m_vecX;
        m_vecX = m_vecY;
        m_vecY = v;
    }
    virtual void disableBack() { m_back = false; }
    virtual void enableBack() { m_back = true; }
    virtual void disableFront() { m_front = false; }
    virtual void enableFront() { m_front = true; }
    virtual void enable() { m_front = true; m_back = true; }
    virtual void disable() { m_front = false; m_back = false; }

    const Vector3 center() const { return m_center; }
    const float spanX() const { return m_spanX; }
    const float spanY() const { return m_spanY; }
    const Vector3 vecX() const { return m_vecX; }
    const Vector3 vecY() const { return m_vecY; }
    virtual Vector3 normal() const { return cross(m_vecX, m_vecY).normalize(); }
    virtual Vector3 normal(const Vector3& v) const { return normal(); }
    //virtual float area() const { return 4.0*m_spanX*m_spanY*cross(m_vecX, m_vecY).length(); }

    virtual void renderGL();
    virtual bool intersect(HitInfo& result, const Ray& ray,
                           float tMin = 0.0f, float tMax = MIRO_TMAX);

protected:
    Vector3 m_center;
    float m_spanX;
    float m_spanY;
    Vector3 m_vecX; // note these not be perpendicular unless one desires a rectangle
    Vector3 m_vecY;
    bool m_back;
    bool m_front;
};

#endif // CSE168_PARALLELOGRAM_H_INCLUDED
