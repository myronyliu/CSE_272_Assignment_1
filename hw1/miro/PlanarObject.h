#ifndef PLANAROBJECT_H_INCLUDED
#define PlANAROBJECT_H_INCLUDED

#include <vector>
#include "Object.h"

class PlanarObject: public virtual Object
{
public:
    PlanarObject() : Object() {}

    void disableBack() { m_back = false; }
    void enableBack() { m_back = true; }
    void disableFront() { m_front = false; }
    void enableFront() { m_front = true; }
    void enable() { m_front = true; m_back = true; }
    void disable() { m_front = false; m_back = false; }
protected:
    bool m_back;
    bool m_front;
};

#endif // PLANAROBJECT_H_INCLUDED
