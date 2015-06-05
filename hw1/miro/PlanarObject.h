#ifndef CSE168_PLANAROBJECT_H_INCLUDED
#define CSE168_PlANAROBJECT_H_INCLUDED

#include <vector>
#include "Object.h"

class PlanarObject: public virtual Object
{
public:
    PlanarObject() {}
    virtual ~PlanarObject() {}

    virtual void disableBack() { m_back = false; }
    virtual void enableBack() { m_back = true; }
    virtual void disableFront() { m_front = false; }
    virtual void enableFront() { m_front = true; }
    virtual void enable() { m_front = true; m_back = true; }
    virtual void disable() { m_front = false; m_back = false; }
protected:
    bool m_back;
    bool m_front;
};

#endif // CSE168_PLANAROBJECT_H_INCLUDED
