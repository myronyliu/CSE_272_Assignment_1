#ifndef PARTICIPATINGMEDIUM_H_INCLUDED
#define PARTICIPATINGMEDIUM_H_INCLUDED

#include <vector>
#include "Material.h"
#include "Vector3.h"

class ParticipatingMedium : public Material
{
public:
    ParticipatingMedium() {};
    ~ParticipatingMedium() {};

    

    void preCalc() {} // use this if you need to

protected:
    Vector3 m_s; // scattering coefficient
    Vector3 m_t; // extinction coefficient
    
};


#endif // PARTICIPATINGMEDIUM_H_INCLUDED
