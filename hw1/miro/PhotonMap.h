#ifndef CSE168_PHOTONMAP_H_INCLUDED
#define CSE168_PHOTONMAP_H_INCLUDED

#include <vector>
#include "Ray.h"
#include "Light.h"

struct PhotonDeposit {
    PhotonDeposit() : m_power(0), m_location(Vector3(0, 0, 0)) {}
    PhotonDeposit(const Vector3& power, const Vector3& location) : m_power(power), m_location(location) {}
    Vector3 m_power;
    Vector3 m_location;
};

typedef std::vector<PhotonDeposit> PhotonMap;

#endif // CSE168_PHOTONMAP_H_INCLUDED