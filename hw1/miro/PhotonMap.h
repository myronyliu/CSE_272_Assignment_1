#ifndef CSE168_PHOTONMAP_H_INCLUDED
#define CSE168_PHOTONMAP_H_INCLUDED

#include <vector>
#include "Vector3.h"

struct PhotonDeposit {
    PhotonDeposit() : m_power(0), m_location(Vector3(0, 0, 0)) {}
    PhotonDeposit(const Vector3& power, const Vector3& location) : m_power(power), m_location(location) {}
    Vector3 m_power;
    Vector3 m_location;
};

class PhotonMap {
public:
    PhotonMap() : m_photonDeposits(std::vector<PhotonDeposit>(0)), m_partitionedPhotonDeposits(std::vector<std::vector<std::vector<std::vector<PhotonDeposit>>>>(0)), m_partitionDimensions(Vector3(0,0,0)) {}
    ~PhotonMap() {}
    void push_back(const PhotonDeposit& photonDeposit) { m_photonDeposits.push_back(photonDeposit); }
    PhotonDeposit operator[](int i) const { return m_photonDeposits[i]; }
    PhotonDeposit& operator[](int i) { return m_photonDeposits[i]; }
    Vector3 powerDensity(const Vector3& x, const float& r);
protected:
    Vector3 m_partitionDimensions; // cell dimensions for partitioning photons into subvolumes
    std::vector<PhotonDeposit> m_photonDeposits;
    std::vector<std::vector<std::vector<std::vector<PhotonDeposit>>>> m_partitionedPhotonDeposits;
};

#endif // CSE168_PHOTONMAP_H_INCLUDED