#define _USE_MATH_DEFINES
#include "PhotonMap.h"

void PhotonMap::push_back(const PhotonDeposit& photonDeposit) {
    m_photonDeposits.push_back(photonDeposit);
    if (photonDeposit.m_location[0] < m_xMin) m_xMin = photonDeposit.m_location[0];
    if (photonDeposit.m_location[1] < m_yMin) m_yMin = photonDeposit.m_location[1];
    if (photonDeposit.m_location[2] < m_zMin) m_zMin = photonDeposit.m_location[2];
    if (photonDeposit.m_location[0] > m_xMax) m_xMax = photonDeposit.m_location[0];
    if (photonDeposit.m_location[1] > m_yMax) m_yMax = photonDeposit.m_location[1];
    if (photonDeposit.m_location[2] > m_zMax) m_zMax = photonDeposit.m_location[2];
}

Vector3 PhotonMap::powerDensity(const Vector3& x, const float& r) {
    Vector3 rho(0, 0, 0);
    if (m_partitionDimensions.length2() <= 0) {
        for (int i = 0; i < m_photonDeposits.size(); i++) {
            if ((m_photonDeposits[i].m_location - x).length2() < r*r) {
                rho += m_photonDeposits[i].m_power / (M_PI*r*r);
            }
        }
    }
    return rho;
}

float PhotonMap::radius(const Vector3& x, const int& n) {
    std::vector<float> displacement2(m_photonDeposits.size(), 0);
    for (int i = 0; i < displacement2.size(); i++) {
        displacement2[i] = (m_photonDeposits[i].m_location - x).length2();
    }
    std::partial_sort(displacement2.begin(), displacement2.begin() + n, displacement2.end());
    return displacement2[n];
}