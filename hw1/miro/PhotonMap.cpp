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
        for (unsigned int i = 0; i < m_photonDeposits.size(); i++) {
            if ((m_photonDeposits[i].m_location - x).length2() < r*r) {
                rho += m_photonDeposits[i].m_power / static_cast<float>(M_PI*r*r);
            }
        }
    }
    return rho;
}

bool comparePhotons(const std::pair<float,PhotonDeposit>& p1, const std::pair<float,PhotonDeposit>& p2) {
    if (p1.first < p2.first) return true;
    else if (p1.first > p2.first) return false;
    PhotonDeposit d1 = p1.second;
    PhotonDeposit d2 = p2.second;
    if (d1.m_location[0] < d2.m_location[0]) return true;
    else if (d1.m_location[0] > d2.m_location[0]) return false;
    if (d1.m_location[1] < d2.m_location[1]) return true;
    else if (d1.m_location[1] > d2.m_location[1]) return false;
    if (d1.m_location[2] < d2.m_location[2]) return true;
    else if (d1.m_location[2] > d2.m_location[2]) return false;
    if (d1.m_power[0] < d2.m_power[0]) return true;
    else if (d1.m_power[0] > d2.m_power[0]) return false;
    if (d1.m_power[1] < d2.m_power[1]) return true;
    else if (d1.m_power[1] > d2.m_power[1]) return false;
    if (d1.m_power[2] < d2.m_power[2]) return true;
    else if (d1.m_power[2] > d2.m_power[2]) return false;
    return false;
}

RadiusDensityPhotons PhotonMap::radiusDensityPhotons(const Vector3& x, const int& n) {
    std::vector<std::pair<float,PhotonDeposit>> displacement2(m_photonDeposits.size());
    for (unsigned int i = 0; i < displacement2.size(); i++) {
        displacement2[i] = std::pair<float, PhotonDeposit>((m_photonDeposits[i].m_location - x).length2(), m_photonDeposits[i]);
    }
    std::partial_sort(displacement2.begin(), displacement2.begin() + n, displacement2.end(),comparePhotons);
    RadiusDensityPhotons rdp;
    rdp.m_radius = sqrt(displacement2[n - 1].first);
    for (int i = 0; i < n; i++) {
        rdp.m_photons.push_back(displacement2[i].second);
        rdp.m_density += displacement2[i].second.m_power;
    }
    rdp.m_density /= M_PI*displacement2[n - 1].first*n;
    return rdp;
}