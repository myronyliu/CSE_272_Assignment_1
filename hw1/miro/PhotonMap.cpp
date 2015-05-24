#define _USE_MATH_DEFINES
#include "PhotonMap.h"

Vector3 PhotonMap::powerDensity(const Vector3& x, const float& r) {
    Vector3 rho(0, 0, 0);
    for (int i = 0; i < m_photonDeposits.size(); i++) {
        if ((m_photonDeposits[i].m_location - x).length2() < r*r) {
            rho += m_photonDeposits[i].m_power/(M_PI*r*r);
        }
    }
    return rho;
}