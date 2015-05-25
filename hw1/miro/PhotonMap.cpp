#define _USE_MATH_DEFINES
#include "PhotonMap.h"


int PhotonMap::getOctant(const Vector3& point) const {
    int oct = 0;
    if (point.x >= m_origin.x) oct |= 4;
    if (point.y >= m_origin.y) oct |= 2;
    if (point.z >= m_origin.z) oct |= 1;
    return oct;
}

void PhotonMap::addPhoton(PhotonDeposit* photon) {
    // If this node doesn't have a data point yet assigned and it is a leaf, then we're done!
    if (isLeafNode()) {
        if (m_photon == NULL) {
            m_photon = photon;
            return;
        }
        else {
            // We're at a leaf, but there's already something here. Split this node so that it has 8 child octants
            // and then insert the old data that was here, along with this new data point
            PhotonDeposit *oldPhoton = m_photon; // Save this data point that was here for a later re-insert
            m_photon = NULL;
            // Split the current node and create new empty trees for each child octant.
            for (int i = 0; i < 8; i++) {
                // Compute new bounding box for this child
                Vector3 newOrigin = m_origin;
                newOrigin.x += m_halfDimensions.x * (i & 4 ? .5f : -.5f);
                newOrigin.y += m_halfDimensions.y * (i & 2 ? .5f : -.5f);
                newOrigin.z += m_halfDimensions.z * (i & 1 ? .5f : -.5f);
                m_children[i] = new PhotonMap(newOrigin, m_halfDimensions*.5f);
            }
            // Re-insert the old point, and insert this new point
            // (We wouldn't need to insert from the root, because we already
            // know it's guaranteed to be in this section of the tree)
            m_children[getOctant(oldPhoton->m_location)]->addPhoton(oldPhoton);
            m_children[getOctant(photon->m_location)]->addPhoton(photon);
        }
    }
    else {
        // We are at an interior node. Insert recursively into the appropriate child octant
        int oct = getOctant(photon->m_location);
        m_children[oct]->addPhoton(photon);
    }
}

// Results holds points within a bounding box defined by min/max points (bmin, bmax)
void PhotonMap::getPhotons(const Vector3& bmin, const Vector3& bmax, std::vector<PhotonDeposit*>& results) {
    // If we're at a leaf node, just see if the current data point is inside the query bounding box
    if (isLeafNode()) {
        if (m_photon != NULL) {
            const Vector3& p = m_photon->m_location;
            if (p.x > bmax.x || p.y > bmax.y || p.z > bmax.z) return;
            if (p.x < bmin.x || p.y < bmin.y || p.z < bmin.z) return;
            results.push_back(m_photon);
        }
    }
    else {
        // We're at an interior node of the tree. Check to see if the query bounding box lies outside the octants of this node.
        for (int i = 0; i < 8; ++i) {
            // Compute the min/max corners of this child octant
            Vector3 cmax = m_children[i]->m_origin + m_children[i]->m_halfDimensions;
            Vector3 cmin = m_children[i]->m_origin - m_children[i]->m_halfDimensions;
            // If the query rectangle is outside the child's bounding box, 
            // then continue
            if (cmax.x < bmin.x || cmax.y < bmin.y || cmax.z<bmin.z) continue;
            if (cmin.x>bmax.x || cmin.y > bmax.y || cmin.z > bmax.z) continue;
            // At this point, we've determined that this child is intersecting 
            // the query bounding box
            m_children[i]->getPhotons(bmin, bmax, results);
        }
    }
}

void PhotonMap::buildOctree(SequentialPhotonMap spm) {
    m_origin = Vector3(spm.xMax() + spm.xMin(), spm.yMax() + spm.yMin(), spm.zMax() + spm.zMin()) / 2;
    m_halfDimensions = Vector3(spm.xMax() - spm.xMin(), spm.yMax() - spm.yMin(), spm.zMax() - spm.zMin()) / 2;
    for (int i = 0; i < spm.nPhotons(); i++) {
        addPhoton(new PhotonDeposit(spm[i]));
    }
}


bool comparePhotons(const std::pair<float, PhotonDeposit>& p1, const std::pair<float, PhotonDeposit>& p2) {
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
    float t = 1.5f*(4.0f / M_PI); // the thing in paranthesis is the ratio of areas between square and enclosed circle
    float r = pow(m_halfDimensions[0] * m_halfDimensions[1] * m_halfDimensions[2], 1.0 / 3.0) / 64;
    std::vector<PhotonDeposit*> photons(0);
    bool keepTrying = true;
    while (keepTrying) {
        if (r > m_halfDimensions.length()) keepTrying = false; // bail on the next iteration
        photons.clear();
        getPhotons(x - Vector3(r, r, r), x + Vector3(r, r, r), photons);
        int m = photons.size();
        if ((float)m / n < t) {
            r *= sqrt(t*n / m);
            continue;
        }
        else { // square is sufficiently large to try our luck
            std::vector<std::pair<float, PhotonDeposit>> displacement2(m);
            for (int i = 0; i < m; i++) {
                displacement2[i] = std::pair<float, PhotonDeposit>((photons[i]->m_location - x).length2(), *photons[i]);
            }
            std::partial_sort(displacement2.begin(), displacement2.begin() + n, displacement2.end(), comparePhotons);
            if (displacement2[n - 1].first > r) {
                r *= sqrt(t*n / m); // maybe can make this tighter, but cannot find an efficient way to do so
                continue;
            }
            else {
                RadiusDensityPhotons rdp;
                rdp.m_radius = sqrt(displacement2[n - 1].first);
                for (int i = 0; i < n; i++) {
                    rdp.m_photons.push_back(displacement2[i].second);
                    rdp.m_density += displacement2[i].second.m_power;
                }
                rdp.m_density /= M_PI*displacement2[n - 1].first*n;
                return rdp;
            }
        }
    }
    return RadiusDensityPhotons();
}


void SequentialPhotonMap::addPhoton(const PhotonDeposit& photonDeposit) {
    m_photonDeposits.push_back(photonDeposit);
    if (photonDeposit.m_location[0] < m_xMin) m_xMin = photonDeposit.m_location[0];
    if (photonDeposit.m_location[1] < m_yMin) m_yMin = photonDeposit.m_location[1];
    if (photonDeposit.m_location[2] < m_zMin) m_zMin = photonDeposit.m_location[2];
    if (photonDeposit.m_location[0] > m_xMax) m_xMax = photonDeposit.m_location[0];
    if (photonDeposit.m_location[1] > m_yMax) m_yMax = photonDeposit.m_location[1];
    if (photonDeposit.m_location[2] > m_zMax) m_zMax = photonDeposit.m_location[2];
}

Vector3 SequentialPhotonMap::powerDensity(const Vector3& x, const float& r) {
    Vector3 rho(0, 0, 0);
    for (int i = 0; i < m_photonDeposits.size(); i++) {
        if ((m_photonDeposits[i].m_location - x).length2() < r*r) {
            rho += m_photonDeposits[i].m_power / (M_PI*r*r);
        }
    }
    return rho;
}


RadiusDensityPhotons SequentialPhotonMap::radiusDensityPhotons(const Vector3& x, const int& n) {
    std::vector<std::pair<float, PhotonDeposit>> displacement2(m_photonDeposits.size());
    for (int i = 0; i < displacement2.size(); i++) {
        displacement2[i] = std::pair<float, PhotonDeposit>((m_photonDeposits[i].m_location - x).length2(), m_photonDeposits[i]);
    }
    std::partial_sort(displacement2.begin(), displacement2.begin() + n, displacement2.end(), comparePhotons);
    RadiusDensityPhotons rdp;
    rdp.m_radius = sqrt(displacement2[n - 1].first);
    for (int i = 0; i < n; i++) {
        rdp.m_photons.push_back(displacement2[i].second);
        rdp.m_density += displacement2[i].second.m_power;
    }
    rdp.m_density /= M_PI*displacement2[n - 1].first*n;
    return rdp;
}
