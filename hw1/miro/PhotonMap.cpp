#define _USE_MATH_DEFINES
#include "PhotonMap.h"

PhotonMap* PhotonMap::getDeepestNode(const Vector3& x) {
    if (isLeafNode()) return this;
    else {
        int splitAxis = m_depth % 3;
        float splitPoint = m_photon->m_location[splitAxis];
        if (x[splitAxis] > splitPoint) {
            if (m_child0->m_XYZ[splitAxis] > m_child1->m_XYZ[splitAxis]) return m_child0->getDeepestNode(x);
            else return m_child1->getDeepestNode(x);
        }
        else {
            if (m_child0->m_xyz[splitAxis] <= m_child1->m_xyz[splitAxis]) return m_child0->getDeepestNode(x);
            else return m_child1->getDeepestNode(x);
        }
    }
}

void PhotonMap::addPhoton(PhotonDeposit newPhotonReference) {
    PhotonDeposit* newPhoton = new PhotonDeposit(newPhotonReference);
    PhotonMap* leafNode = getDeepestNode(newPhoton->m_location);
    if (leafNode->m_photon == NULL) leafNode->m_photon = newPhoton;
    else {
        int splitAxis = leafNode->m_depth % 3;
        Vector3 xyz = leafNode->m_xyz;
        Vector3 XYZ = leafNode->m_XYZ;
        float splitPoint = leafNode->m_photon->m_location[splitAxis];
        XYZ[splitAxis] = splitPoint;
        xyz[splitAxis] = splitPoint;
        if (newPhoton->m_location[splitAxis] > splitPoint) {
            leafNode->m_child0 = new PhotonMap(xyz, leafNode->m_XYZ, newPhoton, leafNode); // LEFT balanced tree (m_child0 gets filled first always)
            leafNode->m_child1 = new PhotonMap(leafNode->m_xyz, XYZ, NULL, leafNode);
        }
        else {
            leafNode->m_child0 = new PhotonMap(leafNode->m_xyz, XYZ, newPhoton, leafNode);
            leafNode->m_child1 = new PhotonMap(xyz, leafNode->m_XYZ, NULL, leafNode);
        }
    }
}
// Results holds points within a bounding box defined by min/max points (bmin, bmax)
void PhotonMap::getPhotons(const Vector3& bmin, const Vector3& bmax, std::vector<PhotonDeposit>& photons) {
    // If we're at a leaf node, just see if the current data point is inside the query bounding box
    if (m_photon != NULL) {
        const Vector3& p = m_photon->m_location;
        if (p.x > bmax.x || p.y > bmax.y || p.z > bmax.z) return;
        if (p.x < bmin.x || p.y < bmin.y || p.z < bmin.z) return;
        photons.push_back(*m_photon);
        if (isLeafNode()) return;
        // We're at an interior node of the tree. Check to see if the query bounding box lies outside the octants of this node.
        bool intersected = true;
        if (m_child0->m_photon == NULL ||
            m_child0->m_XYZ.x < bmin.x || m_child0->m_XYZ.y < bmin.y || m_child0->m_XYZ.z < bmin.z ||
            m_child0->m_xyz.x > bmax.x || m_child0->m_xyz.y > bmax.y || m_child0->m_xyz.z > bmax.z) intersected = false;
        if (intersected == true) m_child0->getPhotons(bmin, bmax, photons);
        intersected = true;
        if (m_child1->m_photon == NULL ||
            m_child1->m_XYZ.x < bmin.x || m_child1->m_XYZ.y < bmin.y || m_child1->m_XYZ.z < bmin.z ||
            m_child1->m_xyz.x > bmax.x || m_child1->m_xyz.y > bmax.y || m_child1->m_xyz.z > bmax.z) intersected = false;
        if (intersected == true) m_child1->getPhotons(bmin, bmax, photons);
    }
}
std::vector<PhotonDeposit> PhotonMap::getNearestPhotons(const Vector3& x, const int& n) {
    std::vector<PhotonDeposit> photons(n);
    return photons;
    const PhotonMap* node = getDeepestNode(x);
    const PhotonMap* lastNode = node;
    int photonCount = 0;
    while (true) {
        if (node->m_photon != NULL) {
            photons[photonCount] = *node->m_photon;
            photonCount++;
            if (photonCount == n) return photons;
        }
        node = node->m_parent;
    }
}

void PhotonMap::buildTree(SequentialPhotonMap spm) {
    m_xyz = Vector3(spm.xMin(), spm.yMin(), spm.zMin());
    m_XYZ = Vector3(spm.xMax(), spm.yMax(), spm.zMax());
    for (int i = 0; i < spm.nPhotons(); i++) addPhoton(spm[i]);
}
bool compareX(const PhotonDeposit& lhs, const PhotonDeposit& rhs) {
    if (lhs.m_location[0] < rhs.m_location[0]) return true;
    else if (lhs.m_location[0] > rhs.m_location[0]) return false;
    else if (lhs.m_location[1] < rhs.m_location[1]) return true;
    else if (lhs.m_location[1] > rhs.m_location[1]) return false;
    else if (lhs.m_location[2] < rhs.m_location[2]) return true;
    else if (lhs.m_location[2] > rhs.m_location[2]) return false;
    else if (lhs.m_power[0] < rhs.m_power[0]) return true;
    else if (lhs.m_power[0] > rhs.m_power[0]) return false;
    else if (lhs.m_power[1] < rhs.m_power[1]) return true;
    else if (lhs.m_power[1] > rhs.m_power[1]) return false;
    else if (lhs.m_power[2] < rhs.m_power[2]) return true;
    else if (lhs.m_power[2] > rhs.m_power[2]) return false;
    else return false;
}
bool compareY(const PhotonDeposit& lhs, const PhotonDeposit& rhs) {
    PhotonDeposit L = lhs;
    PhotonDeposit R = rhs;
    L.m_location = Vector3(L.m_location[1], L.m_location[2], L.m_location[0]);
    R.m_location = Vector3(R.m_location[1], R.m_location[2], R.m_location[0]);
    return compareX(L, R);
}
bool compareZ(const PhotonDeposit& lhs, const PhotonDeposit& rhs) {
    PhotonDeposit L = lhs;
    PhotonDeposit R = rhs;
    L.m_location = Vector3(L.m_location[2], L.m_location[0], L.m_location[1]);
    R.m_location = Vector3(R.m_location[2], R.m_location[0], R.m_location[1]);
    return compareX(L, R);
}
void PhotonMap::buildBalancedTree(SequentialPhotonMap spm) {
    m_xyz = Vector3(spm.xMin(), spm.yMin(), spm.zMin());
    m_XYZ = Vector3(spm.xMax(), spm.yMax(), spm.zMax());
    buildBalancedTree(spm.getPhotons(), 0);
}
void PhotonMap::buildBalancedTree(std::vector<PhotonDeposit>photons, int depth) {
    if (photons.size() == 0) return;
    int medianIndex = photons.size() / 2;
    if (depth % 3 == 0) std::nth_element(photons.begin(), photons.begin() + medianIndex, photons.end(), compareX);
    else if (depth % 3 == 1) std::nth_element(photons.begin(), photons.begin() + medianIndex, photons.end(), compareY);
    else std::nth_element(photons.begin(), photons.begin() + medianIndex, photons.end(), compareZ);
    addPhoton(photons[medianIndex]);
    if (medianIndex > 0) {
        std::vector<PhotonDeposit>photonsL(photons.begin(), photons.begin() + medianIndex - 1);
        buildBalancedTree(photonsL, depth + 1);
    }
    std::vector<PhotonDeposit>photonsR(photons.begin() + medianIndex + 1, photons.end());
    buildBalancedTree(photonsR, depth + 1);
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
    Vector3 dimensions = m_XYZ - m_xyz;
    float V = dimensions[0] * dimensions[1] * dimensions[2];
    float r = pow(V, 1.0 / 3.0) / 64;
    std::vector<PhotonDeposit> photons(0);
    bool keepTrying = true;
    while (keepTrying) {
        if (r > dimensions.length()/2) keepTrying = false; // bail on the next iteration
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
                displacement2[i] = std::pair<float, PhotonDeposit>((photons[i].m_location - x).length2(), photons[i]);
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




////////////////////////////////////////////////////////////////////////////////////////////////////////////





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
