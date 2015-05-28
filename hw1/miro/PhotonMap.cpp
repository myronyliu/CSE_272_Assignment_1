#define _USE_MATH_DEFINES
#include "PhotonMap.h"

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

PhotonMap* PhotonMap::getLeafNode(const Vector3& x) {
    if (isLeafNode()) return this;
    else {
        int splitAxis = m_depth % 3;
        float splitPoint = m_photon->m_location[splitAxis];
        if (x[splitAxis] > splitPoint) {
            if (m_child0->m_XYZ[splitAxis] > m_child1->m_XYZ[splitAxis]) return m_child0->getLeafNode(x);
            else return m_child1->getLeafNode(x);
        }
        else {
            if (m_child0->m_xyz[splitAxis] <= m_child1->m_xyz[splitAxis]) return m_child0->getLeafNode(x);
            else return m_child1->getLeafNode(x);
        }
    }
}

void PhotonMap::addPhoton(PhotonDeposit newPhotonReference) {
    PhotonDeposit* newPhoton = new PhotonDeposit(newPhotonReference);
    PhotonMap* leafNode = getLeafNode(newPhoton->m_location);
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

void PhotonMap::getNearestPhotons(const Vector3& x, const int& k, PhotonMap* inputNode, const int& rootDepth, std::priority_queue<RsqrPhoton>& photons) {
    PhotonMap* node = inputNode->getLeafNode(x);
    PhotonMap* parent = node->m_parent;
    if (parent == NULL) { // this will only ever happen if there is only a single photon, since (in practice) we start at leaf nodes
        photons.push(RsqrPhoton((node->m_photon->m_location - x).length2(), *node->m_photon));
        return;
    }
    PhotonMap* sibling = parent->m_child0;
    if (node == parent->m_child0) sibling = parent->m_child1;

    PhotonDeposit photon;
    float r2;
    while (true) {
        bool photonAdded = false;
        bool siblingIntersected = false;
        if (node->m_photon != NULL) {
            photon = *node->m_photon;
            photonAdded = true;
            r2 = (node->m_photon->m_location - x).length2();
            if (photons.size() < k) photons.push(RsqrPhoton(r2, photon));
            else if (r2 < photons.top().m_r2) {
                photons.pop();
                photons.push(RsqrPhoton(r2, *node->m_photon));
            }
            else photonAdded = false;
        }
        if (node->m_depth == rootDepth) return;
        if (photonAdded == true) {
            int splitAxis = parent->m_depth % 3;
            float d2 = photon.m_location[splitAxis] - parent->m_photon->m_location[splitAxis];
            d2 *= d2;
            if (r2 >= d2) { // sibling division was intersected, so recurse in the sibling
                siblingIntersected = true;
                getNearestPhotons(x, k, sibling, node->m_depth, photons);
            }
        }
        if (parent->m_parent == NULL) { // i.e. if parent is the root of the "entire tree"
            float R2 = (parent->m_photon->m_location - x).length2();
            PhotonDeposit PHOTON = *parent->m_photon;
            if (photons.size() < k) photons.push(RsqrPhoton(R2, PHOTON));
            else if (R2 < photons.top().m_r2) {
                photons.pop();
                photons.push(RsqrPhoton(R2, PHOTON));
            }
            if (siblingIntersected == false && photons.size() < k) {
                getNearestPhotons(x, k, sibling, node->m_depth, photons);
            }
            return;
        }
        node = parent;
        parent = node->m_parent;
        sibling = parent->m_child0;
        if (node == parent->m_child0) sibling = parent->m_child1;
    }
}

std::vector<PhotonDeposit> PhotonMap::getNearestPhotons(const Vector3& x, const int& n) {
    std::priority_queue<RsqrPhoton> photonQueue;
    getNearestPhotons(x, n, getLeafNode(x), m_depth, photonQueue);
    std::vector<PhotonDeposit> photons(photonQueue.size());
    for (int i = 0; i < photonQueue.size(); i++) {
        photons[photonQueue.size() - i - 1] = photonQueue.top().m_photon;
        photonQueue.pop();
    }
    return photons;
}

void PhotonMap::buildTree(SequentialPhotonMap spm) {
    m_xyz = Vector3(spm.xMin(), spm.yMin(), spm.zMin());
    m_XYZ = Vector3(spm.xMax(), spm.yMax(), spm.zMax());
    for (int i = 0; i < spm.nPhotons(); i++) addPhoton(spm[i]);
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
        std::vector<PhotonDeposit>photonsL(photons.begin(), photons.begin() + medianIndex);
        buildBalancedTree(photonsL, depth + 1);
    }
    if (medianIndex + 1 < photons.size()) {
        std::vector<PhotonDeposit>photonsR(photons.begin() + medianIndex + 1, photons.end());
        buildBalancedTree(photonsR, depth + 1);
    }
}


RadiusDensityPhotons PhotonMap::radiusDensityPhotons(const Vector3& x, const int& k) {
    RadiusDensityPhotons rdp;
    rdp.m_photons = getNearestPhotons(x, k);
    float r2 = (rdp.m_photons.back().m_location - x).length2();
    rdp.m_radius = sqrt(r2);
    for (int i = 0; i < rdp.m_photons.size(); i++) rdp.m_density += rdp.m_photons[i].m_power;
    rdp.m_density /= (M_PI*r2);
    return rdp;
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
    rdp.m_density /= M_PI*displacement2[n - 1].first;
    return rdp;
}
