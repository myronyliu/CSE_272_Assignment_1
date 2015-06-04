#define _USE_MATH_DEFINES
#include "PhotonMap.h"

bool compareX(const PhotonDeposit& lhs, const PhotonDeposit& rhs) {
    if (lhs.location()[0] < rhs.location()[0]) return true;
    else if (lhs.location()[0] > rhs.location()[0]) return false;
    else if (lhs.location()[1] < rhs.location()[1]) return true;
    else if (lhs.location()[1] > rhs.location()[1]) return false;
    else if (lhs.location()[2] < rhs.location()[2]) return true;
    else if (lhs.location()[2] > rhs.location()[2]) return false;
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
    L.location() = Vector3(L.location()[1], L.location()[2], L.location()[0]);
    R.location() = Vector3(R.location()[1], R.location()[2], R.location()[0]);
    return compareX(L, R);
}
bool compareZ(const PhotonDeposit& lhs, const PhotonDeposit& rhs) {
    PhotonDeposit L = lhs;
    PhotonDeposit R = rhs;
    L.location() = Vector3(L.location()[2], L.location()[0], L.location()[1]);
    R.location() = Vector3(R.location()[2], R.location()[0], R.location()[1]);
    return compareX(L, R);
}
bool comparePhotons(const std::pair<float, PhotonDeposit>& p1, const std::pair<float, PhotonDeposit>& p2) {
    if (p1.first < p2.first) return true;
    else if (p1.first > p2.first) return false;
    PhotonDeposit d1 = p1.second;
    PhotonDeposit d2 = p2.second;
    if (d1.location()[0] < d2.location()[0]) return true;
    else if (d1.location()[0] > d2.location()[0]) return false;
    if (d1.location()[1] < d2.location()[1]) return true;
    else if (d1.location()[1] > d2.location()[1]) return false;
    if (d1.location()[2] < d2.location()[2]) return true;
    else if (d1.location()[2] > d2.location()[2]) return false;
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
        float splitPoint = m_photon->location()[m_axis];
        if (x[m_axis] < splitPoint) {
            if (m_child0->m_xyz[m_axis] < m_child1->m_XYZ[m_axis]) return m_child0->getLeafNode(x);
            else if (m_child1->m_xyz[m_axis] < m_child0->m_XYZ[m_axis]) return m_child1->getLeafNode(x);
            else return m_child0->getLeafNode(x); // otherwise both children are degenerate (story of my parent's life...)
        }
        else if (x[m_axis] > splitPoint) {
            if (m_child0->m_XYZ[m_axis] > m_child1->m_xyz[m_axis]) return m_child0->getLeafNode(x);
            else if (m_child1->m_XYZ[m_axis] > m_child0->m_xyz[m_axis]) return m_child1->getLeafNode(x);
            else return m_child0->getLeafNode(x);
        }
        else return m_child0->getLeafNode(x);
    }
}

void PhotonMap::addPhoton(PhotonDeposit newPhotonReference) {
    PhotonDeposit* newPhoton = new PhotonDeposit(newPhotonReference);
    PhotonMap* leafNode = getLeafNode(newPhoton->location());
    if (leafNode->m_photon == NULL) leafNode->m_photon = newPhoton;
    else {
        Vector3 xyz = leafNode->m_xyz;
        Vector3 XYZ = leafNode->m_XYZ;
        float splitPoint = newPhoton->location()[m_axis];
        XYZ[m_axis] = splitPoint;
        xyz[m_axis] = splitPoint;
        if (newPhoton->location()[m_axis] <= splitPoint) {
            leafNode->m_child0 = new PhotonMap(xyz, leafNode->m_XYZ, newPhoton, leafNode); // LEFT balanced tree (m_child0 gets filled first always)
            leafNode->m_child1 = new PhotonMap(leafNode->m_xyz, XYZ, NULL, leafNode);
        }
        else {
            leafNode->m_child0 = new PhotonMap(leafNode->m_xyz, XYZ, newPhoton, leafNode);
            leafNode->m_child1 = new PhotonMap(xyz, leafNode->m_XYZ, NULL, leafNode);
        }
        PhotonMap* parent = leafNode->m_parent;
        if (parent != NULL) {
            if (parent->m_child0 == leafNode) return; // sibling is child1 so we are done
            if (parent->m_child0->isLeafNode() == false) return; // sibling is child0, but it already has children, so we are done
            parent->m_child1 = parent->m_child0; // otherwise swap to make tree LEFT balanced
            parent->m_child0 = leafNode;
    }
}
}
// Results holds points within a bounding box defined by min/max points (bmin, bmax)
void PhotonMap::getPhotons(const Vector3& bmin, const Vector3& bmax, std::vector<PhotonDeposit>& photons) {
    // If we're at a leaf node, just see if the current data point is inside the query bounding box
    if (m_photon != NULL) {
        const Vector3& p = m_photon->location();
        if (p.x > bmax.x || p.y > bmax.y || p.z > bmax.z || p.x < bmin.x || p.y < bmin.y || p.z < bmin.z) {}
        else photons.push_back(*m_photon);
        if (isLeafNode()) return;
        // We're at an interior node of the tree. Check to see if the query bounding box lies outside the octants of this node.
        if (m_child0->m_photon == NULL ||
            m_child0->m_XYZ.x < bmin.x || m_child0->m_XYZ.y < bmin.y || m_child0->m_XYZ.z < bmin.z ||
            m_child0->m_xyz.x > bmax.x || m_child0->m_xyz.y > bmax.y || m_child0->m_xyz.z > bmax.z) {}
        else m_child0->getPhotons(bmin, bmax, photons);
        if (m_child1->m_photon == NULL ||
            m_child1->m_XYZ.x < bmin.x || m_child1->m_XYZ.y < bmin.y || m_child1->m_XYZ.z < bmin.z ||
            m_child1->m_xyz.x > bmax.x || m_child1->m_xyz.y > bmax.y || m_child1->m_xyz.z > bmax.z) {}
        else m_child1->getPhotons(bmin, bmax, photons);
    }
}


void PhotonMap::getPhotons(const Vector3& x, const float& r, std::vector<PhotonDeposit>& photons) {
    if (m_photon != NULL) {
        const Vector3& p = m_photon->location();
        if ((p - x).length2() <= r*r) photons.push_back(*m_photon);
        if (isLeafNode()) return;
        // We're at an interior node of the tree. Check to see if the query bounding box lies outside the octants of this node.
        Sphere ball(x, r);
        HitInfo hit;
        for (int i = 0; i < 2; i++) {
            PhotonMap* child = m_child0;
            if (i == 1) child = m_child1;
            if (child->m_photon == NULL ||
                child->m_XYZ[0] < x[0] - r || child->m_XYZ[1] < x[1] - r || child->m_XYZ[2] < x[2] - r ||
                child->m_xyz[0] > x[0] + r || child->m_xyz[1] > x[1] + r || child->m_xyz[2] > x[2] + r) continue;
            else {
                Vector3 corner100 = child->m_xyz;
                Vector3 corner010 = child->m_xyz;
                Vector3 corner001 = child->m_xyz;
                corner100[0] = m_XYZ[0];
                corner010[1] = m_XYZ[1];
                corner001[2] = m_XYZ[2];
                Vector3 boxDims = child->m_XYZ - child->m_xyz;
                if (ball.intersect(hit, Ray(m_xyz, Vector3(1, 0, 0)), 0, boxDims[0]) ||
                    ball.intersect(hit, Ray(m_xyz, Vector3(0, 1, 0)), 0, boxDims[1]) ||
                    ball.intersect(hit, Ray(m_xyz, Vector3(0, 0, 1)), 0, boxDims[2]) ||
                    ball.intersect(hit, Ray(m_XYZ, -Vector3(1, 0, 0)), 0, boxDims[0]) ||
                    ball.intersect(hit, Ray(m_XYZ, -Vector3(0, 1, 0)), 0, boxDims[1]) ||
                    ball.intersect(hit, Ray(m_XYZ, -Vector3(0, 0, 1)), 0, boxDims[2]) ||
                    ball.intersect(hit, Ray(corner100, Vector3(0, 1, 0)), 0, boxDims[1]) ||
                    ball.intersect(hit, Ray(corner100, Vector3(0, 0, 1)), 0, boxDims[2]) ||
                    ball.intersect(hit, Ray(corner010, Vector3(0, 0, 1)), 0, boxDims[2]) ||
                    ball.intersect(hit, Ray(corner010, Vector3(1, 0, 0)), 0, boxDims[0]) ||
                    ball.intersect(hit, Ray(corner001, Vector3(1, 0, 0)), 0, boxDims[0]) ||
                    ball.intersect(hit, Ray(corner001, Vector3(0, 1, 0)), 0, boxDims[1]))
                {
                    child->getPhotons(x, r, photons);
                }
            }
        }
    }
}


void PhotonMap::getNearestPhotons(const Vector3& x, const int& k, std::priority_queue<RsqrPhoton>& photons) {
    PhotonMap* node = getLeafNode(x);
    while (node != this) {
        PhotonMap* parent = node->m_parent;
        PhotonMap* sibling = parent->m_child0;
        if (node == sibling) sibling = parent->m_child1;
        if (node->m_photon == NULL) {
            sibling->getNearestPhotons(x, k, photons);
            node = parent;
            continue;
        }
        PhotonDeposit photon = *node->m_photon;
        float r2 = (photon.location() - x).length2();
        if (photons.size() < k) photons.push(RsqrPhoton(r2, photon));
        else if (r2 < photons.top().m_r2) {
            photons.pop();
            photons.push(RsqrPhoton(r2, photon));
        }
        if (photons.size() < k) sibling->getNearestPhotons(x, k, photons);
        else {
            float d2 = parent->m_photon->location()[parent->m_axis] - x[parent->m_axis];
            d2 *= d2;
            if (d2 < photons.top().m_r2) sibling->getNearestPhotons(x, k, photons);
        }
        node = parent;
    }
    // we've reached this node
    if (node->m_photon == NULL) return;
    PhotonDeposit photon = *node->m_photon;
    float r2 = (photon.location() - x).length2();
    if (photons.size() < k) photons.push(RsqrPhoton(r2, photon));
    else if (r2 < photons.top().m_r2) {
        photons.pop();
        photons.push(RsqrPhoton(r2, photon));
    }
}
std::vector<PhotonDeposit> PhotonMap::getNearestPhotons(const Vector3& x, const int& k) {
    std::priority_queue<RsqrPhoton> photonQueue;
    getNearestPhotons(x, k, photonQueue);
    int n = photonQueue.size();
    std::vector<PhotonDeposit> photons(n);
    for (int i = 0; i < n; i++) {
        photons[n - i - 1] = photonQueue.top().m_photon;
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
    float r = (rdp.m_photons.back().location() - x).length();
    float r2 = (rdp.m_photons.back().location() - x).length2();
    rdp.m_radius = r;
    for (int i = 0; i < rdp.m_photons.size(); i++) rdp.m_density += rdp.m_photons[i].m_power;
    rdp.m_density /= (M_PI*r2);
    return rdp;
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////





void SequentialPhotonMap::addPhoton(const PhotonDeposit& photon) {
    m_photons.push_back(photon);
    if (m_photons.size() == 1) {
        m_xyz = photon.location();
        m_XYZ = photon.location();
    }
    if (photon.location()[0] < m_xyz.x) m_xyz.x = photon.location()[0];
    if (photon.location()[1] < m_xyz.y) m_xyz.y = photon.location()[1];
    if (photon.location()[2] < m_xyz.z) m_xyz.z = photon.location()[2];
    if (photon.location()[0] > m_XYZ.x) m_XYZ.x = photon.location()[0];
    if (photon.location()[1] > m_XYZ.y) m_XYZ.y = photon.location()[1];
    if (photon.location()[2] > m_XYZ.z) m_XYZ.z = photon.location()[2];
}

Vector3 SequentialPhotonMap::powerDensity(const Vector3& x, const float& r) {
    Vector3 rho(0, 0, 0);
    for (int i = 0; i < m_photons.size(); i++) {
        if ((m_photons[i].location() - x).length2() < r*r) {
            rho += m_photons[i].m_power / (M_PI*r*r);
        }
    }
    return rho;
}


RadiusDensityPhotons SequentialPhotonMap::radiusDensityPhotons(const Vector3& x, const int& n) {
    std::vector<std::pair<float, PhotonDeposit>> displacement2(m_photons.size());
    for (int i = 0; i < displacement2.size(); i++) {
        displacement2[i] = std::pair<float, PhotonDeposit>((m_photons[i].location() - x).length2(), m_photons[i]);
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

PhotonMap* SequentialPhotonMap::buildTree() {
    PhotonMap* photonMap = new PhotonMap(m_xyz, m_XYZ);
    for (int i = 0; i < m_photons.size(); i++) photonMap->addPhoton(m_photons[i]);
    return photonMap;
}
void SequentialPhotonMap::buildBalancedTree(int& nPhotons, PhotonMap*& photonMap, std::vector<PhotonDeposit> photons, int depth, bool verbose) {
    if (photons.size() == 0) return;
    int medianIndex = photons.size() / 2;
    if (depth % 3 == 0) std::nth_element(photons.begin(), photons.begin() + medianIndex, photons.end(), compareX);
    else if (depth % 3 == 1) std::nth_element(photons.begin(), photons.begin() + medianIndex, photons.end(), compareY);
    else std::nth_element(photons.begin(), photons.begin() + medianIndex, photons.end(), compareZ);
    photonMap->addPhoton(photons[medianIndex]);
    if (verbose == true && nPhotons % 100 == 0) printf("adding photon %i __________\r", nPhotons);
    nPhotons++;
    if (medianIndex > 0) {
        std::vector<PhotonDeposit>photonsL(photons.begin(), photons.begin() + medianIndex);
        buildBalancedTree(nPhotons, photonMap, photonsL, depth + 1, verbose);
    }
    if (medianIndex + 1 < photons.size()) {
        std::vector<PhotonDeposit>photonsR(photons.begin() + medianIndex + 1, photons.end());
        buildBalancedTree(nPhotons, photonMap, photonsR, depth + 1, verbose);
    }
}
PhotonMap* SequentialPhotonMap::buildBalancedTree(int depth, bool verbose) {
    PhotonMap* photonMap = new PhotonMap(m_xyz, m_XYZ);
    std::vector<PhotonDeposit> photons = m_photons;
    int nPhotons = 0;
    if (verbose == true) buildBalancedTree(nPhotons, photonMap, photons, 0, true);
    else buildBalancedTree(nPhotons, photonMap);
    return photonMap;
}