#ifndef CSE168_PHOTONMAP_H_INCLUDED
#define CSE168_PHOTONMAP_H_INCLUDED

#include <vector>
#include <queue>
#include <algorithm>
#include <functional>
#include "Vector3.h"

struct PhotonDeposit {
    PhotonDeposit() : m_power(Vector3(0, 0, 0)), m_location(Vector3(0, 0, 0)), m_prob(0) {}
    PhotonDeposit(const Vector3& power, const Vector3& location, const float& prob) : m_power(power), m_location(location), m_prob(prob) {}
    PhotonDeposit(const PhotonDeposit& copy) : m_power(copy.m_power), m_location(copy.m_location), m_prob(copy.m_prob) {}
    Vector3 m_location;
    Vector3 m_power;
    float m_prob;
};

struct RadiusDensityPhotons {
    float m_radius;
    Vector3 m_density;
    std::vector<PhotonDeposit> m_photons;

    RadiusDensityPhotons() : m_radius(0), m_density(0), m_photons(std::vector<PhotonDeposit>(0)) {}
};


struct RsqrPhoton {
    float m_r2;
    PhotonDeposit m_photon;
    RsqrPhoton() : m_r2(0), m_photon(PhotonDeposit()) {}
    RsqrPhoton(const float& r2, const PhotonDeposit& photon) : m_r2(r2), m_photon(photon) {}
    bool operator<(const RsqrPhoton& rhs) const {
        if (m_r2 < rhs.m_r2) return true;
        else if (m_r2 > rhs.m_r2) return false;
        else if (m_photon.m_location[0] < rhs.m_photon.m_location[0]) return true;
        else if (m_photon.m_location[0] > rhs.m_photon.m_location[0]) return false;
        else if (m_photon.m_location[1] < rhs.m_photon.m_location[1]) return true;
        else if (m_photon.m_location[1] > rhs.m_photon.m_location[1]) return false;
        else if (m_photon.m_location[2] < rhs.m_photon.m_location[2]) return true;
        else if (m_photon.m_location[2] > rhs.m_photon.m_location[2]) return false;
        else if (m_photon.m_power[0] < rhs.m_photon.m_power[0]) return true;
        else if (m_photon.m_power[0] > rhs.m_photon.m_power[0]) return false;
        else if (m_photon.m_power[1] < rhs.m_photon.m_power[1]) return true;
        else if (m_photon.m_power[1] > rhs.m_photon.m_power[1]) return false;
        else if (m_photon.m_power[2] < rhs.m_photon.m_power[2]) return true;
        else if (m_photon.m_power[2] > rhs.m_photon.m_power[2]) return false;
        else return false;
    }
};


class SequentialPhotonMap {
protected:
    float m_xMin;
    float m_yMin;
    float m_zMin;
    float m_xMax;
    float m_yMax;
    float m_zMax;
    std::vector<PhotonDeposit> m_photonDeposits;
public:
    SequentialPhotonMap() : m_photonDeposits(std::vector<PhotonDeposit>(0)), m_xMin(0), m_yMin(0), m_zMin(0), m_xMax(0), m_yMax(0), m_zMax(0) {}
    ~SequentialPhotonMap() {}

    float xMin() { return m_xMin; }
    float yMin() { return m_yMin; }
    float zMin() { return m_zMin; }
    float xMax() { return m_xMax; }
    float yMax() { return m_yMax; }
    float zMax() { return m_zMax; }
    int nPhotons() { return m_photonDeposits.size(); }
    PhotonDeposit operator[](int i) const { return m_photonDeposits[i]; }
    PhotonDeposit& operator[](int i) { return m_photonDeposits[i]; }
    Vector3 powerDensity(const Vector3& x, const float& r);
    void addPhoton(const PhotonDeposit& photonDeposit);
    RadiusDensityPhotons radiusDensityPhotons(const Vector3& x, const int& n);
    //void buildTree();

    std::vector<PhotonDeposit> getPhotons() { return m_photonDeposits; }
};

// adopted from http://www.brandonpelfrey.com/blog/coding-a-simple-octree/
class PhotonMap {
protected:
    int m_depth;
    int m_axis;
    Vector3 m_xyz;
    Vector3 m_XYZ;
    PhotonMap* m_parent;
    PhotonMap* m_child0;
    PhotonMap* m_child1;
    PhotonDeposit* m_photon; // the photon associated with this octant if this is a leaf node

    void setAxis() { m_axis = m_depth % 3; } // this is just here in case one wishes to define a different splitting convention
    void getPhotons(const Vector3& bmin, const Vector3& bmax, std::vector<PhotonDeposit>& photons);
    void getNearestPhotons(const Vector3& x, const int& k, std::priority_queue<RsqrPhoton>& photons);
public:
    PhotonMap() {};
    PhotonMap(
        const Vector3& xyz, const Vector3& XYZ,
        PhotonDeposit* photon = NULL,
        PhotonMap* parent = NULL,
        PhotonMap* child0 = NULL, PhotonMap* child1 = NULL) :
        m_xyz(xyz), m_XYZ(XYZ), m_photon(photon), m_parent(parent), m_child0(child0), m_child1(child1)
    {
        if (m_parent == NULL) m_depth = 0;
        else m_depth = m_parent->m_depth + 1;
        setAxis();
    }
    PhotonMap(const PhotonMap& copy) :
        m_xyz(copy.m_xyz), m_XYZ(copy.m_XYZ), m_photon(copy.m_photon), m_parent(copy.m_parent), m_child0(copy.m_child0), m_child1(copy.m_child1) {}
    ~PhotonMap() { delete m_child0; delete m_child1; } // recursively delete children

    PhotonMap* getLeafNode(const Vector3& x);
    inline bool isLeafNode() const { return m_child0 == NULL; }
    void addPhoton(PhotonDeposit photon);
    std::vector<PhotonDeposit> getPhotons(const Vector3& bmin, const Vector3& bmax) {
        std::vector<PhotonDeposit> photons;
        getPhotons(bmin, bmax, photons);
        return photons;
    }
    std::vector<PhotonDeposit> getPhotons() { return getPhotons(m_xyz, m_XYZ); }
    std::vector<PhotonDeposit> getNearestPhotons(const Vector3& x, const int& n); // returns the n nearest neighbors of input location x
    void buildTree(SequentialPhotonMap spm);
    void buildBalancedTree(std::vector<PhotonDeposit> spm, int depth = 0);
    void buildBalancedTree(SequentialPhotonMap spm);
    RadiusDensityPhotons radiusDensityPhotons(const Vector3& x, const int& n);

    inline PhotonMap* getSibling() {
        if (m_parent == NULL) return NULL;
        else if (m_parent->m_child0 == this) return m_parent->m_child1;
        else return m_parent->m_child0;
    }
};

#endif // CSE168_PHOTONMAP_H_INCLUDED