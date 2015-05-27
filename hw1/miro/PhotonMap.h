#ifndef CSE168_PHOTONMAP_H_INCLUDED
#define CSE168_PHOTONMAP_H_INCLUDED

#include <list>
#include <vector>
#include <algorithm>
#include "Vector3.h"


struct PhotonDeposit {
    PhotonDeposit() : m_power(0), m_location(Vector3(0, 0, 0)) {}
    PhotonDeposit(const Vector3& power, const Vector3& location) : m_power(power), m_location(location) {}
    PhotonDeposit(const PhotonDeposit& copy) : m_power(copy.m_power), m_location(copy.m_location) {}
    Vector3 m_power;
    Vector3 m_location;
};

struct RadiusDensityPhotons {
    float m_radius;
    Vector3 m_density;
    std::vector<PhotonDeposit> m_photons;

    RadiusDensityPhotons() : m_radius(0), m_density(0), m_photons(std::vector<PhotonDeposit>(0)) {}
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
    Vector3 m_xyz;
    Vector3 m_XYZ;
    PhotonMap* m_parent;
    PhotonMap* m_child0;
    PhotonMap* m_child1;
    PhotonDeposit* m_photon; // the photon associated with this octant if this is a leaf node
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
    }
    PhotonMap(const PhotonMap& copy) :
        m_xyz(copy.m_xyz), m_XYZ(copy.m_XYZ), m_photon(copy.m_photon), m_parent(copy.m_parent), m_child0(copy.m_child0), m_child1(copy.m_child1) {}
    ~PhotonMap() { delete m_child0; delete m_child1; } // recursively delete children
    
    PhotonMap* getDeepestNode(const Vector3& x);
    bool isLeafNode() const { return m_child0 == NULL; }
    void addPhoton(PhotonDeposit photon);
    void getPhotons(const Vector3& bmin, const Vector3& bmax, std::vector<PhotonDeposit>& photons);
    std::vector<PhotonDeposit> getPhotons(const Vector3& bmin, const Vector3& bmax) {
        std::vector<PhotonDeposit> photons; getPhotons(bmin, bmax, photons); return photons;
    }
    std::vector<PhotonDeposit> getPhotons() { return getPhotons(m_xyz, m_XYZ); }
    std::vector<PhotonDeposit> getNearestPhotons(const Vector3& x, const int& n); // returns the n nearest neighbors of input location x
    void buildTree(SequentialPhotonMap spm);
    void buildBalancedTree(std::vector<PhotonDeposit> spm, int depth = 0);
    void buildBalancedTree(SequentialPhotonMap spm);
    RadiusDensityPhotons radiusDensityPhotons(const Vector3& x, const int& n);
};

#endif // CSE168_PHOTONMAP_H_INCLUDED