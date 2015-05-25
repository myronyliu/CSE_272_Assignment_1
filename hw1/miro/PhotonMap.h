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
    void buildOctree();
};


// adopted from http://www.brandonpelfrey.com/blog/coding-a-simple-octree/
class PhotonMap {
protected:
    Vector3 m_origin;
    Vector3 m_halfDimensions;
    PhotonMap* m_children[8];
    PhotonDeposit* m_photon; // the photon associated with this octant if this is a leaf node
    // Children follow a predictable pattern to make accesses simple. - means less than 'origin' in that dimension, + means greater than.
    // child:	0 1 2 3 4 5 6 7
    // x:       - - - - + + + +
    // y:       - - + + - - + +
    // z:       - + - + - + - +
public:
    PhotonMap() : m_origin(Vector3(0, 0, 0)), m_halfDimensions(Vector3(0, 0, 0)), m_photon(NULL) { for (int i = 0; i < 8; i++) m_children[i] = NULL; }
    PhotonMap(const Vector3& origin, const Vector3& halfDimensions) : m_origin(origin), m_halfDimensions(halfDimensions), m_photon(NULL) { for (int i = 0; i < 8; i++) m_children[i] = NULL; }
    PhotonMap(const PhotonMap& copy) : m_origin(copy.m_origin), m_halfDimensions(copy.m_halfDimensions), m_photon(copy.m_photon) {}
    ~PhotonMap() { for (int i = 0; i < 8; i++) delete m_children[i]; } // recursively delete children
    int getOctant(const Vector3& point) const; // Determine which octant of the tree would contain 'point'
    bool isLeafNode() const { return m_children[0] == NULL; }
    void addPhoton(PhotonDeposit* photon);
    void getPhotons(const Vector3& bmin, const Vector3& bmax, std::vector<PhotonDeposit*>& results);
    std::vector<PhotonDeposit> getAllPhotons() {
        std::vector<PhotonDeposit*> photonPointers;
        std::vector<PhotonDeposit> photons;
        float s = 0.05;
        //getPhotons(m_origin - s*m_halfDimensions, m_origin + s*m_halfDimensions, photonPointers);
        Vector3 x(0, 0, 2);
        getPhotons(x - s*m_halfDimensions, x + s*m_halfDimensions, photonPointers);
        for (int i = 0; i < photonPointers.size(); i++) photons.push_back(*photonPointers[i]);
        return photons;
    }
    void buildOctree(SequentialPhotonMap spm);
    RadiusDensityPhotons radiusDensityPhotons(const Vector3& x, const int& n);

    Vector3 origin() { return m_origin; }
    Vector3 halfDimensions() { return m_halfDimensions; }
};

#endif // CSE168_PHOTONMAP_H_INCLUDED