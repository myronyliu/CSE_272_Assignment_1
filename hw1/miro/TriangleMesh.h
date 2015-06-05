#ifndef CSE168_TRIANGLE_MESH_H_INCLUDED
#define CSE168_TRIANGLE_MESH_H_INCLUDED

#include "Matrix4x4.h"
#include <vector>

class TriangleMesh
{
public:
    struct TupleI3
    {
        unsigned int m_x, m_y, m_z;
        TupleI3() : m_x(0), m_y(0), m_z(0) {}
        TupleI3(int x, int y, int z) : m_x(x), m_y(y), m_z(z) {}
    };

    struct VectorR2
    {
        float x, y;
    };

    TriangleMesh();
    TriangleMesh(const std::vector<Vector3>& vertices, const std::vector<TupleI3>& vertexIndices);
    ~TriangleMesh();

    // load from an OBJ file
    bool load(char* file, const Matrix4x4& ctm = Matrix4x4());

    // for single triangles
    void createSingleTriangle();
    inline void setV1(const Vector3& v) {m_vertices[0] = v;}
    inline void setV2(const Vector3& v) {m_vertices[1] = v;}
    inline void setV3(const Vector3& v) {m_vertices[2] = v;}
    inline void setN1(const Vector3& n) {m_normals[0] = n;}
    inline void setN2(const Vector3& n) {m_normals[1] = n;}
    inline void setN3(const Vector3& n) {m_normals[2] = n;}

    

    std::vector<Vector3> vertices()     {return m_vertices;}
    std::vector<Vector3> normals()      {return m_normals;}
    std::vector<TupleI3> vIndices()     {return m_vertexIndices;}
    std::vector<TupleI3> nIndices()     {return m_normalIndices;}
    int numTris()           {return m_numTris;}

    Vector3 vertex(const int& i) { return m_vertices[i]; }
    Vector3 normal(const int& i) { return m_normals[i]; }
    TupleI3 vIndex(const int& i) { return m_vertexIndices[i]; }
    TupleI3 nIndex(const int& i) { return m_normalIndices[i]; }

    void addVertex(const Vector3& v) { m_vertices.push_back(v); }
    void addVertices(const std::vector<Vector3>& v) {
        for (int i = 0; i < v.size(); i++) m_vertices.push_back(v[i]);
    }
    void addTriangle(const TupleI3& vertexIndices, const TupleI3& normalIndices) {
        m_vertexIndices.push_back(vertexIndices);
        m_normalIndices.push_back(normalIndices);
    }
    void addTriangle(const std::vector<TupleI3>& vertexIndices, const std::vector<TupleI3>& normalIndices) {
        for (int i = 0; i < vertexIndices.size(); i++) {
            m_vertexIndices.push_back(vertexIndices[i]);
            m_normalIndices.push_back(normalIndices[i]);
        }
    }

protected:
    void loadObj(FILE* fp, const Matrix4x4& ctm);

    std::vector<Vector3> m_normals;
    std::vector<Vector3> m_vertices;
    std::vector<VectorR2> m_texCoords;

    std::vector<TupleI3> m_normalIndices;
    std::vector<TupleI3> m_vertexIndices;
    std::vector<TupleI3> m_texCoordIndices;
    unsigned int m_numTris;
};


#endif // CSE168_TRIANGLE_MESH_H_INCLUDED
