#ifndef CSE168_POLYGON_MESH_H_INCLUDED
#define CSE168_POLYGON_MESH_H_INCLUDED

#include "Matrix4x4.h"
#include "Scene.h"
#include <vector>

class PolygonMesh
{
public:
    struct TupleI3
    {
        unsigned int m_x, m_y, m_z;
        TupleI3() : m_x(0), m_y(0), m_z(0) {}
        TupleI3(int x, int y, int z) : m_x(x), m_y(y), m_z(z) {}
    };
    struct TupleI4
    {
        unsigned int m_w, m_x, m_y, m_z;
        TupleI4() : m_w(0), m_x(0), m_y(0), m_z(0) {}
        TupleI4(int w, int x, int y, int z) : m_w(w), m_x(x), m_y(y), m_z(z) {}
    };

    struct VectorR2
    {
        float x, y;
    };

    PolygonMesh();
    PolygonMesh(const std::vector<Vector3>& vertices, const std::vector<TupleI3>& vertexIndices);
    ~PolygonMesh();

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
    std::vector<TupleI3> triangleVertexIndices()     {return m_triangleVertexIndices;}
    std::vector<TupleI3> triangleNormalIndices()     {return m_triangleNormalIndices;}
    int numTris()           {return m_numTris;}

    Vector3 vertex(const int& i) { return m_vertices[i]; }
    Vector3 normal(const int& i) { return m_normals[i]; }
    TupleI3 triangleVertexIndex(const int& i) { return m_triangleVertexIndices[i]; }
    TupleI3 triangleNormalIndex(const int& i) { return m_triangleNormalIndices[i]; }
    TupleI4 quadVertexIndex(const int& i) { return m_quadVertexIndices[i]; }
    TupleI4 quadNormalIndex(const int& i) { return m_quadNormalIndices[i]; }

    void addVertex(const Vector3& v) { m_vertices.push_back(v); }
    void addVertices(const std::vector<Vector3>& v) {
        for (int i = 0; i < v.size(); i++) m_vertices.push_back(v[i]);
    }
    void addTriangle(const TupleI3& vertexIndices, const TupleI3& normalIndices) {
        m_triangleVertexIndices.push_back(vertexIndices);
        m_triangleNormalIndices.push_back(normalIndices);
    }
    void addTriangle(const std::vector<TupleI3>& vertexIndices, const std::vector<TupleI3>& normalIndices) {
        for (int i = 0; i < vertexIndices.size(); i++) {
            m_triangleVertexIndices.push_back(vertexIndices[i]);
            m_triangleNormalIndices.push_back(normalIndices[i]);
        }
    }
    void addMeshToScene(Scene* scene);

protected:
    void loadObj(FILE* fp, const Matrix4x4& ctm);

    std::vector<Vector3> m_normals;
    std::vector<Vector3> m_vertices;
    std::vector<VectorR2> m_texCoords;

    std::vector<TupleI3> m_triangleNormalIndices;
    std::vector<TupleI3> m_triangleVertexIndices;
    std::vector<TupleI3> m_triangleTexCoordIndices;

    std::vector<TupleI4> m_quadNormalIndices;
    std::vector<TupleI4> m_quadVertexIndices;
    std::vector<TupleI4> m_quadTexCoordIndices;

    unsigned int m_numTris;
};


#endif // CSE168_POLYGON_MESH_H_INCLUDED
