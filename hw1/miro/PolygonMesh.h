#ifndef CSE168_POLYGON_MESH_H_INCLUDED
#define CSE168_POLYGON_MESH_H_INCLUDED

#include "Matrix4x4.h"
#include "Scene.h"
#include "Lambert.h"
#include "Mirror.h"
#include <vector>

class PolygonMesh
{
public:
    struct TupleI4
    {
        int m_a, m_b, m_c, m_d;
        TupleI4() : m_a(0), m_b(0), m_c(0), m_d(0) {}
        TupleI4(int a, int b, int c, int d) : m_a(a), m_b(b), m_c(c), m_d(d) {}
    };
    struct TupleI3
    {
        int m_a, m_b, m_c;
        TupleI3() : m_a(0), m_b(0), m_c(0) {}
        TupleI3(int a, int b, int c) : m_a(a), m_b(b), m_c(c) {}
        TupleI3(TupleI4 tuple4) : m_a(tuple4.m_a), m_b(tuple4.m_b), m_c(tuple4.m_c) {}
    };
    struct VectorR2
    {
        float m_x, m_y;
        VectorR2() : m_x(0), m_y(0) {}
        VectorR2(float x, float y) : m_x(x), m_y(y) {}
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
    std::vector<TupleI3> triangleNormalIndices()     { return m_triangleNormalIndices; }
    int numTriangles()           { return m_triangleVertexIndices.size(); }
    int numQuads()           { return m_quadVertexIndices.size(); }

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

    std::vector<Material*> m_triangleMaterials;
    std::vector<Material*> m_quadMaterials;
    std::vector<bool> m_triangleIsLight;
    std::vector<bool> m_quadIsLight;

    Material* m_diffuse50 = new Lambert(Vector3(0.5, 0.5, 0.5), Vector3(0, 0, 0));
    Material* m_mirror = new Mirror(Vector3(1, 1, 1), Vector3(0, 0, 0));
    Material* m_absorbing = new Lambert(Vector3(0, 0, 0), Vector3(0, 0, 0));

};


#endif // CSE168_POLYGON_MESH_H_INCLUDED
