#include "PolygonMesh.h"
#include "Triangle.h"


PolygonMesh::PolygonMesh() :
m_normals(0),
m_vertices(0),
m_texCoords(0),
m_triangleNormalIndices(0),
m_triangleVertexIndices(0),
m_triangleTexCoordIndices(0)
{}

PolygonMesh::PolygonMesh(const std::vector<Vector3>& vertices, const std::vector<TupleI3>& vertexIndices) :
m_normals(std::vector<Vector3>(vertexIndices.size())),
m_vertices(vertices),
m_texCoords(0),
m_triangleNormalIndices(std::vector<TupleI3>(vertexIndices.size())),
m_triangleVertexIndices(vertexIndices),
m_triangleTexCoordIndices(0)
{
    for (int i = 0; i < vertexIndices.size(); i++) {
        TupleI3 ti3 = vertexIndices[i];
        Vector3 A = vertices[ti3.m_x];
        Vector3 B = vertices[ti3.m_y];
        Vector3 C = vertices[ti3.m_z];
        m_normals[i] = cross(B - A, C - A).normalize();
        m_triangleNormalIndices[i] = TupleI3(i, i, i);
    }
}

PolygonMesh::~PolygonMesh()
{
    /*delete [] m_normals;
    delete [] m_vertices;
    delete [] m_texCoords;
    delete [] m_triangleNormalIndices;
    delete [] m_triangleVertexIndices;
    delete [] m_triangleTexCoordIndices;*/
}

void PolygonMesh::addMeshToScene(Scene* scene) {
    for (int i = 0; i < m_triangleVertexIndices.size(); i++) {
        Triangle* triangle = new Triangle(this, i);
        scene->addObject(triangle);
    }
}