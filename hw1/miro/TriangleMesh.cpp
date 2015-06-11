#include "TriangleMesh.h"
#include "Triangle.h"


TriangleMesh::TriangleMesh() :
m_normals(0),
m_vertices(0),
m_texCoords(0),
m_normalIndices(0),
m_vertexIndices(0),
m_texCoordIndices(0)
{}

TriangleMesh::TriangleMesh(const std::vector<Vector3>& vertices, const std::vector<TupleI3>& vertexIndices) :
m_normals(std::vector<Vector3>(vertexIndices.size())),
m_vertices(vertices),
m_texCoords(0),
m_normalIndices(std::vector<TupleI3>(vertexIndices.size())),
m_vertexIndices(vertexIndices),
m_texCoordIndices(0)
{
    for (int i = 0; i < vertexIndices.size(); i++) {
        TupleI3 ti3 = vertexIndices[i];
        Vector3 A = vertices[ti3.m_x];
        Vector3 B = vertices[ti3.m_y];
        Vector3 C = vertices[ti3.m_z];
        m_normals[i] = cross(B - A, C - A).normalize();
        m_normalIndices[i] = TupleI3(i, i, i);
    }
}

TriangleMesh::~TriangleMesh()
{
    /*delete [] m_normals;
    delete [] m_vertices;
    delete [] m_texCoords;
    delete [] m_normalIndices;
    delete [] m_vertexIndices;
    delete [] m_texCoordIndices;*/
}

void TriangleMesh::addMeshToScene(Scene* scene) {
    for (int i = 0; i < m_vertexIndices.size(); i++) {
        Triangle* triangle = new Triangle(this, i);
        scene->addObject(triangle);
    }
}