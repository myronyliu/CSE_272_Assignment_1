#include "TriangleMesh.h"
#include "Triangle.h"
#include "Scene.h"


TriangleMesh::TriangleMesh() :
m_normals(0),
m_vertices(0),
m_texCoords(0),
m_normalIndices(0),
m_vertexIndices(0),
m_texCoordIndices(0)
{}

TriangleMesh::TriangleMesh(const std::vector<Vector3>& vertices, const std::vector<TupleI3>& vertexIndices) :
m_normals(0),
m_vertices(vertices),
m_texCoords(0),
m_normalIndices(0),
m_vertexIndices(vertexIndices),
m_texCoordIndices(0)
{}

TriangleMesh::~TriangleMesh()
{
    /*delete [] m_normals;
    delete [] m_vertices;
    delete [] m_texCoords;
    delete [] m_normalIndices;
    delete [] m_vertexIndices;
    delete [] m_texCoordIndices;*/
}
