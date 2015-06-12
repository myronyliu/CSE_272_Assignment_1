#include "PolygonMesh.h"
#include "Triangle.h"
#include "Quad.h"
#include "TriangleLight.h"
#include "QuadLight.h"


PolygonMesh::PolygonMesh() :
m_normals(0),
m_vertices(0),
m_texCoords(0),
m_triangleNormalIndices(0),
m_triangleVertexIndices(0),
m_triangleTexCoordIndices(0)
{
}

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
        Vector3 A = vertices[ti3.m_a];
        Vector3 B = vertices[ti3.m_b];
        Vector3 C = vertices[ti3.m_c];
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
        if (m_triangleIsLight[i] == false) {
            Triangle* triangle = new Triangle(this, i);
            triangle->setMaterial(m_triangleMaterials[i]);
            scene->addObject(triangle);
        }
        else {
            TriangleLight* triangleLight = new TriangleLight(this, i);
            triangleLight->setMaterial(m_triangleMaterials[i]);
            triangleLight->setWattage(triangleLight->area());
            scene->addObject(triangleLight);
        }
    }
    for (int i = 0; i < m_quadVertexIndices.size(); i++) {
        if (m_quadIsLight[i] == false) {
            Quad* quad = new Quad(this, i);
            quad->setMaterial(m_quadMaterials[i]);
            scene->addObject(quad);
        }
        else {
            QuadLight* quadLight = new QuadLight(this, i);
            quadLight->setMaterial(m_quadMaterials[i]);
            quadLight->setWattage(quadLight->area());
            scene->addObject(quadLight);
        }
    }
}