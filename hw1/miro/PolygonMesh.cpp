#include "PolygonMesh.h"
#include "Triangle.h"
#include "Quad.h"
#include "TriangleLight.h"
#include "QuadLight.h"
#include "ParallelogramLight.h"


PolygonMesh::PolygonMesh() :
m_normals(0),
m_vertices(0),
m_texCoords(0),
m_triangleNormalIndices(0),
m_triangleVertexIndices(0),
m_triangleTexCoordIndices(0)
{
}

PolygonMesh::PolygonMesh(const std::vector<Vector3>& vertices, const std::vector<TupleI3>& triangleVertexIndices, const std::vector<TupleI4>& quadVertexIndices) :
m_normals(std::vector<Vector3>(triangleVertexIndices.size()+quadVertexIndices.size())),
m_vertices(vertices),
m_texCoords(0),
m_triangleNormalIndices(std::vector<TupleI3>(triangleVertexIndices.size())),
m_triangleVertexIndices(triangleVertexIndices),
m_triangleTexCoordIndices(0),
m_quadNormalIndices(std::vector<TupleI4>(quadVertexIndices.size())),
m_quadVertexIndices(quadVertexIndices),
m_quadTexCoordIndices(0)
{
    for (int i = 0; i < triangleVertexIndices.size(); i++) {
        TupleI3 ti3 = triangleVertexIndices[i];
        Vector3 A = vertices[ti3.m_a];
        Vector3 B = vertices[ti3.m_b];
        Vector3 C = vertices[ti3.m_c];
        m_normals[i] = cross(B - A, C - A).normalize();
        m_triangleNormalIndices[i] = TupleI3(i, i, i);
    }
    for (int i = 0; i < quadVertexIndices.size(); i++) {
        TupleI4 ti4 = quadVertexIndices[i];
        Vector3 A = vertices[ti4.m_a];
        Vector3 B = vertices[ti4.m_b];
        Vector3 D = vertices[ti4.m_d];
        m_normals[i] = cross(B - A, D - A).normalize();
        m_quadNormalIndices[i] = TupleI4(i, i, i, i);
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
            //QuadLight* quadLight = new QuadLight(this, i);
            //quadLight->setMaterial(m_quadMaterials[i]);
            //quadLight->setWattage(quadLight->area());
            //scene->addObject(quadLight);

            Vector3 A = m_vertices[m_quadVertexIndices[i].m_a];
            Vector3 B = m_vertices[m_quadVertexIndices[i].m_b];
            Vector3 C = m_vertices[m_quadVertexIndices[i].m_c];
            Vector3 D = m_vertices[m_quadVertexIndices[i].m_d];
            Vector3 center = (A + C) / 2;
            Vector3 vecX = (B - A).normalize();
            Vector3 vecY = (D - A).normalize();
            float spanX = (B - A).length() / 2;
            float spanY = (D - A).length() / 2;
            ParallelogramLight* light = new ParallelogramLight(center,vecX,vecY,spanX,spanY);
            light->setMaterial(m_quadMaterials[i]);
            light->setWattage(8*light->area());
            scene->addAreaLight(light,20000);
        }
    }
}