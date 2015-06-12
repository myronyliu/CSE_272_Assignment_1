#include "PolygonMesh.h"
#include "Console.h"
#include <algorithm>

#ifdef WIN32
// disable useless warnings
#pragma warning(disable:4996)
#endif

void
PolygonMesh::createSingleTriangle()
{
    m_normals.resize(3);
    m_vertices.resize(3);
    m_texCoords.resize(3);

    m_texCoords[0].m_x = 0.0f;
    m_texCoords[0].m_y = 0.0f;
    m_texCoords[1].m_x = 1.0f;
    m_texCoords[1].m_y = 0.0f;
    m_texCoords[2].m_x = 0.0f;
    m_texCoords[2].m_y = 1.0f;

    m_triangleNormalIndices.resize(1);
    m_triangleVertexIndices.resize(1);
    m_triangleTexCoordIndices.resize(1);

    m_triangleVertexIndices[0].m_a = 0;
    m_triangleVertexIndices[0].m_b = 1;
    m_triangleVertexIndices[0].m_c = 2;

    m_triangleNormalIndices[0].m_a = 0;
    m_triangleNormalIndices[0].m_b = 1;
    m_triangleNormalIndices[0].m_c = 2;

    m_triangleTexCoordIndices[0].m_a = 0;
    m_triangleTexCoordIndices[0].m_b = 1;
    m_triangleTexCoordIndices[0].m_c = 2;
}

//************************************************************************
// You probably don't want to modify the following functions
// They are for loading .obj files
//************************************************************************

bool
PolygonMesh::load(char* file, const Matrix4x4& ctm)
{
    FILE* fp;
    fopen_s(&fp, file, "rb");
    if (!fp)
    {
        error("Cannot open \"%s\" for reading\n", file);
        return false;
    }
    debug("Loading \"%s\"...\n", file);

    loadObj(fp, ctm);
    debug("Loaded \"%s\" with %d triangles\n", file, numTriangles());
    fclose(fp);

    return true;
}

void
getIndices(char *word, int *vindex, int *tindex, int *nindex)
{
    char *null = " ";
    char *ptr;
    char *tp;
    char *np;

    // by default, the texture and normal pointers are set to the null string
    tp = null;
    np = null;

    // replace slashes with null characters and cause tp and np to point
    // to character immediately following the first or second slash
    for (ptr = word; *ptr != '\0'; ptr++)
    {
        if (*ptr == '/')
        {
            if (tp == null)
                tp = ptr + 1;
            else
                np = ptr + 1;

            *ptr = '\0';
        }
    }

    *vindex = atoi(word);
    *tindex = atoi(tp);
    *nindex = atoi(np);
}


void
PolygonMesh::loadObj(FILE* fp, const Matrix4x4& ctm)
{
    int nv = 0, nt = 0, nn = 0, nf = 0;
    char line[81];

    Matrix4x4 nctm = ctm;
    nctm.invert();
    nctm.transpose();
    nctm.invert();

    int material = -1;

    while (fgets(line, 80, fp) != 0)
    {
        if (line[0] == 'v')
        {
            if (line[1] == 'n')
            {
                float x, y, z;
                sscanf_s(&line[2], "%f %f %f\n", &x, &y, &z);
                Vector3 n(x, y, z);
                m_normals.push_back((nctm*n).normalize());
            }
            else if (line[1] == 't')
            {
                float x, y;
                sscanf_s(&line[2], "%f %f\n", &x, &y);
                m_texCoords.push_back(VectorR2(x, y));
            }
            else
            {
                float x, y, z;
                sscanf_s(&line[1], "%f %f %f\n", &x, &y, &z);
                Vector3 v(x, y, z);
                m_vertices.push_back(ctm*v);
            }
        }
        else if (line[0] == 'f')
        {
            char s1[32], s2[32], s3[32], s4[32];
            int v, t, n;

            int nPolygonVertices = sscanf_s(&line[2], "%s %s %s %s\n", s1, sizeof(s1), s2, sizeof(s2), s3, sizeof(s3), s4, sizeof(s4));

            TupleI4 vIndices(-1, -1, -1, -1);
            TupleI4 nIndices(-1, -1, -1, -1);
            TupleI4 tIndices(-1, -1, -1, -1);
            getIndices(s1, &v, &t, &n);
            if (v) vIndices.m_a = v - 1;
            if (n) nIndices.m_a = n - 1;
            if (t) tIndices.m_a = t - 1;
            getIndices(s2, &v, &t, &n);
            if (v) vIndices.m_b = v - 1;
            if (n) nIndices.m_b = n - 1;
            if (t) tIndices.m_b = t - 1;
            getIndices(s3, &v, &t, &n);
            if (v) vIndices.m_c = v - 1;
            if (n) nIndices.m_c = n - 1;
            if (t) tIndices.m_c = t - 1;
            if (nPolygonVertices > 3) {
                getIndices(s4, &v, &t, &n);
                if (v) vIndices.m_d = v - 1;
                if (n) nIndices.m_d = n - 1;
                if (t) tIndices.m_d = t - 1;
                if (vIndices.m_a >= 0) m_quadVertexIndices.push_back(vIndices);
                if (tIndices.m_a >= 0) m_quadTexCoordIndices.push_back(tIndices);
                if (nIndices.m_a >= 0) m_quadNormalIndices.push_back(nIndices);
                else {   // if no normal was supplied
                    Vector3 e1 = m_vertices[vIndices.m_d] - m_vertices[vIndices.m_a];
                    Vector3 e2 = m_vertices[vIndices.m_b] - m_vertices[vIndices.m_a];
                    m_normals.push_back(cross(e1, e2).normalize());
                    m_quadNormalIndices.push_back(TupleI4(m_normals.size() - 1, m_normals.size() - 1, m_normals.size() - 1, m_normals.size() - 1));
                }
                if (material == 0) {
                    m_quadMaterials.push_back(m_diffuse50);
                    m_quadIsLight.push_back(false);
                }
                else if (material == 1) {
                    m_quadMaterials.push_back(m_mirror);
                    m_quadIsLight.push_back(false);
                }
                else if (material == 2) {
                    m_quadMaterials.push_back(m_absorbing);
                    m_quadIsLight.push_back(true);
                }
            }
            else {
                if (vIndices.m_a >= 0) m_triangleVertexIndices.push_back(TupleI3(vIndices));
                if (tIndices.m_a >= 0) m_triangleTexCoordIndices.push_back(TupleI3(tIndices));
                if (nIndices.m_a >= 0) m_triangleNormalIndices.push_back(TupleI3(nIndices));
                else {   // if no normal was supplied
                    Vector3 e1 = m_vertices[vIndices.m_b] - m_vertices[vIndices.m_a];
                    Vector3 e2 = m_vertices[vIndices.m_c] - m_vertices[vIndices.m_a];
                    m_normals.push_back(cross(e1, e2).normalize());
                    m_triangleNormalIndices.push_back(TupleI3(m_normals.size() - 1, m_normals.size() - 1, m_normals.size() - 1));
                }
                if (material == 0) {
                    m_triangleMaterials.push_back(m_diffuse50);
                    m_triangleIsLight.push_back(false);
                }
                else if (material == 1) {
                    m_triangleMaterials.push_back(m_mirror);
                    m_triangleIsLight.push_back(false);
                }
                else if (material == 2) {
                    m_triangleMaterials.push_back(m_absorbing);
                    m_triangleIsLight.push_back(true);
                }
            }
        }
        else if (line[0] == 'u') {
            char s1[32], s2[32];
            sscanf_s(&line[0], "%s %s\n", s1, sizeof(s1), s2, sizeof(s2));
            if (s2[0] == 'm') material = 1;
            else if (s2[7] == '5') material = 0;
            else if (s2[7] == 'L') material = 2;
        }
    }
}

