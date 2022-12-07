/******************************************************************************
**
**  This file is part of Met.3D -- a research environment for the
**  three-dimensional visual exploration of numerical ensemble weather
**  prediction data.
**
**  Copyright 2015-2020 Marc Rautenhaus
**  Copyright 2022      Julian Münsterberg
**
**  Regional Computing Center, Visualization
**  Universitaet Hamburg, Hamburg, Germany
**
**  Met.3D is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  Met.3D is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with Met.3D.  If not, see <http://www.gnu.org/licenses/>.
**
*******************************************************************************/

#include "tropopausesurfacemesh.h"

// related third party imports
#include "gxfw/gl/typedvertexbuffer.h"
#include "gxfw/gl/indexbuffer.h"

#include <iostream>
#include <fstream>

using namespace Met3D;


/******************************************************************************
***                             MTriangleMesh                               ***
*******************************************************************************/

MTropopauseTriangleMeshSelection::MTropopauseTriangleMeshSelection(uint32_t numVertices,
                             uint32_t numTriangles)
        : MAbstractDataItem(),
          vertices(new Geometry::TropopauseMeshVertex[numVertices]),
          numVertices(numVertices),
          triangles(new Geometry::MTriangle[numTriangles]),
          numTriangles(numTriangles)
{

}

MTropopauseTriangleMeshSelection::~MTropopauseTriangleMeshSelection()
{
    delete[] vertices;
    delete[] triangles;
}

unsigned int MTropopauseTriangleMeshSelection::getMemorySize_kb()
{
    const uint32_t classByteSize = sizeof(MTropopauseTriangleMeshSelection);

    const uint32_t vertexByteSize =
            numVertices * sizeof(Geometry::TropopauseMeshVertex);

    const uint32_t trianglesByteSize =
            sizeof(Geometry::MTriangle) * numTriangles;


    int totalSize_kb = (classByteSize + vertexByteSize + trianglesByteSize) / 1024;

    return totalSize_kb;
}


GL::MVertexBuffer* MTropopauseTriangleMeshSelection::getVertexBuffer(QGLWidget *currentGLContext)
{
    const QString vboID = QString("vbo_triangle_mesh_vertices_%1").arg(getID());

    MGLResourcesManager *glRM = MGLResourcesManager::getInstance();
    glRM->makeCurrent();

    const uint32_t numFloatsPerVertex = sizeof(Geometry::TropopauseMeshVertex)
                                        / sizeof(float);
    const uint32_t numFloats = numVertices * numFloatsPerVertex;

    auto vb = dynamic_cast<GL::MVertexBuffer *>(glRM->getGPUItem(vboID));
    if (vb) { return vb; }

    auto newVb = new GL::MFloatVertexBuffer(vboID, numFloats);

    if (glRM->tryStoreGPUItem(newVb))
    {
        newVb->upload(reinterpret_cast<float*>(vertices),
                      numFloats, currentGLContext);
    }
    else
    {
        LOG4CPLUS_WARN(mlog, "WARNING: cannot store buffer for Tropopause"
                             " vertices in GPU memory.");
        delete newVb;
    }

    return dynamic_cast<GL::MVertexBuffer *>(glRM->getGPUItem(vboID));
}


void MTropopauseTriangleMeshSelection::releaseVertexBuffer()
{
    MGLResourcesManager::getInstance()->releaseGPUItem(getID());
}


GL::MIndexBuffer* MTropopauseTriangleMeshSelection::getIndexBuffer(QGLWidget *currentGLContext)
{
    const QString iboID = QString("ibo_triangle_mesh_indices_%1").arg(getID());//ßTODO:Hier muss auch der richtige indexbuffer rein.

    MGLResourcesManager *glRM = MGLResourcesManager::getInstance();
    glRM->makeCurrent();

    const uint32_t numIntsPerTriangle = sizeof(Geometry::MTriangle)
                                        / sizeof(u_int32_t);

    const uint32_t numIndices = numTriangles * numIntsPerTriangle;

    auto vb = dynamic_cast<GL::MIndexBuffer*>(glRM->getGPUItem(iboID));
    if (vb) { return vb; }

    auto newVb = new GL::MIndexBuffer(iboID, numIndices);

    if (glRM->tryStoreGPUItem(newVb))
    {
        newVb->upload(reinterpret_cast<uint32_t*>(triangles),
                      numIndices, currentGLContext);
    }
    else
    {
        LOG4CPLUS_WARN(mlog, "WARNING: cannot store buffer for Tropopause"
                             " indices in GPU memory.");
        delete newVb;
    }

    return dynamic_cast<GL::MIndexBuffer *>(glRM->getGPUItem(iboID));
}


void MTropopauseTriangleMeshSelection::releaseIndexBuffer()
{
    MGLResourcesManager::getInstance()->releaseGPUItem(getID());
}
