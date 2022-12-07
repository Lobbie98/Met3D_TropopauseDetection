/******************************************************************************
**
**  This file is part of Met.3D -- a research environment for the
**  three-dimensional visual exploration of numerical ensemble weather
**  prediction data.
**
**  Copyright 2017 Marc Rautenhaus
**  Copyright 2017 Michael Kern
**
**  Computer Graphics and Visualization Group
**  Technische Universitaet Muenchen, Garching, Germany
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

#include "frontsurfacemesh.h"

// related third party imports
#include "gxfw/gl/typedvertexbuffer.h"
#include "gxfw/gl/indexbuffer.h"

#include <iostream>
#include <fstream>

using namespace Met3D;


/******************************************************************************
***                             MTriangleMesh                               ***
*******************************************************************************/

MTriangleMeshSelection::MTriangleMeshSelection(uint32_t numVertices,
                             uint32_t numTriangles)
        : MAbstractDataItem(),
          vertices(new Geometry::FrontMeshVertex[numVertices]),
          numVertices(numVertices),
          triangles(new Geometry::MTriangle[numTriangles]),
          numTriangles(numTriangles)
{

}

MTriangleMeshSelection::~MTriangleMeshSelection()
{
    delete[] vertices;
    delete[] triangles;
}

unsigned int MTriangleMeshSelection::getMemorySize_kb()
{
    const uint32_t classByteSize = sizeof(MTriangleMeshSelection);

    const uint32_t vertexByteSize =
            numVertices * sizeof(Geometry::FrontMeshVertex);

    const uint32_t trianglesByteSize =
            sizeof(Geometry::MTriangle) * numTriangles;


    int totalSize_kb = (classByteSize + vertexByteSize + trianglesByteSize) / 1024;

    return totalSize_kb;
}


GL::MVertexBuffer* MTriangleMeshSelection::getVertexBuffer(QGLWidget *currentGLContext)
{
    const QString vboID = QString("vbo_triangle_mesh_vertices_%1").arg(getID());

    MGLResourcesManager *glRM = MGLResourcesManager::getInstance();
    glRM->makeCurrent();

    const uint32_t numFloatsPerVertex = sizeof(Geometry::FrontMeshVertex)
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
        LOG4CPLUS_WARN(mlog, "WARNING: cannot store buffer for 3D front"
                             " vertices in GPU memory.");
        delete newVb;
    }

    return dynamic_cast<GL::MVertexBuffer *>(glRM->getGPUItem(vboID));
}


void MTriangleMeshSelection::releaseVertexBuffer()
{
    MGLResourcesManager::getInstance()->releaseGPUItem(getID());
}


GL::MIndexBuffer* MTriangleMeshSelection::getIndexBuffer(QGLWidget *currentGLContext)
{
    const QString iboID = QString("ibo_triangle_mesh_indices_%1").arg(getID());

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
        LOG4CPLUS_WARN(mlog, "WARNING: cannot store buffer for 3D front"
                             " indices in GPU memory.");
        delete newVb;
    }

    return dynamic_cast<GL::MIndexBuffer *>(glRM->getGPUItem(iboID));
}


void MTriangleMeshSelection::releaseIndexBuffer()
{
    MGLResourcesManager::getInstance()->releaseGPUItem(getID());
}



/******************************************************************************
***                             MNormalCurves                               ***
*******************************************************************************/

MNormalCurvesSelection::MNormalCurvesSelection(uint32_t numNormalCurves)
        : MAbstractDataItem(),
          normalCurves(new Geometry::NormalCurve[numNormalCurves]),
          numNormalCurves(numNormalCurves),
          normalCurveSegment(nullptr),
          numNormalCurveSegments(0),
          restartIndex(nullptr)
{

}

MNormalCurvesSelection::~MNormalCurvesSelection()
{
    delete[] normalCurves;
    if (normalCurveSegment) {delete [] normalCurveSegment;}
    if (restartIndex) {delete [] restartIndex; }
}

unsigned int MNormalCurvesSelection::getMemorySize_kb()
{
    const uint32_t classByteSize = sizeof(MNormalCurvesSelection);

    // get normal curves size
    u_int32_t positionsSize = 0;
    for (u_int32_t i = 0; i < numNormalCurves; i++)
    {
        positionsSize += normalCurves[i].positions.size() * sizeof(QVector3D);
    }

    u_int32_t filteredPositions =
            sizeof (Geometry::NormalCurveVertex) * numNormalCurveSegments;

    const uint32_t normalCurvesByteSize = sizeof(Geometry::NormalCurve)
            * numNormalCurves;

    int totalSize_kb = (classByteSize + normalCurvesByteSize +
                        positionsSize + filteredPositions) / 1024;

    return totalSize_kb;
}


void MNormalCurvesSelection::setNormalCurveSegments(int numSegments,
                            Geometry::NormalCurveVertex *segments,
                            uint32_t numRestart,
                            uint32_t *restartIdx)
{
    numNormalCurveSegments = numSegments;
    normalCurveSegment = segments;
    numRestartIndex = numRestart;
    restartIndex = restartIdx;
}


GL::MVertexBuffer* MNormalCurvesSelection::getVertexBuffer(QGLWidget *currentGLContext)
{
    const QString vboID = QString("vbo_normal_curves_vertices_%1").arg(getID());

    MGLResourcesManager *glRM = MGLResourcesManager::getInstance();
    glRM->makeCurrent();

    auto vb = dynamic_cast<GL::MFloat15VertexBuffer *>(glRM->getGPUItem(vboID));
    if (vb)
    {
        vb->reallocate(reinterpret_cast<float*>(normalCurveSegment),
                       numNormalCurveSegments);
    }
    else
    {

        auto newVb = new GL::MFloat15VertexBuffer(vboID, numNormalCurveSegments);

        if (glRM->tryStoreGPUItem(newVb))
        {
            newVb->upload(reinterpret_cast<float*>(normalCurveSegment),
                          numNormalCurveSegments, currentGLContext);
        }
        else
        {
            LOG4CPLUS_WARN(mlog, "WARNING: cannot store buffer for 3D normal"
                                 " curve vertices in GPU memory.");
            delete newVb;
        }
    }
    return dynamic_cast<GL::MVertexBuffer *>(glRM->getGPUItem(vboID));
}


void MNormalCurvesSelection::releaseVertexBuffer()
{
    MGLResourcesManager::getInstance()->releaseGPUItem(getID());
}


GL::MIndexBuffer* MNormalCurvesSelection::getIndexBuffer(QGLWidget *currentGLContext)
{
    const QString iboID = QString("ibo_normal_curves_mesh_indices_%1").arg(getID());

    MGLResourcesManager *glRM = MGLResourcesManager::getInstance();
    glRM->makeCurrent();

    auto ib = dynamic_cast<GL::MIndexBuffer*>(glRM->getGPUItem(iboID));

    if (ib)
    {
        ib->upload(reinterpret_cast<uint32_t*>(restartIndex),
                   numRestartIndex, currentGLContext);
    }
    else
    {
        auto newVb = new GL::MIndexBuffer(iboID, numRestartIndex);

        if (glRM->tryStoreGPUItem(newVb))
        {
            newVb->upload(reinterpret_cast<uint32_t*>(restartIndex),
                          numRestartIndex, currentGLContext);
        }
        else
        {
            LOG4CPLUS_WARN(mlog, "WARNING: cannot store buffer for 3D normal "
                                 " curve indices in GPU memory.");
            delete newVb;
        }
    }
    return dynamic_cast<GL::MIndexBuffer *>(glRM->getGPUItem(iboID));
}


void MNormalCurvesSelection::releaseIndexBuffer()
{
    MGLResourcesManager::getInstance()->releaseGPUItem(getID());
}


M3DFrontSelection::M3DFrontSelection() : MAbstractDataItem()
{

}

M3DFrontSelection::~M3DFrontSelection()
{
    delete normalCurvesSelection;
    delete triangleMeshSelection;
}


unsigned int M3DFrontSelection::getMemorySize_kb()
{
    const uint32_t classByteSize = sizeof(M3DFrontSelection);

    int totalSize_kb = classByteSize / 1024
            + triangleMeshSelection->getMemorySize_kb()
            + normalCurvesSelection->getMemorySize_kb();
    return totalSize_kb;
}

M3DFrontFilterSource::~M3DFrontFilterSource()
{
    if (inputMeshSource){ delete inputMeshSource; }
}

M3DFrontSelection* M3DFrontSelectionSource::getData(MDataRequest request)
{
    return dynamic_cast<M3DFrontSelection*>(MScheduledDataSource::getData(request));
}


void M3DFrontFilterSource::setInputMeshSource(
        M3DFrontSelectionSource* dataSource)
{
    inputMeshSource = dataSource;
    registerInputSource(inputMeshSource);
    enablePassThrough(inputMeshSource);
}

