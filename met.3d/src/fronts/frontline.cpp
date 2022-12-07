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

#include "frontline.h"

#include "gxfw/gl/typedvertexbuffer.h"

using namespace Met3D;

MLineSelection::MLineSelection(uint32_t numVertices)
    : MAbstractDataItem(),
      vertices(new Geometry::FrontLineVertex[numVertices]),
      numVertices(numVertices)
{

}

MLineSelection::~MLineSelection()
{
    delete[] vertices;
    delete[] indices;
}

unsigned int MLineSelection::getMemorySize_kb()
{
    const uint32_t classByteSize = sizeof(MLineSelection);

    const uint32_t vertexByteSize =
            numVertices * sizeof(Geometry::FrontLineVertex);

    const u_int32_t indexByteSize = numIndices * sizeof (u_int32_t);

    return (classByteSize + vertexByteSize + indexByteSize) / 1024;
}


void MLineSelection::setIndexArray(QVector<u_int32_t> indexArray)
{
    numIndices = indexArray.size();
    indices = new u_int32_t[numIndices];
    for (u_int32_t i = 0; i < numIndices; i++)
    {
        indices[i] = indexArray.at(i);
    }
}

GL::MVertexBuffer* MLineSelection::getVertexBuffer(QGLWidget *currentGLContext)
{
    const QString vboID = QString("vbo_line_vertices_%1").arg(getID());

    MGLResourcesManager *glRM = MGLResourcesManager::getInstance();
    glRM->makeCurrent();

    const uint32_t numFloatsPerVertex = sizeof(Geometry::FrontLineVertex)
                                        / sizeof(float);
    const uint32_t numFloats = numVertices * numFloatsPerVertex;

    auto vb = dynamic_cast<GL::MVertexBuffer *>(glRM->getGPUItem(vboID));
    if (vb) { return vb; }

    auto newVb = new GL::MFloatVertexBuffer(vboID, numFloats);

    if (glRM->tryStoreGPUItem(newVb))
    {
        newVb->upload(reinterpret_cast<float*>(vertices),
                      numVertices, currentGLContext);
    }
    else
    {
        LOG4CPLUS_WARN(mlog, "WARNING: cannot store buffer for volume"
                             " bbox in GPU memory.");
        delete newVb;
    }

    return dynamic_cast<GL::MVertexBuffer *>(glRM->getGPUItem(vboID));
}


void MLineSelection::releaseVertexBuffer()
{
    MGLResourcesManager::getInstance()->releaseGPUItem(getID());
}


GL::MIndexBuffer* MLineSelection::getIndexBuffer(QGLWidget *currentGLContext)
{
    const QString iboID = QString("ibo_line_indices_%1").arg(getID());

    MGLResourcesManager *glRM = MGLResourcesManager::getInstance();
    glRM->makeCurrent();

    auto vb = dynamic_cast<GL::MIndexBuffer*>(glRM->getGPUItem(iboID));
    if (vb) { return vb; }

    auto newVb = new GL::MIndexBuffer(iboID, numIndices);

    if (glRM->tryStoreGPUItem(newVb))
    {
        newVb->upload(reinterpret_cast<uint32_t*>(indices),
                      numIndices, currentGLContext);
    }
    else
    {
        LOG4CPLUS_WARN(mlog, "WARNING: cannot store buffer for volume"
                             " bbox in GPU memory.");
        delete newVb;
    }

    return dynamic_cast<GL::MIndexBuffer *>(glRM->getGPUItem(iboID));
}

void MLineSelection::releaseIndexBuffer()
{
    MGLResourcesManager::getInstance()->releaseGPUItem(getID());
}


M2DFrontSelection::M2DFrontSelection() : MAbstractDataItem()
{

}

M2DFrontSelection::~M2DFrontSelection()
{
    delete lineSelection;
    delete normalCurvesSelection;
}


unsigned int M2DFrontSelection::getMemorySize_kb()
{
    const uint32_t classByteSize = sizeof(M3DFrontSelection);

    return classByteSize / 1024
            + lineSelection->getMemorySize_kb()
            + normalCurvesSelection->getMemorySize_kb();
}

M2DFrontFilterSource::~M2DFrontFilterSource()
{
    if (inputLineSource){ delete inputLineSource; }
}

M2DFrontSelection* M2DFrontSelectionSource::getData(MDataRequest request)
{
    return dynamic_cast<M2DFrontSelection*>(MScheduledDataSource::getData(request));
}


void M2DFrontFilterSource::setInputLineSource(
        M2DFrontSelectionSource* dataSource)
{
    inputLineSource = dataSource;
    registerInputSource(inputLineSource);
    enablePassThrough(inputLineSource);
}


