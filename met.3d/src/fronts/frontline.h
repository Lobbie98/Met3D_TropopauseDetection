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

#ifndef MET_3D_FRONTLINE_H
#define MET_3D_FRONTLINE_H

// standard library imports

// related third party imports

// local application imports
#include "data/structuredgrid.h"
#include "data/datarequest.h"
#include "data/scheduleddatasource.h"
#include "util/geometry.h"
#include "gxfw/gl/indexbuffer.h"
#include "fronts/frontsurfacemesh.h"

namespace Met3D
{

class MLineSelection : public MAbstractDataItem
{
public:
    /**
      The constructor allocates the data array

     * @param numVertices
     * @param numIndices
     */
    explicit MLineSelection(uint32_t numVertices);

    /** Destructur frees memory field **/
    ~MLineSelection();

    unsigned int getMemorySize_kb() override;

    // vertices
    inline Geometry::FrontLineVertex getVertex(uint index) const
    { return vertices[index]; }

    const Geometry::FrontLineVertex* getVertices() const
    { return vertices; }

    Geometry::FrontLineVertex* getVerticesData()
    { return vertices; }

    uint32_t getNumVertices() { return numVertices; }

    inline void setVertex(uint i, Geometry::FrontLineVertex vertex)
    {vertices[i] = vertex ; }

    GL::MVertexBuffer* getVertexBuffer(QGLWidget *currentGLContext = 0);


    // Triangles
    u_int32_t getIndex(const uint index) const
    { return indices[index]; }

    const u_int32_t* getIndeces() const
    { return indices; }

    u_int32_t* getIndicesData()
    { return indices; }

    uint32_t getNumIndices() const { return numIndices; }

    inline void setIndex(uint i, u_int32_t index)
    {indices[i] = index; }

    void setIndexArray(QVector<u_int32_t> indexArray);

    GL::MIndexBuffer*  getIndexBuffer(QGLWidget *currentGLContext = 0);

    void releaseVertexBuffer();
    void releaseIndexBuffer();

protected:

    Geometry::FrontLineVertex *vertices;
    uint32_t                   numVertices;
    u_int32_t                 *indices;
    uint32_t                   numIndices;

private:

};


class M2DFrontSelection : public MAbstractDataItem
{
public:
    M2DFrontSelection();
    ~M2DFrontSelection();

    unsigned int getMemorySize_kb() override;

    void setLineSelection(MLineSelection *inLine)
    {lineSelection = inLine;}

    void setNormalCurvesSelection(MNormalCurvesSelection *inNormalCurves)
    {normalCurvesSelection = inNormalCurves;}

    MLineSelection* getLineSelection() {return lineSelection; }
    MNormalCurvesSelection* getNormalCurvesSelection() {return normalCurvesSelection; }

protected:
    MLineSelection         *lineSelection;
    MNormalCurvesSelection *normalCurvesSelection;
};


class M2DFrontSelectionSource : public MScheduledDataSource
{
public:
    M2DFrontSelection* getData(MDataRequest request) override;
};


class M2DFrontFilterSource : public M2DFrontSelectionSource
{
public:
    ~M2DFrontFilterSource();
    M2DFrontSelection* produceData(MDataRequest request) override = 0;

    void setInputLineSource(M2DFrontSelectionSource* dataSource);

protected:
    M2DFrontSelectionSource*  inputLineSource;
};


}

#endif // MET_3D_FRONTSURFACEMESH_H
