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

#ifndef MET_3D_FRONTSURFACEMESH_H
#define MET_3D_FRONTSURFACEMESH_H

// standard library imports

// related third party imports

// local application imports
#include "data/structuredgrid.h"
#include "data/datarequest.h"
#include "data/scheduleddatasource.h"
#include "gxfw/gl/indexbuffer.h"
#include "util/geometry.h"

namespace Met3D
{


class MTriangleMeshSelection : public MAbstractDataItem
{
public:
    /**
      The constructor allocates the data array

     * @param numVertices
     * @param numTriangles
     */
    explicit MTriangleMeshSelection(uint32_t numVertices,
                                    uint32_t numTriangles);

    /** Destructor frees memory fields. */
    ~MTriangleMeshSelection();

    unsigned int getMemorySize_kb() override;

    // vertices
    inline Geometry::FrontMeshVertex getVertex(uint index) const
    { return vertices[index]; }

    const Geometry::FrontMeshVertex* getVertices() const
    { return vertices; }

    Geometry::FrontMeshVertex* getVerticesData()
    { return vertices; }

    uint32_t getNumVertices() { return numVertices; }

    inline void setVertex(uint i, Geometry::FrontMeshVertex vertex)
    {vertices[i] = vertex ; }

    GL::MVertexBuffer* getVertexBuffer(QGLWidget *currentGLContext = 0);

    // Triangles
    const Geometry::MTriangle getTriangle(const uint index) const
    { return triangles[index]; }

    const Geometry::MTriangle* getTriangles() const
    { return triangles; }

    Geometry::MTriangle* getTrianglesData()
    { return triangles; }

    uint32_t getNumTriangles() const { return numTriangles; }

    inline void setTriangle(uint i, Geometry::MTriangle triangle)
    {triangles[i] = triangle; }

    uint32_t getNumTriangleIndices() const { return numTriangles * 3; }

    GL::MIndexBuffer*  getIndexBuffer(QGLWidget *currentGLContext = 0);

    void releaseVertexBuffer();
    void releaseIndexBuffer();

protected:
    Geometry::FrontMeshVertex *vertices;
    uint32_t                  numVertices;

    Geometry::MTriangle *triangles;
    uint32_t            numTriangles;

private:

};


class MNormalCurvesSelection : public MAbstractDataItem
{
public:
    /**
      The constructor allocates the data array

     * @param numNormalCurves
     */
    explicit MNormalCurvesSelection(uint32_t numNormalCurves);

    /** Destructor frees memory fields. */
    ~MNormalCurvesSelection();

    unsigned int getMemorySize_kb() override;


    // normal curves
    Geometry::NormalCurve getNormalCurve(const uint index) const
    { return normalCurves[index]; }

    const Geometry::NormalCurve* getNormalCurves() const
    { return normalCurves; }

    Geometry::NormalCurve* getNormalCurvesData()
    { return normalCurves; }

    uint32_t getNumNormalCurves() const { return numNormalCurves; }

    void setNormalCurve(uint i, Geometry::NormalCurve normalCurve)
    {normalCurves[i] = normalCurve ; }

    void setNormalCurveSegments(int numSegments,
                                Geometry::NormalCurveVertex *segments,
                                u_int32_t numRestart,
                                uint32_t *restartIdx);

    uint32_t getNumNormalCurveSegments() {return numNormalCurveSegments; }
    Geometry::NormalCurveVertex* getNormalCurveSegments() {return normalCurveSegment; }
    Geometry::NormalCurveVertex getNormalCurveSegment(int i) {return normalCurveSegment[i]; }
    uint32_t getNumRestartIndex() {return numRestartIndex; }
    uint32_t* getRestartIndex() {return restartIndex; }

    GL::MVertexBuffer*  getVertexBuffer(QGLWidget *currentGLContext = 0);
    GL::MIndexBuffer*  getIndexBuffer(QGLWidget *currentGLContext = 0);

    void releaseVertexBuffer();
    void releaseIndexBuffer();

protected:

    Geometry::NormalCurve *normalCurves;
    uint32_t               numNormalCurves;

    Geometry::NormalCurveVertex *normalCurveSegment;
    uint32_t numNormalCurveSegments;
    uint32_t numRestartIndex;
    uint32_t *restartIndex;

private:

};


class M3DFrontSelection : public MAbstractDataItem
{
public:
    M3DFrontSelection();
    ~M3DFrontSelection();

    unsigned int getMemorySize_kb() override;

    void setTriangleMeshSelection(MTriangleMeshSelection *inTriangleMesh)
    {triangleMeshSelection = inTriangleMesh;}

    void setNormalCurvesSelection(MNormalCurvesSelection *inNormalCurves)
    {normalCurvesSelection = inNormalCurves;}

    MTriangleMeshSelection* getTriangleMeshSelection() {return triangleMeshSelection; }
    MNormalCurvesSelection* getNormalCurvesSelection() {return normalCurvesSelection; }

protected:
    MTriangleMeshSelection *triangleMeshSelection;
    MNormalCurvesSelection *normalCurvesSelection;
};


class M3DFrontSelectionSource : public MScheduledDataSource
{
public:
    M3DFrontSelection* getData(MDataRequest request) override;
};


class M3DFrontFilterSource : public M3DFrontSelectionSource
{
public:
    ~M3DFrontFilterSource();
    M3DFrontSelection* produceData(MDataRequest request) override = 0;

    void setInputMeshSource(M3DFrontSelectionSource* dataSource);

protected:
    M3DFrontSelectionSource*  inputMeshSource;
};


}

#endif // MET_3D_FRONTSURFACEMESH_H
