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

#ifndef MET_3D_TRIANGLEMESHSELECTIONTROPOPAUSE_H
#define MET_3D_TRIANGLEMESHSELECTIONTROPOPAUSE_H

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


class MTropopauseTriangleMeshSelection : public MAbstractDataItem
{
public:
    /**
      The constructor allocates the data array

     * @param numVertices
     * @param numTriangles
     */
    explicit MTropopauseTriangleMeshSelection(uint32_t numVertices,
                                    uint32_t numTriangles);

    /** Destructor frees memory fields. */
    ~MTropopauseTriangleMeshSelection();

    unsigned int getMemorySize_kb() override;

    // vertices
    inline Geometry::TropopauseMeshVertex getVertex(uint index) const
    { return vertices[index]; }

    const Geometry::TropopauseMeshVertex* getVertices() const
    { return vertices; }

    Geometry::TropopauseMeshVertex* getVerticesData()
    { return vertices; }

    uint32_t getNumVertices() { return numVertices; }

    inline void setVertex(uint i, Geometry::TropopauseMeshVertex vertex)
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
    Geometry::TropopauseMeshVertex *vertices;
    uint32_t                  numVertices;

    Geometry::MTriangle *triangles;
    uint32_t            numTriangles;

private:

};


}

#endif // MET_3D_TRIANGLEMESHSELECTIONTROPOPAUSE_H
