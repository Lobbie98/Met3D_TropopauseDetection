/******************************************************************************
**
**  This file is part of Met.3D -- a research environment for the
**  three-dimensional visual exploration of numerical ensemble weather
**  prediction data.
**
**  Copyright 2018 Marc Rautenhaus
**  Copyright 2018 Michael Kern
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

#ifndef MMARCHINGCUBES_H
#define MMARCHINGCUBES_H

// standard library imports
#include <QtCore>
#include <array>

// related third party imports

// local imports
#include "data/structuredgrid.h"
#include "geometry.h"

namespace Met3D
{

class MMarchingCubes
{
public:
    explicit MMarchingCubes(MStructuredGrid* grid, MStructuredGrid* zGrid = nullptr);
//    explicit MMarchingCubes();

    void computeMeshOnCPU(const float isovalue);

//    QVector<QVector3D>* getIntersectionPoints() { return &intersectionPoints; }
//    QVector<QVector3D>* getIntersectionNormals() { return &intersectionNormals; }
//    QVector<MTriangle>* getTriangles() { return &triangles; }
    QVector<QVector3D>* getFlattenInterPoints() { return &flattenPoints; }
    QVector<QVector3D>* getFlattenInterNormals() { return &flattenNormals; }
    QVector<QVector3D>* getFlattenInterNormalsZ() { return &flattenNormalsZ; }
    QVector<Geometry::MTriangle>* getFlattenTriangles() { return &flattenTriangles; }

private:

    void precompute();

    void computeVoxelIndices(const double isovalue);
    void computeIntersectionPoints(const double isovalue);
    void generateTriangles();

    //void generateMesh();
    //void flatten();

    void initializeVoxel(const uint32_t k, const uint32_t j, const uint32_t i,
                         Geometry::MVoxel& voxel) const;

    MStructuredGrid*    inputGrid;
    MStructuredGrid*    heightGrid;
    uint32_t            nx;
    uint32_t            ny;
    uint32_t            nz;
    uint32_t            numEdges;
    uint32_t            numVoxels;
    QVector<QVector3D>  gridPoints;
    QVector<QVector3D>  normals;
    QVector<QVector3D>  normalsZ;
    QVector<uint8_t>    voxelIndices;
    QVector<uint8_t>    voxelNumTriangles;
    QVector<uint32_t>   voxelMinIndexTriangle;
    QVector<uint32_t>   intersectionIndices;
    QVector<QVector3D>  flattenPoints;
    QVector<QVector3D>  flattenNormals;
    QVector<QVector3D>  flattenNormalsZ;
    QVector<Geometry::MTriangle>  flattenTriangles;
    uint32_t            atomicEdgeCount;
    uint32_t            atomicTriangleCount;

    void lerpAtVoxel(const double isovalue,
                     const uint32_t k,
                     const uint32_t j,
                     const uint32_t i,
                     const uint8_t edge);

    inline QVector3D vec3Lerp(const float isovalue,
                              const float valP, const float valN,
                              const QVector3D& vP, const QVector3D& vN) const;

    inline void getEdgeVertices(const uint8_t edge, uint8_t& vP, uint8_t& vN) const;

    inline float getVoxelValue(const uint32_t k, const uint32_t j,
                               const uint32_t i, const uint8_t v) const;

    inline QVector3D getPosition(const uint32_t k, const uint32_t j,
                                 const uint32_t i, const uint8_t v) const;

    inline QVector3D getNormal(const uint32_t k, const uint32_t j,
                               const uint32_t i, const uint8_t v) const;

    inline QVector3D getNormalZ(const uint32_t k, const uint32_t j,
                                const uint32_t i, const uint8_t v) const;

    inline uint32_t getPointIndex(const uint32_t k, const uint32_t j,
                                  const uint32_t i, const uint32_t edgeIndex);
};

}

#endif //MMARCHINGCUBES_H
