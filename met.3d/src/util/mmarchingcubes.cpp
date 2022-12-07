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
#include "mmarchingcubes.h"

#include <omp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <chrono>

// related third party imports
#include <log4cplus/loggingmacros.h>

#define MEASURE_CPU_TIME
#define COMPUTE_PARALLEL

using namespace Met3D;

#define qNaN (std::numeric_limits<float>::quiet_NaN())
#define iNaN (std::numeric_limits<unsigned int>::max())

MMarchingCubes::MMarchingCubes(MStructuredGrid* grid, MStructuredGrid* zGrid)
    : inputGrid(grid),
      heightGrid(zGrid),
      nx(grid->getNumLons() - 1),
      ny(grid->getNumLats() - 1),
      nz(grid->getNumLevels() - 1),
      numEdges(3 * grid->getNumValues()
               - grid->getNumLons() * grid->getNumLevels()
               - grid->getNumLons() * grid->getNumLats()
               - grid->getNumLats() * grid->getNumLevels()),
      numVoxels(nx * ny * nz),
      gridPoints(grid->getNumValues()),
      normals(grid->getNumValues()),
      normalsZ((zGrid) ? grid->getNumValues() : 0),
      voxelIndices(numVoxels, 0),
      voxelNumTriangles(numVoxels, 0),
      voxelMinIndexTriangle(numVoxels, 0),
      intersectionIndices(numEdges, iNaN),
      flattenPoints(numEdges, QVector3D(qNaN, qNaN, qNaN)),
      flattenNormals(numEdges, QVector3D(qNaN, qNaN, qNaN)),
      flattenNormalsZ((zGrid) ? numEdges : 0, QVector3D(qNaN, qNaN, qNaN)),
      atomicEdgeCount(0),
      atomicTriangleCount(0)
{
    precompute();
}


//! TODO this matrix is symmetric! Reduce the size!
int edgeTable[256] = {
        0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
        0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
        0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
        0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
        0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
        0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
        0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
        0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
        0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
        0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
        0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
        0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
        0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
        0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
        0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
        0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
        0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
        0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
        0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
        0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
        0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
        0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
        0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
        0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
        0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
        0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
        0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
        0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
        0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
        0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
        0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
        0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };

//! TODO this matrix is symmetric! Reduce the size!
int triTable[256][16] = {
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
        {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
        {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
        {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
        {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
        {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
        {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
        {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
        {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
        {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
        {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
        {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
        {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
        {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
        {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
        {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
        {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
        {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
        {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
        {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
        {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
        {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
        {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
        {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
        {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
        {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
        {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
        {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
        {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
        {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
        {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
        {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
        {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
        {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
        {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
        {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
        {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
        {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
        {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
        {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
        {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
        {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
        {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
        {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
        {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
        {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
        {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
        {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
        {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
        {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
        {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
        {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
        {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
        {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
        {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
        {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
        {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
        {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
        {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
        {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
        {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
        {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
        {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
        {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
        {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
        {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
        {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
        {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
        {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
        {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
        {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
        {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
        {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
        {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
        {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
        {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
        {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
        {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
        {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
        {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
        {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
        {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
        {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
        {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
        {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
        {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
        {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
        {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
        {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
        {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
        {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
        {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
        {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
        {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
        {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
        {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
        {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
        {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
        {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
        {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
        {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
        {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
        {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
        {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
        {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
        {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
        {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
        {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
        {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
        {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
        {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
        {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
        {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
        {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
        {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
        {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
        {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
        {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
        {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
        {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
        {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
        {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
        {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
        {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
        {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
        {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
        {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
        {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
        {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
        {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
        {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
        {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
        {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
        {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
        {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
        {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
        {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
        {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
        {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
        {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
        {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
        {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
        {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
        {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
        {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
        {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
        {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
        {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
        {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
        {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
        {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
        {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
        {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
        {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
        {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
        {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
        {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
        {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
        {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
        {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
        {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
        {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
        {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
        {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
        {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
        {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
        {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
        {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
        {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1} };


void MMarchingCubes::computeMeshOnCPU(const float isovalue)
{
    // 1)
    computeVoxelIndices(isovalue);

    // 2)
    computeIntersectionPoints(isovalue);

    // 3)
    generateTriangles();

    // 4) not necessary
    //flatten();
}


void MMarchingCubes::precompute()
{
    printf("\t -> compute grid points and normals...");
#ifdef MEASURE_CPU_TIME
    auto start = std::chrono::system_clock::now();
#endif

    const float DELTA_LAT_M = 1.112E5; // ~111.2km

    for (uint32_t k = 0; k <= nz; ++k)
    {
        auto kN = std::min(int(k + 1), int(nz));
        auto kP = std::max(int(k - 1), 0);

#ifdef COMPUTE_PARALLEL
#pragma omp parallel for collapse(2)
#endif
        for (uint32_t j = 0; j <= ny; ++j)
        {
            for (uint32_t i = 0; i <= nx; ++i)
            {
                const float lon = inputGrid->getLons()[i];
                const float lat = inputGrid->getLats()[j];
                const float DELTA_LON_M = std::cos(lat / 180.0f * M_PI) * DELTA_LAT_M;

                const float pressure = inputGrid->getPressure(k, j, i);

                const uint32_t pIndex = INDEX3zyx(k, j, i, ny + 1, nx + 1);

                gridPoints[pIndex] = QVector3D(lon, lat, pressure);

                auto jN = std::min(int(j + 1), int(ny));
                auto jP = std::max(int(j - 1), 0);

                auto iN = std::min(int(i + 1), int(nx));
                auto iP = std::max(int(i - 1), 0);

                const float dx = inputGrid->getLons()[1] - inputGrid->getLons()[0];
                const float dy = inputGrid->getLats()[1] - inputGrid->getLats()[0];

                const float deltaDegreeX = (iN - iP) * dx;
                const float deltaDegreeY = (jN - jP) * dy;

                const float deltahPa = inputGrid->getPressure(kN, j, i) - inputGrid->getPressure(kP, j, i);
                const float deltaLon = DELTA_LON_M * deltaDegreeX;
                const float deltaLat = DELTA_LAT_M * deltaDegreeY;

                float diffLon = inputGrid->getValue(k, j, iN) - inputGrid->getValue(k, j, iP);
                float diffLat = inputGrid->getValue(k, jN, i) - inputGrid->getValue(k, jP, i);
                float diffZ = inputGrid->getValue(kN, j, i) - inputGrid->getValue(kP, j, i);

                QVector3D normal(diffLon / deltaDegreeX,
                                 diffLat / deltaDegreeY,
                                 diffZ / deltahPa);

                normals[pIndex] = -normal;

                if (heightGrid)
                {
                    const float deltaZ = heightGrid->getValue(kN, j, i) - heightGrid->getValue(kP, j, i);

                    QVector3D normalZ(diffLon / deltaLon,
                                      diffLat / deltaLat,
                                      diffZ / deltaZ);

                    normalsZ[pIndex] = normalZ;

                }


            }
        }
    }

#ifdef MEASURE_CPU_TIME
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    LOG4CPLUS_DEBUG(mlog, " done in " << elapsed.count() << "ms.\n");
#endif
}



void MMarchingCubes::computeVoxelIndices(const double isovalue)
{
    printf("\t -> compute voxel edge indices...");
#ifdef MEASURE_CPU_TIME
    auto start = std::chrono::system_clock::now();
#endif

    atomicTriangleCount = 0;

    for (uint32_t k = 0; k < nz; ++k)
    {
#ifdef COMPUTE_PARALLEL
#pragma omp parallel for collapse(2)
#endif
        for (uint32_t j = 0; j < ny; ++j)
        {
            for (uint32_t i = 0; i < nx; ++i)
            {
                uint8_t voxelIndex = 0;

                voxelIndex |= (inputGrid->getValue(k, j, i) < isovalue)              ? 0x1 : 0;
                voxelIndex |= (inputGrid->getValue(k, j, i + 1) < isovalue)          ? 0x2 : 0;
                voxelIndex |= (inputGrid->getValue(k, j + 1, i + 1) < isovalue)      ? 0x4 : 0;
                voxelIndex |= (inputGrid->getValue(k, j + 1, i) < isovalue)          ? 0x8 : 0;
                voxelIndex |= (inputGrid->getValue(k + 1, j, i) < isovalue)          ? 0x10 : 0;
                voxelIndex |= (inputGrid->getValue(k + 1, j, i + 1) < isovalue)      ? 0x20 : 0;
                voxelIndex |= (inputGrid->getValue(k + 1, j + 1, i + 1) < isovalue)  ? 0x40 : 0;
                voxelIndex |= (inputGrid->getValue(k + 1, j + 1, i) < isovalue)      ? 0x80 : 0;

                const uint32_t pIndex = INDEX3zyx(k, j, i, ny, nx);
                voxelIndices[pIndex] = voxelIndex;

                if (voxelIndex == 0 || voxelIndex == 0xFF) { continue; }

                uint8_t numVoxelTris = 0;

                for (uint32_t ii = 0; triTable[voxelIndex][ii] != -1; ii += 3)
                {
                    numVoxelTris++;
                }

                uint32_t triIndex = 0;
#ifdef COMPUTE_PARALLEL
#pragma omp atomic capture
#endif
                {
                    triIndex = atomicTriangleCount;
                    atomicTriangleCount += numVoxelTris;
                }
                voxelNumTriangles[pIndex] = numVoxelTris;
                voxelMinIndexTriangle[pIndex] = triIndex;
            }
        }
    }


    flattenTriangles.resize(atomicTriangleCount);
//    flattenPoints.resize(numEdges);
//    flattenNormals.resize(numEdges);

#ifdef MEASURE_CPU_TIME
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    LOG4CPLUS_DEBUG(mlog, " done in " << elapsed.count() << "ms.\n");
#endif
}

void MMarchingCubes::computeIntersectionPoints(const double isovalue)
{
    printf("\t -> compute intersection points at edges...");
#ifdef MEASURE_CPU_TIME
    auto start = std::chrono::system_clock::now();
#endif

    atomicEdgeCount = 0;

    for (uint32_t k = 0; k < nz; ++k)
    {
#ifdef COMPUTE_PARALLEL
#pragma omp parallel for collapse(2)
#endif
        for (uint32_t j = 0; j < ny; ++j)
        {
            for (uint32_t i = 0; i < nx; ++i)
            {
                const uint32_t index = INDEX3zyx(k, j, i, ny, nx);
                const uint8_t voxelIndex = voxelIndices[index];

                if (voxelIndex == 0 || voxelIndex == 0xFF) { continue; }

                const int edgeIndex = edgeTable[voxelIndex];

                // last slice boundary --> collect upper edges
                if (k == nz - 1)
                {
                    // Loop over all edges in x-direction
                    // ___ ___ ___ ___
                    if (j == ny - 1 && edgeIndex & 0x40)
                    {
                        lerpAtVoxel(isovalue, k, j, i, 6);
                    }
                    if (edgeIndex & 0x10)
                    {
                        lerpAtVoxel(isovalue, k, j, i, 4);
                    }

                    // Loop over all edges in y-direction
                    // |   |   |   |
                    if (i == nx - 1 && edgeIndex & 0x20)
                    {
                        lerpAtVoxel(isovalue, k, j, i, 5);
                    }
                    if (edgeIndex & 0x80)
                    {
                        lerpAtVoxel(isovalue, k, j, i, 7);
                    }
                }
                // collect lower edges of voxels

                // Loop over all edges in x-direction
                // ___ ___ ___ ___
                if (j == ny - 1 && edgeIndex & 0x4)
                {
                    lerpAtVoxel(isovalue, k, j, i, 2);
                }
                if (edgeIndex & 0x1)
                {
                    lerpAtVoxel(isovalue, k, j, i, 0);
                }

                // Loop over all edges in y-direction
                // |   |   |   |
                if (i == nx - 1 && edgeIndex & 0x2)
                {
                    lerpAtVoxel(isovalue, k, j, i, 1);
                }
                if (edgeIndex & 0x8)
                {
                    lerpAtVoxel(isovalue, k, j, i, 3);
                }

                // Loop over all edges in z-direction
                // |/_   _|/_   _|/_   _|/
                if (edgeIndex & 0x100 )
                {
                    lerpAtVoxel(isovalue, k, j, i, 8);
                }
                if (edgeIndex & 0x200 && i == nx - 1)
                {
                    lerpAtVoxel(isovalue, k, j, i, 9);
                }
                if (edgeIndex & 0x400 && j == ny - 1 && i == nx - 1)
                {
                    lerpAtVoxel(isovalue, k, j, i, 10);
                }
                if (edgeIndex & 0x800 && j == ny - 1)
                {
                    lerpAtVoxel(isovalue, k, j, i, 11);
                }

            }
        }
    }

    flattenPoints.resize(atomicEdgeCount);
    flattenNormals.resize(atomicEdgeCount);
    if (heightGrid) { flattenNormalsZ.resize(atomicEdgeCount); }

#ifdef MEASURE_CPU_TIME
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    LOG4CPLUS_DEBUG(mlog, " done in " << elapsed.count() << "ms.\n");
#endif
}


void MMarchingCubes::generateTriangles()
{
    printf("\t -> generate triangles...");
    auto start = std::chrono::system_clock::now();

    for (uint32_t k = 0; k < nz; ++k)
    {
#ifdef COMPUTE_PARALLEL
#pragma omp parallel for collapse(2)
#endif
        for (uint32_t j = 0; j < ny; ++j)
        {
            for (uint32_t i = 0; i < nx; ++i)
            {
                const uint32_t index = INDEX3zyx(k, j, i, ny, nx);
                const uint8_t voxelIndex = voxelIndices[index];

                if (voxelIndex == 0 || voxelIndex == 0xFF) { continue; }

                uint32_t minIndexTriangle = voxelMinIndexTriangle[index];
                uint32_t numTriangles = voxelNumTriangles[index];

                bool degeneratedTriangle = false;

                //for (uint32_t ii = 0; triTable[voxelIndex][ii] != -1; ii += 3)
                for (uint32_t kk = 0; kk < numTriangles; ++kk)
                {
                    uint32_t ii = kk * 3;
                    Geometry::MTriangle triangle = { 0, 0, 0 };

                    for (uint32_t jj = 0; jj < 3; ++jj)
                    {
                        const auto edgeIndex = triTable[voxelIndex][ii + jj];
                        const uint32_t pIndex = getPointIndex(k, j, i, edgeIndex);
                        const uint32_t fIndex = intersectionIndices[pIndex];

                        degeneratedTriangle = (fIndex == iNaN);
                        if (degeneratedTriangle) { break; }

                        triangle.indices[jj] = intersectionIndices[pIndex];
                    }

                    if (!degeneratedTriangle) { flattenTriangles[minIndexTriangle + kk] = triangle; }
                    else { flattenTriangles[minIndexTriangle + kk] = { 0, 0, 0 }; }
                }
            }
        }
    }

#ifdef MEASURE_CPU_TIME
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    LOG4CPLUS_DEBUG(mlog, " done in " << elapsed.count() << "ms.\n");
#endif
}


void MMarchingCubes::initializeVoxel(const uint32_t k, const uint32_t j,
                                     const uint32_t i,
                                     Geometry::MVoxel& voxel) const
{
    voxel.i = i; voxel.j = j; voxel.k = k;

    for (uint8_t cc = 0; cc < 8; ++cc)
    {
        voxel.values[cc] = getVoxelValue(k, j, i, cc);
        voxel.positions[cc] = getPosition(k, j, i, cc);
        voxel.normals[cc] = getNormal(k, j, i, cc);
    }
}


void MMarchingCubes::lerpAtVoxel(const double isovalue,
                                 const uint32_t k,
                                 const uint32_t j,
                                 const uint32_t i,
                                 const uint8_t edge)
{
    uint8_t vP = 0;
    uint8_t vN = 0;
    getEdgeVertices(edge, vP, vN);

    const float valP = getVoxelValue(k, j, i, vP);
    const QVector3D posP = getPosition(k, j, i, vP);
    const QVector3D normalP = getNormal(k, j, i, vP);

    const float valN = getVoxelValue(k, j, i, vN);
    const QVector3D posN = getPosition(k, j, i, vN);
    const QVector3D normalN = getNormal(k, j, i, vN);

    const uint32_t pIndex = getPointIndex(k, j, i, edge);

    uint32_t fIndex = 0;

#ifdef COMPUTE_PARALLEL
#pragma omp atomic capture
#endif
    {
        fIndex = atomicEdgeCount;
        atomicEdgeCount++;
    }

    if (std::isnan(valN) || std::isnan(valP))
    {
        flattenPoints[fIndex] = QVector3D(qNaN, qNaN, qNaN);
        flattenNormals[fIndex] = QVector3D(qNaN, qNaN, qNaN);;

        if (heightGrid)
        {
            flattenNormalsZ[fIndex] = QVector3D(qNaN, qNaN, qNaN);;
        }

        return;
    }

    intersectionIndices[pIndex] = fIndex;
    flattenPoints[fIndex] = vec3Lerp(isovalue, valP, valN, posP, posN);
    flattenNormals[fIndex] = vec3Lerp(isovalue, valP, valN, normalP, normalN);

    if (heightGrid)
    {
        const QVector3D normalZP = getNormalZ(k, j, i, vP);
        const QVector3D normalZN = getNormalZ(k, j, i, vN);

        flattenNormalsZ[fIndex] = vec3Lerp(isovalue, valP, valN, normalZP, normalZN);
    }
}


QVector3D MMarchingCubes::vec3Lerp(const float isovalue,
                                   const float valP, const float valN,
                                   const QVector3D& vP, const QVector3D& vN) const
{
    const double PRECISION = 1E-30;

    if (std::abs(isovalue - valP) < PRECISION) { return vP; }
    if (std::abs(isovalue - valN) < PRECISION) { return vN; }
    if (std::abs(valP - valN) < PRECISION)     { return vP; }
//    if (std::isnan(valP)) { return vN; }
//    if (std::isnan(valN)) { return vP; }

    const double mu = (isovalue - valP) / (valN - valP);

    return vP + mu * (vN - vP);
}


void MMarchingCubes::getEdgeVertices(const uint8_t edge,
                                     uint8_t& vP, uint8_t& vN) const
{
    if (edge == 0) { vP = 0; vN = 1; }
    if (edge == 1) { vP = 1; vN = 2; }
    if (edge == 2) { vP = 2; vN = 3; }
    if (edge == 3) { vP = 3; vN = 0; }
    if (edge == 4) { vP = 4; vN = 5; }
    if (edge == 5) { vP = 5; vN = 6; }
    if (edge == 6) { vP = 6; vN = 7; }
    if (edge == 7) { vP = 7; vN = 4; }
    if (edge == 8) { vP = 0; vN = 4; }
    if (edge == 9) { vP = 1; vN = 5; }
    if (edge == 10) { vP = 2; vN = 6; }
    if (edge == 11) { vP = 3; vN = 7; }
}


float MMarchingCubes::getVoxelValue(const uint32_t k, const uint32_t j,
                                    const uint32_t i, const uint8_t v) const
{
    if (v == 0) { return inputGrid->getValue(k, j, i); }
    if (v == 1) { return inputGrid->getValue(k, j, i + 1); }
    if (v == 2) { return inputGrid->getValue(k, j + 1, i + 1); }
    if (v == 3) { return inputGrid->getValue(k, j + 1, i); }
    if (v == 4) { return inputGrid->getValue(k + 1, j, i); }
    if (v == 5) { return inputGrid->getValue(k + 1, j, i + 1); }
    if (v == 6) { return inputGrid->getValue(k + 1, j + 1, i + 1); }
    if (v == 7) { return inputGrid->getValue(k + 1, j + 1, i); }
    return 0;
}


QVector3D MMarchingCubes::getPosition(const uint32_t k, const uint32_t j,
                                      const uint32_t i, const uint8_t v) const
{
    if (v == 0) { return gridPoints[INDEX3zyx(k, j, i, ny + 1, nx + 1)]; }
    if (v == 1) { return gridPoints[INDEX3zyx(k, j, i + 1, ny + 1, nx + 1)]; }
    if (v == 2) { return gridPoints[INDEX3zyx(k, j + 1, i + 1, ny + 1, nx + 1)]; }
    if (v == 3) { return gridPoints[INDEX3zyx(k, j + 1, i, ny + 1, nx + 1)]; }
    if (v == 4) { return gridPoints[INDEX3zyx(k + 1, j, i, ny + 1, nx + 1)]; }
    if (v == 5) { return gridPoints[INDEX3zyx(k + 1, j, i + 1, ny + 1, nx + 1)]; }
    if (v == 6) { return gridPoints[INDEX3zyx(k + 1, j + 1, i + 1, ny + 1, nx + 1)]; }
    if (v == 7) { return gridPoints[INDEX3zyx(k + 1, j + 1, i, ny + 1, nx + 1)]; }

    return QVector3D(0, 0, 0);
}


QVector3D MMarchingCubes::getNormal(const uint32_t k, const uint32_t j,
                                      const uint32_t i, const uint8_t v) const
{
    if (v == 0) { return normals[INDEX3zyx(k, j, i, ny + 1, nx + 1)]; }
    if (v == 1) { return normals[INDEX3zyx(k, j, i + 1, ny + 1, nx + 1)]; }
    if (v == 2) { return normals[INDEX3zyx(k, j + 1, i + 1, ny + 1, nx + 1)]; }
    if (v == 3) { return normals[INDEX3zyx(k, j + 1, i, ny + 1, nx + 1)]; }
    if (v == 4) { return normals[INDEX3zyx(k + 1, j, i, ny + 1, nx + 1)]; }
    if (v == 5) { return normals[INDEX3zyx(k + 1, j, i + 1, ny + 1, nx + 1)]; }
    if (v == 6) { return normals[INDEX3zyx(k + 1, j + 1, i + 1, ny + 1, nx + 1)]; }
    if (v == 7) { return normals[INDEX3zyx(k + 1, j + 1, i, ny + 1, nx + 1)]; }

    return QVector3D(0, 0, 0);
}


QVector3D MMarchingCubes::getNormalZ(const uint32_t k, const uint32_t j,
                                     const uint32_t i, const uint8_t v) const
{
    if (v == 0) { return normalsZ[INDEX3zyx(k, j, i, ny + 1, nx + 1)]; }
    if (v == 1) { return normalsZ[INDEX3zyx(k, j, i + 1, ny + 1, nx + 1)]; }
    if (v == 2) { return normalsZ[INDEX3zyx(k, j + 1, i + 1, ny + 1, nx + 1)]; }
    if (v == 3) { return normalsZ[INDEX3zyx(k, j + 1, i, ny + 1, nx + 1)]; }
    if (v == 4) { return normalsZ[INDEX3zyx(k + 1, j, i, ny + 1, nx + 1)]; }
    if (v == 5) { return normalsZ[INDEX3zyx(k + 1, j, i + 1, ny + 1, nx + 1)]; }
    if (v == 6) { return normalsZ[INDEX3zyx(k + 1, j + 1, i + 1, ny + 1, nx + 1)]; }
    if (v == 7) { return normalsZ[INDEX3zyx(k + 1, j + 1, i, ny + 1, nx + 1)]; }

    return QVector3D(0, 0, 0);
}


uint32_t MMarchingCubes::getPointIndex( const uint32_t k, const uint32_t j,
                                        const uint32_t i, const uint32_t edgeIndex)
{
    const uint32_t dimX = inputGrid->getNumLons();
    const uint32_t dimY = inputGrid->getNumLats();
    const uint32_t dimZ = inputGrid->getNumLevels();

    const uint32_t numXEdges = dimZ * nx * dimY;
    const uint32_t numYEdges = dimZ * ny * dimX;

    // x-edges
    if (edgeIndex == 0) { return INDEX3zyx(k, j, i, dimY, nx); }
    if (j == ny - 1 && edgeIndex == 2) { return INDEX3zyx(k, j + 1, i, dimY, nx); }
    if (edgeIndex == 2) { return getPointIndex(k, j + 1, i, 0); }
    // y-edges
    if (edgeIndex == 3) { return INDEX3zyx(k, j, i, ny, dimX) + numXEdges; }
    if (i == nx - 1 && edgeIndex == 1) { return INDEX3zyx(k, j, i + 1, ny, dimX) + numXEdges; }
    if (edgeIndex == 1) { return getPointIndex(k, j, i + 1, 3); }

    if (k == nz - 1)
    {
        // x-edges
        if (edgeIndex == 4) { return INDEX3zyx(k + 1, j, i, dimY, nx); }
        if (j == ny - 1 && edgeIndex == 6) { return INDEX3zyx(k + 1, j + 1, i, dimY, nx); }
        if (edgeIndex == 6) { return getPointIndex(k, j + 1, i, 0); }
        // y-edges
        if (edgeIndex == 7) { return INDEX3zyx(k + 1, j, i, ny, dimX) + numXEdges; }
        if (i == nx - 1 && edgeIndex == 5) { return INDEX3zyx(k + 1, j, i + 1, ny, dimX) + numXEdges; }
        if (edgeIndex == 5) { return getPointIndex(k, j, i + 1, 7); }
    }
    else
    {
        // x-edges
        if (edgeIndex == 4) { return getPointIndex(k + 1, j, i, 0); }
        if (edgeIndex == 6) { return getPointIndex(k + 1, j, i, 2); }
        // y-edges
        if (edgeIndex == 5) { return getPointIndex(k + 1, j, i, 1); }
        if (edgeIndex == 7) { return getPointIndex(k + 1, j, i, 3); }
    }

    const uint32_t numXYEdges = numXEdges + numYEdges;

    // z-edges
    if (edgeIndex == 8) {  return INDEX3zyx(k, j, i, dimY, dimX) + numXYEdges; }
    if (edgeIndex == 9) {  return INDEX3zyx(k, j, i + 1, dimY, dimX) + numXYEdges; }
    if (edgeIndex == 10) { return INDEX3zyx(k, j + 1, i + 1, dimY, dimX) + numXYEdges; }
    if (edgeIndex == 11) { return INDEX3zyx(k, j + 1, i, dimY, dimX) + numXYEdges; }

    return 0;
}
