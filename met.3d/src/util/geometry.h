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

#ifndef MET_3D_GEOMETRY_H
#define MET_3D_GEOMETRY_H

// related third party imports
#include <QVector>
#include <QVector3D>
#include <QVector4D>

namespace Met3D
{
namespace Geometry
{

struct MTriangle
{
    uint32_t indices[3];
};


struct MVoxel
{
    uint32_t i, j, k;
    float values[8];
    QVector3D positions[8];
    QVector3D normals[8];
};

//Julian M: Tropopausen Vertex.
struct TropopauseMeshVertex
{
    QVector3D position;
    QVector3D normal;
    float firstDeriv;
};

struct FrontMeshVertex
{
    QVector3D position;
    QVector3D normal; //
    QVector3D nCEnd;
    float     tfp;
    float     abz;
    float     strength;
    float     type; //warm or cold
    float     breadth;
    float     slope;
};

struct FrontLineVertex
{
    QVector3D position;
    QVector3D nCEnd;
    float     tfp;
    float     abz;
    float     strength;
    float     type; //warm or cold
    float     breadth;
};


struct NormalCurve
{
    QVector<QVector3D> positions;
    float              tfp;
    float              abz;
    float              strength;
    float              type;
    float              breadth;
};


struct NormalCurveVertex
{
    QVector3D   position;
    QVector3D   start;
    QVector3D   end;
    float       tfp;
    float       abz;
    float       strength;
    float       type;
    float       breadth;
};

}
}

#endif //MET_3D_GEOMETRY_H
