/******************************************************************************
**
**  This file is part of Met.3D -- a research environment for the
**  three-dimensional visual exploration of numerical ensemble weather
**  prediction data.
**
**  Copyright 2015-2018 Marc Rautenhaus
**  Copyright 2015      Michael Kern
**  Copyright 2017-2018 Bianca Tost
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
******************************************************************************/

// Compute shader that computes normal curve line vertices from given
// initial points.

/*****************************************************************************
 ***                             CONSTANTS
 *****************************************************************************/

const int MAX_ISOSURFACES = 10;

const int VALUE_MODE_LENGTH = 0;
const int VALUE_MODE_VALUE = 1;
const int VALUE_MODE_AVG_DIFFERENCE = 2;
const int VALUE_MODE_TOTAL_DIFFERENCE = 3;

// Vertical level type; see structuredgrid.h.
const int SURFACE_2D = 0;
const int PRESSURE_LEVELS_3D = 1;
const int HYBRID_SIGMA_PRESSURE_3D = 2;
const int POTENTIAL_VORTICITY_2D = 3;
const int LOG_PRESSURE_LEVELS_3D = 4;
const int AUXILIARY_PRESSURE_3D = 5;

//const int SHADOWS_OFF = 0;
//const int SHADOWS_MAP = 1;
//const int SHADOWS_VOLUME_AND_RAY = 2;

const float EARTH_RADIUS_KM = 6371;
const float DELTA_LAT_M = 1.112E5;
const float M_PI = 3.1415926535897932384626433832795;

/*****************************************************************************
 ***                             UNIFORMS
 *****************************************************************************/

// Normal curve vertex element.
struct InitPointNC
{
    vec4 position;
};

struct EndPosLengthNC
{
    vec4 position;
};


// Shader buffer storage object of inti points for normal curves.
layout (std430, binding=0) buffer NormalCurveIntegrationBuffer
{
    InitPointNC startVertex[];
};

// Shader buffer storage object of computed normal curve vertices.
layout (std430, binding=1) writeonly buffer NormalCurveLineBuffer
{
    EndPosLengthNC endVertex[];
};


// matrices
// ========
uniform sampler3D dataVolume; // SIGMA_HYBRID_PRESSURE | PRESSURE_LEVELS // fle grid
uniform sampler2D surfacePressure; // SIGMA_HYBRID_PRESSURE
uniform sampler1D hybridCoefficients; // SIGMA_HYBRID_PRESSURE
uniform sampler1D pressureTable; // FOR PRESSURE_LEVELS only
uniform sampler2D pressureTexCoordTable2D; // HYBRID_SIGMA
uniform sampler3D auxPressureField3D_hPa; // AUXILIARY_PRESSURE_3D
uniform sampler1D transferFunction;
uniform sampler1D lonLatLevAxes;

uniform sampler3D dataVolumeShV; // SIGMA_HYBRID_PRESSURE | PRESSURE_LEVELS // integration / detection var grid
uniform sampler2D surfacePressureShV; // SIGMA_HYBRID_PRESSURE
uniform sampler1D hybridCoefficientsShV; // SIGMA_HYBRID_PRESSURE
uniform sampler1D pressureTableShV; // FOR PRESSURE_LEVELS only
uniform sampler2D pressureTexCoordTable2DShV; // HYBRID_SIGMA
uniform sampler3D auxPressureField3DShV_hPa; // AUXILIARY_PRESSURE_3D
uniform sampler1D transferFunctionShV;
uniform sampler1D lonLatLevAxesShV;

uniform bool    isoEnables[MAX_ISOSURFACES];
uniform float   isoValues[MAX_ISOSURFACES];
uniform int     numIsoValues;

// vectors
// =======

uniform vec3    volumeBottomSECrnr;
uniform vec3    volumeTopNWCrnr;

// scalars
// =======

uniform float   integrationStepSize; // Integration step size
uniform float   minValue;
uniform uint    bisectionSteps;
uniform int     maxNumIterations;
uniform int     numNormalCurves;

uniform int     integrationMode; // integration direction (backwards: -1, forwards: 1)
uniform int     valueMode; // Determine which values are stored at normal curve vertices

// ECMWF-specific
// ==============
uniform vec2    pToWorldZParams;
uniform mat4    mvpMatrix;

// Normal Curve specific
// =====================

uniform float   isoValueStop;

// BBox
uniform vec2 minBBox;
uniform vec2 maxBBox;


/*****************************************************************************
 ***                             INCLUDES
 *****************************************************************************/
// include global definitions
#include "volume_defines.glsl"
// include global structs
#include "volume_global_structs_utils.glsl"
// include hybrid model volume sampling methods
#include "volume_hybrid_utils.glsl"
// include model level volume with auxiliary pressure field sampling methods
#include "volume_auxiliarypressure_utils.glsl"
// include pressure levels volume sampling methods
#include "volume_pressure_utils.glsl"
// defines subroutines and auxiliary ray-casting functions
#include "volume_sample_utils.glsl"
// subroutines for second grid
#include "volume_sample_shv_utils.glsl"


/*****************************************************************************
 ***                              UTILS
 *****************************************************************************/

// Correct the position of any detected iso-surface by using the
// bisection algorithm.
void bisectionCorrection(inout vec3 position,
                         inout float locator,
                         in vec3 prevPosition,
                         in float threshold)
{
    vec3 centerPosition;
    const uint NUM_BISECTION_STEPS = 5;
    float locatorCenter = locator;
    float prevLocator = locator;

    for (int i = 0; i < NUM_BISECTION_STEPS; ++i)
    {
        centerPosition = (position + prevPosition) / 2.0;
        locatorCenter = sampleDataAtPos(centerPosition);

        //bool condition = checkThreshold(scalarCenter, centerPosition);

        if (locatorCenter <= 0)
        {
            position = centerPosition;
            locator = locatorCenter;
        } else {
            prevPosition = centerPosition;
            prevLocator = locatorCenter;
        }
    }
}


float degreesToRadians(float angle)
{
    return angle / 180. * M_PI;
}


// Compute great circle distance between two points on the earth surface
float gcDistanceEarth(const vec2 pointA, const vec2 pointB)
{
    float dlon = degreesToRadians(pointB.x) - degreesToRadians(pointA.x);
    float dlat = degreesToRadians(pointB.y) - degreesToRadians(pointA.y);

    float sin_dlat = sin(dlat/2.);
    float sin_dlon = sin(dlon/2.);
    float a = sin_dlat * sin_dlat
            + cos(pointA.y) * cos(pointB.y)
            * sin_dlon * sin_dlon;
    float c = 2. * asin(min(1.,sqrt(a)));
    return c * EARTH_RADIUS_KM;
}


vec3 computeHorizontalGradient(in vec3 pos, bool geometricLength)
{
    //float PI = 3.1415926535897932384626433832795;
    float deltaLat = (geometricLength) ? 1.12E5 : 1.0;

    vec3 pos_east = vec3(min(pos.x + dataExtentShV.deltaLon, dataExtentShV.dataSECrnr.x), pos.yz);
    vec3 pos_west = vec3(max(pos.x - dataExtentShV.deltaLon, dataExtentShV.dataNWCrnr.x), pos.yz);

    vec3 pos_north = vec3(pos.x, min(pos.y + dataExtentShV.deltaLat, dataExtentShV.dataNWCrnr.y), pos.z);
    vec3 pos_south = vec3(pos.x, max(pos.y - dataExtentShV.deltaLat, dataExtentShV.dataSECrnr.y), pos.z);

    float deltaLon = deltaLat * cos(degreesToRadians(pos.y));

    float hx = (pos_east.x - pos_west.x) * deltaLon;
    float hy = (pos_north.y - pos_south.y) * deltaLat;

    float x1 = sampleShadingDataAtPos(pos_east);
    float x2 = sampleShadingDataAtPos(pos_west);

    float y1 = sampleShadingDataAtPos(pos_north);
    float y2 = sampleShadingDataAtPos(pos_south);

    return normalize(vec3((x1 - x2) / hx, (y1 - y2) / hy, 0.0));
}


void computeLine(uint pointIndex)
{

    float deltaLatKM = 111.2; //km
    float deltaLatM = 1.12E5; //m

    // Initialize gradient step size.
    vec3 h_gradient = initGradientSamplingStep(dataExtentShV);

    int index = int(pointIndex);

    // Measure line length in KM
    float lineLengthKM = 0;
    float segmentLength = 0.;

    // Determine the entry position of the normal curve
    vec3 currentPos = startVertex[index].position.xyz;

    // save height in hPa
    float height = currentPos.z;

    // Compute z-world coordinate
    float worldZ = (log(currentPos.z) - pToWorldZParams.x) * pToWorldZParams.y;
    currentPos.z = worldZ;

    vec3 prevPos = currentPos;

    // Start integrating process and computation of normal curve.
    float locator = sampleDataAtPos(currentPos);
    float tfp = startVertex[index].position.w;

    int direction;
    if (tfp < 0)
    {
       direction = -1;
    }
    else
    {
        direction = 1;
    }

    for (int i = 0; i < maxNumIterations; i++)
    {
        vec3 gradient = computeHorizontalGradient(currentPos, false);
        prevPos = currentPos;
        currentPos -= gradient * integrationStepSize * direction;

        if(sampleShadingDataAtPos(currentPos) < minValue
           || currentPos.x < minBBox.x || currentPos.x > maxBBox.x
           || currentPos.y < minBBox.y || currentPos.y > maxBBox.y)
        {
            endVertex[index].position.xy = prevPos.xy;
            endVertex[index].position.z = height;
            endVertex[index].position.w = lineLengthKM;
            break;
        }

        locator = sampleDataAtPos(currentPos);


        if(locator <= 0.)
        {
            bisectionCorrection(currentPos, locator, prevPos, 0.0);
            endVertex[index].position.xy = currentPos.xy;
            endVertex[index].position.z = height;
            segmentLength = abs(gcDistanceEarth(prevPos.xy, currentPos.xy));
            lineLengthKM += segmentLength;
            endVertex[index].position.w = lineLengthKM;
            break;
        }

        segmentLength = abs(gcDistanceEarth(prevPos.xy, currentPos.xy));
        lineLengthKM += segmentLength;
    }
}


/*****************************************************************************
 ***                           COMPUTE SHADER
 *****************************************************************************/

// Shader that integrates only in one certain direction
shader CSsingleIntegration()
{
    uint pointIndex = gl_GlobalInvocationID.x;

    if (pointIndex >= numNormalCurves) { return; }

    computeLine(pointIndex);
}

/*****************************************************************************
 ***                             PROGRAMS
 *****************************************************************************/

program SingleIntegration
{
    cs(430)=CSsingleIntegration() : in(local_size_x = 32);
};
