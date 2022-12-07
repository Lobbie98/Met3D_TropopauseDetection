/******************************************************************************
**
**  This file is part of Met.3D -- a research environment for the
**  three-dimensional visual exploration of numerical ensemble weather
**  prediction data.
**
**  Copyright 2015 Marc Rautenhaus
**  Copyright 2015 Michael Kern
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

// requires volume_global_structs_utils.glsl
// requires volume_<datatype>_utils.glsl
// requires volume_sample_utils.glsl

uniform DataVolumeExtent dataExtentSource;

// sample shading data at given world position
float sampleSourceDataAtPos(in vec3 pos)
{
    // case PRESSURE_LEVEL_3D
    if (dataExtentSource.levelType == 0)
    {
        return samplePressureLevelVolumeAtPos(dataVolumeSource, dataExtentSource,
                                              pressureTableSource, pos);
    }
    // case HYBRID_SIGMA_PRESSURE_3D
    else if (dataExtentSource.levelType == 1)
    {
        return sampleHybridSigmaVolumeAtPos(dataVolumeSource, dataExtentSource,
                                            surfacePressureSource, hybridCoefficientsSource,
                                            pos);
    }

    return 0;
}


// sample shading data at given world position
float sampleSourceDataAtPos_maxNeighbour(in vec3 pos)
{
    // case PRESSURE_LEVEL_3D
    if (dataExtentSource.levelType == 0)
    {
        return samplePressureLevelVolumeAtPos_maxNeighbour(
                                        dataVolumeSource, dataExtentSource, pos);
    }
    // case HYBRID_SIGMA_PRESSURE_3D
    else if (dataExtentSource.levelType == 1)
    {
        return sampleHybridSigmaVolumeAtPos_maxNeighbour(
                                        dataVolumeSource, dataExtentSource,
                                        surfacePressureSource, hybridCoefficientsSource,
                                        pos);
    }

    return 0;
}



vec3 sourceGradientAtPos(in vec3 pos, in vec3 h, float ztop, float zbot)
{
    vec3 gradient;

    // 1. Sample with horizontal displacement. Check if the grid is cyclic in
    // longitude.

    float deltaLat = 1.112E2; // km

    vec3 pos_east;
    vec3 pos_west;
    if (dataExtentSource.gridIsCyclicInLongitude)
    {
        pos_east = vec3(pos.x + h.x, pos.yz);
        pos_west = vec3(pos.x - h.x, pos.yz);
    }
    else
    {
        // Non-cyclic grids: clamp to data boundary.
        pos_east = vec3(min(pos.x + h.x, dataExtentSource.dataSECrnr.x), pos.yz);
        pos_west = vec3(max(pos.x - h.x, dataExtentSource.dataNWCrnr.x), pos.yz);
    }

    float x1 = sampleSourceDataAtPos(pos_east);
    float x2 = sampleSourceDataAtPos(pos_west);
    float hx = (pos_east.x - pos_west.x) * cos(pos.y / 180.0 * 3.14) * deltaLat;

    vec3 pos_north = vec3(pos.x, min(pos.y + h.y, dataExtentSource.dataNWCrnr.y), pos.z);
    vec3 pos_south = vec3(pos.x, max(pos.y - h.y, dataExtentSource.dataSECrnr.y), pos.z);
    float y1 = sampleSourceDataAtPos(pos_north);
    float y2 = sampleSourceDataAtPos(pos_south);
    float hy = (pos_north.y - pos_south.y) * deltaLat;

    // 2. Sample with vertical displacement, considering data volume
    // boundaries.

    vec3 pos_top = vec3(pos.xy, min(pos.z + h.z, ztop));
    vec3 pos_bot = vec3(pos.xy, max(pos.z - h.z, zbot));
    float z1 = sampleSourceDataAtPos(pos_top);
    float z2 = sampleSourceDataAtPos(pos_bot);
    float hz = pos_top.z - pos_bot.z;

    // 3. Compute gradient.

    gradient = vec3((x1 - x2) / hx, (y1 - y2) / hy, (z1 - z2) / abs(hz));

    return normalize(gradient);
}

vec3 sourceGradientPressureLevel(in vec3 pos, in vec3 h)
{
    // Find the pressue levels enclosing pos.z to estimate hz for the
    // vertical difference.

    // Compute pressure from world z coordinate.
    float p = exp(pos.z / pToWorldZParams.y + pToWorldZParams.x);

    // Binary search to find model levels k, k1 that enclose pressure level p.
    int k = 0;
    int k1 = dataExtentSource.nLev - 1;
    int vertOffset = dataExtentSource.nLon + dataExtentSource.nLat;

    while (abs(k1 - k) > 1)
    {
        int kMid = (k + k1) / 2;
        float pMid = texelFetch(lonLatLevAxesSource, kMid + vertOffset, 0).a;
        if (p >= pMid) k = kMid; else k1 = kMid;
    }

    float lnPk  = log(texelFetch(lonLatLevAxesSource, k + vertOffset, 0).a);
    float lnPk1 = log(texelFetch(lonLatLevAxesSource, k1 + vertOffset, 0).a);

    float hz1 = (lnPk - pToWorldZParams.x) * pToWorldZParams.y;
    float hz2 = (lnPk1 - pToWorldZParams.x) * pToWorldZParams.y;
    vec3 h_ = vec3(h.x, h.y, abs(hz1 - hz2));

    float ztop = dataExtentSource.dataNWCrnr.z;
    float zbot = dataExtentSource.dataSECrnr.z;

    return sourceGradientAtPos(pos, h_, ztop, zbot);
}

vec3 sourceGradientModelLevels(in vec3 pos, in vec3 h)
{
    // Determine the approximate distance between the model levels above and
    // below the current position ("hz").
    float zbot, ztop;
    getHybridSigmaBotTopLevelAtPos(dataExtentSource, surfacePressureSource,
                                   hybridCoefficientsSource, pos, zbot, ztop);

    float hz = getHybridApproxWorldZLevelDistanceAtPos(dataExtentSource,
                                                       surfacePressureSource,
                                                       hybridCoefficientsSource, pos);
    vec3 h_ = vec3(h.x, h.y, hz);

    return sourceGradientAtPos(pos, h_, ztop, zbot);
}

vec3 computeSourceGradient(in vec3 pos, in vec3 h)
{
    // case PRESSURE_LEVEL_3D
    if (dataExtentSource.levelType == 0)
    {
        return sourceGradientPressureLevel(pos, h);
    }
    // case HYBRID_SIGMA_PRESSURE_3D
    else if (dataExtentSource.levelType == 1)
    {
        return sourceGradientModelLevels(pos, h);
    }
    else
    {
        return vec3(0, 0, 0);
    }
}
