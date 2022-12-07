/******************************************************************************
**
**  This file is part of Met.3D -- a research environment for the
**  three-dimensional visual exploration of numerical ensemble weather
**  prediction data.
**
**  Copyright 2021 Andreas Beckert
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

/*****************************************************************************
 ***                             CONSTANTS
 *****************************************************************************/

/* - use const datatype instead of #define */
const uint MAX_FRAGMENTS_PER_PIXEL = 128;
const float MIN_OPACITY = 0.05;

const int SURFACE_2D = 0;
const int PRESSURE_LEVELS_3D = 1;
const int HYBRID_SIGMA_PRESSURE_3D = 2;
const int POTENTIAL_VORTICITY_2D = 3;
const int LOG_PRESSURE_LEVELS_3D = 4;
const int AUXILIARY_PRESSURE_3D = 5;

const int TRANSFER_COLORMODE = 0;
const int PRESSURE_COLORMODE = 1;
const int TFP_COLORMODE = 2;
const int STRENGTH_THETAW_COLORMODE = 3;
const int STRENGTH_THETA_COLORMODE = 4;
const int SLOPE_COLORMODE = 5;
const int BREADTH_COLORMODE = 6;

/*****************************************************************************
 ***                             INTERFACES
 *****************************************************************************/

/* - connection between different stages.
   - defines custom input and output variables
   - acts like a struct
*/

interface VStoGS
{
    smooth vec3 pos;
    smooth vec3 normal;
    smooth float tfp;
    smooth float abz;
    smooth float strength;
    smooth float type;
    smooth float breadth;
    smooth float slope;
    smooth float alpha;
    smooth float shading;
    flat   bool isShadow;
};

interface GStoFS
{
    smooth vec3 worldPos;
    smooth vec3 normal;
    smooth float tfp;
    smooth float abz;
    smooth float strength;
    smooth float type;
    smooth float breadth;
    smooth float slope;
    smooth float alpha;
    smooth float shading;
};



/*****************************************************************************
 ***                             UNIFORMS
 *****************************************************************************/

/* - variables common to all shader stages */
uniform mat4 mvpMatrix;
uniform vec2 pToWorldZParams;
uniform vec3 lightDirection;
uniform vec3 cameraPosition;

uniform sampler1D tfShading;
uniform vec2 tfShadingMinMax;

uniform sampler1D tfTFP;
uniform vec2 tfTFPMinMax;

uniform sampler1D tfABZ;
uniform vec2 tfABZMinMax;

uniform sampler1D tfSlope;
uniform vec2 tfSlopeMinMax;

uniform sampler1D tfBreadth;
uniform vec2 tfBreadthMinMax;

uniform bool useTFPFilter;
uniform bool useABZFilter;
uniform bool useSlopeFilter;
uniform bool useBreadthFilter;

uniform vec3 bboxMin;
uniform vec3 bboxMax;

uniform bool displayFrontTypes;
uniform bool displayColdSideFront;

uniform float shadowHeight;
uniform vec4 shadowColor;
//uniform int shadingMode;

uniform int lightingMode;

// GPU based filter for frontal strength
uniform bool useFSFilter;
uniform sampler1D tfFS;
uniform vec2 tfFSMinMax;

// OIT
#include "oit_utils.glsl"
// functions:
    // insertNewFragmentOIT


// Shadow mapping
#include "shadow_mapping.glsl"
// functions:
    // computeShadowMappingFactor
    //


// contains precomputed z-coordinates
//uniform sampler1D pressureTable; // PRESSURE_LEVEL
//// contains hybrid coefficients a,b
//uniform sampler1D hybridCoefficients; // HYBRID_SIGMA
//// contains surface pressure at grid point (i, j)
//uniform sampler2D surfacePressure; // HYBRID_SIGMA
//// contains pressure field
//uniform sampler3D auxPressureField3D_hPa; // AUXILIARY_PRESSURE_3D
//// contains precomputed z-coordinates
//uniform sampler2D pressureTexCoordTable2D; // HYBRID_SIGMA

//#ifdef ENABLE_MINMAX_ACCELERATION
//// Regular grid that uniformly divides space and stores min/max for each voxel.
//// Is used to skip regions in which an isosurface cannot be located. See
//// MStructuredGrid::getMinMaxAccelTexture3D().
//uniform sampler3D minMaxAccel3D;
//#endif

uniform sampler3D dataVolume;
//uniform sampler1D lonLatLevAxes;

//uniform int     numIsoValues;
//uniform uint    bisectionSteps;
//uniform float   isoValues[8];


/*****************************************************************************
 ***                             INCLUDES
 *****************************************************************************/

/*
  - you can include several files via: #include "filename.glsl"
  - #include is simply replaced by glfx with included code
  - error messages carry an index, indicating which file caused the error
  */

//#include "volume_defines.glsl"

// include global structs
//#include "volume_global_structs_utils.glsl"
// include hybrid model volume sampling methods
//#include "volume_hybrid_utils.glsl"
// include model level volume with auxiliary pressure field sampling methods
//#include "volume_auxiliarypressure_utils.glsl"
// include pressure levels volume sampling methods
//#include "volume_pressure_utils.glsl"
// defines subroutines and auxiliary ray-casting functions
//#include "volume_sample_utils.glsl"


/*****************************************************************************
 ***                           VERTEX SHADER
 *****************************************************************************/

shader VStriSimple(in vec3 pos : 0,
                   in vec3 normal : 1,
                   in float tfp : 2,
                   in float abz : 3,
                   in float strength : 4,
                   in float type : 5,
                   in float breadth : 6,
                   in float slope : 7,
                   in float alpha : 8,
                   in float shading : 9,
                   out VStoGS Output)
{
    float worldZ = (log(pos.z) - pToWorldZParams.x) * pToWorldZParams.y;
    vec3 worldPos = vec3(pos.x, pos.y, worldZ);
    Output.normal = normal;
    Output.pos = worldPos;
    Output.tfp = tfp;
    Output.abz = abz;
    Output.strength = strength;
    Output.type = type;
    Output.breadth = breadth;
    Output.slope = slope;
    Output.alpha = alpha;
    Output.shading = shading;
    Output.isShadow = false;
}


shader VStriSimpleShadow(in vec3 pos : 0,
                         in vec3 normal : 1,
                         in float tfp : 2,
                         in float abz : 3,
                         in float strength : 4,
                         in float type : 5,
                         in float breadth : 6,
                         in float slope : 7,
                         in float alpha : 8,
                         in float shading : 9,
                         out VStoGS Output)
{
    float realWorldZ = (log(pos.z) - pToWorldZParams.x) * pToWorldZParams.y;
    float worldZ = shadowHeight;
    vec3 worldPos = vec3(pos.x, pos.y, realWorldZ);
    Output.normal = normal;
    Output.pos = worldPos;
    Output.tfp = tfp;
    Output.abz = abz;
    Output.strength = strength;
    Output.type = type;
    Output.breadth = breadth;
    Output.slope = slope;
    Output.alpha = alpha;
    Output.shading = shading;
    Output.isShadow = true;
}

/*****************************************************************************
 ***                          GEOMETRY SHADER
 *****************************************************************************/

shader GStri(in VStoGS Input[], out GStoFS Output)
{
    for (int i = 0; i < 3; ++i)
    {
        vec3 pos = Input[i].pos;
        vec3 normal = Input[i].normal;
        float tfp = Input[i].tfp;
        bool isShadow = Input[i].isShadow;

        // filter out cold-side frontal surfaces
        if (!displayColdSideFront && tfp < 0) { break; }

        // discard all triangles that are located outside the bounding box
        if (pos.x < bboxMin.x || pos.x > bboxMax.x
            || pos.y < bboxMin.y || pos.y > bboxMax.y) { break; }

        if (isShadow)
        {
            if (pos.z > bboxMax.z || pos.z < bboxMin.z) { break; }

            pos.z = shadowHeight;
        }
        else
        {
            if (pos.z > bboxMax.z) { break; }
            if (pos.z < bboxMin.z) { break; }
        }

        // correct the vertical component of the normal
        float pressure = exp(pos.z / pToWorldZParams.y + pToWorldZParams.x);
        float nextP = pressure + 1;
        float nextWorldZ = (log(nextP) - pToWorldZParams.x) * pToWorldZParams.y;
        float deltaZ = nextWorldZ - pos.z;
        normal.z /= deltaZ;
        normalize(normal);

        // project vertex onto screen and output vars to fragment shader
        gl_Position = mvpMatrix * vec4(pos, 1);
        Output.worldPos = pos;
        Output.normal = normal;
        Output.tfp = Input[i].tfp;
        Output.abz = Input[i].abz;
        Output.strength = Input[i].strength;
        Output.type = Input[i].type;
        Output.breadth = Input[i].breadth;
        Output.slope = Input[i].slope;
        Output.alpha = Input[i].alpha;
        Output.shading = Input[i].shading;

        EmitVertex();
    }
}

/*****************************************************************************
 ***                          FRAGMENT SHADER
 *****************************************************************************/


shader FStriOIT(in GStoFS Input, out vec4 fragColour)
{
    float tfp = Input.tfp;
    float abz = Input.abz;
    float strength = Input.strength;
    float slope = Input.slope;
    float breadth = Input.breadth;
    float aFilter = Input.alpha;
    float shading = Input.shading;

    if (displayColdSideFront && tfp < 0)
    {
        tfp *= -1;
    }

    float aTFP = 1;
    float aFS = 1;
    float aABZ = 1;
    float aSlope = 1;
    float aBreadth = 1;

    if (useTFPFilter)
    {
        float tTFP = min(max((tfp - tfTFPMinMax.x) /
                             (tfTFPMinMax.y - tfTFPMinMax.x), 0), 1);
        aTFP = texture(tfTFP, tTFP).a;
    }

    if (useABZFilter)
    {
        float tABZ = min(max((abz - tfABZMinMax.x) /
                             (tfABZMinMax.y - tfABZMinMax.x), 0), 1);
        aABZ = texture(tfABZ, tABZ).a;
    }

    if (useSlopeFilter)
    {
        float tSlope = min(max((slope - tfSlopeMinMax.x) /
                             (tfSlopeMinMax.y - tfSlopeMinMax.x), 0), 1);
        aSlope = texture(tfSlope, tSlope).a;
    }

    if (useBreadthFilter)
    {
        float tBreadth = min(max((breadth - tfBreadthMinMax.x) /
                             (tfBreadthMinMax.y - tfBreadthMinMax.x), 0), 1);
        aBreadth = texture(tfBreadth, tBreadth).a;
    }

    if (useFSFilter)
    {
        float tFS = 0;
        if (breadth > 0.1)
        {
            tFS = min(max(((strength / breadth * 100) - tfFSMinMax.x) /
                                 (tfFSMinMax.y - tfFSMinMax.x), 0), 1);
        }
        aFS = texture(tfFS, tFS).a;
    }

    vec3 surfaceColor = vec3(1);

    if (displayFrontTypes)
    {
        bool type = bool(Input.type);

        // Colors from HCL editor http://hclwizard.org/hclwizard/
        surfaceColor = (type) ? vec3(0.541, 0.098, 0.137) : vec3(0.008, 0.247, 0.647);
    }
    else
    {
        float t = 1;
        t = min(max((shading - tfShadingMinMax.x)
                / (tfShadingMinMax.y - tfShadingMinMax.x), 0), 1);
        surfaceColor = texture(tfShading, t).rgb;
    }

//    vec3 surfaceColor = texture(tfStrength, tStrength2).rgb;
    float surfaceAlpha = aFilter * aTFP * aABZ * aSlope * aBreadth * aFS;

    if (surfaceAlpha < MIN_OPACITY) { discard; return;  }

    // compute lighting

    vec3 lightColor = vec3(1, 1, 1);

    vec3 n = normalize(Input.normal);
    vec3 v = normalize(cameraPosition - Input.worldPos);
    vec3 l = normalize(-lightDirection);
    vec3 h = normalize(v + l);

    vec3 Ka = 0.2 * surfaceColor;
    vec3 Kd = 0.6 * surfaceColor;
    float Ks = 0.2;
    float m = 10;

    vec3 shadingColor = Ka;

    float shadowFactor = computeShadowMappingFactor(Input.worldPos);

    if (!inShadowMappingMode)
    {
        lightColor *= (0.5 + shadowFactor * 0.5);
    }

    if (lightingMode <= 1)
    {
        float diffuse = (lightingMode == 0) ? clamp(abs(dot(n, l)), 0, 1) : clamp(dot(n, l), 0, 1);
        vec3 diffuseColor = Kd * lightColor * diffuse;

        float specular = (lightingMode == 0) ? max(0., pow(clamp(abs(dot(n, h)), 0., 1.), m)) : max(0., pow(clamp(dot(n, h), 0., 1.), m));
        vec3 specularColor = Ks * lightColor * specular;

        shadingColor += diffuseColor + specularColor;
    }
    else
    {
        int numLightSources = lightingMode;

        vec3 lights[2];
        vec3 lightColors[2];

        // First light from light direction, second from head light.
        lights[0] = normalize(-lightDirection);
        lights[1] = v;

        lightColors[0] = lightColor;
        // Scale the head light according to cos( light dir, view dir).
        lightColors[1] = lightColors[0] * 1.0 - clamp(dot(lights[0], lights[1]), 0., 1.);

        // For each light source add diffuse and specular terms:
        // =====================================================

        for (int i = 0; i < numLightSources; i++)
        {
            // Half vector, mid-way between view and light directions.
            vec3 h = normalize(v + lights[i]);

            float diffuseLight = (lightingMode == 2) ? clamp(abs(dot(n, lights[i])), 0., 1.)
                                                        : clamp(dot(n, lights[i]), 0., 1.);
            vec3 diffuse = Kd * lightColors[i] * diffuseLight;

            float specularLight = (lightingMode == 2) ?  max(0., pow(clamp(abs(dot(n, h)), 0., 1.), m))
                                                      : max(0., pow(clamp(dot(n, h), 0., 1.), m));
            vec3 specular = Ks * lightColors[i] * specularLight;

            shadingColor += diffuse + specular;
        }
    }

    vec4 color = vec4(shadingColor, surfaceAlpha);

    if (!inShadowMappingMode)
    {
        insertNewFragmentOIT(color, gl_FragCoord.xyz);
        discard;
    }
    else
    {
        if (surfaceAlpha < 0.3) { discard; }
    }
    discard;
}

shader FStriShadow(in GStoFS Input, out vec4 fragColour)
{
    float tfp = Input.tfp;
    float abz = Input.abz;
    float strength = Input.strength;
    float slope = Input.slope;
    float breadth = Input.breadth;
    float aFilter = Input.alpha;

    bool type = bool(Input.type);

    if (displayColdSideFront && tfp < 0)
    {
        tfp *= -1;
    }

    float aTFP = 1;
    float aFS = 1;
    float aABZ = 1;
    float aSlope = 1;
    float aBreadth = 1;

    if (useTFPFilter)
    {
        float tTFP = min(max((tfp - tfTFPMinMax.x) /
                             (tfTFPMinMax.y - tfTFPMinMax.x), 0), 1);
        aTFP = texture(tfTFP, tTFP).a;
    }

    if (useABZFilter)
    {
        float tABZ = min(max((abz - tfABZMinMax.x) /
                             (tfABZMinMax.y - tfABZMinMax.x), 0), 1);
        aABZ = texture(tfABZ, tABZ).a;
    }

    if (useSlopeFilter)
    {
        float tSlope = min(max((slope - tfSlopeMinMax.x) /
                             (tfSlopeMinMax.y - tfSlopeMinMax.x), 0), 1);
        aSlope = texture(tfSlope, tSlope).a;
    }

    if (useBreadthFilter)
    {
        float tBreadth = min(max((breadth - tfBreadthMinMax.x) /
                             (tfBreadthMinMax.y - tfBreadthMinMax.x), 0), 1);
        aBreadth = texture(tfBreadth, tBreadth).a;
    }

    if (useFSFilter)
    {
        float tFS = 0;
        if (breadth > 0.1)
        {
            tFS = min(max(((strength / breadth * 100) - tfFSMinMax.x) /
                                 (tfFSMinMax.y - tfFSMinMax.x), 0), 1);
        }
        aFS = texture(tfFS, tFS).a;
    }

    float surfaceAlpha = aFilter * aTFP * aABZ * aSlope * aBreadth * aFS;

    if (surfaceAlpha < MIN_OPACITY) { discard; return; }

    vec4 surfaceColor = vec4(shadowColor.rgb, surfaceAlpha * shadowColor.a);

    insertNewFragmentOIT(surfaceColor, gl_FragCoord.xyz, true);
    discard;
}




/*****************************************************************************
 ***                             PROGRAMS
 *****************************************************************************/

program TriangleFilteringOIT
{
    vs(430)=VStriSimple();
    gs(430)=GStri() : in(triangles), out(triangle_strip, max_vertices = 6);
    fs(430)=FStriOIT() : in(early_fragment_tests);
};

program TriangleFilteringShadowMap
{
    vs(430)=VStriSimple();
    gs(430)=GStri() : in(triangles), out(triangle_strip, max_vertices = 6);
    fs(430)=FStriOIT();
};

program TriangleFilteringShadow
{
    vs(430)=VStriSimpleShadow();
    gs(430)=GStri() : in(triangles), out(triangle_strip, max_vertices = 6);
    fs(430)=FStriShadow() : in(early_fragment_tests);
};
