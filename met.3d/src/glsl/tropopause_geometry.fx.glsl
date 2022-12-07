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
    smooth float firstDeriv;
    flat   bool isShadow;
};

interface GStoFS
{
    smooth vec3 worldPos;
    smooth vec3 normal;
    smooth float firstDeriv;
};



/*****************************************************************************
 ***                             UNIFORMS
 *****************************************************************************/

/* - variables common to all shader stages */
uniform mat4 mvpMatrix;
uniform vec2 pToWorldZParams;
uniform vec3 lightDirection;
uniform vec3 cameraPosition;

uniform vec4   colour;

uniform sampler1D tfFD;
uniform vec2 tfFDMinMax;

uniform bool useFDFilter;

uniform vec3 bboxMin;
uniform vec3 bboxMax;

uniform bool displayNegative;

uniform float shadowHeight;
uniform vec4 shadowColor;

uniform int lightingMode;

// OIT
#include "oit_utils.glsl"
// functions:
    // insertNewFragmentOIT


// Shadow mapping
#include "shadow_mapping.glsl"
// functions:
    // computeShadowMappingFactor
    //



uniform sampler3D dataVolume;


/*****************************************************************************
 ***                           VERTEX SHADER
 *****************************************************************************/

shader VStriSimple(in vec3 pos : 0,
                   in vec3 normal : 1,
                   in float firstDeriv: 2,
                   out VStoGS Output)
{
    float worldZ = (log(pos.z) - pToWorldZParams.x) * pToWorldZParams.y;
    vec3 worldPos = vec3(pos.x, pos.y, worldZ);
    Output.normal = normal;
    Output.pos = worldPos;
    Output.firstDeriv = firstDeriv;
    Output.isShadow = false;
}


shader VStriSimpleShadow(in vec3 pos : 0,
                         in vec3 normal : 1,
                         in float firstDeriv : 2,
                         out VStoGS Output)
{
    float realWorldZ = (log(pos.z) - pToWorldZParams.x) * pToWorldZParams.y;
    float worldZ = shadowHeight;
    vec3 worldPos = vec3(pos.x, pos.y, realWorldZ);
    Output.normal = normal;
    Output.pos = worldPos;
    Output.firstDeriv = firstDeriv;
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

        bool isShadow = Input[i].isShadow;

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
        Output.firstDeriv = Input[i].firstDeriv;

        EmitVertex();
    }
}

/*****************************************************************************
 ***                          FRAGMENT SHADER
 *****************************************************************************/


shader FStriOIT(in GStoFS Input, out vec4 fragColour)
{
    float firstDeriv = Input.firstDeriv;

    float surfaceAlpha = 1;
    if(firstDeriv < 0.0)
    {
        discard;
        return;
    }
    if (useFDFilter)
    {
        float tFD = min(max((firstDeriv - tfFDMinMax.x) /
                             (tfFDMinMax.y - tfFDMinMax.x), 0), 1);
        surfaceAlpha = texture(tfFD, tFD).a;
    }

    if (surfaceAlpha < MIN_OPACITY) { discard; return;  }

    // compute lighting
    vec3 lightColor = vec3(1, 1, 1);

    vec3 n = normalize(Input.normal);
    vec3 v = normalize(cameraPosition - Input.worldPos);
    vec3 l = normalize(-lightDirection);
    vec3 h = normalize(v + l);

    vec3 colourvec = colour.xyz;
    vec3 Ka = 0.2 * colourvec;
    vec3 Kd = 0.6 * colourvec;
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
    float firstDeriv = Input.firstDeriv;

    float surfaceAlpha = 1;

    if (useFDFilter)
    {
        float tfirstDeriv = min(max((firstDeriv - tfFDMinMax.x) /
                             (tfFDMinMax.y - tfFDMinMax.x), 0), 1);
        surfaceAlpha = texture(tfFD, tfirstDeriv).a;
    }

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
