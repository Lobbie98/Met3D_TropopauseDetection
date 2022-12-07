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
// OIT
#include "oit_utils.glsl"
// functions:
    // insertNewFragmentOIT


// Shadow mapping
#include "shadow_mapping.glsl"
// functions:
    // computeShadowMappingFactor


/*****************************************************************************
 ***                             CONSTANTS
 *****************************************************************************/

/* - use const datatype instead of #define */
const float MIN_OPACITY = 0.05;

/*****************************************************************************
 ***                             INTERFACES
 *****************************************************************************/

interface VStoGS
{
    smooth vec3 pos;
    smooth float tfp;
    smooth float abz;
    smooth float strength;
    smooth float type;
    smooth float breadth;
    smooth float alpha;
    smooth float shading;
    flat bool isShadow;
};

interface GStoFS
{
    smooth vec3 worldPos;
    smooth vec3 normalPS;
    smooth float tfp;
    smooth float abz;
    smooth float strength;
    smooth float type;
    smooth float breadth;
    smooth float alpha;
    smooth float shading;
};



/*****************************************************************************
 ***                             UNIFORMS
 *****************************************************************************/

uniform float   tubeRadius;
uniform float   worldZPos;            // parameters to convert p[hPa] to world z
uniform mat4    mvpMatrix;

uniform bool    normalized;

uniform vec3    cameraPosition;
uniform vec3    lightDirection;

uniform bool toDepth;
uniform sampler2D depthTex;

uniform sampler1D transferFunction;

uniform float   tfMinimum;
uniform float   tfMaximum;

uniform sampler1D tfTFP;
uniform vec2 tfTFPMinMax;

uniform sampler1D tfABZ;
uniform vec2 tfABZMinMax;

uniform sampler1D tfBreadth;
uniform vec2 tfBreadthMinMax;

uniform bool useTFPFilter;
uniform bool useABZFilter;
uniform bool useBreadthFilter;

uniform vec3 bboxMin;
uniform vec3 bboxMax;

uniform vec4    shadowColor;
uniform float shadowHeight;

// GPU based filter for frontal strength
uniform bool useFSFilter;
uniform sampler1D tfFS;
uniform vec2 tfFSMinMax;


/*****************************************************************************
 ***                           VERTEX SHADER
 *****************************************************************************/

shader VStubeSimple(in vec2 pos : 0,
                    in float tfp : 1,
                    in float abz : 2,
                    in float strength : 3,
                    in float type : 4,
                    in float breadth : 5,
                    in float alpha : 6,
                    in float shading : 7,
                    out VStoGS Output)
{
    if (pos.x ==-1 && pos.y == -1)
    {
        Output.pos = vec3(-1, -1, -1);
    }
    else
    {
        Output.pos = vec3(pos, worldZPos);
    }
    Output.alpha = alpha;
    Output.shading = shading;
    Output.tfp = tfp;
    Output.abz = abz;
    Output.strength = strength;
    Output.breadth = breadth;
    Output.type = type;
    Output.isShadow = false;
}


shader VStubeSimpleShadow(in vec2 pos : 0,
                          in float tfp : 1,
                          in float abz : 2,
                          in float strength : 3,
                          in float type : 4,
                          in float breadth : 5,
                          in float alpha : 6,
                          in float shading : 7,
                          out VStoGS Output)
{
    if (pos.x ==-1 && pos.y == -1)
    {
        Output.pos = vec3(-1, -1, -1);
    }
    else
    {
        Output.pos = vec3(pos, worldZPos);
    }
    Output.shading = shading;
    Output.alpha = alpha;
    Output.tfp = tfp;
    Output.abz = abz;
    Output.strength = strength;
    Output.type = type;
    Output.breadth = breadth;
    Output.isShadow = true;
}


/*****************************************************************************
 ***                          GEOMETRY SHADER
 *****************************************************************************/
// calculate new basis
void calculateRayBasis( in vec3 prevPos,
                        in vec3 pos,
                        out vec3 normal,
                        out vec3 binormal)
{
    vec3 lineDir = normalize(pos - prevPos);

    normal = cross(lineDir,vec3(0,0,1));
    //normal = normalize(normal);
    if (length(normal) <= 0.01)
        {
            normal = cross(lineDir, vec3(1, 0, 0));

            if (length(normal) <= 0.01)
            {
                normal = cross(lineDir, vec3(0, 0, 1));
            }
        }


    normal = normalize(normal);

    binormal = cross(lineDir,normal);
    binormal = normalize(binormal);

}


/*****************************************************************************/
shader GStube(in VStoGS Input[], out GStoFS Output)
{

    Output.alpha = 0;
    Output.shading = 0;

    if (Input[1].pos.z == -1 || Input[2].pos.z == -1)
    {
        return;
    }

    //    o ----o______o-----o
    vec3 pos0, pos1, pos2, pos3;

    float value1, value2, transp1, transp2;

    value1 = Input[1].shading;
    value2 = Input[2].shading;

    transp1 = Input[1].alpha;
    transp2 = Input[2].alpha;

    if (Input[0].pos.z == -1 && Input[3].pos.z == -1)
    {
        pos0 = pos1 = Input[1].pos;
        pos2 = pos3 = Input[2].pos;
    }
    else if (Input[0].pos.z == -1)
    {
        pos0 = pos1 = Input[1].pos;
        pos2 = Input[2].pos;
        pos3 = Input[3].pos;
    }
    else if (Input[3].pos.z == -1)
    {
        pos0 = Input[0].pos;
        pos1 = Input[1].pos;
        pos2 = pos3 = Input[2].pos;
    }
    else
    {
        pos0 = Input[0].pos;
        pos1 = Input[1].pos;
        pos2 = Input[2].pos;
        pos3 = Input[3].pos;
    }

    // set tubes that are located outside the bounding box 100% transparent
    if (pos1.x < bboxMin.x || pos1.x > bboxMax.x
        || pos1.y < bboxMin.y || pos1.y > bboxMax.y
        || pos2.x < bboxMin.x || pos2.x > bboxMax.x
        || pos2.y < bboxMin.y || pos2.y > bboxMax.y)
    {
        transp1 = 0;
        transp2 = 0;
    }

    bool isShadow = Input[0].isShadow;
    if (isShadow)
    {
        pos0.z = pos1.z = pos2.z = pos3.z = shadowHeight;
    }


    // generate tube

    vec3 normalPrev, binormalPrev, normalNext, binormalNext;

    calculateRayBasis(pos0, pos2, normalPrev, binormalPrev);

    calculateRayBasis(pos1, pos3, normalNext, binormalNext);

    int tube_segments = 8;
    float angle_t = 360. / float(tube_segments);

    vec3 tube_positions[8];

    for (int t = 0; t <= tube_segments; t++)
    {
        float angle = radians(angle_t * t);
        float cosi = cos(angle) * tubeRadius;
        float sini = sin(angle) * tubeRadius;

        vec3 prev_world_pos = pos1 + cosi * normalPrev + sini * binormalPrev;
        vec3 next_world_pos = pos2 + cosi * normalNext + sini * binormalNext;

        Output.shading = value1;
        Output.alpha = transp1;
        Output.tfp = Input[1].tfp;
        Output.abz = Input[1].abz;
        Output.strength = Input[1].strength;
        Output.breadth = Input[1].breadth;
        gl_Position = mvpMatrix * vec4(prev_world_pos,1);
        Output.normalPS = normalize(prev_world_pos - pos1);
        Output.worldPos = prev_world_pos;

        EmitVertex();

        Output.shading = value2;
        Output.alpha = transp2;
        Output.tfp = Input[2].tfp;
        Output.abz = Input[2].abz;
        Output.strength = Input[2].strength;
        Output.breadth = Input[2].breadth;
        gl_Position = mvpMatrix * vec4(next_world_pos,1);
        Output.normalPS = normalize(next_world_pos - pos2);
        Output.worldPos = next_world_pos;

        EmitVertex();

        tube_positions[t] = next_world_pos;
    }

    EndPrimitive();
}

/*****************************************************************************
 ***                          FRAGMENT SHADER
 *****************************************************************************/

void getBlinnPhongColor(in vec3 worldPos, in vec3 normalDir, in float t,
                        out vec3 color)
{
    const vec3 lightColor = vec3(1,1,1);

    if (!normalized)
    {
        t = (t - tfMinimum) / (tfMaximum - tfMinimum);
    }

    const vec3 ambientColor = texture(transferFunction, t).rgb;

    const vec3 kA = 0.3 * ambientColor;
    const vec3 kD = 0.5 * ambientColor;
    const float kS = 0.2;
    const float s = 10;

    const vec3 n = normalize(normalDir);
    const vec3 v = normalize(cameraPosition - worldPos);
    const vec3 l = normalize(-lightDirection); // specialCase
    const vec3 h = normalize(v + l);

    vec3 diffuse = kD * clamp(abs(dot(n,l)),0.0,1.0) * lightColor;
    vec3 specular = kS * pow(clamp(abs(dot(n,h)),0.0,1.0),s) * lightColor;

    color = kA + diffuse + specular;
}


shader FSGeom(in GStoFS Input, out vec4 fragColor)
{
    if (Input.valuePS == -1) { discard; }

    vec3 color;
    getBlinnPhongColor(Input.worldPos, Input.normalPS, Input.valuePS, color);
    fragColor = vec4(color,1);

}


shader FStubeOIT(in GStoFS Input, out vec4 fragColor)
{
    if (Input.shading == -1) { discard; }

    float tfp = Input.tfp;
    float abz = Input.abz;
    float strength = Input.strength;
    float breadth = Input.breadth;

    float aTFP = 1;
    float aABZ = 1;
    float aFS = 1;
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

    float alpha = Input.alpha * aTFP * aABZ * aBreadth * aFS;


    if (alpha < MIN_OPACITY) { discard; return;  }

    float shadowFactor = computeShadowMappingFactor(Input.worldPos);

    vec3 shadingColor;
    getBlinnPhongColor(Input.worldPos, Input.normalPS, Input.shading, shadingColor);

    vec4 color = vec4(shadingColor, alpha);
//    fragColor = vec4(color, alpha);

    if (!inShadowMappingMode)
    {
        insertNewFragmentOIT(color, gl_FragCoord.xyz);
        discard;
    }
    else
    {
        if (alpha < 0.3) { discard; }
    }

}


shader FStubeShadow(in GStoFS Input, out vec4 fragColor)
{
    if (Input.shading == -1) { discard; }

    float tfp = Input.tfp;
    float abz = Input.abz;
    float strength = Input.strength;
    float breadth = Input.breadth;

    float aTFP = 1;
    float aABZ = 1;
    float aBreadth = 1;
    float aFS = 1;

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

    float alpha = Input.alpha * aTFP * aABZ * aBreadth * aFS;

    if (alpha < MIN_OPACITY) { discard; return;  }

    vec4 tubeColor = vec4(shadowColor.rgb, alpha * shadowColor.a);

    insertNewFragmentOIT(tubeColor, gl_FragCoord.xyz, true);
    discard;

}


/*****************************************************************************
 ***                             PROGRAMS
 *****************************************************************************/

program TubeFilteringOIT
{
    vs(430)=VStubeSimple();
    gs(430)=GStube() : in(lines_adjacency), out(triangle_strip, max_vertices = 128);
    fs(430)=FStubeOIT() : in(early_fragment_tests);
};

program TubeFilteringShadowMap
{
    vs(430)=VStubeSimple();
    gs(430)=GStube() : in(lines_adjacency), out(triangle_strip, max_vertices = 128);
    fs(430)=FStubeOIT();
};


program TubeFilteringShadow
{
    vs(430)=VStubeSimpleShadow();
    gs(430)=GStube() : in(lines_adjacency), out(triangle_strip, max_vertices = 128);
    fs(430)=FStubeShadow() : in(early_fragment_tests);
};
