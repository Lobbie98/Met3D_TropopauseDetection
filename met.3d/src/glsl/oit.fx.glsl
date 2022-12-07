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
******************************************************************************/

/*****************************************************************************
 ***                             CONSTANTS
 *****************************************************************************/

/* - use const datatype instead of #define */
const uint MAX_FRAGMENTS_PER_PIXEL = 128;
const float MIN_OPACITY = 0.1;

//#define DEBUG_FRAGMENT_COUNT

/*****************************************************************************
 ***                             INTERFACES
 *****************************************************************************/

/* - connection between different stages.
   - defines custom input and output variables
   - acts like a struct

    interface VStoFS {
        smooth vec3 pos;
    };
*/


/*****************************************************************************
 ***                             UNIFORMS
 *****************************************************************************/

/* - variables common to all shader stages */

// OIT
layout(binding=0, offset=0)   uniform atomic_uint  atomicFragCount;
layout(binding=0, r32ui)      uniform uimage2D     headBuffer;
layout(binding=1, r32ui)      uniform uimage2D     headBufferShadow;
//layout(binding=1, rgba32ui)   uniform uimageBuffer linkedList;

// Shader buffer storage object that represents the fragment linked list
struct Fragment
{
    uvec3 data;
};


layout (std430, binding=0) buffer LinkedListBuffer
{
    Fragment linkedList[];
};


/*****************************************************************************
 ***                             INCLUDES
 *****************************************************************************/

/*
  - you can include several files via: #include "filename.glsl"
  - #include is simply replaced by glfx with included code
  - error messages carry an index, indicating which file caused the error
  */


/*****************************************************************************
 ***                           VERTEX SHADER
 *****************************************************************************/


shader VSOIT()
{
    vec4 screenPos = vec4(0, 0, 0, 1);

    if (gl_VertexID == 0) {         screenPos.xy = vec2(-1, -1); }
    else if (gl_VertexID == 1) {    screenPos.xy = vec2( 1, -1); }
    else if (gl_VertexID == 2) {    screenPos.xy = vec2(-1,  1); }
    else if (gl_VertexID == 3) {    screenPos.xy = vec2( 1,  1); }

    gl_Position = screenPos;
}

/*****************************************************************************
 ***                          GEOMETRY SHADER
 *****************************************************************************/

/*****************************************************************************
 ***                          FRAGMENT SHADER
 *****************************************************************************/

// Buffer to sort fragments
uint colorList[MAX_FRAGMENTS_PER_PIXEL];
float depthList[MAX_FRAGMENTS_PER_PIXEL];

shader FSOIT(out vec4 fragColour)
{
    uint headIDX = imageLoad(headBuffer, ivec2(gl_FragCoord.xy)).r;
    uint headShadowIDX = imageLoad(headBufferShadow, ivec2(gl_FragCoord.xy)).r;
    //uvec3 shadowData = imageLoad(headBufferShadow, ivec2(gl_FragCoord.xy)).rgb;

    uint fragCount = 0;

    if (headShadowIDX != 0)
    {
        vec4 shadowColor = vec4(1,1,1,0);
        float shadowDepth = 1;

        while (headShadowIDX > 0 && fragCount < MAX_FRAGMENTS_PER_PIXEL)
        {
            Fragment frag = linkedList[headShadowIDX];

            vec4 color = unpackUnorm4x8(frag.data.x);
            float depth = uintBitsToFloat(frag.data.y);

            shadowColor.rgb = min(shadowColor.rgb, color.rgb);
            shadowColor.a = max(shadowColor.a, color.a);

            shadowDepth = min(shadowDepth, depth);

            headShadowIDX = uint(frag.data.z);
        }

        Fragment frag = linkedList[headIDX];

        depthList[0] = shadowDepth;
        colorList[0] = packUnorm4x8(shadowColor);
        fragCount = 1;
    }

    if (headIDX != 0)
    {
        Fragment frag = linkedList[headIDX];

        uint color = frag.data.x;
        float depth = uintBitsToFloat(frag.data.y);

        // Init first dummy element;
        if (fragCount == 0 && headIDX != 0)
        {
            colorList[0] = color;
            depthList[0] = depth;

            headIDX = uint(frag.data.z);

            fragCount = 1;
        }

        // Loop through fragment list
        while (headIDX > 0 && fragCount < MAX_FRAGMENTS_PER_PIXEL)
        {
            frag = linkedList[headIDX];

            uint color = frag.data.x;
            float depth = uintBitsToFloat(frag.data.y);

            uint cc = fragCount;

            while (cc > 0 && depthList[cc - 1] > depth)
            {
                depthList[cc] = depthList[cc - 1];
                colorList[cc] = colorList[cc - 1];
                cc--;
            }

            if (cc < MAX_FRAGMENTS_PER_PIXEL)
            {
                depthList[cc] = depth;
                colorList[cc] = color;
            }

            fragCount++;
            headIDX = uint(frag.data.z);
        }
    }

    if (fragCount > 0)
    {
        // 2) Blend front-to-back
        fragColour = vec4(0);

        for (int i = 0; i < fragCount; i++)
        {
            vec4 color = unpackUnorm4x8(colorList[i]);

            fragColour.rgb += (1. - fragColour.a) * color.a * color.rgb;
            fragColour.a += (1. - fragColour.a) * color.a;

        }

        if (fragColour.a != 0) { fragColour.rgb = fragColour.rgb / fragColour.a; }

#ifdef DEBUG_FRAGMENT_COUNT
        if (fragCount == 1) fragColour = vec4(0, 0, 1, 1);
        else if (fragCount <= 5) fragColour = vec4(0, 1, 0.5, 1);
        else if (fragCount <= 10) fragColour = vec4(0, 1, 0, 1);
        else if (fragCount <= 20) fragColour = vec4(0.5, 1, 0, 1);
        else if (fragCount <= 30) fragColour = vec4(1, 1, 0, 1);
        else if (fragCount <= 40) fragColour = vec4(1, 0.5, 0, 1);
        else if (fragCount >= 50) fragColour = vec4(1, 0, 0, 1);
#endif
    }
    else { discard; }
}


/*****************************************************************************
 ***                             PROGRAMS
 *****************************************************************************/

program OIT
{
    vs(430)=VSOIT();
    fs(430)=FSOIT();
};
