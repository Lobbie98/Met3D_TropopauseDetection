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

// Uniforms needed to gather fragments in head buffer

// Atomic fragment counter
layout(binding=0, offset=0)   uniform atomic_uint  atomicFragCount;
// Head buffer with pointer to first elements in linked list
layout(binding=0, r32ui)      uniform uimage2D     headBuffer;
layout(binding=1, r32ui)      uniform uimage2D     headBufferShadow;

// Shader buffer storage object that represents the fragment linked list
struct Fragment
{
    uvec3 data;
};

layout (std430, binding=0) buffer LinkedListBuffer
{
    Fragment linkedList[];
};

// Maximum number of fragments for the global linked list
uniform uint maxNumFragments;


// Insert new fragment into linked list and head buffer
void insertNewFragmentOIT(in vec4 color, in vec3 fragCoords, in bool shadowMode = false)
{
    //
//    uint headIDX = imageLoad(headBuffer, ivec2(fragCoords.xy)).r;
//    Fragment curFrag = linkedList[headIDX];
//    vec4 col = unpackUnorm4x8(curFrag.data.x);
//
//    if (shadowMode && headIDX > 1) { return; }

    // Increment number of total gathered fragments and obtain the last index
    uint nextFragIDX = atomicCounterIncrement(atomicFragCount);

    if (shadowMode)
    {
        // If index is smaller than the maximum number of fragments,
        // then store the fragment in the linked list and head buffer
        if (nextFragIDX < maxNumFragments)
        {
            // Get pointer to old first element and set head pointer to this (currently inserted) element
            uint oldHeadIDX = imageAtomicExchange(headBufferShadow, ivec2(fragCoords.xy), nextFragIDX);
            Fragment frag;
            // Write fragment data (color, depth, pointer to next element)
            frag.data = uvec3(packUnorm4x8(color), floatBitsToUint(fragCoords.z), oldHeadIDX);
            // Write fragment to linked list at the global index
            linkedList[nextFragIDX] = frag;
        }

//        memoryBarrier();
//        uvec3 prevData = imageLoad(headBufferShadow, ivec2(fragCoords.xy)).rgb;
//
//        vec4 shadeColor = unpackUnorm4x8(prevData.x);
//        shadeColor.rgb = min(shadeColor.rgb, color.rgb);
//        shadeColor.a = max(shadeColor.a, color.a);
//
//        memoryBarrier();
//        uvec3 data = uvec3(packUnorm4x8(shadeColor), floatBitsToUint(fragCoords.z), 1);
//        imageStore(headBufferShadow, ivec2(fragCoords.xy), uvec4(data, 0));

        return;
    }



    // If index is smaller than the maximum number of fragments,
    // then store the fragment in the linked list and head buffer
    if (nextFragIDX < maxNumFragments)
    {
        // Get pointer to old first element and set head pointer to this (currently inserted) element
        uint oldHeadIDX = imageAtomicExchange(headBuffer, ivec2(fragCoords.xy), nextFragIDX);
        Fragment frag;
        // Write fragment data (color, depth, pointer to next element)
        frag.data = uvec3(packUnorm4x8(color), floatBitsToUint(fragCoords.z), oldHeadIDX);
        // Write fragment to linked list at the global index
        linkedList[nextFragIDX] = frag;
    }
}



