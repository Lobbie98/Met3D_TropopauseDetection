
/*****************************************************************************
 ***                             CONSTANTS
 *****************************************************************************/

const float EARTH_RADIUS = 6371;
const float M_PI = 3.1415926535897932384626433832795;

 /*****************************************************************************
  ***                             UNIFORMS
  *****************************************************************************/

layout(r32f, binding = 0) uniform image3D inputGrid;
layout(r32f, binding = 1) uniform image3D resultGrid;

//!TODO temporary workaround to avoid imageSize bug
uniform ivec3 dimensions;

// For default Gaussian / Box blur.
uniform int kernelSize;
uniform int stdDeviation;

uniform vec2 boundaryNW;
uniform vec2 gridResolution;

/*****************************************************************************
 ***                             INCLUDES
 *****************************************************************************/

/*****************************************************************************
 ***                             FUNCTIONS
 *****************************************************************************/

float gcDistanceEarth(const vec2 pointA, const vec2 pointB)
{
    float dlon = pointB.x - pointA.x;
    float dlat = pointB.y - pointA.y;

    float sin_dlat = sin(dlat/2.);
    float sin_dlon = sin(dlon/2.);
    float a = sin_dlat * sin_dlat
            + cos(pointA.y) * cos(pointB.y)
            * sin_dlon * sin_dlon;
    float c = 2. * asin(min(1.,sqrt(a)));
    return c * EARTH_RADIUS;
}


double degreesToRadians(double angle)
{
    return angle / 180. * M_PI;
}


// A := amplitude
float gaussian2D(const vec2 point, const vec2 mean, const vec2 stdDev)
{
    float A = 1.0f / (2 * M_PI * stdDev.x * stdDev.y);

    float distX = point.x - mean.x;
    float distY = point.y - mean.y;

    float exp1 = (distX * distX) / (2 * stdDev.x * stdDev.x);
    float exp2 = (distY * distY) / (2 * stdDev.y * stdDev.y);

    return A * exp(-exp1 - exp2);
}


float gaussian1D(const float x, const float stdDev)
{
    float denom = 2 * stdDev * stdDev;

    //float A = 1.0f / (stdDev * sqrt(2 * M_PI));
    float A = 1.0f / (stdDev * stdDev * 2 * M_PI);

    return A * exp(- (x * x) / denom);
}


/*****************************************************************************
 ***                           COMPUTE SHADER
 *****************************************************************************/

shader computeBoxBlur()
{
    //! TODO: (mk, Apr2018): Obviously imageSize returns wrong z dimension (here 1)
    //! TODO: I had no idea why this did not work, probably a bug?
    //! TODO: Investigate what causes this problem and how to fix it.
    //ivec3 dims = imageSize(inputGrid);

    int x = int(gl_GlobalInvocationID.x);
    //equivalent to int(gl_LocalInvocationID.x + gl_WorkGroupID.x * gl_WorkGroupSize.x);
    int y = int(gl_GlobalInvocationID.y);
    //equivalent to int(gl_LocalInvocationID.y + gl_WorkGroupID.y * gl_WorkGroupSize.y);
    int z = int(gl_GlobalInvocationID.z);
    //equivalent to int(gl_LocalInvocationID.z + gl_WorkGroupID.z * gl_WorkGroupSize.z);

    // Only process data if we are in the data domain. This is important as
    // the number of threads can exceed the actual size of our data.
    if (x < dimensions.x && y < dimensions.y && z < dimensions.z)
    {
        int offset = kernelSize / 2;
        int numWeights = 0;
        float value = 0;

        for (int j = int(y - offset); j <= int(y + offset); ++j)
        {
            if (j < 0 || j >= dimensions.y) { continue; }

            for (int i = int(x - offset); i <= int(x + offset); ++i)
            {
                if (i < 0 || i >= dimensions.x) { continue; }

                value += imageLoad(inputGrid, ivec3(i, j, z)).r;
                numWeights++;
            }
        }

        imageStore(resultGrid, ivec3(x, y, z), vec4(value / float(numWeights), 0, 0, 0));
    }
}


shader computeBoxBlurECMWF()
{
    //! TODO: (mk, Apr2018): Obviously imageSize returns wrong z dimension (here 1)
    //! TODO: I had no idea why this did not work, probably a bug?
    //! TODO: Investigate what causes this problem and how to fix it.
    //ivec3 dims = imageSize(inputGrid);

    int x = int(gl_GlobalInvocationID.x);
    int y = int(gl_GlobalInvocationID.y);
    int z = int(gl_GlobalInvocationID.z);

    // Only process data if we are in the data domain. This is important as
    // the number of threads can exceed the actual size of our data.
    if (x < dimensions.x && y < dimensions.y && z < dimensions.z)
    {
        float sumWeights = 8.0;

        int xN = min(x + 1, dimensions.x - 1);
        int xP = max(x - 1, 0);

        int yN = min(y + 1, dimensions.y - 1);
        int yP = max(y - 1, 0);

        float A = imageLoad(inputGrid, ivec3(x, yP, z)).r;
        float B = imageLoad(inputGrid, ivec3(x, yN, z)).r;
        float C = imageLoad(inputGrid, ivec3(xN, y, z)).r;
        float D = imageLoad(inputGrid, ivec3(xP, y, z)).r;
        float X = imageLoad(inputGrid, ivec3(x, y, z)).r;

        float result = (4 * X + A + B + C + D) / sumWeights;

        imageStore(resultGrid, ivec3(x, y, z), vec4(result, 0, 0, 0));
    }
}


shader computeBoxGaussianDistanceBlur()
{
    //! TODO: (mk, Apr2018): Obviously imageSize returns wrong z dimension (here 1)
    //! TODO: I had no idea why this did not work, probably a bug?
    //! TODO: Investigate what causes this problem and how to fix it.
    //ivec3 dims = imageSize(inputGrid);

    int x = int(gl_GlobalInvocationID.x);
    int y = int(gl_GlobalInvocationID.y);
    int z = int(gl_GlobalInvocationID.z);

    // Only process data if we are in the data domain. This is important as
    // the number of threads can exceed the actual size of our data.
    if (x < dimensions.x && y < dimensions.y && z < dimensions.z)
    {
        float lon = boundaryNW.x + x * gridResolution.x;
        float lat = boundaryNW.y + y * gridResolution.y;

        vec2 mean = vec2(degreesToRadians(lon), degreesToRadians(lat));
        vec2 stdDev = vec2(stdDeviation, stdDeviation);

        int offset = kernelSize / 2;
        float sumWeights = 0;
        float value = 0;

        for (int j = int(y - offset); j <= int(y + offset); ++j)
        {
            if (j < 0 || j >= dimensions.y) { continue; }

            for (int i = int(x - offset); i <= int(x + offset); ++i)
            {
                if (i < 0 || i >= dimensions.x) { continue; }

                float curLon = boundaryNW.x + i * gridResolution.x;
                float curLat = boundaryNW.y + j * gridResolution.y;

                vec2 curPoint = vec2(degreesToRadians(curLon), degreesToRadians(curLat));

                float distance = gcDistanceEarth(curPoint, mean);
                //float weight = gaussian2D(curPoint, mean, stdDev, 1);
                float weight = gaussian1D(distance, stdDev.x);

                value += weight * imageLoad(inputGrid, ivec3(i, j, z)).r;
                sumWeights += weight;
            }
        }

        imageStore(resultGrid, ivec3(x, y, z), vec4(value / sumWeights, 0, 0, 0));
    }
}

 /*****************************************************************************
  ***                             PROGRAMS
  *****************************************************************************/

program BoxBlur
{
    cs(430)=computeBoxBlur() : in(local_size_x = 16, local_size_y = 12, local_size_z = 8);
};

program BoxBlurECMWF
{
    cs(430)=computeBoxBlurECMWF() : in(local_size_x = 16, local_size_y = 12, local_size_z = 8);
};

program GaussianBlur
{
    cs(430)=computeBoxGaussianDistanceBlur() : in(local_size_x = 16, local_size_y = 12, local_size_z = 8);
};
