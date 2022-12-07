
uniform mat4        lightMatrix;
uniform sampler2D   shadowMap;
uniform bool        inShadowMappingMode;

float computeShadowMappingFactor(in vec3 worldPos)
{
    float shadowFactor = 0;

    vec4 lightSpacePos = lightMatrix * vec4(worldPos, 1);
    vec3 lightPosNDC = lightSpacePos.xyz / lightSpacePos.w;
    lightPosNDC = lightPosNDC * 0.5 + 0.5;

    float currentDepth = lightPosNDC.z;
    //float closestDepth = texture(shadowMap, lightPosNDC.xy).r;

    const float bias = 0.001;

    vec2 texelSize = 1.0 / textureSize(shadowMap, 0);
    const int kernelSize = 2;

    for(int x = -kernelSize; x <= kernelSize; ++x)
    {
        for(int y = -kernelSize; y <= kernelSize; ++y)
        {
            float depth = texture(shadowMap,
                                  lightPosNDC.xy + vec2(x, y) * texelSize).r;
            shadowFactor += (currentDepth - bias < depth) ? 1.0 : 0.0;
        }
    }

    shadowFactor /= (kernelSize * 2.0 + 1) * (kernelSize * 2.0 + 1);

    return shadowFactor;
}