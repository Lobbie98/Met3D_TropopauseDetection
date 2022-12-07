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

#include "frontdetectionactor.h"

#include <QFileDialog>

// related third party imports
#include <netcdf>

using namespace Met3D;
using namespace netCDF;
using namespace netCDF::exceptions;

#include <log4cplus/loggingmacros.h>

#define qNaN (std::numeric_limits<float>::quiet_NaN())

/******************************************************************************
***                             PUBLIC METHODS                              ***
*******************************************************************************/


MFrontDetectionActor::MFrontDetectionActor()
        : MNWPMultiVarActor(),
          MBoundingBoxInterface(this, MBoundingBoxConnectionType::VOLUME),
          suppressUpdates(false),
          isComputing(false),
          applySettingsClickProperty(nullptr),
          autoComputeProperty(nullptr),
          autoCompute(false),
          computeOnGPUProperty(nullptr),
          computeOnGPU(false),
          inputVarGroupProperty(nullptr),
          detectionVariableIndexProperty(nullptr),
          detectionVariableIndex(0),
          detectionVariable(nullptr),
          windUVarIndexProp(nullptr),
          windUVariableIndex(0),
          windUVar(nullptr),
          windVVarIndexProp(nullptr),
          windVVariableIndex(0),
          windVVar(nullptr),
          zVarIndexProp(nullptr),
          zVariableIndex(0),
          zVar(nullptr),
          displayOptionsGroupProperty(nullptr),
          render3DFrontProperty(nullptr),
          render3DFront(false),
          render3DNCProperty(nullptr),
          render3DNC(false),
          render2DFrontProperty(nullptr),
          render2DFront(false),
          render2DNCProperty(nullptr),
          render2DNC(false),
          filterGroupProperty(nullptr),
          fpmaDistanceProperty(nullptr),
          fpmaDistanceValue_km(75.),
          transferFunctionTFPProperty(nullptr),
          transferFunctionTFP(nullptr),
          transferFunctionTFPTexUnit(-1),
          transferFunctionFSProperty(nullptr),
          transferFunctionFS(nullptr),
          transferFunctionFSTexUnit(-1),
          genericFilterGroupProperty(nullptr),
          addNormalCurveFilterProperty(nullptr),
          optionalFilterProperty(nullptr),
          transferFunctionABZProperty(nullptr),
          transferFunctionABZ(nullptr),
          transferFunctionABZTexUnit(-1),
          transferFunctionBreadthProperty(nullptr),
          transferFunctionBreadth(nullptr),
          transferFunctionBreadthTexUnit(-1),
          transferFunctionSlopeProperty(nullptr),
          transferFunctionSlope(nullptr),
          transferFunctionSlopeTexUnit(-1),
          ncOvertracingProperty(nullptr),
          ncOvertracing(0),
          fleIsoProperty(nullptr),
          fleIsoValue(0.0),
          appearanceGroupProperty(nullptr),
          shadingGroupProperty(nullptr),
          shadingModeProperty(nullptr),
          shadingMode(0),
          shadingVariableIndexProperty(nullptr),
          shadingVariableIndex(0),
          shadingVariable(nullptr),
          shadingVarModeProperty(nullptr),
          shadingVarMode(PER_100KM),
          transferFunctionShadingProperty(nullptr),
          transferFunctionShading(nullptr),
          transferFunctionShadingTexUnit(-1),
          shading2dGroupProperty(nullptr),
          frontElevationProperty(nullptr),
          frontElevation2d_hPa(925.0),
          frontsTubeRadiusProp(nullptr),
          frontsTubeRadius(0.1),
          showFrontTypesProperty(nullptr),
          showFrontTypes(false),
          showColdSideFrontProperty(nullptr),
          showColdSideFront(false),
          lightShadowGroupProperty(nullptr),
          shadowHeightProperty(nullptr),
          shadowHeight(1045),
          shadowColorProperty(nullptr),
          shadowColor(QColor(70, 70, 70, 150)),
          lightingModeProperty(nullptr),
          lightingMode(0),
          normalCurvesGroupProp(nullptr),
          trigger3DNormalCurveFilter(false),
          trigger2DNormalCurveFilter(false),
          normalCurvesTubeRadiusProp(nullptr),
          normalCurvesTubeRadius(0.05),
          ncShadingVarIndexProp(nullptr),
          ncShadingVariableIndex(0),
          ncShadingVar(nullptr),
          seedPointSpacingProp(nullptr),
          seedPointSpacing(2),
          seedPointSpacingZProp(nullptr),
          seedPointSpacingZ(100),
          ncShadingModeProperty(nullptr),
          ncShadingMode(0),
          transferFunctionNCShadingProperty(nullptr),
          transferFunctionNCShading(nullptr),
          transferFunctionNCShadingTexUnit(-1),
          useTFPTransferFunction(false),
          useFSTransferFunction(false),
          useABZTransferFunction(false),
          useSlopeTransferFunction(false),
          useBreadthTransferFunction(false),
          fleSource(nullptr),
          tfpSource(nullptr),
          abzSource(nullptr),
          partialDerivFilter(nullptr),
          frontDetection2DSource(nullptr),
          frontDetection3DSource(nullptr),
          frontDetection3DSourceGPU(nullptr),
          frontSurfaceSelection(nullptr),
          frontLineSelection(nullptr),
          frontMesh3D(nullptr),
          frontLines2D(nullptr),
          normalCurves3D(nullptr),
          normalCurves2D(nullptr),
          fronts3DSurfaceShader(nullptr),
          frontTubeShaderOIT(nullptr),
          normalCurvesTubeShaderOIT(nullptr),

          vbo3DFronts(nullptr),
          ibo3DFronts(nullptr),
          vbo3DFrontsShading(nullptr),

          vbo3DNormalCurves(nullptr),
          ibo3DNormalCurves(nullptr),
          vbo3DNormalCurvesShading(nullptr),

          vbo2DFronts(nullptr),
          ibo2DFronts(nullptr),
          vbo2DFrontsShading(nullptr),

          vbo2DNormalCurves(nullptr),
          ibo2DNormalCurves(nullptr),
          vbo2DNormalCurvesShading(nullptr),

          alphaShading3DFronts(0),
          alphaShading3DNormalCurves(0),

          alphaShading2DFronts(0),
          alphaShading2DNormalCurves(0),

          prevComputedOnGPU(false),
          recomputeAlpha3DFronts(true),
          recomputeShading3DFronts(true),
          recomputeAlpha3DNormalCurves(true),
          recomputeShading3DNormalCurves(true),
          recomputeAlpha2DFronts(true),
          recomputeShading2DFronts(true),
          recomputeAlpha2DNormalCurves(true),
          recomputeShading2DNormalCurves(true),
          restartIndex(std::numeric_limits<u_int32_t>::max())
{
    bBoxConnection =
            new MBoundingBoxConnection(this,
                                       MBoundingBoxConnectionType::VOLUME);

    // Transfer function.
    // Scan currently available actors for transfer functions. Add TFs to
    // the list displayed in the combo box of the transferFunctionProperty.
    QStringList availableTFs;
    availableTFs << "None";
    MGLResourcesManager *glRM = MGLResourcesManager::getInstance();
    foreach (MActor *mactor, glRM->getActors())
    {
        if (MTransferFunction1D *tf = dynamic_cast<MTransferFunction1D*>(mactor))
        {
            availableTFs << tf->transferFunctionName();
        }
    }

    beginInitialiseQtProperties();

    setActorType("Front Detection Actor");
    setName(getActorType());


    /**************************************************************************
                            MAIN PROPERTIES
    **************************************************************************/

    // auto compute active fronts, a front is active when its render[2/3]DFront
    // property is activated
    autoComputeProperty = addProperty(BOOL_PROPERTY, "automatic update",
                                      actorPropertiesSupGroup);
    properties->mBool()->setValue(autoComputeProperty, autoCompute);

    computeOnGPUProperty = addProperty(BOOL_PROPERTY, "compute normal curve integration on GPU",
                                      actorPropertiesSupGroup);
    properties->mBool()->setValue(computeOnGPUProperty, computeOnGPU);

    // compute active fronts, a front is active when its render[2/3]DFront
    // property is activated
    applySettingsClickProperty = addProperty(CLICK_PROPERTY, "update",
                                        actorPropertiesSupGroup);

    /**************************************************************************
                            INPUT VAR PROPERTIES
    **************************************************************************/
    inputVarGroupProperty = addProperty(GROUP_PROPERTY,
                                      "input variables",
                                      actorPropertiesSupGroup);

    detectionVariableIndexProperty = addProperty(
            ENUM_PROPERTY, "detection variable", inputVarGroupProperty);
    detectionVariableIndexProperty->setToolTip(
                "example: wet-bulb potential temperature");

    windUVarIndexProp = addProperty(
            ENUM_PROPERTY, "eastward wind (optional)", inputVarGroupProperty);
    windUVarIndexProp->setToolTip("u-wind component");

    windVVarIndexProp = addProperty(
            ENUM_PROPERTY, "northward wind (optional)", inputVarGroupProperty);
    windVVarIndexProp->setToolTip("v-wind component");

    zVarIndexProp = addProperty(
            ENUM_PROPERTY, "geopotential height (optional)", inputVarGroupProperty);
    zVarIndexProp->setToolTip("variable to calculating the slope of 3D front surfaces");

    /**************************************************************************
                            DISPLAY OPTIONS
    **************************************************************************/

    displayOptionsGroupProperty = addProperty(GROUP_PROPERTY,
                                              "display options",
                                              actorPropertiesSupGroup);

    render3DFrontProperty = addProperty(BOOL_PROPERTY, "3D fronts",
                                        displayOptionsGroupProperty);
    properties->mBool()->setValue(render3DFrontProperty, render3DFront);

    render3DNCProperty = addProperty(BOOL_PROPERTY, "3D normal curves",
                                        displayOptionsGroupProperty);
    properties->mBool()->setValue(render3DNCProperty, render3DNC);

    render2DFrontProperty = addProperty(BOOL_PROPERTY, "2D fronts",
                                        displayOptionsGroupProperty);
    properties->mBool()->setValue(render2DFrontProperty, render2DFront);

    render2DNCProperty = addProperty(BOOL_PROPERTY, "2D normal curves",
                                        displayOptionsGroupProperty);
    properties->mBool()->setValue(render2DNCProperty, render2DNC);

    /**************************************************************************
                            FILTER OPTIONS
    **************************************************************************/
    filterGroupProperty = addProperty(GROUP_PROPERTY,
                                      "filter options",
                                      actorPropertiesSupGroup);

    fpmaDistanceProperty = addProperty(DECORATEDDOUBLE_PROPERTY,
                                        "5-point-mean-axis distance",
                                         filterGroupProperty);
    properties->setDDouble(fpmaDistanceProperty, fpmaDistanceValue_km,
                           0., 1000., 1, 5, " km");
    fpmaDistanceProperty->setToolTip(
                "Set the distance between the 5-point-mean-axis. \n"
                "The distance ideally corresponds to the smoothing radius.");

    transferFunctionTFPProperty = addProperty(ENUM_PROPERTY,
                                              "transfer function TFP filter",
                                              filterGroupProperty);
    properties->mEnum()->setEnumNames(transferFunctionTFPProperty, availableTFs);
    transferFunctionTFPProperty->setToolTip(
                "Set alpha value of transfer function for fuzzy filtering.");

    transferFunctionFSProperty = addProperty(ENUM_PROPERTY,
                                              "transfer function frontal strength filter",
                                              filterGroupProperty);
    properties->mEnum()->setEnumNames(transferFunctionFSProperty, availableTFs);
    transferFunctionFSProperty->setToolTip(
                "Set alpha values of transfer function for fuzzy filtering.");

    genericFilterGroupProperty = addProperty(GROUP_PROPERTY,
                                             "generic filter",
                                             filterGroupProperty);

    addNormalCurveFilterProperty = addProperty(CLICK_PROPERTY,
                                               "add normal curve filter",
                                               genericFilterGroupProperty);
    addNormalCurveFilterProperty->setToolTip(
                "At least one generic filter is recommended, \n"
                "which filters according to the frontal strength \n"
                "(change per 100 km of detection variable).\n"
                "Replacement of ABZ filter of Hewson (1998)");

    optionalFilterProperty = addProperty(GROUP_PROPERTY,
                                         "optional filter",
                                         filterGroupProperty);

    transferFunctionABZProperty = addProperty(ENUM_PROPERTY,
                                              "transfer function abz filter",
                                              optionalFilterProperty);
    properties->mEnum()->setEnumNames(transferFunctionABZProperty, availableTFs);
    transferFunctionABZProperty->setToolTip(
                "Set alpha value of transfer function for fuzzy filtering. \n"
                "ABZ filter after Hewson (1998), estimates frontal strength");

    transferFunctionBreadthProperty = addProperty(ENUM_PROPERTY,
                                              "transfer function breadth filter",
                                              optionalFilterProperty);
    properties->mEnum()->setEnumNames(transferFunctionBreadthProperty, availableTFs);
    transferFunctionBreadthProperty->setToolTip(
                "Set alpha value of transfer function for fuzzy filtering. \n"
                "Breadth of frontal zone");

    transferFunctionSlopeProperty = addProperty(ENUM_PROPERTY,
                                              "transfer function filter slope (3d)",
                                              optionalFilterProperty);
    properties->mEnum()->setEnumNames(transferFunctionSlopeProperty, availableTFs);
    transferFunctionSlopeProperty->setToolTip(
                "Set alpha value of transfer function for fuzzy filtering. \n"
                "Slope of frontal frontal surface (only for 3D Fronts)");

    ncOvertracingProperty = addProperty(INT_PROPERTY,
                                        "normal curve over-tracing (km)",
                                        optionalFilterProperty);
    properties->setInt(ncOvertracingProperty, ncOvertracing, 0, 1000, 10);
    ncOvertracingProperty->setToolTip(
                "If normal curves reached the cold side of front, \n"
                "normal curve integration will stop. However, sometimes the \n"
                "front location equation detects some islands within the \n"
                "frontal zone and normal curve tracing will stop. Here it is \n"
                "possible to trace normal curves over these islands. Set the \n"
                "length in km for how long you would like to continue the normal \n"
                "curve tracing over such islands. If the frontal zone \n"
                "continues after these islands, normal curve tracing will \n"
                "ignored these islands.");

    fleIsoProperty = addProperty(SCIENTIFICDOUBLE_PROPERTY, "iso value",
                                 optionalFilterProperty);
    properties->setSciDouble(fleIsoProperty, fleIsoValue,
                             -1.00, 1.00, 2, 0.01);
    fleIsoProperty->setToolTip("Iso value of front location equation. \n"
                               "Default and recommendation value = 0.0");

    /**************************************************************************
                            FRONT APPEARANCE PROPERTIES
    **************************************************************************/
    appearanceGroupProperty = addProperty(GROUP_PROPERTY,
                                                 "front appearance",
                                                 actorPropertiesSupGroup);
    shadingGroupProperty = addProperty(GROUP_PROPERTY,
                                     "shading",
                                     appearanceGroupProperty);

    QStringList shadingModes;
    shadingModes << "Shading Var" << "Pressure" << "TFP" << "ABZ"
                 << "Slope" << "Breadth of Frontal Zone";

    shadingModeProperty = addProperty(ENUM_PROPERTY, "shading mode", shadingGroupProperty);
    properties->mEnum()->setEnumNames(shadingModeProperty, shadingModes);
    properties->mEnum()->setValue(shadingModeProperty, shadingMode);

    shadingVariableIndexProperty = addProperty(
            ENUM_PROPERTY, "shading variable", shadingGroupProperty);

    QStringList shadingVarModes;
    shadingVarModes << filterTypeToString(ABSOLUTE_CHANGE)
                    << filterTypeToString(PER_100KM)
                    << filterTypeToString(VALUE_AT_FRONT);

    shadingVarModeProperty = addProperty(ENUM_PROPERTY, "shading variable mode",
                                         shadingGroupProperty);
    properties->mEnum()->setEnumNames(shadingVarModeProperty, shadingVarModes);
    properties->mEnum()->setValue(shadingVarModeProperty, shadingVarMode);


    transferFunctionShadingProperty = addProperty(ENUM_PROPERTY,
                                                  "transfer function shading",
                                                  shadingGroupProperty);
    properties->mEnum()->setEnumNames(transferFunctionShadingProperty,
                                      availableTFs);


    shading2dGroupProperty = addProperty(GROUP_PROPERTY,
                                         "appearance 2D",
                                         appearanceGroupProperty);

    frontElevationProperty = addProperty(DECORATEDDOUBLE_PROPERTY,
                                        "elevation fronts (2d)",
                                         shading2dGroupProperty);
    properties->setDDouble(frontElevationProperty, frontElevation2d_hPa,
                           1., 1050., 2, 0.1, " hPa");

    frontsTubeRadiusProp = addProperty(
            DOUBLE_PROPERTY, "fronts tube radius", shading2dGroupProperty);
    properties->setDouble(frontsTubeRadiusProp, frontsTubeRadius, 0.01, 1.0, 2, 0.01);

    appearanceGroupProperty->addSubProperty(bBoxConnection->getProperty());

    showFrontTypesProperty = addProperty(BOOL_PROPERTY, "show warm/cold fronts (3d)",
                                         appearanceGroupProperty);
    properties->mBool()->setValue(showFrontTypesProperty, showFrontTypes);


    showColdSideFrontProperty = addProperty(BOOL_PROPERTY, "show cold-side fronts (3d)",
                                            appearanceGroupProperty);
    properties->mBool()->setValue(showColdSideFrontProperty, showColdSideFront);

    lightShadowGroupProperty = addProperty(GROUP_PROPERTY,
                                           "light and shadow",
                                           appearanceGroupProperty);


    shadowColorProperty = addProperty(COLOR_PROPERTY, "shadow color", lightShadowGroupProperty);
    properties->mColor()->setValue(shadowColorProperty, shadowColor);

    shadowHeightProperty = addProperty(DOUBLE_PROPERTY, "shadow elevation", lightShadowGroupProperty);
    properties->setDouble(shadowHeightProperty, shadowHeight, 0, 1050, 0, 5);

    QStringList modesLst;
    modesLst << "double-sided" << "single-sided" << "double-sided + headlight" << "single-sided + headlight";
    lightingModeProperty = addProperty(ENUM_PROPERTY, "lighting mode", lightShadowGroupProperty);
    properties->mEnum()->setEnumNames(lightingModeProperty, modesLst);
    properties->mEnum()->setValue(lightingModeProperty, lightingMode);


    /**************************************************************************
                            Normal curves
    **************************************************************************/

    normalCurvesGroupProp = addProperty(GROUP_PROPERTY, "normal curve settings",
                                        actorPropertiesSupGroup);

    normalCurvesTubeRadiusProp = addProperty(
            DOUBLE_PROPERTY, "normal curves tube radius", normalCurvesGroupProp);
    properties->setDouble(normalCurvesTubeRadiusProp, normalCurvesTubeRadius,
                          0.01, 1.0, 2, 0.01);

    seedPointSpacingProp = addProperty(DOUBLE_PROPERTY,
                                          "seed point spacing",
                                       normalCurvesGroupProp);
    properties->setDouble(seedPointSpacingProp, seedPointSpacing,
                          0.1, 100, 1, 0.1);

    seedPointSpacingZProp = addProperty(INT_PROPERTY,
                                       "seed point spacing (z) [hPa]",
                                        normalCurvesGroupProp);
    properties->setInt(seedPointSpacingZProp, seedPointSpacingZ,
                       10, 500, 10);

    QStringList ncShadingList;
    ncShadingList << "curve length (km)" << "shading value" << "avg. difference"
                  << "total difference";
    ncShadingModeProperty = addProperty(ENUM_PROPERTY, "nc shading mode",
                                        normalCurvesGroupProp);
    properties->mEnum()->setEnumNames(ncShadingModeProperty, ncShadingList);
    properties->mEnum()->setValue(ncShadingModeProperty, ncShadingMode);

    ncShadingVarIndexProp = addProperty(
            ENUM_PROPERTY, "nc shading variable", normalCurvesGroupProp);

    transferFunctionNCShadingProperty = addProperty(ENUM_PROPERTY,
                                                "transfer function nc shading",
                                                    normalCurvesGroupProp);
    properties->mEnum()->setEnumNames(transferFunctionNCShadingProperty,
                                      availableTFs);


    connect(glRM, SIGNAL(actorCreated(MActor*)),
            SLOT(onActorCreated(MActor*)));
    connect(glRM, SIGNAL(actorDeleted(MActor*)),
            SLOT(onActorDeleted(MActor*)));
    connect(glRM, SIGNAL(actorRenamed(MActor*,QString)),
            SLOT(onActorRenamed(MActor*,QString)));

    endInitialiseQtProperties();
}


MFrontDetectionActor::~MFrontDetectionActor()
{
//    if (frontSurfaceSelection)
//    {
//        curFilter->releaseData(frontSurfaceSelection);
//    }
}

/******************************************************************************
***                            PUBLIC METHODS                               ***
*******************************************************************************/

#define SHADER_VERTEX_ATTRIBUTE 0

void MFrontDetectionActor::reloadShaderEffects()
{
    LOG4CPLUS_DEBUG(mlog, "loading shader programs");

    beginCompileShaders(4);

    compileShadersFromFileWithProgressDialog(
                fronts3DSurfaceShader, "src/glsl/fronts3d_geometry.fx.glsl");
    compileShadersFromFileWithProgressDialog(
                frontTubeShaderOIT,
                "src/glsl/tube_frontlines_OIT.fx.glsl");
    compileShadersFromFileWithProgressDialog(
                normalCurvesTubeShaderOIT,
                "src/glsl/tube_normalcurves_OIT.fx.glsl");

    endCompileShaders();
}


void MFrontDetectionActor::saveConfiguration(QSettings *settings)
{
    MNWPMultiVarActor::saveConfiguration(settings);

    settings->beginGroup(getSettingsID());

    MBoundingBoxInterface::saveConfiguration(settings);


    settings->setValue("render3DFront", render3DFront);
    settings->setValue("render2DFront", render2DFront);
    settings->setValue("render3DNC", render3DNC);
    settings->setValue("render2DNC", render2DNC);
    settings->setValue("autoCompute", autoCompute);
    settings->setValue("computeOnGPU", computeOnGPU);
    settings->setValue("frontElevation2d_hPa", frontElevation2d_hPa);
    settings->setValue("detectionVariableIndex", detectionVariableIndex);
    settings->setValue("fleIsoValue", fleIsoValue);
    settings->setValue("fpmaDistanceValue", fpmaDistanceValue_km);
    settings->setValue("useTFPTransferFunction", useTFPTransferFunction);
    settings->setValue("useFSTransferFunction", useFSTransferFunction);
    settings->setValue("useABZTransferFunction", useABZTransferFunction);
    settings->setValue("useSlopeTransferFunction", useSlopeTransferFunction);
    settings->setValue("useBreadthTransferFunction", useBreadthTransferFunction);
    settings->setValue("ncOvertracing", ncOvertracing);

    settings->beginWriteArray("filterTfsFrontDetection",
                              normalCurvesFilterSetList.size());
    int index = 0;
    foreach (auto filter, normalCurvesFilterSetList)
    {
        settings->setArrayIndex(index);
        settings->setValue("enabled", filter->enabled);
        settings->setValue("variable", filter->variableIndex);
        settings->setValue("transferFunctionProperty", properties->getEnumItem(
                           filter->transferFunctionProperty));
        settings->setValue("type", filterTypeToString(filter->type));
        index++;
    }
    settings->endArray(); // "filterTfsFrontDetection"

    settings->setValue("showFrontTypes", showFrontTypes);
    settings->setValue("showColdSideFront", showColdSideFront);
    settings->setValue("shadingMode", shadingMode);
    settings->setValue("shadingVarMode", filterTypeToString(shadingVarMode));
    settings->setValue("frontsTubeRadius", frontsTubeRadius);

    settings->setValue("shadingVariableIndex", shadingVariableIndex);
    settings->setValue("windUVariableIndex", windUVariableIndex);
    settings->setValue("windVVariableIndex", windVVariableIndex);
    settings->setValue("zVariableIndex", zVariableIndex);

    settings->setValue("normalCurvesTubeRadius", normalCurvesTubeRadius);
    settings->setValue("ncShadingVariableIndex", ncShadingVariableIndex);
    settings->setValue("seedPointSpacing", seedPointSpacing);
    settings->setValue("seedPointSpacingZ", seedPointSpacingZ);
    settings->setValue("ncShadingMode", ncShadingMode);

    settings->setValue("transferFunctionShading",
                       properties->getEnumItem(
                               transferFunctionShadingProperty));
    settings->setValue("transferFunctionTFPProperty",
                       properties->getEnumItem(
                               transferFunctionTFPProperty));
    settings->setValue("transferFunctionFSProperty",
                       properties->getEnumItem(
                               transferFunctionFSProperty));
    settings->setValue("transferFunctionABZProperty",
                       properties->getEnumItem(
                               transferFunctionABZProperty));
    settings->setValue("transferFunctionSlopeProperty",
                       properties->getEnumItem(
                               transferFunctionSlopeProperty));
    settings->setValue("transferFunctionBreadthProperty",
                       properties->getEnumItem(
                               transferFunctionBreadthProperty));
    settings->setValue("transferFunctionNCShadingProperty",
                       properties->getEnumItem(
                               transferFunctionNCShadingProperty));

    settings->setValue("shadowColor", shadowColor);
    settings->setValue("shadowHeight", shadowHeight);
    settings->setValue("lightingMode", lightingMode);
    settings->endGroup();
}


void MFrontDetectionActor::loadConfiguration(QSettings *settings)
{
    MNWPMultiVarActor::loadConfiguration(settings);

    suppressUpdates = true;

    settings->beginGroup(getSettingsID());

    MBoundingBoxInterface::loadConfiguration(settings);


    // ********************* MAIN PROPERTIES *********************
    render3DFront = settings->value("render3DFront").toBool();
    properties->mBool()->setValue(render3DFrontProperty, render3DFront);

    render2DFront = settings->value("render2DFront").toBool();
    properties->mBool()->setValue(render2DFrontProperty, render2DFront);

    render3DNC = settings->value("render3DNC").toBool();
    properties->mBool()->setValue(render3DNCProperty, render3DNC);

    render2DNC = settings->value("render2DNC").toBool();
    properties->mBool()->setValue(render2DNCProperty, render2DNC);

    autoCompute = settings->value("autoCompute").toBool();
    properties->mBool()->setValue(autoComputeProperty, autoCompute);

    computeOnGPU = settings->value("computeOnGPU").toBool();
    properties->mBool()->setValue(computeOnGPUProperty, computeOnGPU);

    frontElevation2d_hPa = settings->value("frontElevation2d_hPa").toDouble();
    properties->mDDouble()->setValue(frontElevationProperty,
                                     frontElevation2d_hPa);

    detectionVariableIndex = settings->value("detectionVariableIndex").toInt();
    properties->mInt()->setValue(detectionVariableIndexProperty,
                                 detectionVariableIndex);

    fpmaDistanceValue_km = settings->value("fpmaDistanceValue").toDouble();
    properties->mDDouble()->setValue(fpmaDistanceProperty,
                                     fpmaDistanceValue_km);


    // ********************* APPEARANCE PROPERTIES *********************

    fleIsoValue = settings->value("fleIsoValue").toDouble();
    properties->mDouble()->setValue(fleIsoProperty, fleIsoValue);

    useTFPTransferFunction = settings->value("useTFPTransferFunction").toBool();

    useFSTransferFunction = settings->value("useFSTransferFunction").toBool();

    useABZTransferFunction = settings->value("useABZTransferFunction").toBool();

    useSlopeTransferFunction =
            settings->value("useSlopeTransferFunction").toBool();

    useBreadthTransferFunction =
            settings->value("useBreadthTransferFunction").toBool();

    ncOvertracing = settings->value("ncOvertracingProperty").toInt();
    properties->mInt()->setValue(ncOvertracingProperty, ncOvertracing);

    int numFilterSet = settings->beginReadArray("filterTfsFrontDetection");
    for (int i = 0; i < numFilterSet; i++)
    {
        settings->setArrayIndex(i);
        bool enabled = settings->value("enabled", false).toBool();
        int varIndex = settings->value("variable", 0).toInt();
        QString tfName = settings->value("transferFunctionProperty", "").toString();
        filterType type = stringToFilterType(settings->value("type").toString());
        loadNormalCurvesFilter(i, enabled, varIndex, tfName, type);
    }
    settings->endArray();

    showFrontTypes = settings->value("showFrontTypes").toBool();
    properties->mBool()->setValue(showColdSideFrontProperty, showFrontTypes);

    shadingMode = settings->value("shadingMode").toInt();
    properties->mInt()->setValue(shadingModeProperty, shadingMode);

    shadingVarMode = stringToFilterType(
                settings->value("shadingVarMode").toString());
    properties->mInt()->setValue(shadingVarModeProperty, shadingVarMode);

    frontsTubeRadius = settings->value("frontsTubeRadius", 0.3).toFloat();
    properties->mDouble()->setValue(frontsTubeRadiusProp, frontsTubeRadius);

    shadingVariableIndex = settings->value("shadingVariableIndex").toInt();
    properties->mInt()->setValue(shadingVariableIndexProperty,
                                 shadingVariableIndex);

    windUVariableIndex = settings->value("windUVariableIndex").toInt();
    properties->mInt()->setValue(windUVarIndexProp,
                                 windUVariableIndex);

    windVVariableIndex = settings->value("windVVariableIndex").toInt();
    properties->mInt()->setValue(windVVarIndexProp,
                                 windVVariableIndex);

    zVariableIndex = settings->value("zVariableIndex").toInt();
    properties->mInt()->setValue(zVarIndexProp,
                                 zVariableIndex);

//    normalCurvesEnabled = settings->value("normalCurvesEnabled").toBool();
//    properties->mBool()->setValue(normalCurvesEnabledProp, normalCurvesEnabled);

    normalCurvesTubeRadius = settings->value("normalCurvesTubeRadius", 0.2).toFloat();
    properties->mDouble()->setValue(normalCurvesTubeRadiusProp, normalCurvesTubeRadius);

    ncShadingVariableIndex = settings->value("ncShadingVariableIndex").toInt();
    properties->mInt()->setValue(ncShadingVarIndexProp,
                                 ncShadingVariableIndex);

    seedPointSpacing = settings->value("seedPointSpacing").toFloat();
    properties->mDouble()->setValue(seedPointSpacingProp, seedPointSpacing);

    seedPointSpacingZ = settings->value("seedPointSpacingZ", 100).toInt();
    properties->mInt()->setValue(seedPointSpacingZProp, seedPointSpacingZ);

    ncShadingMode = settings->value("ncShadingMode").toInt();
    properties->mEnum()->setValue(ncShadingModeProperty, ncShadingMode);

    QString tfName = settings->value("transferFunctionShading", "None").toString();
    if (!setTransferFunction(transferFunctionShadingProperty, tfName))
    {
    }

    tfName = settings->value("transferFunctionTFPProperty", "None").toString();
    if (!setTransferFunction(transferFunctionTFPProperty, tfName))
    {
    }

    tfName = settings->value("transferFunctionFSProperty", "None").toString();
    if (!setTransferFunction(transferFunctionFSProperty, tfName))
    {
    }

    tfName = settings->value("transferFunctionABZProperty", "None").toString();
    if (!setTransferFunction(transferFunctionABZProperty, tfName))
    {
    }

    tfName = settings->value("transferFunctionSlopeProperty", "None").toString();
    if (!setTransferFunction(transferFunctionSlopeProperty, tfName))
    {
    }

    tfName = settings->value("transferFunctionBreadthProperty", "None").toString();
    if (!setTransferFunction(transferFunctionBreadthProperty, tfName))
    {
    }

    tfName = settings->value("transferFunctionNCShadingProperty", "None").toString();
    if (!setTransferFunction(transferFunctionNCShadingProperty, tfName))
    {
    }

    setTransferFunctionFromProperty(transferFunctionShadingProperty,
                                    &transferFunctionShading);
    setTransferFunctionFromProperty(transferFunctionTFPProperty,
                                    &transferFunctionTFP);
    setTransferFunctionFromProperty(transferFunctionFSProperty,
                                    &transferFunctionFS);
    setTransferFunctionFromProperty(transferFunctionABZProperty,
                                    &transferFunctionABZ);
    setTransferFunctionFromProperty(transferFunctionSlopeProperty,
                                    &transferFunctionSlope);
    setTransferFunctionFromProperty(transferFunctionBreadthProperty,
                                    &transferFunctionBreadth);
    setTransferFunctionFromProperty(transferFunctionNCShadingProperty,
                                    &transferFunctionNCShading);

    shadowColor = settings->value("shadowColor",
                                  QColor(70, 70, 70)).value<QColor>();
    properties->mColor()->setValue(shadowColorProperty, shadowColor);

    shadowHeight = settings->value("shadowHeight", 1045.0).toFloat();
    properties->mDouble()->setValue(shadowHeightProperty, shadowHeight);

    lightingMode = settings->value("lightingMode").toInt();
    properties->mEnum()->setValue(lightingModeProperty, lightingMode);

    if (!variables.empty())
    {
        detectionVariable = static_cast<MNWP3DVolumeActorVariable*>(
                variables.at(detectionVariableIndex));
        //potTemperatureVar = static_cast<MNWP3DVolumeActorVariable*>(
        //        variables.at(potTemperatureVariableIndex));
        shadingVariable = static_cast<MNWP3DVolumeActorVariable*>(
                variables.at(shadingVariableIndex));
        windUVar = static_cast<MNWP3DVolumeActorVariable*>(
                variables.at(windUVariableIndex));
        windVVar = static_cast<MNWP3DVolumeActorVariable*>(
                variables.at(windVVariableIndex));
        zVar = static_cast<MNWP3DVolumeActorVariable*>(
                variables.at(zVariableIndex));
        ncShadingVar = static_cast<MNWP3DVolumeActorVariable*>(
                variables.at(ncShadingVariableIndex));
    }


    settings->endGroup();

    suppressUpdates = false;
    //requestFrontLines();
    emitActorChangedSignal();
}


bool MFrontDetectionActor::setTransferFunction(QtProperty* tfProp, const QString &tfName)
{
    QStringList tfNames = properties->mEnum()->enumNames(tfProp);
    int tfIndex = tfNames.indexOf(tfName);

    if (tfIndex >= 0)
    {
        properties->mEnum()->setValue(tfProp, tfIndex);
        return true;
    }

    // Set transfer function property to "None".
    properties->mEnum()->setValue(tfProp, 0);

    return false; // The given tf name could not be found.
}


void MFrontDetectionActor::setTransferFunctionFromProperty(
        QtProperty* tfProp, MTransferFunction1D** transferFunc)
{
    MGLResourcesManager *glRM = MGLResourcesManager::getInstance();

    QString tfName = properties->getEnumItem(tfProp);

    if (tfName == "None")
    {
        *transferFunc = nullptr;
        return;
    }

    // Find the selected transfer function in the list of actors from the
    // resources manager. Not very efficient, but works well enough for the
    // small number of actors at the moment..
    foreach (MActor *actor, glRM->getActors())
    {
        if (MTransferFunction1D *tf = dynamic_cast<MTransferFunction1D *>(actor))
        {
            if (tf->transferFunctionName() == tfName)
            {
                *transferFunc = tf;
                return;
            }
        }
    }
}


MFrontDetectionActor::NormalCurvesFilterSettings::
NormalCurvesFilterSettings(MFrontDetectionActor* a,
                           uint8_t filterNumber,
                           QList<QString> variableNames,
                           bool enabled,
                           int variableIndex,
                           MTransferFunction1D *transferFunction,
                           filterType type,
                           int textureUnitTransferFunction)
    : enabled(enabled),
      variableIndex(variableIndex),
      transferFunction(transferFunction),
      textureUnitTransferFunction(textureUnitTransferFunction),
      varNameList(variableNames),
      type(type)
{
    //MActor *a = this;
    MQtProperties *properties = a->getQtProperties();

    a->beginInitialiseQtProperties();

    QString propertyTitle = QString("nc filter #%1").arg(filterNumber);

    groupProperty = a->addProperty(GROUP_PROPERTY, propertyTitle);

    enabledProperty = a->addProperty(BOOL_PROPERTY, "enabled", groupProperty);
    properties->mBool()->setValue(enabledProperty, enabled);

    variableProperty = a->addProperty(ENUM_PROPERTY,
                              "filter variable",
                              groupProperty);

    int tmpVarIndex = variableIndex;

    properties->mEnum()->setEnumNames(
                variableProperty,
                varNameList);
    properties->mEnum()->setValue(
                variableProperty,
                tmpVarIndex);

    // Transfer function.
    // Scan currently available actors for transfer functions. Add TFs to
    // the list displayed in the combo box of the transferFunctionProperty.
    QStringList availableTFs;
    availableTFs << "None";
    MGLResourcesManager *glRM = MGLResourcesManager::getInstance();
    foreach (MActor *mactor, glRM->getActors())
    {
        if (MTransferFunction1D *tf = dynamic_cast<MTransferFunction1D*>(mactor))
        {
            availableTFs << tf->transferFunctionName();
        }
    }

    int tfIndex = -1;
    if (transferFunction != nullptr)
    {
        for (int i = 0; i < availableTFs.size(); i++)
        {
            if (availableTFs[i] == transferFunction->transferFunctionName())
            {
                tfIndex = i;
                break;
            }
        }
    }
    transferFunctionProperty = a->addProperty(ENUM_PROPERTY,
                                                  "transfer function",
                                                  groupProperty);
    properties->mEnum()->setEnumNames(transferFunctionProperty, availableTFs);
    properties->mEnum()->setValue(transferFunctionProperty, tfIndex);
    transferFunctionProperty->setToolTip("This transfer function is used "
                                         "for mapping frontal strength to "
                                         "frontal's colour.");
    QStringList typeModes;

    typeModes << filterTypeToString(ABSOLUTE_CHANGE)
              << filterTypeToString(PER_100KM)
              << filterTypeToString(VALUE_AT_FRONT);

    typeProperty = a->addProperty(ENUM_PROPERTY, "filterType", groupProperty);
    properties->mEnum()->setEnumNames(typeProperty, typeModes);
    properties->mEnum()->setValue(typeProperty, type);

    removeProperty = a->addProperty(CLICK_PROPERTY, "remove", groupProperty);

    a->endInitialiseQtProperties();
}


void MFrontDetectionActor::addNormalCurvesFilter(
        int8_t filterNumber)
{
    QList<QString> variableList;
    for (int vi = 0; vi < variables.size(); vi++)
    {
        MNWPActorVariable* var = variables.at(vi);

        variableList << var->variableName;
    }
    int varIndex = -1;
    if (varNameList.size() > 0)
    {
        varIndex = 0;
    }
    NormalCurvesFilterSettings *ncFilterSettings =
            new NormalCurvesFilterSettings(this, filterNumber, variableList,
                                           false, varIndex);
    normalCurvesFilterSetList.append(ncFilterSettings);
    genericFilterGroupProperty->addSubProperty(
    normalCurvesFilterSetList.back()->groupProperty);

    if (ncFilterSettings->enabled)
    {
        bool isConnected;
        isConnected = connect(ncFilterSettings->transferFunction,
                              SIGNAL(actorChanged()),
                              this, SLOT(recomputeFrontsAlpha));
        assert(isConnected);
    }
}


void MFrontDetectionActor::loadNormalCurvesFilter(
        int8_t filterNumber, bool enabled, int varIndex, QString tfName,
        filterType type)
{
    QList<QString> variableList;
    for (int vi = 0; vi < variables.size(); vi++)
    {
        MNWPActorVariable* var = variables.at(vi);

        variableList << var->variableName;
    }
    if (varIndex > varNameList.size())
    {
        varIndex = -1;
    }

    MTransferFunction1D *tfFilter = nullptr;
    MGLResourcesManager *glRM = MGLResourcesManager::getInstance();
    foreach (MActor *actor, glRM->getActors())
    {
        if (MTransferFunction1D *tf = dynamic_cast<MTransferFunction1D*>(actor))
        {
            if (tf->transferFunctionName() == tfName)
            {
                tfFilter = tf;

                NormalCurvesFilterSettings *ncFilterSettings =
                        new NormalCurvesFilterSettings(this, filterNumber, variableList,
                                                       enabled, varIndex, tfFilter, type);
                normalCurvesFilterSetList.append(ncFilterSettings);
                genericFilterGroupProperty->addSubProperty(
                normalCurvesFilterSetList.back()->groupProperty);

                if (ncFilterSettings->enabled)
                {
                    bool isConnected;
                    isConnected = connect(ncFilterSettings->transferFunction,
                                        SIGNAL(actorChanged()),
                                        this, SLOT(recomputeFrontsAlpha()));
                    assert(isConnected);
                }
            }
        }
    }
}


bool MFrontDetectionActor::removeNormalCurvesFilter(
        NormalCurvesFilterSettings *filter)
{
    int index = normalCurvesFilterSetList.indexOf(filter);
    // Don't ask for confirmation during application start.
    if (MSystemManagerAndControl::getInstance()->applicationIsInitialized())
    {
        QMessageBox yesNoBox;
        yesNoBox.setWindowTitle("Delete filter");
        yesNoBox.setText(QString("Do you really want to delete "
                                 "filter #%1").arg(index + 1));
        yesNoBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        yesNoBox.setDefaultButton(QMessageBox::No);

        if (yesNoBox.exec() != QMessageBox::Yes)
        {
            return false;
        }
    }

    // disconnect filter connection to this actor
    disconnect(filter->transferFunction, 0, this, 0);

    // Redraw is only needed if filter sets to delete are enabled.
    bool needsRedraw = filter->enabled;
    MQtProperties *properties = this->getQtProperties();
    properties->mBool()->setValue(filter->enabledProperty,
                                  false);
    genericFilterGroupProperty->removeSubProperty(
                filter->groupProperty);
    normalCurvesFilterSetList.remove(index);

    QString text = filter->groupProperty->propertyName();
    text = text.mid(0, text.indexOf("#") + 1);
    // Rename filter sets after deleted filter to
    // close gap left by deleted filter sets.
    for (; index < normalCurvesFilterSetList.size(); index++)
    {
        filter->groupProperty->setPropertyName(
                    text + QString::number(index + 1));
    }
    return needsRedraw;
}


MNWPActorVariable* MFrontDetectionActor::createActorVariable(
        const MSelectableDataSource& dataSource)
{
    MNWP3DVolumeActorVariable *newVar = new MNWP3DVolumeActorVariable(this);

    newVar->dataSourceID = dataSource.dataSourceID;
    newVar->levelType = dataSource.levelType;
    newVar->variableName = dataSource.variableName;

    return newVar;
}


const QList<MVerticalLevelType> MFrontDetectionActor::supportedLevelTypes()
{
    return (QList<MVerticalLevelType>()
            << HYBRID_SIGMA_PRESSURE_3D << PRESSURE_LEVELS_3D << AUXILIARY_PRESSURE_3D);
}


void MFrontDetectionActor::onBoundingBoxChanged()
{
    labels.clear();
    if (suppressActorUpdates())
    {
        return;
    }
    // Switching to no bounding box only needs a redraw, but no recomputation
    // because it disables rendering of the actor.
    if (bBoxConnection->getBoundingBox() == nullptr)
    {
        emitActorChangedSignal();
        return;
    }

    emitActorChangedSignal();
}


/******************************************************************************
***                               PUBLIC SLOTS                              ***
*******************************************************************************/


void MFrontDetectionActor::asynchronous3DFrontsAvailable(MDataRequest request)
{
    LOG4CPLUS_DEBUG(mlog, "MFrontDetectionActor::asynchronous3DFrontsAvailable");
    if (frontMesh3D && !prevComputedOnGPU)
    {
        frontMesh3D->releaseIndexBuffer();
        frontMesh3D->releaseVertexBuffer();

        normalCurves3D->releaseVertexBuffer();
        normalCurves3D->releaseIndexBuffer();

        frontDetection3DSource->releaseData(frontSurfaceSelection);
    }
    if (frontMesh3D && prevComputedOnGPU)
    {
        frontMesh3D->releaseIndexBuffer();
        frontMesh3D->releaseVertexBuffer();

        frontDetection3DSourceGPU->releaseData(frontSurfaceSelection);
    }

    if (computeOnGPU)
    {
        frontSurfaceSelection = frontDetection3DSourceGPU->getData(request);
        prevComputedOnGPU = true;
    }
    else
    {
        frontSurfaceSelection = frontDetection3DSource->getData(request);
        prevComputedOnGPU = false;
    }

    frontMesh3D = frontSurfaceSelection->getTriangleMeshSelection();
    normalCurves3D = frontSurfaceSelection->getNormalCurvesSelection();

    num3DFrontVertices = frontMesh3D->getNumVertices();

    vbo3DFronts = frontMesh3D->getVertexBuffer();
    ibo3DFronts = frontMesh3D->getIndexBuffer();

    applySettingsClickProperty->setEnabled(true);
    actorPropertiesSupGroup->setEnabled(true);

    suppressUpdates = false;
    isComputing = false;

    trigger3DNormalCurveFilter = true;

    recomputeAlpha3DFronts = true;
    recomputeShading3DFronts = true;
    recomputeAlpha3DNormalCurves = true;
    recomputeShading3DNormalCurves = true;

    emitActorChangedSignal();
}


void MFrontDetectionActor::asynchronous2DFrontsAvailable(
        MDataRequest request)
{
    if (frontLines2D)
    {
        frontLines2D->releaseIndexBuffer();
        frontLines2D->releaseVertexBuffer();

        normalCurves2D->releaseVertexBuffer();
        normalCurves2D->releaseIndexBuffer();

        frontDetection2DSource->releaseData(frontLineSelection);
    }
    frontLineSelection = frontDetection2DSource->getData(request);

    frontLines2D = frontLineSelection->getLineSelection();
    normalCurves2D = frontLineSelection->getNormalCurvesSelection();

    vbo2DFronts = frontLines2D->getVertexBuffer();
    ibo2DFronts = frontLines2D->getIndexBuffer();

    num2DFrontVertices = frontLines2D->getNumVertices();

    applySettingsClickProperty->setEnabled(true);
    actorPropertiesSupGroup->setEnabled(true);

    isComputing = false;
    suppressUpdates = false;

    trigger2DNormalCurveFilter = true;

    recomputeAlpha2DFronts = true;
    recomputeShading2DFronts = true;
    recomputeAlpha2DNormalCurves = true;
    recomputeShading2DNormalCurves = true;

    emitActorChangedSignal();

}


void MFrontDetectionActor::onAddActorVariable(MNWPActorVariable *var)
{
    varNameList << var->variableName;

    // Temporarily save variable indices.
    int tmpVarIndex = detectionVariableIndex;
    int tmpShadingVarIndex = shadingVariableIndex;
    int tmpWindUVarIndex = windUVariableIndex;
    int tmpWindVVarIndex = windVVariableIndex;
    int tmpZVarIndex = zVariableIndex;
    int tmpNCShadingVarIndex = ncShadingVariableIndex;

    // Update enum lists.
    properties->mEnum()->setEnumNames(detectionVariableIndexProperty, varNameList);
    properties->mEnum()->setValue(detectionVariableIndexProperty, tmpVarIndex);
    properties->mEnum()->setEnumNames(shadingVariableIndexProperty, varNameList);
    properties->mEnum()->setValue(shadingVariableIndexProperty, tmpShadingVarIndex);
    properties->mEnum()->setEnumNames(windUVarIndexProp, varNameList);
    properties->mEnum()->setValue(windUVarIndexProp, tmpWindUVarIndex);
    properties->mEnum()->setEnumNames(windVVarIndexProp, varNameList);
    properties->mEnum()->setValue(windVVarIndexProp, tmpWindVVarIndex);
    properties->mEnum()->setEnumNames(zVarIndexProp, varNameList);
    properties->mEnum()->setValue(zVarIndexProp, tmpZVarIndex);
    properties->mEnum()->setEnumNames(ncShadingVarIndexProp, varNameList);
    properties->mEnum()->setValue(ncShadingVarIndexProp, tmpNCShadingVarIndex);

    foreach (auto filter, normalCurvesFilterSetList)
    {
        filter->varNameList << var->variableName;
        int tmpVarIndex = filter->variableIndex;

        properties->mEnum()->setEnumNames(
                    filter->variableProperty,
                    filter->varNameList);
        properties->mEnum()->setValue(
                    filter->variableProperty,
                    tmpVarIndex);
    }

    refreshEnumsProperties(nullptr);
}


void MFrontDetectionActor::onDeleteActorVariable(MNWPActorVariable *var)
{
    int i = variables.indexOf(var);

    // Update variableIndex and shadingVariableIndex if these point to
    // the removed variable or to one with a lower index.
    if (i <= detectionVariableIndex)
    {
        detectionVariableIndex = std::max(-1, detectionVariableIndex - 1);
    }
    if (i <= shadingVariableIndex)
    {
        shadingVariableIndex = std::max(-1, shadingVariableIndex - 1);
    }
    if (i <= windUVariableIndex)
    {
        windUVariableIndex = std::max(-1, windUVariableIndex - 1);
    }
    if (i <= windVVariableIndex)
    {
        windVVariableIndex = std::max(-1, windVVariableIndex - 1);
    }
    if (i <= zVariableIndex)
    {
        zVariableIndex = std::max(-1, zVariableIndex - 1);
    }
    if (i <= ncShadingVariableIndex)
    {
        ncShadingVariableIndex = std::max(-1, ncShadingVariableIndex - 1);
    }

    // Temporarily save variable indices.
    int tmpVarIndex = detectionVariableIndex;
    int tmpShadingVarIndex = shadingVariableIndex;
    int tmpWindUVarIndex = windUVariableIndex;
    int tmpWindVVarIndex = windVVariableIndex;
    int tmpZVarIndex = zVariableIndex;
    int tmpNCShadingVarIndex = ncShadingVariableIndex;

    // Remove the variable name from the enum lists.
    varNameList.removeAt(i);

    // Update enum lists.
    properties->mEnum()->setEnumNames(detectionVariableIndexProperty, varNameList);
    properties->mEnum()->setValue(detectionVariableIndexProperty, tmpVarIndex);
    properties->mEnum()->setEnumNames(shadingVariableIndexProperty, varNameList);
    properties->mEnum()->setValue(shadingVariableIndexProperty, tmpShadingVarIndex);
    properties->mEnum()->setEnumNames(windUVarIndexProp, varNameList);
    properties->mEnum()->setValue(windUVarIndexProp, tmpWindUVarIndex);
    properties->mEnum()->setEnumNames(windVVarIndexProp, varNameList);
    properties->mEnum()->setValue(windVVarIndexProp, tmpWindVVarIndex);
    properties->mEnum()->setEnumNames(zVarIndexProp, varNameList);
    properties->mEnum()->setValue(zVarIndexProp, tmpZVarIndex);
    properties->mEnum()->setEnumNames(ncShadingVarIndexProp, varNameList);
    properties->mEnum()->setValue(ncShadingVarIndexProp, tmpNCShadingVarIndex);

    foreach (auto* filter, normalCurvesFilterSetList)
    {
        if(i <= filter->variableIndex)
        {
            filter->variableIndex =
                    std::max(-1, filter->variableIndex - 1);
        }
        int tempVarIndex = filter->variableIndex;

        // Remove the variable name from the enum lists.
        filter->varNameList.removeAt(i);

        // Update enum lists.
        properties->mEnum()->setEnumNames(
                    filter->variableProperty,
                    filter->varNameList);

        properties->mEnum()->setValue(
                    filter->variableProperty,
                    tempVarIndex);

    }

    refreshEnumsProperties(var);
}


void MFrontDetectionActor::onChangeActorVariable(MNWPActorVariable *var)
{
    int varIndex = variables.indexOf(var);
    foreach (auto filter, normalCurvesFilterSetList)
    {
        filter->varNameList.replace(varIndex, var->variableName);
        int tmpVarIndex = filter->variableIndex;

        properties->mEnum()->setEnumNames(
                    filter->variableProperty,
                    filter->varNameList);
        properties->mEnum()->setValue(
                    filter->variableProperty,
                    tmpVarIndex);
    }
    refreshEnumsProperties(nullptr);
}


void MFrontDetectionActor::onActorCreated(MActor *actor)
{
    // If the new actor is a transfer function, add it to the list of
    // available transfer functions.
    if (MTransferFunction1D *tf = dynamic_cast<MTransferFunction1D *>(actor))
    {
        // Don't render while the properties are being updated.
        enableEmissionOfActorChangedSignal(false);

        int indexS = properties->mEnum()->value(transferFunctionShadingProperty);
        int indexTFP = properties->mEnum()->value(transferFunctionTFPProperty);
        int indexFS = properties->mEnum()->value(transferFunctionFSProperty);
        int indexABZ = properties->mEnum()->value(transferFunctionABZProperty);
        int indexSlope = properties->mEnum()->value(transferFunctionSlopeProperty);
        int indexBreadth = properties->mEnum()->value(transferFunctionBreadthProperty);
        int indexNCShading = properties->mEnum()->value(transferFunctionNCShadingProperty);

        QStringList availableTFs = properties->mEnum()->enumNames(transferFunctionTFPProperty);
        availableTFs << tf->transferFunctionName();

        properties->mEnum()->setEnumNames(transferFunctionShadingProperty, availableTFs);
        properties->mEnum()->setEnumNames(transferFunctionTFPProperty, availableTFs);
        properties->mEnum()->setEnumNames(transferFunctionFSProperty, availableTFs);
        properties->mEnum()->setEnumNames(transferFunctionABZProperty, availableTFs);
        properties->mEnum()->setEnumNames(transferFunctionSlopeProperty, availableTFs);
        properties->mEnum()->setEnumNames(transferFunctionBreadthProperty, availableTFs);
        properties->mEnum()->setEnumNames(transferFunctionNCShadingProperty, availableTFs);

        properties->mEnum()->setValue(transferFunctionShadingProperty, indexS);
        properties->mEnum()->setValue(transferFunctionTFPProperty, indexTFP);
        properties->mEnum()->setValue(transferFunctionFSProperty, indexFS);
        properties->mEnum()->setValue(transferFunctionABZProperty, indexABZ);
        properties->mEnum()->setValue(transferFunctionSlopeProperty, indexSlope);
        properties->mEnum()->setValue(transferFunctionBreadthProperty, indexBreadth);
        properties->mEnum()->setValue(transferFunctionNCShadingProperty, indexNCShading);

        foreach (auto* filter, normalCurvesFilterSetList)
        {
            availableTFs.clear();
            int index = properties->mEnum()->value(
                        filter->transferFunctionProperty);
            availableTFs = properties->mEnum()->enumNames(
                        filter->transferFunctionProperty);
            availableTFs << tf->transferFunctionName();
            properties->mEnum()->setEnumNames(filter->transferFunctionProperty,
                        availableTFs);
            properties->mEnum()->setValue(
                        filter->transferFunctionProperty, index);
        }

        enableEmissionOfActorChangedSignal(true);
    }
}


void MFrontDetectionActor::onActorDeleted(MActor *actor)
{
    // If the deleted actor is a transfer function, remove it from the list of
    // available transfer functions.
    if (MTransferFunction1D *tf = dynamic_cast<MTransferFunction1D *>(actor))
    {
        enableEmissionOfActorChangedSignal(false);

        int indexS = properties->mEnum()->value(transferFunctionShadingProperty);
        int indexTFP = properties->mEnum()->value(transferFunctionTFPProperty);
        int indexFS = properties->mEnum()->value(transferFunctionFSProperty);
        int indexABZ = properties->mEnum()->value(transferFunctionABZProperty);
        int indexSlope = properties->mEnum()->value(transferFunctionSlopeProperty);
        int indexBreadth = properties->mEnum()->value(transferFunctionBreadthProperty);
        int indexNCShading = properties->mEnum()->value(transferFunctionNCShadingProperty);

        QStringList availableTFs =
                properties->mEnum()->enumNames(transferFunctionTFPProperty);
        availableTFs << tf->transferFunctionName();

        // If the deleted transfer function is currently connected to this
        // variable, set current transfer function to "None" (index 0).
        if (availableTFs.at(indexS) == tf->getName()) { indexS = 0; }
        if (availableTFs.at(indexTFP) == tf->getName()) { indexTFP = 0; }
        if (availableTFs.at(indexFS) == tf->getName()) { indexFS = 0; }
        if (availableTFs.at(indexSlope) == tf->getName()) { indexSlope = 0; }
        if (availableTFs.at(indexBreadth) == tf->getName()) { indexBreadth = 0; }
        if (availableTFs.at(indexNCShading) == tf->getName()) { indexNCShading = 0; }

        availableTFs.removeOne(tf->getName());

        properties->mEnum()->setEnumNames(transferFunctionShadingProperty, availableTFs);
        properties->mEnum()->setEnumNames(transferFunctionTFPProperty, availableTFs);
        properties->mEnum()->setEnumNames(transferFunctionFSProperty, availableTFs);
        properties->mEnum()->setEnumNames(transferFunctionABZProperty, availableTFs);
        properties->mEnum()->setEnumNames(transferFunctionSlopeProperty, availableTFs);
        properties->mEnum()->setEnumNames(transferFunctionBreadthProperty, availableTFs);
        properties->mEnum()->setEnumNames(transferFunctionNCShadingProperty, availableTFs);

        properties->mEnum()->setValue(transferFunctionShadingProperty, indexS);
        properties->mEnum()->setEnumNames(transferFunctionTFPProperty, availableTFs);
        properties->mEnum()->setValue(transferFunctionFSProperty, indexFS);
        properties->mEnum()->setValue(transferFunctionABZProperty, indexABZ);
        properties->mEnum()->setValue(transferFunctionSlopeProperty, indexSlope);
        properties->mEnum()->setValue(transferFunctionBreadthProperty, indexBreadth);
        properties->mEnum()->setValue(transferFunctionNCShadingProperty, indexNCShading);

        foreach (auto* filter, normalCurvesFilterSetList)
        {
            QString tFName = properties->getEnumItem(filter->transferFunctionProperty);
            availableTFs = properties->mEnum()->enumNames(
                        filter->transferFunctionProperty);
            availableTFs.removeOne(tf->getName());
            // Get the current index of the transfer function selected. If the
            // transfer function is the one to be deleted, the selection is set to
            // 'None'.
            int index = availableTFs.indexOf(tFName);

            properties->mEnum()->setEnumNames(filter->transferFunctionProperty,
                        availableTFs);
            properties->mEnum()->setValue(
                        filter->transferFunctionProperty, index);
        }

        enableEmissionOfActorChangedSignal(true);
    }
}


void MFrontDetectionActor::onActorRenamed(MActor *actor,
                                            QString oldName)
{
    // If the renamed actor is a transfer function, change its name in the list
    // of available transfer functions.
    if (MTransferFunction1D *tf = dynamic_cast<MTransferFunction1D *>(actor))
    {
        // Don't render while the properties are being updated.
        enableEmissionOfActorChangedSignal(false);

        int indexS = properties->mEnum()->value(transferFunctionShadingProperty);
        int indexTFP = properties->mEnum()->value(transferFunctionTFPProperty);
        int indexFS = properties->mEnum()->value(transferFunctionFSProperty);
        int indexABZ = properties->mEnum()->value(transferFunctionABZProperty);
        int indexSlope = properties->mEnum()->value(transferFunctionSlopeProperty);
        int indexBreadth = properties->mEnum()->value(transferFunctionBreadthProperty);
        int indexNCShading = properties->mEnum()->value(transferFunctionNCShadingProperty);

        QStringList availableTFs = properties->mEnum()->enumNames(transferFunctionTFPProperty);

        // Replace affected entry.
        availableTFs[availableTFs.indexOf(oldName)] = tf->getName();

        properties->mEnum()->setEnumNames(transferFunctionShadingProperty, availableTFs);
        properties->mEnum()->setEnumNames(transferFunctionTFPProperty, availableTFs);
        properties->mEnum()->setEnumNames(transferFunctionFSProperty, availableTFs);
        properties->mEnum()->setEnumNames(transferFunctionABZProperty, availableTFs);
        properties->mEnum()->setEnumNames(transferFunctionSlopeProperty, availableTFs);
        properties->mEnum()->setEnumNames(transferFunctionBreadthProperty, availableTFs);
        properties->mEnum()->setEnumNames(transferFunctionNCShadingProperty, availableTFs);

        properties->mEnum()->setValue(transferFunctionShadingProperty, indexS);
        properties->mEnum()->setValue(transferFunctionTFPProperty, indexTFP);
        properties->mEnum()->setValue(transferFunctionFSProperty, indexFS);
        properties->mEnum()->setValue(transferFunctionABZProperty, indexABZ);
        properties->mEnum()->setValue(transferFunctionSlopeProperty, indexSlope);
        properties->mEnum()->setValue(transferFunctionBreadthProperty, indexBreadth);
        properties->mEnum()->setValue(transferFunctionNCShadingProperty, indexNCShading);

        foreach (auto* filter, normalCurvesFilterSetList)
        {
            int index = properties->mEnum()->value(
                        filter->transferFunctionProperty);
            availableTFs = properties->mEnum()->enumNames(
                        filter->transferFunctionProperty);

            // Replace affected entry.
            availableTFs[availableTFs.indexOf(oldName)] = tf->getName();

            properties->mEnum()->setEnumNames(
                        filter->transferFunctionProperty,
                                              availableTFs);
            properties->mEnum()->setValue(
                        filter->transferFunctionProperty, index);
        }

        enableEmissionOfActorChangedSignal(true);
    }
}


void MFrontDetectionActor::dataFieldChangedEvent()
{
    if (autoCompute)
    {
        if (render2DFront) triggerAsynchronous2DFrontsRequest();
        if (render3DFront) triggerAsynchronous3DFrontsRequest();
    }

}


void MFrontDetectionActor::updateShadow()
{
    emitActorChangedSignal();
}


/******************************************************************************
***                             PRIVATE METHODS                             ***
*******************************************************************************/


void MFrontDetectionActor::initializeActorResources()
{
    // Parent initialisation.
    MNWPMultiVarActor::initializeActorResources();

    MGLResourcesManager *glRM = MGLResourcesManager::getInstance();

    if (!variables.empty())
    {
        detectionVariable = static_cast<MNWP3DVolumeActorVariable*>(
                variables.at(detectionVariableIndex));
        shadingVariable = static_cast<MNWP3DVolumeActorVariable*>(
                variables.at(shadingVariableIndex));
        windUVar = static_cast<MNWP3DVolumeActorVariable*>(
                variables.at(windUVariableIndex));
        windVVar = static_cast<MNWP3DVolumeActorVariable*>(
                variables.at(windVVariableIndex));
        zVar = static_cast<MNWP3DVolumeActorVariable*>(
                variables.at(zVariableIndex));
        ncShadingVar = static_cast<MNWP3DVolumeActorVariable*>(
                variables.at(ncShadingVariableIndex));
    }

    varNameList.clear();
    for (int vi = 0; vi < variables.size(); vi++)
    {
        MNWPActorVariable* var = variables.at(vi);
        varNameList << var->variableName;
    }

    properties->mEnum()->setEnumNames(detectionVariableIndexProperty, varNameList);
    properties->mEnum()->setValue(detectionVariableIndexProperty, detectionVariableIndex);
    properties->mEnum()->setEnumNames(shadingVariableIndexProperty, varNameList);
    properties->mEnum()->setValue(shadingVariableIndexProperty, shadingVariableIndex);
    properties->mEnum()->setEnumNames(windUVarIndexProp, varNameList);
    properties->mEnum()->setValue(windUVarIndexProp, windUVariableIndex);
    properties->mEnum()->setEnumNames(windVVarIndexProp, varNameList);
    properties->mEnum()->setValue(windVVarIndexProp, windVVariableIndex);
    properties->mEnum()->setEnumNames(zVarIndexProp, varNameList);
    properties->mEnum()->setValue(zVarIndexProp, zVariableIndex);
    properties->mEnum()->setEnumNames(ncShadingVarIndexProp, varNameList);
    properties->mEnum()->setValue(ncShadingVarIndexProp, ncShadingVariableIndex);

    bool loadShaders = false;

    loadShaders |= glRM->generateEffectProgram("fronts3d_geometry",
                                               fronts3DSurfaceShader);
    loadShaders |= glRM->generateEffectProgram("tube_frontlines_OIT",
                                               frontTubeShaderOIT);
    loadShaders |= glRM->generateEffectProgram("tube_normalcurves_OIT",
                                               normalCurvesTubeShaderOIT);
    if (loadShaders) reloadShaderEffects();

    // Create filter / source pipeline
    MSystemManagerAndControl *sysMC = MSystemManagerAndControl::getInstance();
    MAbstractScheduler *scheduler = sysMC->getScheduler("MultiThread");
    MAbstractMemoryManager *memoryManager = sysMC->getMemoryManager("NWP");

    fleSource = std::make_shared<MFrontLocationEquationSource>();
    fleSource->setScheduler(scheduler);
    fleSource->setMemoryManager(memoryManager);

    tfpSource = std::make_shared<MThermalFrontParameterSource>();
    tfpSource->setScheduler(scheduler);
    tfpSource->setMemoryManager(memoryManager);

    abzSource = std::make_shared<MAdjacentBaroclinicZoneSource>();
    abzSource->setScheduler(scheduler);
    abzSource->setMemoryManager(memoryManager);

    partialDerivFilter = std::make_shared<MPartialDerivativeFilter>();
    partialDerivFilter->setScheduler(scheduler);
    partialDerivFilter->setMemoryManager(memoryManager);

    frontDetection2DSource = std::make_shared<MFrontDetection2DSource>();
    frontDetection2DSource->setScheduler(scheduler);
    frontDetection2DSource->setMemoryManager(memoryManager);

    frontDetection3DSource = std::make_shared<MFrontDetection3DSource>();
    frontDetection3DSource->setScheduler(scheduler);
    frontDetection3DSource->setMemoryManager(memoryManager);

    frontDetection3DSourceGPU = std::make_shared<MFrontDetection3DSourceGPU>();
    frontDetection3DSourceGPU->setScheduler(scheduler);
    frontDetection3DSourceGPU->setMemoryManager(memoryManager);

    transferFunctionShadingTexUnit = assignTextureUnit();
    transferFunctionTFPTexUnit = assignTextureUnit();
    transferFunctionFSTexUnit = assignTextureUnit();
    transferFunctionABZTexUnit = assignTextureUnit();
    transferFunctionSlopeTexUnit = assignTextureUnit();
    transferFunctionBreadthTexUnit = assignTextureUnit();
    transferFunctionNCShadingTexUnit = assignTextureUnit();

    sampleSubroutines.resize(MVerticalLevelType::SIZE_LEVELTYPES);
    normalCompSubroutines.resize(MVerticalLevelType::SIZE_LEVELTYPES);

    bool connected = connect(frontDetection3DSource.get(),
            SIGNAL(dataRequestCompleted(MDataRequest)),
            this, SLOT(asynchronous3DFrontsAvailable(MDataRequest)));
    assert(connected);

    connected = connect(frontDetection3DSourceGPU.get(),
            SIGNAL(dataRequestCompleted(MDataRequest)),
            this, SLOT(asynchronous3DFrontsAvailable(MDataRequest)));
    assert(connected);

    connected = connect(frontDetection2DSource.get(),
            SIGNAL(dataRequestCompleted(MDataRequest)),
            this, SLOT(asynchronous2DFrontsAvailable(MDataRequest)));

    assert(connected);


    sampleSubroutines[PRESSURE_LEVELS_3D]
            << "samplePressureLevel"
            << "pressureLevelGradient";

    sampleSubroutines[HYBRID_SIGMA_PRESSURE_3D]
            << "sampleHybridLevel"
            << "hybridLevelGradient";

    sampleSubroutines[AUXILIARY_PRESSURE_3D]
            << "sampleAuxiliaryPressure"
            << "auxiliaryPressureGradient";

    normalCompSubroutines[PRESSURE_LEVELS_3D]
            << "samplePressureLevel"
            << "pressureLevelGradient";

    normalCompSubroutines[HYBRID_SIGMA_PRESSURE_3D]
            << "sampleHybridLevel"
            << "hybridLevelGradient";

    normalCompSubroutines[AUXILIARY_PRESSURE_3D]
            << "sampleAuxiliaryPressure"
            << "auxiliaryPressureGradient";
}


void MFrontDetectionActor::onQtPropertyChanged(QtProperty *property)
{
    if (suppressUpdates) { return; }

    // TODO: how to handle auto computation?
    // Parent signal processing.
    MNWPMultiVarActor::onQtPropertyChanged(property);

    // ********************* MAIN PROPERTIES *********************
    if (property == autoComputeProperty)
    {
        autoCompute = properties->mBool()->value(autoComputeProperty);
        emitActorChangedSignal();
    }

    if (property == computeOnGPUProperty)
    {
        computeOnGPU = properties->mBool()->value(computeOnGPUProperty);
        emitActorChangedSignal();
    }

    if (property == render3DFrontProperty)
    {
        render3DFront = properties->mBool()->value(render3DFrontProperty);
        emitActorChangedSignal();
    }

    if (property == render3DNCProperty)
    {
        render3DNC = properties->mBool()->value(render3DNCProperty);
        trigger3DNormalCurveFilter = true;
        emitActorChangedSignal();
    }

    else if (property == render2DFrontProperty)
    {
        render2DFront = properties->mBool()->value(render2DFrontProperty);
        emitActorChangedSignal();
    }

    if (property == render2DNCProperty)
    {
        render2DNC = properties->mBool()->value(render2DNCProperty);
        trigger2DNormalCurveFilter = true;
        emitActorChangedSignal();
    }

    else if (property == frontElevationProperty)
    {
        frontElevation2d_hPa = properties->mDDouble()
                ->value(frontElevationProperty);
        emitActorChangedSignal();
    }

    else if (property == applySettingsClickProperty)
    {
        if (render3DFront) triggerAsynchronous3DFrontsRequest();
        if (render2DFront) triggerAsynchronous2DFrontsRequest();
        emitActorChangedSignal();
    }

    else if (property == detectionVariableIndexProperty)
    {
        detectionVariableIndex = properties->mEnum()
                ->value(detectionVariableIndexProperty);
        if (detectionVariableIndex < 0) return;

        if (detectionVariableIndex >= variables.size())
        {
            detectionVariableIndex = variables.size() - 1;
            properties->mEnum()->setValue(detectionVariableIndexProperty,
                                          detectionVariableIndex);
        }

        detectionVariable = static_cast<MNWP3DVolumeActorVariable*>(
                variables[detectionVariableIndex]);
        emitActorChangedSignal();
    }

    else if (property == fpmaDistanceProperty)
    {
        fpmaDistanceValue_km = properties->mDDouble()
                ->value(fpmaDistanceProperty);
        emitActorChangedSignal();
    }

    // ********************* APPEARANCE PROPERTIES ********************

    else if (property == fleIsoProperty)
    {
        fleIsoValue = properties->mSciDouble()->value(fleIsoProperty);
        emitActorChangedSignal();
    }

    else if (property == ncOvertracingProperty)
    {
        ncOvertracing = properties->mInt()->value(ncOvertracingProperty);
        emitActorChangedSignal();
    }

    else if (property == frontsTubeRadiusProp)
    {
        frontsTubeRadius = properties->mDouble()->value(frontsTubeRadiusProp);
        emitActorChangedSignal();
    }

    else if (property == showColdSideFrontProperty
             || property == showFrontTypesProperty)
    {
        showColdSideFront = properties
                ->mBool()->value(showColdSideFrontProperty);
        showFrontTypes = properties->mBool()->value(showFrontTypesProperty);
        emitActorChangedSignal();
    }

    else if (property == shadingModeProperty)
    {
        shadingMode = properties->mEnum()->value(shadingModeProperty);
        if (shadingMode != 0)
        {
            shadingVarModeProperty->setEnabled(false);
            shadingVariableIndexProperty->setEnabled(false);
        }
        else
        {
            shadingVarModeProperty->setEnabled(true);
            shadingVariableIndexProperty->setEnabled(true);
        }
        recomputeShading3DFronts = true;
        recomputeShading3DNormalCurves = true;
        recomputeShading2DFronts = true;
        recomputeShading2DNormalCurves = true;
        emitActorChangedSignal();
    }

    else if (property == shadingVarModeProperty)
    {
        shadingVarMode = stringToFilterType(
                    properties->getEnumItem(shadingVarModeProperty));
        recomputeShading3DFronts = true;
        recomputeShading3DNormalCurves = true;
        recomputeShading2DFronts = true;
        recomputeShading2DNormalCurves = true;
        emitActorChangedSignal();
    }

    else if (property == shadingVariableIndexProperty)
    {
        shadingVariableIndex = properties->mEnum()->value(shadingVariableIndexProperty);
        if (shadingVariableIndex < 0) return;

        if (shadingVariableIndex >= variables.size())
        {
            shadingVariableIndex = variables.size() - 1;
            properties->mEnum()->setValue(shadingVariableIndexProperty, shadingVariableIndex);
        }

        shadingVariable = static_cast<MNWP3DVolumeActorVariable*>(
                variables[shadingVariableIndex]);
        recomputeShading3DFronts = true;
        recomputeShading2DFronts = true;
        emitActorChangedSignal();
    }

    else if (property == windUVarIndexProp)
    {
        windUVariableIndex = properties->mEnum()->value(windUVarIndexProp);
        if (windUVariableIndex < 0) return;

        if (windUVariableIndex >= variables.size())
        {
            windUVariableIndex = variables.size() - 1;
            properties->mEnum()->setValue(windUVarIndexProp, windUVariableIndex);
        }
        windUVar = static_cast<MNWP3DVolumeActorVariable*>(
                variables.at(windUVariableIndex));
        emitActorChangedSignal();
    }

    else if (property == windVVarIndexProp)
    {
        windVVariableIndex = properties->mEnum()->value(windVVarIndexProp);
        if (windVVariableIndex < 0) return;

        if (windVVariableIndex >= variables.size())
        {
            windVVariableIndex = variables.size() - 1;
            properties->mEnum()->setValue(windVVarIndexProp, windVVariableIndex);
        }
        windVVar = static_cast<MNWP3DVolumeActorVariable*>(
                variables.at(windVVariableIndex));
        emitActorChangedSignal();
    }

    else if (property == zVarIndexProp)
    {
        zVariableIndex = properties->mEnum()->value(zVarIndexProp);
        if (zVariableIndex < 0) return;

        if (zVariableIndex >= variables.size())
        {
            zVariableIndex = variables.size() - 1;
            properties->mEnum()->setValue(zVarIndexProp, zVariableIndex);
        }
        zVar = static_cast<MNWP3DVolumeActorVariable*>(
                variables.at(zVariableIndex));
        emitActorChangedSignal();
    }

    else if (property == normalCurvesTubeRadiusProp)
    {
        normalCurvesTubeRadius = properties->mDouble()->value(normalCurvesTubeRadiusProp);
        emitActorChangedSignal();
    }


    else if (property == ncShadingVarIndexProp)
    {
        ncShadingVariableIndex = properties->mEnum()->value(ncShadingVarIndexProp);
        if (ncShadingVariableIndex < 0) return;

        if (ncShadingVariableIndex >= variables.size())
        {
            ncShadingVariableIndex = variables.size() - 1;
            properties->mEnum()->setValue(ncShadingVarIndexProp, ncShadingVariableIndex);
        }
        ncShadingVar = static_cast<MNWP3DVolumeActorVariable*>(
                variables.at(ncShadingVariableIndex));

        trigger2DNormalCurveFilter = true;
        trigger3DNormalCurveFilter = true;

        recomputeShading2DNormalCurves = true;
        recomputeShading3DNormalCurves = true;
        emitActorChangedSignal();
    }

    else if (property == seedPointSpacingProp)
    {
        seedPointSpacing = properties->mDouble()->value(seedPointSpacingProp);
        trigger2DNormalCurveFilter = true;
        trigger3DNormalCurveFilter = true;
        emitActorChangedSignal();
    }

    else if (property == seedPointSpacingZProp)
    {
        seedPointSpacingZ = properties->mInt()->value(seedPointSpacingZProp);
        trigger3DNormalCurveFilter = true;
        emitActorChangedSignal();
    }

    else if (property == ncShadingModeProperty)
    {
        ncShadingMode = properties->mEnum()->value(ncShadingModeProperty);
        recomputeShading2DNormalCurves = true;
        recomputeShading3DNormalCurves = true;
        emitActorChangedSignal();
    }

    else if (property == transferFunctionShadingProperty
             || property == transferFunctionTFPProperty
             || property == transferFunctionFSProperty
             || property == transferFunctionABZProperty
             || property == transferFunctionSlopeProperty
             || property == transferFunctionBreadthProperty
             || property == transferFunctionNCShadingProperty)
    {
        setTransferFunctionFromProperty(transferFunctionShadingProperty,
                                        &transferFunctionShading);
        setTransferFunctionFromProperty(transferFunctionTFPProperty,
                                        &transferFunctionTFP);
        setTransferFunctionFromProperty(transferFunctionFSProperty,
                                        &transferFunctionFS);
        setTransferFunctionFromProperty(transferFunctionABZProperty,
                                        &transferFunctionABZ);
        setTransferFunctionFromProperty(transferFunctionSlopeProperty,
                                        &transferFunctionSlope);
        setTransferFunctionFromProperty(transferFunctionBreadthProperty,
                                        &transferFunctionBreadth);
        setTransferFunctionFromProperty(transferFunctionNCShadingProperty,
                                        &transferFunctionNCShading);
        useTFPTransferFunction = (transferFunctionTFP != nullptr);
        useFSTransferFunction = (transferFunctionFS != nullptr);
        useABZTransferFunction = (transferFunctionABZ != nullptr);
        useSlopeTransferFunction = (transferFunctionSlope != nullptr);
        useBreadthTransferFunction = (transferFunctionBreadth != nullptr);

        if (transferFunctionShading != nullptr)

        if (suppressActorUpdates())
        {
            return;
        }

        emitActorChangedSignal();
    }

    else if (property == shadowColorProperty
             || property == shadowHeightProperty
             || property == lightingModeProperty)
    {
        shadowColor = properties->mColor()->value(shadowColorProperty);
        shadowHeight = properties->mDouble()->value(shadowHeightProperty);
        lightingMode = properties->mEnum()->value(lightingModeProperty);
        emitActorChangedSignal();
    }

    else if (property == addNormalCurveFilterProperty)
    {
        enableEmissionOfActorChangedSignal(false);
        int filterNumber =
                normalCurvesFilterSetList.size();
        addNormalCurvesFilter(filterNumber + 1);
        enableEmissionOfActorChangedSignal(true);
        emitActorChangedSignal();
    }

    //NormalCurvesFilterSettings *ncFilterSet = nullptr;
    foreach (NormalCurvesFilterSettings* filter,
             normalCurvesFilterSetList)
    {
        if (property == filter->enabledProperty)
        {
            filter->enabled = properties->mBool()->value(
                        filter->enabledProperty);
            if (filter->enabled)
            {
                bool isConnected;
                isConnected = connect(filter->transferFunction,
                                    SIGNAL(actorChanged()),
                                    this, SLOT(recomputeFrontsAlpha()));
                //assert(isConnected);
            }
            else
            {
                disconnect(filter->transferFunction, 0, this, 0);
            }
            recomputeAlpha3DFronts = true;
            recomputeAlpha3DNormalCurves = true;
            recomputeAlpha2DFronts = true;
            recomputeAlpha2DNormalCurves = true;
            emitActorChangedSignal();
        }
        else if (property == filter->variableProperty)
        {
            filter->variableIndex = properties->mEnum()
                    ->value(filter->variableProperty);
            recomputeAlpha3DFronts = true;
            recomputeAlpha3DNormalCurves = true;
            recomputeAlpha2DFronts = true;
            recomputeAlpha2DNormalCurves = true;
            emitActorChangedSignal();
        }
        else if (property == filter->transferFunctionProperty)
        {
            if (filter->enabled)
            {
                disconnect(filter->transferFunction, 0, this, 0);
            }
            setNCFilterTransferFunctionFromProperty(filter);
            if (filter->enabled)
            {
                bool isConnected;
                isConnected = connect(filter->transferFunction,
                                    SIGNAL(actorChanged()),
                                    this, SLOT(recomputeFrontsAlpha()));
                //assert(isConnected);
            }
            recomputeAlpha3DFronts = true;
            recomputeAlpha2DFronts = true;
            recomputeAlpha2DNormalCurves = true;
            recomputeAlpha3DNormalCurves = true;
            emitActorChangedSignal();
        }
        else if (property == filter->typeProperty)
        {
            filter->type = stringToFilterType(
                        properties->getEnumItem(filter->typeProperty));
            recomputeAlpha3DFronts = true;
            recomputeAlpha3DNormalCurves = true;
            recomputeAlpha2DFronts = true;
            recomputeAlpha2DNormalCurves = true;
            emitActorChangedSignal();
        }
        else if (property == filter->removeProperty)
        {
            bool redraw = removeNormalCurvesFilter(filter);
            if (redraw)
            {
                recomputeAlpha3DFronts = true;
                recomputeAlpha3DNormalCurves = true;
                recomputeAlpha2DFronts = true;
                recomputeAlpha2DNormalCurves = true;
            }
            emitActorChangedSignal();
        }
    }
}


void MFrontDetectionActor::recomputeFrontsAlpha()
{
    recomputeAlpha3DFronts = true;
    recomputeAlpha3DNormalCurves = true;
    recomputeAlpha2DFronts = true;
    recomputeAlpha2DNormalCurves = true;
    emitActorChangedSignal();
}

// unused, since use OIT render methods
void MFrontDetectionActor::renderToCurrentContext(
        MSceneViewGLWidget *sceneView)
{
    Q_UNUSED(sceneView);
}


// for translucent objects, this method is called twice
void MFrontDetectionActor::renderTransparencyToCurrentContext(
        MSceneViewGLWidget *sceneView)
{
    // calling render methods with OIT
    render3DFrontSurfaceOIT(sceneView);
    render3DNormalCurvesOIT(sceneView);
    render2DFrontTubesOIT(sceneView);
    render2DNormalCurvesOIT(sceneView);
}


// for translucent objects, this method is called twice
void MFrontDetectionActor::render3DFrontSurfaceOIT(
        MSceneViewGLWidget *sceneView)
{
    LOG4CPLUS_DEBUG(mlog, "MFrontDetectionActor::render3DFrontSurfaceOIT");
    if (transferFunctionShading == nullptr
        || !frontSurfaceSelection
        || !render3DFront
        || num3DFrontVertices < 1
            )
    {
        return;
    }

    if (alphaShading3DFronts.size() < 1
            || alphaShading3DFronts.size() != num3DFrontVertices)
    {
        alphaShading3DFronts = QVector<QVector2D>(num3DFrontVertices, QVector2D(1.0, 1.0));
        recomputeAlpha3DFronts = true;
        recomputeShading3DFronts = true;
    }

    if (recomputeAlpha3DFronts)
    {
        for (int i = 0; i < num3DFrontVertices; i++)
        {
            alphaShading3DFronts[i][0] = 1;
        }
        foreach (auto filter, normalCurvesFilterSetList)
        {
            if (!filter->enabled || filter->transferFunction == nullptr) continue;
            MStructuredGrid* f = variables.at(filter->variableIndex)->grid;
            switch(filter->type)
            {
                case ABSOLUTE_CHANGE:
                {
#pragma omp parallel for
                    for (int p = 0; p < num3DFrontVertices; p++)
                    {
                        QVector3D startPos = frontMesh3D->getVertex(p).position;
                        QVector3D endPos = frontMesh3D->getVertex(p).nCEnd;
                        float vStart = f->interpolateValue(startPos);
                        float vEnd = f->interpolateValue(endPos);
                        float filterValue = vStart - vEnd;
                        alphaShading3DFronts[p][0] *=
                                filter->transferFunction->getColorValue(
                                    filterValue).alphaF();
                    }
                    break;
                }
                case PER_100KM:
                {
#pragma omp parallel for
                    for (int p = 0; p < num3DFrontVertices; p++)
                    {
                        QVector3D startPos = frontMesh3D->getVertex(p).position;
                        QVector3D endPos = frontMesh3D->getVertex(p).nCEnd;
                        float vStart = f->interpolateValue(startPos);
                        float vEnd = f->interpolateValue(endPos);
                        float filterValue = (vStart - vEnd)/frontMesh3D->getVertex(p).breadth * 100;
                        alphaShading3DFronts[p][0] *=
                                filter->transferFunction->getColorValue(
                                    filterValue).alphaF();
                    }
                    break;
                }
                case VALUE_AT_FRONT:
                {
#pragma omp parallel for
                    for (int p = 0; p < num3DFrontVertices; p++)
                    {
                        QVector3D startPos = frontMesh3D->getVertex(p).position;
                        float vStart = f->interpolateValue(startPos);
                        alphaShading3DFronts[p][0] *=
                                filter->transferFunction->getColorValue(
                                    vStart).alphaF();
                    }
                    break;
                }
            }
        }
        recomputeAlpha3DFronts = false;
    }

    if (recomputeShading3DFronts)
    {
        if (shadingMode == 0)
        {
            if (shadingVariableIndex < 0) return;
            MStructuredGrid* f = variables.at(shadingVariableIndex)->grid;
            switch (shadingVarMode)
            {
                case ABSOLUTE_CHANGE:
                {
#pragma omp parallel for
                    for (int p = 0; p < num3DFrontVertices; p++)
                    {
                        QVector3D startPos = frontMesh3D->getVertex(p).position;
                        QVector3D endPos = frontMesh3D->getVertex(p).nCEnd;
                        float vStart = f->interpolateValue(startPos);
                        float vEnd = f->interpolateValue(endPos);
                        alphaShading3DFronts[p][1] = vStart - vEnd;
                    }
                    break;
                }
                case PER_100KM:
                {
#pragma omp parallel for
                    for (int p = 0; p < num3DFrontVertices; p++)
                    {
                        QVector3D startPos = frontMesh3D->getVertex(p).position;
                        QVector3D endPos = frontMesh3D->getVertex(p).nCEnd;
                        float vStart = f->interpolateValue(startPos);
                        float vEnd = f->interpolateValue(endPos);
                        alphaShading3DFronts[p][1] = (vStart - vEnd)
                                /frontMesh3D->getVertex(p).breadth * 100;
                    }
                    break;
                }
                case VALUE_AT_FRONT:
                {
#pragma omp parallel for
                    for (int p = 0; p < num3DFrontVertices; p++)
                    {
                        QVector3D startPos = frontMesh3D->getVertex(p).position;
                        alphaShading3DFronts[p][1] = f->interpolateValue(startPos);
                    }
                    break;
                }
            }
        }
        else if (shadingMode == 1)
        {
#pragma omp parallel for
            for (int p = 0; p < num3DFrontVertices; p++)
            {
                alphaShading3DFronts[p][1] = frontMesh3D->getVertex(p).position.z();
            }
        }
        else if (shadingMode == 2)
        {
#pragma omp parallel for
            for (int p = 0; p < num3DFrontVertices; p++)
            {
                alphaShading3DFronts[p][1] = frontMesh3D->getVertex(p).tfp;
            }
        }
        else if (shadingMode == 3)
        {
#pragma omp parallel for
            for (int p = 0; p < num3DFrontVertices; p++)
            {
                alphaShading3DFronts[p][1] = frontMesh3D->getVertex(p).abz;
            }
        }
        else if (shadingMode == 4)
        {
#pragma omp parallel for
            for (int p = 0; p < num3DFrontVertices; p++)
            {
                alphaShading3DFronts[p][1] = frontMesh3D->getVertex(p).slope;
            }
        }
        else if (shadingMode == 5)
        {
#pragma omp parallel for
            for (int p = 0; p < num3DFrontVertices; p++)
            {
                alphaShading3DFronts[p][1] = frontMesh3D->getVertex(p).breadth;
            }
        }
        recomputeShading3DFronts = false;
    }


    float dataMinZ;
    float dataMaxZ;
    QVector3D minBBox;
    QVector3D maxBBox;

    if (bBoxConnection->getBoundingBox() == nullptr)
    {
        dataMinZ = static_cast<float>(sceneView->worldZfromPressure(
                                          detectionVariable->grid->getBottomDataVolumePressure_hPa()));
        dataMaxZ = static_cast<float>(sceneView->worldZfromPressure(
                                          detectionVariable->grid->getTopDataVolumePressure_hPa()));
        minBBox = QVector3D(detectionVariable->grid->getWestDataVolumeCorner_lon(),
                            detectionVariable->grid->getSouthDataVolumeCorner_lat(),
                            dataMinZ);
        maxBBox = QVector3D(detectionVariable->grid->getEastDataVolumeCorner_lon(),
                            detectionVariable->grid->getNorthDataVolumeCorner_lat(),
                            dataMaxZ);
    }
    else
    {
        dataMinZ = static_cast<float>(
            sceneView->worldZfromPressure(bBoxConnection->bottomPressure_hPa()));
        dataMaxZ = static_cast<float>(
            sceneView->worldZfromPressure(bBoxConnection->topPressure_hPa()));
        minBBox = QVector3D(bBoxConnection->westLon(), bBoxConnection->southLat(), dataMinZ);
        maxBBox = QVector3D(bBoxConnection->eastLon(), bBoxConnection->northLat(), dataMaxZ);
    }
    // 1) Render shadows
    // 2) Render transparent triangle geometry

    for (auto i = 0; i < 2; ++i)
    {
        if (sceneView->shadowMappingInProgress())
        {
            fronts3DSurfaceShader->bindProgram("TriangleFilteringShadowMap");
        } else
        {
            fronts3DSurfaceShader->bindProgram((i == 0) ? "TriangleFilteringShadow"
                                               : "TriangleFilteringOIT");
        }

        fronts3DSurfaceShader->setUniformValue("mvpMatrix",
                                      *(sceneView->getModelViewProjectionMatrix()));
        fronts3DSurfaceShader->setUniformValue("pToWorldZParams",
                                      sceneView->pressureToWorldZParameters());
        fronts3DSurfaceShader->setUniformValue("lightDirection",
                                      sceneView->getLightDirection());
        fronts3DSurfaceShader->setUniformValue("cameraPosition",
                                      sceneView->getCamera()->getOrigin());

        fronts3DSurfaceShader->setUniformValue("bboxMax", maxBBox);
        fronts3DSurfaceShader->setUniformValue("bboxMin", minBBox);

        fronts3DSurfaceShader->setUniformValue("displayFrontTypes", showFrontTypes);
        fronts3DSurfaceShader->setUniformValue("displayColdSideFront", showColdSideFront);

        fronts3DSurfaceShader->setUniformValue("useTFPFilter", useTFPTransferFunction);
        fronts3DSurfaceShader->setUniformValue("useFSFilter", useFSTransferFunction);
        fronts3DSurfaceShader->setUniformValue("useABZFilter", useABZTransferFunction);
        fronts3DSurfaceShader->setUniformValue("useSlopeFilter", useSlopeTransferFunction);
        fronts3DSurfaceShader->setUniformValue("useBreadthFilter",
                                               useBreadthTransferFunction);

        float shadowHeightWorldZ = sceneView->worldZfromPressure(shadowHeight);
        fronts3DSurfaceShader->setUniformValue("shadowColor", shadowColor);
        fronts3DSurfaceShader->setUniformValue("shadowHeight",
                                               shadowHeightWorldZ);

        fronts3DSurfaceShader->setUniformValue("shadingMode", shadingMode);

        fronts3DSurfaceShader->setUniformValue("lightingMode", lightingMode);

        if (transferFunctionShading)
        {
            transferFunctionShading->getTexture()->bindToTextureUnit(transferFunctionShadingTexUnit);CHECK_GL_ERROR;
            fronts3DSurfaceShader->setUniformValue("tfShading", transferFunctionShadingTexUnit); CHECK_GL_ERROR;
            fronts3DSurfaceShader->setUniformValue(
                        "tfShadingMinMax", QVector2D(
                            transferFunctionShading->getMinimumValue(),
                            transferFunctionShading->getMaximumValue()));
        }
        if (transferFunctionTFP != nullptr)
        {
            transferFunctionTFP->getTexture()->bindToTextureUnit(transferFunctionTFPTexUnit);
            fronts3DSurfaceShader->setUniformValue("tfTFP", transferFunctionTFPTexUnit); CHECK_GL_ERROR;

            fronts3DSurfaceShader->setUniformValue("tfTFPMinMax",
                                                   QVector2D(transferFunctionTFP->getMinimumValue(),
                                                             transferFunctionTFP->getMaximumValue()));
        }

        if (transferFunctionFS)
        {
            transferFunctionFS->getTexture()->bindToTextureUnit(transferFunctionFSTexUnit);
            fronts3DSurfaceShader->setUniformValue("tfFS", transferFunctionFSTexUnit); CHECK_GL_ERROR;

            fronts3DSurfaceShader->setUniformValue("tfFSMinMax",
                                                   QVector2D(transferFunctionFS->getMinimumValue(),
                                                             transferFunctionFS->getMaximumValue()));
        }

        if (transferFunctionABZ != nullptr)
        {
            transferFunctionABZ->getTexture()->bindToTextureUnit(transferFunctionABZTexUnit);
            fronts3DSurfaceShader->setUniformValue("tfABZ", transferFunctionABZTexUnit); CHECK_GL_ERROR;

            fronts3DSurfaceShader->setUniformValue("tfABZMinMax",
                    QVector2D(transferFunctionABZ->getMinimumValue(),
                              transferFunctionABZ->getMaximumValue()));
        }

        if (transferFunctionSlope != nullptr)
        {
            transferFunctionSlope->getTexture()->bindToTextureUnit(transferFunctionSlopeTexUnit);
            fronts3DSurfaceShader->setUniformValue("tfSlope", transferFunctionSlopeTexUnit); CHECK_GL_ERROR;

            fronts3DSurfaceShader->setUniformValue("tfSlopeMinMax",
                    QVector2D(transferFunctionSlope->getMinimumValue(),
                              transferFunctionSlope->getMaximumValue()));
        }

        if (transferFunctionBreadth != nullptr)
        {
            transferFunctionBreadth->getTexture()->bindToTextureUnit(transferFunctionBreadthTexUnit);
            fronts3DSurfaceShader->setUniformValue("tfBreadth", transferFunctionBreadthTexUnit); CHECK_GL_ERROR;

            fronts3DSurfaceShader->setUniformValue("tfBreadthMinMax",
                                                   QVector2D(transferFunctionBreadth->getMinimumValue(),
                                                             transferFunctionBreadth->getMaximumValue()));
        }

        fronts3DSurfaceShader->setUniformValue(
                    "inShadowMappingMode", sceneView->shadowMappingInProgress());
        fronts3DSurfaceShader->setUniformValue(
                    "lightMatrix", *sceneView->getLightOrthoProjectionMatrix());

        sceneView->getShadowMap()->bindToTextureUnit(
                    static_cast<GLuint>(sceneView->getShadowMapTexUnit()));
        fronts3DSurfaceShader->setUniformValue(
                    "shadowMap", sceneView->getShadowMapTexUnit());
        CHECK_GL_ERROR;

        sceneView->setOITUniforms(fronts3DSurfaceShader);

        const QString vboIDa = QString("frontsurface_vbo_shading#%1").arg(myID);
        uploadVec2ToVertexBuffer(alphaShading3DFronts, vboIDa, &vbo3DFrontsShading,
                                 sceneView);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,
                     ibo3DFronts->getIndexBufferObject()); CHECK_GL_ERROR;
        glBindBuffer(GL_ARRAY_BUFFER,
                     vbo3DFronts->getVertexBufferObject()); CHECK_GL_ERROR;

        int fontMeshVertexSize = sizeof(Geometry::FrontMeshVertex);
        long bytePosition = 0;

        // vertex
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, fontMeshVertexSize,
                              nullptr); CHECK_GL_ERROR;
        // normals
        bytePosition += 3 * sizeof(float);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, fontMeshVertexSize,
                              (const GLvoid*) bytePosition); CHECK_GL_ERROR;
        // tfp
        bytePosition += 6 * sizeof(float); // skip normal curves end vertex
        glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, fontMeshVertexSize,
                              (const GLvoid*) bytePosition); CHECK_GL_ERROR;
        //abz
        bytePosition += sizeof(float);
        glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, fontMeshVertexSize,
                              (const GLvoid*) bytePosition); CHECK_GL_ERROR;
        //strength
        bytePosition += sizeof(float);
        glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, fontMeshVertexSize,
                              (const GLvoid*) bytePosition); CHECK_GL_ERROR;
        // type (warm o. cold)
        bytePosition += sizeof(float);
        glVertexAttribPointer(5, 1, GL_FLOAT, GL_FALSE, fontMeshVertexSize,
                              (const GLvoid*) bytePosition); CHECK_GL_ERROR;
        // breadth
        bytePosition += sizeof(float);
        glVertexAttribPointer(6, 1, GL_FLOAT, GL_FALSE, fontMeshVertexSize,
                              (const GLvoid*) bytePosition); CHECK_GL_ERROR;
        // slope
        bytePosition += sizeof(float);
        glVertexAttribPointer(7, 1, GL_FLOAT, GL_FALSE, fontMeshVertexSize,
                              (const GLvoid*) bytePosition); CHECK_GL_ERROR;

        //alpha
        glBindBuffer(GL_ARRAY_BUFFER, vbo3DFrontsShading->getVertexBufferObject()); CHECK_GL_ERROR;
        glVertexAttribPointer(8, 1, GL_FLOAT, GL_FALSE, sizeof(QVector2D),
                              nullptr); CHECK_GL_ERROR;
        //shading
        glVertexAttribPointer(9, 1, GL_FLOAT, GL_FALSE, sizeof(QVector2D),
                              (const GLvoid*) sizeof(float)); CHECK_GL_ERROR;

        glEnableVertexAttribArray(0); CHECK_GL_ERROR;
        glEnableVertexAttribArray(1); CHECK_GL_ERROR;
        glEnableVertexAttribArray(2); CHECK_GL_ERROR;
        glEnableVertexAttribArray(3); CHECK_GL_ERROR;
        glEnableVertexAttribArray(4); CHECK_GL_ERROR;
        glEnableVertexAttribArray(5); CHECK_GL_ERROR;
        glEnableVertexAttribArray(6); CHECK_GL_ERROR;
        glEnableVertexAttribArray(7); CHECK_GL_ERROR;
        glEnableVertexAttribArray(8); CHECK_GL_ERROR;
        glEnableVertexAttribArray(9); CHECK_GL_ERROR;

        // Create linked list for fragments
        glPolygonMode(GL_FRONT_AND_BACK, (renderAsWireFrame) ? GL_LINE : GL_FILL); CHECK_GL_ERROR;
        glDrawElements(GL_TRIANGLES, frontMesh3D->getNumTriangleIndices(),
                       GL_UNSIGNED_INT, nullptr); CHECK_GL_ERROR;
    }

    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

}


void MFrontDetectionActor::render3DNormalCurvesOIT(
        MSceneViewGLWidget *sceneView)
{
    // ########################################################################
    // ### CHECK IF ALL VALUES ARE VALID
    // ########################################################################

    if (transferFunctionNCShading == nullptr
        || !frontSurfaceSelection
        || !render3DNC)
    {
        return;
    }

    if (trigger3DNormalCurveFilter)
    {
        filter3DNormalCurves();
        trigger3DNormalCurveFilter = false;
    }

    if (!normalCurves3D
        || normalCurves3D->getNumNormalCurveSegments() < 1)
    {
            return;
    }

    int numNormalCurveSegments = normalCurves3D->getNumNormalCurveSegments();

    if (alphaShading3DNormalCurves.size() < 1
            || alphaShading3DNormalCurves.size() != numNormalCurveSegments)
    {
        alphaShading3DNormalCurves =
                QVector<QVector2D>(numNormalCurveSegments, QVector2D(1.0, 1.0));
        recomputeAlpha3DNormalCurves = true;
        recomputeShading3DNormalCurves = true;
    }

    if (recomputeAlpha3DNormalCurves)
    {
        for (int i = 0; i < numNormalCurveSegments; i++)
        {
            alphaShading3DNormalCurves[i][0] = 1;
        }
        foreach (auto filter, normalCurvesFilterSetList)
        {
            if (!filter->enabled || filter->transferFunction == nullptr) continue;
            MStructuredGrid* aGrid = variables.at(filter->variableIndex)->grid;
            switch(filter->type)
            {
            case ABSOLUTE_CHANGE:
            {
                for (int i = 0; i < numNormalCurveSegments; i++)
                {
                    Geometry::NormalCurveVertex nc = normalCurves3D->getNormalCurveSegment(i);
                    float vStart = aGrid->interpolateValue(nc.start);
                    float vEnd = aGrid->interpolateValue(nc.end);
                    float filterValue = vStart - vEnd;
                    alphaShading3DNormalCurves[i][0] *=
                            filter->transferFunction->getColorValue(
                                filterValue).alphaF();
                }
                break;
            }
            case PER_100KM:
            {
                for (int i = 0; i < numNormalCurveSegments; i++)
                {
                    Geometry::NormalCurveVertex nc = normalCurves3D->getNormalCurveSegment(i);
                    float vStart = aGrid->interpolateValue(nc.start);
                    float vEnd = aGrid->interpolateValue(nc.end);
                    float filterValue = (vStart - vEnd)
                            / nc.breadth * 100.;
                    alphaShading3DNormalCurves[i][0] *=
                            filter->transferFunction->getColorValue(
                                filterValue).alphaF();
                }
                break;
            }
            case VALUE_AT_FRONT:
            {
                for (int i = 0; i < numNormalCurveSegments; i++)
                {
                    Geometry::NormalCurveVertex nc = normalCurves3D->getNormalCurveSegment(i);
                    float vStart = aGrid->interpolateValue(nc.start);
                    alphaShading3DNormalCurves[i][0] *=
                            filter->transferFunction->getColorValue(
                                vStart).alphaF();
                }
                break;
            }
            }
        }
        recomputeAlpha3DNormalCurves = false;
    }

    if (recomputeShading3DNormalCurves)
    {
        if (ncShadingMode == 1
            || ncShadingMode == 2
            || ncShadingMode == 3)
        {
            if (ncShadingVariableIndex < 0) return;
            MStructuredGrid* sGrid = variables.at(ncShadingVariableIndex)->grid;
            if (ncShadingMode == 1) // shading value
            {
                for (int i = 0; i < numNormalCurveSegments; i++)
                {
                    Geometry::NormalCurveVertex nc = normalCurves3D->getNormalCurveSegment(i);
                    float v = sGrid->interpolateValue(nc.start);
                    alphaShading3DNormalCurves[i][1] = v;
                }
            }
            else if (ncShadingMode == 2) // avg. difference
            {
                for (int i = 0; i < numNormalCurveSegments; i++)
                {
                    Geometry::NormalCurveVertex nc = normalCurves3D->getNormalCurveSegment(i);
                    float vStart = sGrid->interpolateValue(nc.start);
                    float vEnd = sGrid->interpolateValue(nc.end);
                    float shadingValue = (vStart - vEnd)
                            / nc.breadth * 100.;
                    alphaShading3DNormalCurves[i][1] = shadingValue;
                }
            }
            else if (ncShadingMode == 3) // total difference
            {
                for (int i = 0; i < numNormalCurveSegments; i++)
                {
                    Geometry::NormalCurveVertex nc = normalCurves3D->getNormalCurveSegment(i);
                    float vStart = sGrid->interpolateValue(
                                nc.start);
                    float vEnd = sGrid->interpolateValue(
                                nc.end);
                    alphaShading3DNormalCurves[i][1] = vStart - vEnd;
                }
            }
        }
        recomputeShading3DNormalCurves = false;
    }

    float dataMinZ;
    float dataMaxZ;
    QVector3D minBBox;
    QVector3D maxBBox;

    if (bBoxConnection->getBoundingBox() == nullptr)
    {
        dataMinZ = static_cast<float>(sceneView->worldZfromPressure(
                                          detectionVariable->grid->getBottomDataVolumePressure_hPa()));
        dataMaxZ = static_cast<float>(sceneView->worldZfromPressure(
                                          detectionVariable->grid->getTopDataVolumePressure_hPa()));
        minBBox = QVector3D(detectionVariable->grid->getWestDataVolumeCorner_lon(),
                            detectionVariable->grid->getSouthDataVolumeCorner_lat(),
                            dataMinZ);
        maxBBox = QVector3D(detectionVariable->grid->getEastDataVolumeCorner_lon(),
                            detectionVariable->grid->getNorthDataVolumeCorner_lat(),
                            dataMaxZ);
    }
    else
    {
        dataMinZ = static_cast<float>(
            sceneView->worldZfromPressure(bBoxConnection->bottomPressure_hPa()));
        dataMaxZ = static_cast<float>(
            sceneView->worldZfromPressure(bBoxConnection->topPressure_hPa()));
        minBBox = QVector3D(bBoxConnection->westLon(), bBoxConnection->southLat(), dataMinZ);
        maxBBox = QVector3D(bBoxConnection->eastLon(), bBoxConnection->northLat(), dataMaxZ);
    }

    // ########################################################################
    // ### BIND BUFFER AND SET UNIFORM VALUES IN GLSL
    // ########################################################################

    for (auto i = 0; i < 2; ++i)
    {
        if (sceneView->shadowMappingInProgress())
        {
            normalCurvesTubeShaderOIT->bindProgram("TubeFilteringShadowMap3D");
        } else
        {
            normalCurvesTubeShaderOIT->bindProgram((i == 0) ? "TubeFilteringShadow3D"
                                               : "TubeFilteringOIT3D");
        }

        normalCurvesTubeShaderOIT->setUniformValue(
                "pToWorldZParams", sceneView->pressureToWorldZParameters());

        normalCurvesTubeShaderOIT->setUniformValue(
                "mvpMatrix", *(sceneView->getModelViewProjectionMatrix())); CHECK_GL_ERROR;
        normalCurvesTubeShaderOIT->setUniformValue(
                "cameraPosition", sceneView->getCamera()->getOrigin()); CHECK_GL_ERROR;

        // Lighting direction from scene view.
        normalCurvesTubeShaderOIT->setUniformValue(
                "lightDirection", sceneView->getLightDirection()); CHECK_GL_ERROR;

        normalCurvesTubeShaderOIT->setUniformValue("bboxMax", maxBBox);
        normalCurvesTubeShaderOIT->setUniformValue("bboxMin", minBBox);

        transferFunctionNCShading->getTexture()->bindToTextureUnit(
                transferFunctionNCShadingTexUnit);
        normalCurvesTubeShaderOIT->setUniformValue(
                "transferFunction", transferFunctionNCShadingTexUnit);

        normalCurvesTubeShaderOIT->setUniformValue(
                "tfMinimum", transferFunctionNCShading->getMinimumValue());
        normalCurvesTubeShaderOIT->setUniformValue(
                "tfMaximum", transferFunctionNCShading->getMaximumValue());

        normalCurvesTubeShaderOIT->setUniformValue(
                "normalized", false);

        normalCurvesTubeShaderOIT->setUniformValue("tubeRadius", normalCurvesTubeRadius); CHECK_GL_ERROR;

        float shadowHeightWorldZ = sceneView->worldZfromPressure(shadowHeight);
        normalCurvesTubeShaderOIT->setUniformValue("shadowColor", shadowColor);
        normalCurvesTubeShaderOIT->setUniformValue("shadowHeight",
                                               shadowHeightWorldZ);

        normalCurvesTubeShaderOIT->setUniformValue("useTFPFilter", useTFPTransferFunction);
        normalCurvesTubeShaderOIT->setUniformValue("useABZFilter", useABZTransferFunction);
        normalCurvesTubeShaderOIT->setUniformValue("useFSFilter", useFSTransferFunction);
        normalCurvesTubeShaderOIT->setUniformValue("useBreadthFilter", useBreadthTransferFunction);

        if (transferFunctionTFP)
        {
            transferFunctionTFP->getTexture()->bindToTextureUnit(transferFunctionTFPTexUnit);
            normalCurvesTubeShaderOIT->setUniformValue("tfTFP", transferFunctionTFPTexUnit); CHECK_GL_ERROR;

            normalCurvesTubeShaderOIT->setUniformValue("tfTFPMinMax",
                                                QVector2D(transferFunctionTFP->getMinimumValue(),
                                                transferFunctionTFP->getMaximumValue()));
        }

        if (transferFunctionABZ)
        {
            transferFunctionABZ->getTexture()->bindToTextureUnit(transferFunctionABZTexUnit);
            normalCurvesTubeShaderOIT->setUniformValue("tfABZ", transferFunctionABZTexUnit); CHECK_GL_ERROR;

            normalCurvesTubeShaderOIT->setUniformValue("tfABZMinMax",
                    QVector2D(transferFunctionABZ->getMinimumValue(),
                              transferFunctionABZ->getMaximumValue()));
        }

        if (transferFunctionFS)
        {
            transferFunctionFS->getTexture()->bindToTextureUnit(transferFunctionFSTexUnit);
            normalCurvesTubeShaderOIT->setUniformValue("tfFS", transferFunctionFSTexUnit); CHECK_GL_ERROR;

            normalCurvesTubeShaderOIT->setUniformValue("tfFSMinMax",
                                                   QVector2D(transferFunctionFS->getMinimumValue(),
                                                             transferFunctionFS->getMaximumValue()));
        }

        if (transferFunctionBreadth)
        {
            transferFunctionBreadth->getTexture()->bindToTextureUnit(transferFunctionBreadthTexUnit);
            normalCurvesTubeShaderOIT->setUniformValue("tfBreadth",
                                             transferFunctionBreadthTexUnit);
            CHECK_GL_ERROR;
            normalCurvesTubeShaderOIT->setUniformValue(
                        "tfBreadthMinMax",
                        QVector2D(transferFunctionBreadth->getMinimumValue(),
                                  transferFunctionBreadth->getMaximumValue()));
        }

        normalCurvesTubeShaderOIT->setUniformValue(
                    "inShadowMappingMode", sceneView->shadowMappingInProgress());

        normalCurvesTubeShaderOIT->setUniformValue(
                    "shadingMode", ncShadingMode);

        sceneView->getShadowMap()->bindToTextureUnit(
                    static_cast<GLuint>(sceneView->getShadowMapTexUnit()));
        normalCurvesTubeShaderOIT->setUniformValue(
                    "shadowMap", sceneView->getShadowMapTexUnit());
        CHECK_GL_ERROR;

        sceneView->setOITUniforms(normalCurvesTubeShaderOIT);
        // ########################################################################
        // ### GET RESOURCE MANAGER AND UPLOAD DATA TO VERTEX BUFFER
        // ########################################################################

        const QString vboIDs = QString("normal_curves_vbo_shading#%1").arg(myID);
        uploadVec2ToVertexBuffer(alphaShading3DNormalCurves, vboIDs,
                                 &vbo3DNormalCurvesShading, sceneView);

        //front vertices
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,
                     ibo3DNormalCurves->getIndexBufferObject()); CHECK_GL_ERROR;
        glBindBuffer(GL_ARRAY_BUFFER,
                     vbo3DNormalCurves->getVertexBufferObject()); CHECK_GL_ERROR;

        int normalCurveVertexSize = sizeof(Geometry::NormalCurveVertex);
        long bytePosition = 0;

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, normalCurveVertexSize,
                              nullptr); CHECK_GL_ERROR;

        bytePosition += 9 * sizeof(float); // tfp, 9 because of start and end points  (2 * QVector3D)
        glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, normalCurveVertexSize,
                              (const GLvoid*) bytePosition); CHECK_GL_ERROR;

        bytePosition += sizeof(float); // abz
        glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, normalCurveVertexSize,
                              (const GLvoid*) bytePosition); CHECK_GL_ERROR;

        bytePosition += sizeof(float); // strength
        glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, normalCurveVertexSize,
                              (const GLvoid*) bytePosition); CHECK_GL_ERROR;

        bytePosition += sizeof(float); // type
        glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, normalCurveVertexSize,
                              (const GLvoid*) bytePosition); CHECK_GL_ERROR;

        bytePosition += sizeof(float); // breadth
        glVertexAttribPointer(5, 1, GL_FLOAT, GL_FALSE, normalCurveVertexSize,
                             (const GLvoid*) bytePosition); CHECK_GL_ERROR;

        glBindBuffer(GL_ARRAY_BUFFER, vbo3DNormalCurvesShading->getVertexBufferObject()); CHECK_GL_ERROR;
        glVertexAttribPointer(6, 1, GL_FLOAT, GL_FALSE, sizeof(QVector2D),
                              nullptr); CHECK_GL_ERROR;

        glVertexAttribPointer(7, 1, GL_FLOAT, GL_FALSE, sizeof(QVector2D),
                              (const GLvoid*) sizeof(float)); CHECK_GL_ERROR;

        glEnableVertexAttribArray(0); CHECK_GL_ERROR;
        glEnableVertexAttribArray(1); CHECK_GL_ERROR;
        glEnableVertexAttribArray(2); CHECK_GL_ERROR;
        glEnableVertexAttribArray(3); CHECK_GL_ERROR;
        glEnableVertexAttribArray(4); CHECK_GL_ERROR;
        glEnableVertexAttribArray(5); CHECK_GL_ERROR;
        glEnableVertexAttribArray(6); CHECK_GL_ERROR;
        glEnableVertexAttribArray(7); CHECK_GL_ERROR;

        glPrimitiveRestartIndex(restartIndex);
        glEnable(GL_PRIMITIVE_RESTART);

        glEnable(GL_CULL_FACE); CHECK_GL_ERROR;
        glCullFace(GL_FRONT);

        glPolygonMode(GL_BACK, GL_FILL); CHECK_GL_ERROR;
        glDrawElements(GL_LINE_STRIP_ADJACENCY, numNormalCurveSegments, GL_UNSIGNED_INT, nullptr);
        CHECK_GL_ERROR;

    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);  CHECK_GL_ERROR;
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0); CHECK_GL_ERROR;
    glDisable(GL_CULL_FACE); CHECK_GL_ERROR;
    glDisable(GL_PRIMITIVE_RESTART);
}


void MFrontDetectionActor::render2DFrontTubesOIT(
        MSceneViewGLWidget *sceneView)
{
    // ########################################################################
    // ### CHECK IF ALL VALUES ARE VALID
    // ########################################################################
    if (transferFunctionShading == nullptr
        || !frontLineSelection
        || !render2DFront
        || num2DFrontVertices < 1)
    {
        return;
    }


    MNWP2DHorizontalActorVariable *varThermal =
            static_cast<MNWP2DHorizontalActorVariable*>(
                variables.at(detectionVariableIndex)
                );

    // Don't render 2d front lines if horizontal slice position is outside the
    // data domain.

    if (varThermal->grid->getLevelType() != SURFACE_2D &&
            (varThermal->grid->getBottomDataVolumePressure_hPa() < frontElevation2d_hPa
             || varThermal->grid->getTopDataVolumePressure_hPa() > frontElevation2d_hPa))
    {
        return;
    }


    if (alphaShading2DFronts.size() < 1
            || alphaShading2DFronts.size() != num2DFrontVertices)
    {
        alphaShading2DFronts = QVector<QVector2D>(num2DFrontVertices, QVector2D(1.0, 1.0));
        recomputeAlpha2DFronts = true;
        recomputeShading2DFronts = true;
    }


    if (recomputeAlpha2DFronts)
    {
        for (int i = 0; i < num2DFrontVertices; i++)
        {
            alphaShading2DFronts[i][0] = 1;
        }
        foreach (auto filter, normalCurvesFilterSetList)
        {
            if (!filter->enabled || filter->transferFunction == nullptr) continue;
            MStructuredGrid* aGrid = variables.at(filter->variableIndex)->grid;
            switch(filter->type)
            {
            case ABSOLUTE_CHANGE:
            {
                for (int i = 0; i < num2DFrontVertices; i++)
                {
                    float vStart = aGrid->interpolateValue(
                                frontLines2D->getVertex(i).position);
                    float vEnd = aGrid->interpolateValue(
                                frontLines2D->getVertex(i).nCEnd);
                    float filterValue = vStart - vEnd;
                    alphaShading2DFronts[i][0] *=
                            filter->transferFunction->getColorValue(
                                filterValue).alphaF();
                }
                break;
            }
            case PER_100KM:
            {
                for (int i = 0; i < num2DFrontVertices; i++)
                {
                    float vStart = aGrid->interpolateValue(
                                    frontLines2D->getVertex(i).position);
                    float vEnd = aGrid->interpolateValue(
                                    frontLines2D->getVertex(i).nCEnd);
                    float filterValue = (vStart - vEnd)
                            / frontLines2D->getVertex(i).breadth * 100.;
                    alphaShading2DFronts[i][0] *=
                            filter->transferFunction->getColorValue(
                                filterValue).alphaF();
                }
                break;
            }
            case VALUE_AT_FRONT:
            {
                for (int i = 0; i < num2DFrontVertices; i++)
                {
                    float vStart = aGrid->interpolateValue(
                                    frontLines2D->getVertex(i).position);
                    alphaShading2DFronts[i][0] *=
                            filter->transferFunction->getColorValue(
                                vStart).alphaF();
                }
                break;
            }
            }
        }
        recomputeAlpha2DFronts = false;
    }

    if (recomputeShading2DFronts)
    {
        if (shadingMode == 0
            || shadingMode == 4)
        {
            if (shadingVariableIndex < 0) return;
            MStructuredGrid* sGrid = variables.at(shadingVariableIndex)->grid;
            switch (shadingVarMode)
            {
            case ABSOLUTE_CHANGE:
            {
                for (int i = 0; i < num2DFrontVertices; i++)
                {
                    float vStart = sGrid->interpolateValue(
                                    frontLines2D->getVertex(i).position);
                    float vEnd = sGrid->interpolateValue(
                                    frontLines2D->getVertex(i).nCEnd);
                    alphaShading2DFronts[i][1] = vStart - vEnd;
                }
                break;
            }
            case PER_100KM:
            {
                for (int i = 0; i < num2DFrontVertices; i++)
                {
                    float vStart = sGrid->interpolateValue(
                                    frontLines2D->getVertex(i).position);
                    float vEnd = sGrid->interpolateValue(
                                    frontLines2D->getVertex(i).nCEnd);
                    alphaShading2DFronts[i][1] = (vStart - vEnd)
                            / frontLines2D->getVertex(i).breadth * 100.;
                }
                break;
            }
            case VALUE_AT_FRONT:
            {
                for (int i = 0; i < num2DFrontVertices; i++)
                {
                    alphaShading2DFronts[i][1] = sGrid->interpolateValue(
                                frontLines2D->getVertex(i).position);
                }
                break;
            }
            }
        }
        else if (shadingMode == 1)
        {
            for (int i = 0; i < num2DFrontVertices; i++)
            {
                alphaShading2DFronts[i][1] = frontElevation2d_hPa;
            }
        }
        else if (shadingMode == 2)
        {
            for (int i = 0; i < num2DFrontVertices; i++)
            {
                alphaShading2DFronts[i][1] = frontLines2D->getVertex(i).tfp;
            }
        }
        else if (shadingMode == 3)
        {
            for (int i = 0; i < num2DFrontVertices; i++)
            {
                alphaShading2DFronts[i][1] = frontLines2D->getVertex(i).abz;
            }
        }
        else if (shadingMode == 5)
        {
            for (int i = 0; i < num2DFrontVertices; i++)
            {
                alphaShading2DFronts[i][1] = frontLines2D->getVertex(i).breadth;
            }
        }
        recomputeShading2DFronts = false;
    }

    float dataMinZ;
    float dataMaxZ;
    QVector3D minBBox;
    QVector3D maxBBox;

    if (bBoxConnection->getBoundingBox() == nullptr)
    {
        dataMinZ = static_cast<float>(sceneView->worldZfromPressure(
                                          detectionVariable->grid->getBottomDataVolumePressure_hPa()));
        dataMaxZ = static_cast<float>(sceneView->worldZfromPressure(
                                          detectionVariable->grid->getTopDataVolumePressure_hPa()));
        minBBox = QVector3D(detectionVariable->grid->getWestDataVolumeCorner_lon(),
                            detectionVariable->grid->getSouthDataVolumeCorner_lat(),
                            dataMinZ);
        maxBBox = QVector3D(detectionVariable->grid->getEastDataVolumeCorner_lon(),
                            detectionVariable->grid->getNorthDataVolumeCorner_lat(),
                            dataMaxZ);
    }
    else
    {
        dataMinZ = static_cast<float>(
            sceneView->worldZfromPressure(bBoxConnection->bottomPressure_hPa()));
        dataMaxZ = static_cast<float>(
            sceneView->worldZfromPressure(bBoxConnection->topPressure_hPa()));
        minBBox = QVector3D(bBoxConnection->westLon(), bBoxConnection->southLat(), dataMinZ);
        maxBBox = QVector3D(bBoxConnection->eastLon(), bBoxConnection->northLat(), dataMaxZ);
    }

    float worldZPos = sceneView->worldZfromPressure(frontElevation2d_hPa);
    // ########################################################################
    // ### BIND BUFFER AND SET UNIFORM VALUES IN GLSL
    // ########################################################################

    for (auto i = 0; i < 2; ++i)
    {
        if (sceneView->shadowMappingInProgress())
        {
            frontTubeShaderOIT->bindProgram("TubeFilteringShadowMap");
        } else
        {
            frontTubeShaderOIT->bindProgram((i == 0) ? "TubeFilteringShadow"
                                               : "TubeFilteringOIT");
        }

        frontTubeShaderOIT->setUniformValue(
                "worldZPos",
                worldZPos);
        frontTubeShaderOIT->setUniformValue(
                "mvpMatrix", *(sceneView->getModelViewProjectionMatrix())); CHECK_GL_ERROR;
        frontTubeShaderOIT->setUniformValue(
                "cameraPosition", sceneView->getCamera()->getOrigin()); CHECK_GL_ERROR;

        // Lighting direction from scene view.
        frontTubeShaderOIT->setUniformValue(
                "lightDirection", sceneView->getLightDirection()); CHECK_GL_ERROR;

        frontTubeShaderOIT->setUniformValue("bboxMax", maxBBox);
        frontTubeShaderOIT->setUniformValue("bboxMin", minBBox);

        transferFunctionShading->getTexture()->bindToTextureUnit(
                transferFunctionShadingTexUnit);
        frontTubeShaderOIT->setUniformValue(
                "transferFunction", transferFunctionShadingTexUnit);

        frontTubeShaderOIT->setUniformValue(
                "tfMinimum", transferFunctionShading->getMinimumValue());
        frontTubeShaderOIT->setUniformValue(
                "tfMaximum", transferFunctionShading->getMaximumValue());

        frontTubeShaderOIT->setUniformValue(
                "normalized", false);

        frontTubeShaderOIT->setUniformValue("tubeRadius", frontsTubeRadius); CHECK_GL_ERROR;

        float shadowHeightWorldZ = sceneView->worldZfromPressure(shadowHeight);
        frontTubeShaderOIT->setUniformValue("shadowColor", shadowColor);
        frontTubeShaderOIT->setUniformValue("shadowHeight",
                                               shadowHeightWorldZ);

        frontTubeShaderOIT->setUniformValue("useTFPFilter", useTFPTransferFunction);
        frontTubeShaderOIT->setUniformValue("useABZFilter", useABZTransferFunction);
        frontTubeShaderOIT->setUniformValue("useBreadthFilter", useBreadthTransferFunction);
        frontTubeShaderOIT->setUniformValue("useFSFilter", useFSTransferFunction);

        if (transferFunctionTFP)
        {
            transferFunctionTFP->getTexture()->bindToTextureUnit(transferFunctionTFPTexUnit);
            frontTubeShaderOIT->setUniformValue("tfTFP", transferFunctionTFPTexUnit); CHECK_GL_ERROR;

            frontTubeShaderOIT->setUniformValue("tfTFPMinMax",
                                                QVector2D(transferFunctionTFP->getMinimumValue(),
                                                transferFunctionTFP->getMaximumValue()));
        }

        if (transferFunctionABZ)
        {
            transferFunctionABZ->getTexture()->bindToTextureUnit(transferFunctionABZTexUnit);
            frontTubeShaderOIT->setUniformValue("tfABZ", transferFunctionABZTexUnit); CHECK_GL_ERROR;

            frontTubeShaderOIT->setUniformValue("tfABZMinMax",
                    QVector2D(transferFunctionABZ->getMinimumValue(),
                              transferFunctionABZ->getMaximumValue()));
        }


        if (transferFunctionFS)
        {
            transferFunctionFS->getTexture()->bindToTextureUnit(transferFunctionFSTexUnit);
            frontTubeShaderOIT->setUniformValue("tfFS", transferFunctionFSTexUnit); CHECK_GL_ERROR;

            frontTubeShaderOIT->setUniformValue("tfFSMinMax",
                                                   QVector2D(transferFunctionFS->getMinimumValue(),
                                                             transferFunctionFS->getMaximumValue()));
        }

        if (transferFunctionBreadth)
        {
            transferFunctionBreadth->getTexture()->bindToTextureUnit(transferFunctionBreadthTexUnit);
            frontTubeShaderOIT->setUniformValue("tfBreadth",
                                             transferFunctionBreadthTexUnit);
            CHECK_GL_ERROR;
            frontTubeShaderOIT->setUniformValue(
                        "tfBreadthMinMax",
                        QVector2D(transferFunctionBreadth->getMinimumValue(),
                                  transferFunctionBreadth->getMaximumValue()));
        }

        frontTubeShaderOIT->setUniformValue(
                    "inShadowMappingMode", sceneView->shadowMappingInProgress());

        sceneView->getShadowMap()->bindToTextureUnit(
                    static_cast<GLuint>(sceneView->getShadowMapTexUnit()));
        frontTubeShaderOIT->setUniformValue(
                    "shadowMap", sceneView->getShadowMapTexUnit());
        CHECK_GL_ERROR;

        sceneView->setOITUniforms(frontTubeShaderOIT);
        // ########################################################################
        // ### GET RESOURCE MANAGER AND UPLOAD DATA TO VERTEX BUFFER
        // ########################################################################

        const QString vboIDs = QString("frontline_vbo_shading#%1").arg(myID);
        uploadVec2ToVertexBuffer(alphaShading2DFronts, vboIDs, &vbo2DFrontsShading,
                                 sceneView);

        //front vertices
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,
                     ibo2DFronts->getIndexBufferObject()); CHECK_GL_ERROR;
        glBindBuffer(GL_ARRAY_BUFFER,
                     vbo2DFronts->getVertexBufferObject()); CHECK_GL_ERROR;

        int fontLineVertexSize = sizeof(Geometry::FrontLineVertex);
        long bytePosition = 0;

        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, fontLineVertexSize,
                              nullptr); CHECK_GL_ERROR;

        bytePosition += 6 * sizeof(float); // tfp 6 because there of ncEnd point (QVector3D)
        glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, fontLineVertexSize,
                              (const GLvoid*) bytePosition); CHECK_GL_ERROR;

        bytePosition += sizeof(float); // abz
        glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, fontLineVertexSize,
                              (const GLvoid*) bytePosition); CHECK_GL_ERROR;

        bytePosition += sizeof(float); // strength
        glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, fontLineVertexSize,
                              (const GLvoid*) bytePosition); CHECK_GL_ERROR;

        bytePosition += sizeof(float); // type
        glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, fontLineVertexSize,
                              (const GLvoid*) bytePosition); CHECK_GL_ERROR;

        bytePosition += sizeof(float); // breadth
        glVertexAttribPointer(5, 1, GL_FLOAT, GL_FALSE, fontLineVertexSize,
                             (const GLvoid*) bytePosition); CHECK_GL_ERROR;

        glBindBuffer(GL_ARRAY_BUFFER, vbo2DFrontsShading->getVertexBufferObject()); CHECK_GL_ERROR;
        glVertexAttribPointer(6, 1, GL_FLOAT, GL_FALSE, sizeof(QVector2D),
                              nullptr); CHECK_GL_ERROR;

        glVertexAttribPointer(7, 1, GL_FLOAT, GL_FALSE, sizeof(QVector2D),
                              (const GLvoid*) sizeof(float)); CHECK_GL_ERROR;

        glEnableVertexAttribArray(0); CHECK_GL_ERROR;
        glEnableVertexAttribArray(1); CHECK_GL_ERROR;
        glEnableVertexAttribArray(2); CHECK_GL_ERROR;
        glEnableVertexAttribArray(3); CHECK_GL_ERROR;
        glEnableVertexAttribArray(4); CHECK_GL_ERROR;
        glEnableVertexAttribArray(5); CHECK_GL_ERROR;
        glEnableVertexAttribArray(6); CHECK_GL_ERROR;
        glEnableVertexAttribArray(7); CHECK_GL_ERROR;

//        u_int32_t restartIndex = std::numeric_limits<u_int32_t>::max();
        glPrimitiveRestartIndex(restartIndex);
        glEnable(GL_PRIMITIVE_RESTART);

        glEnable(GL_CULL_FACE); CHECK_GL_ERROR;
        glCullFace(GL_FRONT);

        glPolygonMode(GL_BACK, GL_FILL); CHECK_GL_ERROR;
        glDrawElements(GL_LINE_STRIP_ADJACENCY, num2DFrontVertices,
                       GL_UNSIGNED_INT, nullptr); CHECK_GL_ERROR;

    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);  CHECK_GL_ERROR;
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0); CHECK_GL_ERROR;
    glDisable(GL_CULL_FACE); CHECK_GL_ERROR;
    glDisable(GL_PRIMITIVE_RESTART);
}


void MFrontDetectionActor::render2DNormalCurvesOIT(
        MSceneViewGLWidget *sceneView)
{
    // ########################################################################
    // ### CHECK IF ALL VALUES ARE VALID
    // ########################################################################
    if (transferFunctionNCShading == nullptr
        || !frontLineSelection
        || !render2DNC)
    {
        return;
    }

    if (trigger2DNormalCurveFilter)
    {
        filter2DNormalCurves();
        trigger2DNormalCurveFilter = false;
    }

    if (!normalCurves2D
        || normalCurves2D->getNumNormalCurveSegments() < 1)
    {
            return;
    }


    MNWP2DHorizontalActorVariable *varThermal =
            static_cast<MNWP2DHorizontalActorVariable*>(
                variables.at(detectionVariableIndex)
                );

    // Don't render if horizontal slice position is outside the
    // data domain.
    if (varThermal->grid->getLevelType() != SURFACE_2D &&
            (varThermal->grid->getBottomDataVolumePressure_hPa() < frontElevation2d_hPa
             || varThermal->grid->getTopDataVolumePressure_hPa() > frontElevation2d_hPa))
    {
        return;
    }

    int numNormalCurveSegments = normalCurves2D->getNumNormalCurveSegments();

    if (alphaShading2DNormalCurves.size() < 1
            || alphaShading2DNormalCurves.size() != numNormalCurveSegments)
    {
        alphaShading2DNormalCurves = QVector<QVector2D>(numNormalCurveSegments, QVector2D(1.0, 1.0));
        recomputeAlpha2DNormalCurves = true;
        recomputeShading2DNormalCurves = true;
    }

    if (recomputeAlpha2DNormalCurves)
    {
        for (int i = 0; i < numNormalCurveSegments; i++)
        {
            alphaShading2DNormalCurves[i][0] = 1;
        }
        foreach (auto filter, normalCurvesFilterSetList)
        {
            if (!filter->enabled || filter->transferFunction == nullptr) continue;
            MStructuredGrid* aGrid = variables.at(filter->variableIndex)->grid;
            switch(filter->type)
            {
            case ABSOLUTE_CHANGE:
            {
                for (int i = 0; i < numNormalCurveSegments; i++)
                {
                    Geometry::NormalCurveVertex nc = normalCurves2D->getNormalCurveSegment(i);
                    float vStart = aGrid->interpolateValue(nc.start);
                    float vEnd = aGrid->interpolateValue(nc.end);
                    float filterValue = vStart - vEnd;
                    alphaShading2DNormalCurves[i][0] *=
                            filter->transferFunction->getColorValue(
                                filterValue).alphaF();
                }
                break;
            }
            case PER_100KM:
            {
                for (int i = 0; i < numNormalCurveSegments; i++)
                {
                    Geometry::NormalCurveVertex nc = normalCurves2D->getNormalCurveSegment(i);
                    float vStart = aGrid->interpolateValue(nc.start);
                    float vEnd = aGrid->interpolateValue(nc.end);
                    float filterValue = (vStart - vEnd) / nc.breadth * 100.;
                    alphaShading2DNormalCurves[i][0] *=
                            filter->transferFunction->getColorValue(
                                filterValue).alphaF();
                }
                break;
            }
            case VALUE_AT_FRONT:
            {
                for (int i = 0; i < numNormalCurveSegments; i++)
                {
                    Geometry::NormalCurveVertex nc = normalCurves2D->getNormalCurveSegment(i);
                    float vStart = aGrid->interpolateValue(nc.start);
                    alphaShading2DNormalCurves[i][0] *=
                            filter->transferFunction->getColorValue(
                                vStart).alphaF();
                }
                break;
            }
            }
        }
        recomputeAlpha2DNormalCurves = false;
    }

    if (recomputeShading2DNormalCurves)
    {
        if (ncShadingMode == 1
            || ncShadingMode == 2
            || ncShadingMode == 3)
        {
            if (ncShadingVariableIndex < 0) return;
            MStructuredGrid* sGrid = variables.at(ncShadingVariableIndex)->grid;
            if (ncShadingMode == 1) // shading value
            {
                for (int i = 0; i < numNormalCurveSegments; i++)
                {
                    Geometry::NormalCurveVertex nc = normalCurves2D->getNormalCurveSegment(i);
                    float v = sGrid->interpolateValue(nc.start);
                    alphaShading2DNormalCurves[i][1] = v;
                }
            }
            else if (ncShadingMode == 2) // avg. difference
            {
                for (int i = 0; i < numNormalCurveSegments; i++)
                {
                    Geometry::NormalCurveVertex nc = normalCurves2D->getNormalCurveSegment(i);
                    float vStart = sGrid->interpolateValue(nc.start);
                    float vEnd = sGrid->interpolateValue(nc.end);
                    float shadingValue = (vStart - vEnd) / nc.breadth * 100.;
                    alphaShading2DNormalCurves[i][1] = shadingValue;
                }
            }
            else if (ncShadingMode == 3) // total difference
            {
                for (int i = 0; i < numNormalCurveSegments; i++)
                {
                    Geometry::NormalCurveVertex nc = normalCurves2D->getNormalCurveSegment(i);
                    float vStart = sGrid->interpolateValue(nc.start);
                    float vEnd = sGrid->interpolateValue(nc.end);
                    alphaShading2DNormalCurves[i][1] = vStart - vEnd;
                }
            }
        }
        recomputeShading2DNormalCurves = false;
    }


    float dataMinZ;
    float dataMaxZ;
    QVector3D minBBox;
    QVector3D maxBBox;

    if (bBoxConnection->getBoundingBox() == nullptr)
    {
        dataMinZ = static_cast<float>(sceneView->worldZfromPressure(
                                          detectionVariable->grid->getBottomDataVolumePressure_hPa()));
        dataMaxZ = static_cast<float>(sceneView->worldZfromPressure(
                                          detectionVariable->grid->getTopDataVolumePressure_hPa()));
        minBBox = QVector3D(detectionVariable->grid->getWestDataVolumeCorner_lon(),
                            detectionVariable->grid->getSouthDataVolumeCorner_lat(),
                            dataMinZ);
        maxBBox = QVector3D(detectionVariable->grid->getEastDataVolumeCorner_lon(),
                            detectionVariable->grid->getNorthDataVolumeCorner_lat(),
                            dataMaxZ);
    }
    else
    {
        dataMinZ = static_cast<float>(
            sceneView->worldZfromPressure(bBoxConnection->bottomPressure_hPa()));
        dataMaxZ = static_cast<float>(
            sceneView->worldZfromPressure(bBoxConnection->topPressure_hPa()));
        minBBox = QVector3D(bBoxConnection->westLon(), bBoxConnection->southLat(), dataMinZ);
        maxBBox = QVector3D(bBoxConnection->eastLon(), bBoxConnection->northLat(), dataMaxZ);
    }

//    float worldZPos = sceneView->worldZfromPressure(frontElevation2d_hPa);
    // ########################################################################
    // ### BIND BUFFER AND SET UNIFORM VALUES IN GLSL
    // ########################################################################

    for (auto i = 0; i < 2; ++i)
    {
        if (sceneView->shadowMappingInProgress())
        {
            normalCurvesTubeShaderOIT->bindProgram("TubeFilteringShadowMap3D");
        } else
        {
            normalCurvesTubeShaderOIT->bindProgram((i == 0) ? "TubeFilteringShadow3D"
                                               : "TubeFilteringOIT3D");
        }

        normalCurvesTubeShaderOIT->setUniformValue(
                "pToWorldZParams", sceneView->pressureToWorldZParameters());

//        normalCurvesTubeShaderOIT->setUniformValue(
//                "worldZPos",
//                worldZPos);
        normalCurvesTubeShaderOIT->setUniformValue(
                "mvpMatrix", *(sceneView->getModelViewProjectionMatrix())); CHECK_GL_ERROR;
        normalCurvesTubeShaderOIT->setUniformValue(
                "cameraPosition", sceneView->getCamera()->getOrigin()); CHECK_GL_ERROR;

        // Lighting direction from scene view.
        normalCurvesTubeShaderOIT->setUniformValue(
                "lightDirection", sceneView->getLightDirection()); CHECK_GL_ERROR;

        normalCurvesTubeShaderOIT->setUniformValue("bboxMax", maxBBox);
        normalCurvesTubeShaderOIT->setUniformValue("bboxMin", minBBox);

        transferFunctionNCShading->getTexture()->bindToTextureUnit(
                transferFunctionNCShadingTexUnit);
        normalCurvesTubeShaderOIT->setUniformValue(
                "transferFunction", transferFunctionNCShadingTexUnit);

        normalCurvesTubeShaderOIT->setUniformValue(
                "tfMinimum", transferFunctionNCShading->getMinimumValue());
        normalCurvesTubeShaderOIT->setUniformValue(
                "tfMaximum", transferFunctionNCShading->getMaximumValue());

        normalCurvesTubeShaderOIT->setUniformValue(
                "normalized", false);

        normalCurvesTubeShaderOIT->setUniformValue("tubeRadius", normalCurvesTubeRadius); CHECK_GL_ERROR;

        float shadowHeightWorldZ = sceneView->worldZfromPressure(shadowHeight);
        normalCurvesTubeShaderOIT->setUniformValue("shadowColor", shadowColor);
        normalCurvesTubeShaderOIT->setUniformValue("shadowHeight",
                                               shadowHeightWorldZ);

        normalCurvesTubeShaderOIT->setUniformValue("useTFPFilter", useTFPTransferFunction);
        normalCurvesTubeShaderOIT->setUniformValue("useABZFilter", useABZTransferFunction);
        normalCurvesTubeShaderOIT->setUniformValue("useFSFilter", useFSTransferFunction);
        normalCurvesTubeShaderOIT->setUniformValue("useBreadthFilter", useBreadthTransferFunction);

        if (transferFunctionTFP)
        {
            transferFunctionTFP->getTexture()->bindToTextureUnit(transferFunctionTFPTexUnit);
            normalCurvesTubeShaderOIT->setUniformValue("tfTFP", transferFunctionTFPTexUnit); CHECK_GL_ERROR;

            normalCurvesTubeShaderOIT->setUniformValue("tfTFPMinMax",
                                                QVector2D(transferFunctionTFP->getMinimumValue(),
                                                transferFunctionTFP->getMaximumValue()));
        }

        if (transferFunctionABZ)
        {
            transferFunctionABZ->getTexture()->bindToTextureUnit(transferFunctionABZTexUnit);
            normalCurvesTubeShaderOIT->setUniformValue("tfABZ", transferFunctionABZTexUnit); CHECK_GL_ERROR;

            normalCurvesTubeShaderOIT->setUniformValue("tfABZMinMax",
                    QVector2D(transferFunctionABZ->getMinimumValue(),
                              transferFunctionABZ->getMaximumValue()));
        }

        if (transferFunctionFS)
        {
            transferFunctionFS->getTexture()->bindToTextureUnit(transferFunctionFSTexUnit);
            normalCurvesTubeShaderOIT->setUniformValue("tfFS", transferFunctionFSTexUnit); CHECK_GL_ERROR;

            normalCurvesTubeShaderOIT->setUniformValue("tfFSMinMax",
                                                   QVector2D(transferFunctionFS->getMinimumValue(),
                                                             transferFunctionFS->getMaximumValue()));
        }

        if (transferFunctionBreadth)
        {
            transferFunctionBreadth->getTexture()->bindToTextureUnit(transferFunctionBreadthTexUnit);
            normalCurvesTubeShaderOIT->setUniformValue("tfBreadth",
                                             transferFunctionBreadthTexUnit);
            CHECK_GL_ERROR;
            normalCurvesTubeShaderOIT->setUniformValue(
                        "tfBreadthMinMax",
                        QVector2D(transferFunctionBreadth->getMinimumValue(),
                                  transferFunctionBreadth->getMaximumValue()));
        }

        normalCurvesTubeShaderOIT->setUniformValue(
                    "inShadowMappingMode", sceneView->shadowMappingInProgress());

        normalCurvesTubeShaderOIT->setUniformValue(
                    "shadingMode", ncShadingMode);

        sceneView->getShadowMap()->bindToTextureUnit(
                    static_cast<GLuint>(sceneView->getShadowMapTexUnit()));
        normalCurvesTubeShaderOIT->setUniformValue(
                    "shadowMap", sceneView->getShadowMapTexUnit());
        CHECK_GL_ERROR;

        sceneView->setOITUniforms(normalCurvesTubeShaderOIT);
        // ########################################################################
        // ### GET RESOURCE MANAGER AND UPLOAD DATA TO VERTEX BUFFER
        // ########################################################################

        const QString vboIDs = QString("normal_curves_vbo_shading#%1").arg(myID);
        uploadVec2ToVertexBuffer(alphaShading2DNormalCurves, vboIDs,
                                 &vbo2DNormalCurvesShading, sceneView);

        //front vertices
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,
                     ibo2DNormalCurves->getIndexBufferObject()); CHECK_GL_ERROR;
        glBindBuffer(GL_ARRAY_BUFFER,
                     vbo2DNormalCurves->getVertexBufferObject()); CHECK_GL_ERROR;

        int normalCurveVertexSize = sizeof(Geometry::NormalCurveVertex);
        long bytePosition = 0;

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, normalCurveVertexSize,
                              nullptr); CHECK_GL_ERROR;

        bytePosition += 9 * sizeof(float); // tfp, 9 because of start and end points  2 * (QVector3D)
        glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, normalCurveVertexSize,
                              (const GLvoid*) bytePosition); CHECK_GL_ERROR;

        bytePosition += sizeof(float); // abz
        glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, normalCurveVertexSize,
                              (const GLvoid*) bytePosition); CHECK_GL_ERROR;

        bytePosition += sizeof(float); // strength
        glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, normalCurveVertexSize,
                              (const GLvoid*) bytePosition); CHECK_GL_ERROR;

        bytePosition += sizeof(float); // type
        glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, normalCurveVertexSize,
                              (const GLvoid*) bytePosition); CHECK_GL_ERROR;

        bytePosition += sizeof(float); // breadth
        glVertexAttribPointer(5, 1, GL_FLOAT, GL_FALSE, normalCurveVertexSize,
                             (const GLvoid*) bytePosition); CHECK_GL_ERROR;

        glBindBuffer(GL_ARRAY_BUFFER, vbo2DNormalCurvesShading->getVertexBufferObject()); CHECK_GL_ERROR;
        glVertexAttribPointer(6, 1, GL_FLOAT, GL_FALSE, sizeof(QVector2D),
                              nullptr); CHECK_GL_ERROR;

        glVertexAttribPointer(7, 1, GL_FLOAT, GL_FALSE, sizeof(QVector2D),
                              (const GLvoid*) sizeof(float)); CHECK_GL_ERROR;

        glEnableVertexAttribArray(0); CHECK_GL_ERROR;
        glEnableVertexAttribArray(1); CHECK_GL_ERROR;
        glEnableVertexAttribArray(2); CHECK_GL_ERROR;
        glEnableVertexAttribArray(3); CHECK_GL_ERROR;
        glEnableVertexAttribArray(4); CHECK_GL_ERROR;
        glEnableVertexAttribArray(5); CHECK_GL_ERROR;
        glEnableVertexAttribArray(6); CHECK_GL_ERROR;
        glEnableVertexAttribArray(7); CHECK_GL_ERROR;

//        u_int32_t restartIndex = std::numeric_limits<u_int32_t>::max();
        glPrimitiveRestartIndex(restartIndex);
        glEnable(GL_PRIMITIVE_RESTART);

        glEnable(GL_CULL_FACE); CHECK_GL_ERROR;
        glCullFace(GL_FRONT);

        glPolygonMode(GL_BACK, GL_FILL); CHECK_GL_ERROR;
        glDrawElements(GL_LINE_STRIP_ADJACENCY, numNormalCurveSegments,
                       GL_UNSIGNED_INT, nullptr); CHECK_GL_ERROR;

    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);  CHECK_GL_ERROR;
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0); CHECK_GL_ERROR;
    glDisable(GL_CULL_FACE); CHECK_GL_ERROR;
    glDisable(GL_PRIMITIVE_RESTART);
}


void MFrontDetectionActor::triggerAsynchronous3DFrontsRequest()
{
    LOG4CPLUS_DEBUG(mlog, "MFrontDetectionActor::triggerAsynchronous3DFrontsRequest");
    if (variables.size() < 1 || getViews().empty() || isComputing
        || !detectionVariable->hasData()
        || !windUVar->hasData()
        || !windVVar->hasData()
        || !zVar->hasData()
//        || !locatorVar->hasData()
            ) { return; }

    detectionVariable = static_cast<MNWP3DVolumeActorVariable*>(variables.at(detectionVariableIndex));
//    potTemperatureVar =
//            static_cast<MNWP3DVolumeActorVariable*>(variables.at(potTemperatureVariableIndex));
    windUVar = static_cast<MNWP3DVolumeActorVariable*>(variables.at(windUVariableIndex));
    windVVar = static_cast<MNWP3DVolumeActorVariable*>(variables.at(windVVariableIndex));
    zVar = static_cast<MNWP3DVolumeActorVariable*>(variables.at(zVariableIndex));

    if (detectionVariable->grid == nullptr) { return; }
//    if (potTemperatureVar->grid == nullptr) { return; }

    if (!detectionVariable->pendingRequests.empty()
//         || !potTemperatureVar->pendingRequests.empty()
           || !windUVar->pendingRequests.empty()
              || !windVVar->pendingRequests.empty()
                 || !zVar->pendingRequests.empty()
//                    || !locatorVar->pendingRequests.empty()
            ) { return; }

    suppressUpdates = true;
    isComputing = true;

    actorPropertiesSupGroup->setEnabled(false);
    applySettingsClickProperty->setEnabled(false);

    if (frontSurfaceSelection)
    {
        MGLResourcesManager *glRM = MGLResourcesManager::getInstance();
        glRM->releaseGPUItem(vbo3DFronts);
        glRM->releaseGPUItem(ibo3DFronts);
    }

    MDataRequestHelper rh = detectionVariable->constructAsynchronousDataRequest();

    rh.insert("FRONTS_ISOVALUE", QString::number(fleIsoValue));
    rh.insert("FRONTS_WINDU_VAR", windUVar->variableName);
    rh.insert("FRONTS_WINDV_VAR", windVVar->variableName);
    rh.insert("FRONTS_Z_VAR", zVar->variableName);
    rh.insert("FRONTS_NC_OVER_TRACING", QString::number(ncOvertracing));
    rh.insert("FRONTS_NC_INTEGRATION_GPU",  QString::number(1));
    rh.insert("FPMA_DISTANCE", QString::number(fpmaDistanceValue_km));


    partialDerivFilter->setInputSource(detectionVariable->dataSource);
    fleSource->setInputSource(detectionVariable->dataSource);
    tfpSource->setInputSource(detectionVariable->dataSource);
    abzSource->setInputSource(detectionVariable->dataSource);

    if(computeOnGPU)
    {
        frontDetection3DSourceGPU->setFLESource(fleSource.get());
        frontDetection3DSourceGPU->setTFPSource(tfpSource.get());
        frontDetection3DSourceGPU->setABZSource(abzSource.get());
        frontDetection3DSourceGPU->setDetectionVariableSource(
                    detectionVariable->dataSource);
        frontDetection3DSourceGPU->setDetectionVariablePartialDeriveSource(
                    partialDerivFilter.get());

        frontDetection3DSourceGPU->setZSource(zVar->dataSource);

        frontDetection3DSourceGPU->setWindUSource(windUVar->dataSource);
        frontDetection3DSourceGPU->setWindVSource(windVVar->dataSource);

        // Request values
        frontDetection3DSourceGPU->requestData(rh.request());
    }
    else
    {
        frontDetection3DSource->setFLESource(fleSource.get());
        frontDetection3DSource->setTFPSource(tfpSource.get());
        frontDetection3DSource->setABZSource(abzSource.get());
        frontDetection3DSource->setDetectionVariableSource(
                    detectionVariable->dataSource);
        frontDetection3DSource->setDetectionVariablePartialDeriveSource(
                    partialDerivFilter.get());

        frontDetection3DSource->setZSource(zVar->dataSource);

        frontDetection3DSource->setWindUSource(windUVar->dataSource);
        frontDetection3DSource->setWindVSource(windVVar->dataSource);

        // Request values
        frontDetection3DSource->requestData(rh.request());
    }
}


void MFrontDetectionActor::triggerAsynchronous2DFrontsRequest()
{
    // Return if fronts are being computed.
    if (variables.size() < 1 || getViews().empty() || isComputing) {return;}

    // Check if detection variable exists.
    MNWPActorVariable *detectionVar =
            variables.at(detectionVariableIndex);
    if (detectionVar->grid == nullptr) { return; }

    // Freeze front detection settings while computing front lines and normal
    // curves. Will be unfreezed when request is executed in
    // asynchronousFrontLinesAvailable.
    actorPropertiesSupGroup->setEnabled(false);
    applySettingsClickProperty->setEnabled(false);

    isComputing = true;
    suppressUpdates = true;

    if (frontLineSelection)
    {
        // ToDo release data
//        frontDetection2dSource->releaseData(frontLineSelection);
        MGLResourcesManager *glRM = MGLResourcesManager::getInstance();
        glRM->releaseGPUItem(vbo2DFronts);
        glRM->releaseGPUItem(ibo2DFronts);
    }


    fleSource->setInputSource(detectionVar->dataSource);

    tfpSource->setInputSource(detectionVar->dataSource);
    abzSource->setInputSource(detectionVar->dataSource);

    partialDerivFilter->setInputSource(detectionVariable->dataSource);

    frontDetection2DSource->setFLESource(fleSource.get());
    frontDetection2DSource->setTFPSource(tfpSource.get());
    frontDetection2DSource->setABZSource(abzSource.get());
    frontDetection2DSource->setDetectionVariableSource(detectionVar->dataSource);
    frontDetection2DSource->setDetectionVariablePartialDeriveSource(partialDerivFilter.get());
    frontDetection2DSource->setWindUSource(windUVar->dataSource);
    frontDetection2DSource->setWindVSource(windVVar->dataSource);

    MDataRequestHelper rh = detectionVar->
            constructAsynchronousDataRequest();

    rh.insert("FRONTS_ISOVALUE",
                             QString::number(fleIsoValue));
    rh.insert("FRONTS_WINDU_VAR", windUVar->variableName);
    rh.insert("FRONTS_WINDV_VAR", windVVar->variableName);
    rh.insert("FRONTS_ELEVATION", QString::number(frontElevation2d_hPa));
    rh.insert("FRONTS_NC_OVER_TRACING", QString::number(ncOvertracing));
    rh.insert("FRONTS_NC_INTEGRATION_GPU",  QString::number(0));
    rh.insert("FPMA_DISTANCE", QString::number(fpmaDistanceValue_km));
    frontDetection2DSource->requestData(rh.request());

}


void MFrontDetectionActor::filter3DNormalCurves()
{
    if (!frontSurfaceSelection) { return; }

    float dataMinZ, dataMaxZ;
    float westLon, eastLon, southLat, northLat;

    if (bBoxConnection->getBoundingBox() == nullptr)
    {
        dataMinZ = detectionVariable->grid->getBottomDataVolumePressure_hPa();
        dataMaxZ = detectionVariable->grid->getTopDataVolumePressure_hPa();
        westLon = detectionVariable->grid->getWestDataVolumeCorner_lon();
        eastLon = detectionVariable->grid->getEastDataVolumeCorner_lon();
        southLat = detectionVariable->grid->getSouthDataVolumeCorner_lat();
        northLat = detectionVariable->grid->getNorthDataVolumeCorner_lat();
    }
    else
    {
        dataMinZ = bBoxConnection->bottomPressure_hPa();
        dataMaxZ = bBoxConnection->topPressure_hPa();
        westLon = bBoxConnection->westLon();
        eastLon = bBoxConnection->eastLon();
        southLat = bBoxConnection->southLat();
        northLat = bBoxConnection->northLat();
    }


    // Determine seed points grid spacing.
    const float gridSpaceLon = seedPointSpacing;
    const float gridSpaceLat = seedPointSpacing;
    const float gridSpaceHeight = seedPointSpacingZ;

    if ((dataMinZ * dataMaxZ * gridSpaceLon * gridSpaceLat * gridSpaceHeight) <= 0)
    { return; }

    // Compute data extent in lon, lat and height domain.
    const float dataExtentLon = std::abs(eastLon - westLon);
    const float dataExtentLat = std::abs(northLat - southLat);
    const float dataExtentHeight = std::abs(dataMaxZ - dataMinZ);

    const uint32_t numCellsLon = dataExtentLon / gridSpaceLon + 1;
    const uint32_t numCellsLat = dataExtentLat / gridSpaceLat + 1;
    const uint32_t numCellsHeight = dataExtentHeight / gridSpaceHeight + 1;

    QVector<GLint> ghostGrid(numCellsLon * numCellsLat * numCellsHeight, 0);

    QVector3D bboxMin(westLon, southLat, dataMinZ);
    QVector3D bboxMax(eastLon, northLat, dataMaxZ);
    QVector3D boxExtent(dataExtentLon, dataExtentLat, dataExtentHeight);

    // normal curve filtering:
    QVector<Geometry::NormalCurveVertex> normalCurveVertexFiltered;
    QVector<u_int32_t> normalCurveIndexFiltered;
    u_int32_t i = 0;

    for (uint32_t k = 0; k < normalCurves3D->getNumNormalCurves(); ++k)
    {
        Geometry::NormalCurve normalCurve = normalCurves3D->getNormalCurve(k);
        if (normalCurve.positions.empty())
        { continue; }
        QVector3D startPosition = normalCurve.positions.first();

        // compute world z from pressure
//        position.setZ(sceneView->worldZfromPressure(position.z()));

        // check if start position of normal curve is valid and within the
        // current bounding box
        if (startPosition.x() == qNaN || std::isnan(startPosition.x())) { continue; }

        if (startPosition.x() < bboxMin.x() || startPosition.x() > bboxMax.x()
                || startPosition.y() < bboxMin.y() || startPosition.y() > bboxMax.y()
                || startPosition.z() < bboxMin.z() || startPosition.z() > bboxMax.z())
        { continue; }

        // compute indices
        QVector3D normTexCoords = (startPosition - bboxMin) / boxExtent;

        int texCoordX = int(normTexCoords.x() * (numCellsLon - 1));
        int texCoordY = int(normTexCoords.y() * (numCellsLat - 1));
        int texCoordZ = int(normTexCoords.z() * (numCellsHeight - 1));

        const int cellIndex = INDEX3zyx(texCoordZ, texCoordY, texCoordX, numCellsLat, numCellsLon);

        // critical section
        int ghostVisited = ghostGrid[cellIndex];
//        const float MIN_OPACITY = 0.1;

        if (!ghostVisited)
        {
            // compute alpha value at current vertex position
            const float tfp = normalCurve.tfp;

            if (tfp < 0 && !showColdSideFront) { continue; }

            QVector3D start = normalCurve.positions.first();
            QVector3D end = normalCurve.positions.last();

            for (int p = 0; p < normalCurve.positions.size(); p++)
            {
                normalCurveIndexFiltered.append(i);
                Geometry::NormalCurveVertex segment;
                segment.position = normalCurve.positions.at(p);
                segment.start = start;
                segment.end = end;
                segment.tfp = normalCurve.tfp;
                segment.abz = normalCurve.abz;
                segment.strength = normalCurve.strength;
                segment.type = normalCurve.type;
                segment.breadth = normalCurve.breadth;
                normalCurveVertexFiltered.append(segment);
                i++;
            }
            normalCurveIndexFiltered.append(restartIndex);
            ghostGrid[cellIndex] = 1;
        }
    }
    // change format from QVector to array pointer
    uint32_t numNormalCurveSegments = normalCurveVertexFiltered.size();
    Geometry::NormalCurveVertex *normalCurveSegment =
            new Geometry::NormalCurveVertex[numNormalCurveSegments];

    for (uint32_t k = 0; k < numNormalCurveSegments; k++)
    {
        normalCurveSegment[k] = normalCurveVertexFiltered.at(k);
    }

    uint32_t numRestartIndex = normalCurveIndexFiltered.size();
    uint32_t *restartIndex = new uint32_t[numRestartIndex];

    for (uint32_t k = 0; k < numRestartIndex; k++)
    {
        restartIndex[k] = normalCurveIndexFiltered.at(k);
    }

    normalCurves3D->setNormalCurveSegments(numNormalCurveSegments,
                                           normalCurveSegment,
                                           numRestartIndex,
                                           restartIndex);

    normalCurves3D->releaseVertexBuffer();
    normalCurves3D->releaseIndexBuffer();

    vbo3DNormalCurves = normalCurves3D->getVertexBuffer();
    ibo3DNormalCurves = normalCurves3D->getIndexBuffer();
    //num3DNormalCurveVertices = normalCurves3D->getNumNormalCurveSegments();
}


void MFrontDetectionActor::filter2DNormalCurves()
{
    if (!frontLineSelection
         || num2DFrontVertices < 1) { return; }

    // Determine seed points grid spacing.
    const float gridSpaceLon = seedPointSpacing;
    const float gridSpaceLat = seedPointSpacing;

    if ((gridSpaceLon * gridSpaceLat) <= 0)
    { return; }

    float westLon, eastLon, southLat, northLat;

    if (bBoxConnection->getBoundingBox() == nullptr)
    {
        westLon = detectionVariable->grid->getWestDataVolumeCorner_lon();
        eastLon = detectionVariable->grid->getEastDataVolumeCorner_lon();
        southLat = detectionVariable->grid->getSouthDataVolumeCorner_lat();
        northLat = detectionVariable->grid->getNorthDataVolumeCorner_lat();
    }
    else
    {
        westLon = bBoxConnection->westLon();
        eastLon = bBoxConnection->eastLon();
        southLat = bBoxConnection->southLat();
        northLat = bBoxConnection->northLat();
    }

    // Compute data extent in lon, lat and height domain.
    const float dataExtentLon = std::abs(eastLon - westLon);
    const float dataExtentLat = std::abs(northLat - southLat);

    const uint32_t numCellsLon = dataExtentLon / gridSpaceLon + 1;
    const uint32_t numCellsLat = dataExtentLat / gridSpaceLat + 1;

    QVector<GLint> ghostGrid(numCellsLon * numCellsLat, 0);

    QVector2D bboxMin(westLon, southLat);
    QVector2D bboxMax(eastLon, northLat);
    QVector2D boxExtent(dataExtentLon, dataExtentLat);

    // normal curve filtering:
    QVector<Geometry::NormalCurveVertex> normalCurveVertexFiltered;
    QVector<u_int32_t> normalCurveIndexFiltered;
    u_int32_t i = 0;

    for (uint32_t k = 0; k < normalCurves2D->getNumNormalCurves(); ++k)
    {
        Geometry::NormalCurve normalCurve = normalCurves2D->getNormalCurve(k);
        if (normalCurve.positions.empty())
        { continue; }
        QVector3D position = normalCurve.positions.first();

        // check if start position of normal curve is valid and within the
        // current bounding box
        if (position.x() == qNaN || std::isnan(position.x())) { continue; }

        if (position.x() < bboxMin.x() || position.x() > bboxMax.x()
                || position.y() < bboxMin.y() || position.y() > bboxMax.y())
        { continue; }

        // compute indices
        QVector2D normTexCoords = (QVector2D(position.x(),  position.y()) - bboxMin) / boxExtent;

        int texCoordX = int(normTexCoords.x() * (numCellsLon - 1));
        int texCoordY = int(normTexCoords.y() * (numCellsLat - 1));

        const int cellIndex = INDEX2yx(texCoordY, texCoordX, numCellsLon);

        // critical section
        int ghostVisited = ghostGrid[cellIndex];
//        const float MIN_OPACITY = 0.1;

        if (!ghostVisited)
        {
            // compute alpha value at current vertex position
            const float tfp = normalCurve.tfp;

            if (tfp < 0 && !showColdSideFront) { continue; }

            QVector3D start = normalCurve.positions.first();
            QVector3D end = normalCurve.positions.last();

            for (int p = 0; p < normalCurve.positions.size(); p++)
            {
                normalCurveIndexFiltered.append(i);
                Geometry::NormalCurveVertex segment;
                segment.position = normalCurve.positions.at(p);
                segment.start = start;
                segment.end = end;
                segment.tfp = normalCurve.tfp;
                segment.abz = normalCurve.abz;
                segment.strength = normalCurve.strength;
                segment.type = normalCurve.type;
                segment.breadth = normalCurve.breadth;
                normalCurveVertexFiltered.append(segment);
                i++;
            }
            normalCurveIndexFiltered.append(restartIndex);
            ghostGrid[cellIndex] = 1;
        }
    }
    // change format from QVector to array pointer
    uint32_t numNormalCurveSegments = normalCurveVertexFiltered.size();
    Geometry::NormalCurveVertex *normalCurveSegment =
            new Geometry::NormalCurveVertex[numNormalCurveSegments];

    for (uint32_t k = 0; k < numNormalCurveSegments; k++)
    {
        normalCurveSegment[k] = normalCurveVertexFiltered.at(k);
    }

    uint32_t numRestartIndex = normalCurveIndexFiltered.size();
    uint32_t *restartIndex = new uint32_t[numRestartIndex];

    for (uint32_t k = 0; k < numRestartIndex; k++)
    {
        restartIndex[k] = normalCurveIndexFiltered.at(k);
    }

    normalCurves2D->setNormalCurveSegments(numNormalCurveSegments,
                                           normalCurveSegment,
                                           numRestartIndex,
                                           restartIndex);

    normalCurves2D->releaseVertexBuffer();
    normalCurves2D->releaseIndexBuffer();

    vbo2DNormalCurves = normalCurves2D->getVertexBuffer();
    ibo2DNormalCurves = normalCurves2D->getIndexBuffer();
}


void MFrontDetectionActor::setNCFilterTransferFunctionFromProperty(
        NormalCurvesFilterSettings* filter)
{
    MGLResourcesManager *glRM = MGLResourcesManager::getInstance();

    QString tfName = properties->getEnumItem(filter->transferFunctionProperty);

    if (tfName == "None")
    {
        filter->transferFunction = nullptr;
        return;
    }

    // Find the selected transfer function in the list of actors from the
    // resources manager. Not very efficient, but works well enough for the
    // small number of actors at the moment..
    foreach (MActor *actor, glRM->getActors())
    {
        if (MTransferFunction1D *tf = dynamic_cast<MTransferFunction1D*>(actor))
        {
            if (tf->transferFunctionName() == tfName)
            {
                filter->transferFunction  = tf;
                return;
            }
        }
    }
}


void MFrontDetectionActor::refreshEnumsProperties(MNWPActorVariable *var)
{
    Q_UNUSED(var);
}


QString MFrontDetectionActor::filterTypeToString(filterType type)
{
    {
        switch (type)
        {
            case ABSOLUTE_CHANGE: return "absoluteChangeWithinFrontalZone";
            case PER_100KM: return "changePer100km";
            case VALUE_AT_FRONT: return "valueAtFrontline";
        }
        return "changePer100km";
    }
}


MFrontDetectionActor::filterType MFrontDetectionActor::stringToFilterType(
        QString type)
{
    if (type == "absoluteChangeWithinFrontalZone")
    {
        return ABSOLUTE_CHANGE;
    }
    else if (type == "changePer100km")
    {
        return PER_100KM;
    }
    else if (type == "valueAtFrontline")
    {
        return VALUE_AT_FRONT;
    }
    else
    {
        return PER_100KM;
    }
}
