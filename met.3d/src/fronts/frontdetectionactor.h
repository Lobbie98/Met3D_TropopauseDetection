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

#ifndef MET_3D_FRONTDETECTIONACTOR_H
#define MET_3D_FRONTDETECTIONACTOR_H

// standard library imports
#include <memory>

#include <QtCore>
#include <src/util/mmarchingcubes.h>

#include "data/datarequest.h"
#include "gxfw/nwpmultivaractor.h"
#include "gxfw/nwpactorvariableproperties.h"

#include "gxfw/gl/shadereffect.h"
#include "gxfw/gl/indexbuffer.h"
#include "gxfw/gl/typedvertexbuffer.h"
#include "gxfw/gl/shaderstoragebufferobject.h"

#include "data/smoothfilter.h"
#include "data/partialderivativefilter.h"

#include "fronts/vectormagnitudefilter.h"
#include "fronts/frontdetection3Dsource.h"
#include "fronts/frontdetection3DsourceGPU.h"
#include "fronts/frontlocationequationsource.h"
#include "fronts/thermalfrontparametersource.h"
#include "fronts/adjacentbarocliniczonesource.h"
#include "fronts/frontdetection2Dsource.h"

namespace Met3D
{
class MFrontDetectionActor : public MNWPMultiVarActor,
                               public MBoundingBoxInterface
{
Q_OBJECT

public:

    explicit MFrontDetectionActor();
    ~MFrontDetectionActor() override;

    QString getSettingsID() override
    { return "FrontDetectionActor3D"; }

    void reloadShaderEffects() override;

    void saveConfiguration(QSettings *settings) override;
    void loadConfiguration(QSettings *settings) override;

    MNWPActorVariable* createActorVariable(
            const MSelectableDataSource& dataSource) override;

    const QList<MVerticalLevelType> supportedLevelTypes() override;

    void onBoundingBoxChanged() override;

public slots:
    void asynchronous3DFrontsAvailable(MDataRequest request);

    void asynchronous2DFrontsAvailable(MDataRequest request);

    void onAddActorVariable(MNWPActorVariable *var) override;

    void onDeleteActorVariable(MNWPActorVariable *var) override;

    void onChangeActorVariable(MNWPActorVariable *var) override;

    void onActorCreated(MActor *actor);

    void onActorDeleted(MActor *actor);

    void onActorRenamed(MActor *actor, QString oldName);

    void dataFieldChangedEvent() override;

    void updateShadow();

    void recomputeFrontsAlpha();

private:
    void initializeActorResources() override;

    void onQtPropertyChanged(QtProperty *property) override;

    void renderToCurrentContext(MSceneViewGLWidget *sceneView) override;

    void render3DFrontSurfaceOIT(MSceneViewGLWidget *sceneView);

    void render3DNormalCurvesOIT(MSceneViewGLWidget *sceneView);

    void render2DFrontTubesOIT(MSceneViewGLWidget *sceneView);

    void render2DNormalCurvesOIT(MSceneViewGLWidget *sceneView);

    void renderTransparencyToCurrentContext(MSceneViewGLWidget *sceneView) override;

    void triggerAsynchronous3DFrontsRequest();

    void triggerAsynchronous2DFrontsRequest();

    void filter3DNormalCurves();

    void filter2DNormalCurves();

    bool setTransferFunction(QtProperty* tfProp, const QString &tfName);

    void setTransferFunctionFromProperty(
            QtProperty* tfProp, MTransferFunction1D** transferFunc);

    void addNormalCurvesFilter(int8_t filterNumber);

    typedef enum {
        ABSOLUTE_CHANGE = 0,
        PER_100KM = 1,
        VALUE_AT_FRONT = 2,
    } filterType;

    void loadNormalCurvesFilter(
            int8_t filterNumber, bool enabled, int varIndex, QString tfName,
            filterType type);
    /**
    * Refresh all variables select boxes if a variable was added or removed.
     */
    void refreshEnumsProperties(MNWPActorVariable *var);

    static QString filterTypeToString(filterType type);

    static filterType stringToFilterType(QString type);

    struct NormalCurvesFilterSettings
    {
        NormalCurvesFilterSettings(MFrontDetectionActor* a,
                                   uint8_t filterNumber,
                                   QList<QString> variableNames,
                                   bool enabled = false,
                                   int variableIndex = -1,
                                   MTransferFunction1D
                                   *transferFunction = nullptr,
                                   filterType = PER_100KM,
                                   int textureUnitTransferFunction = -1);
    public:
        bool                 enabled;
        int8_t               variableIndex;
        MTransferFunction1D *transferFunction;
        int                  textureUnitTransferFunction;
        QList<QString>       varNameList;
        filterType           type;

        QtProperty *groupProperty;
        QtProperty *enabledProperty;
        QtProperty *variableProperty;
        QtProperty *transferFunctionProperty;
        QtProperty *typeProperty;
        QtProperty *removeProperty;
    };

    friend struct   NormalCurvesFilterSettings;

    bool removeNormalCurvesFilter(NormalCurvesFilterSettings *filter);

    void setNCFilterTransferFunctionFromProperty(
            NormalCurvesFilterSettings* filter);

    bool suppressUpdates;
    bool isComputing;

    // List that stores the names of all registered actor variables (for the
    // enum properties that allow the user to select observed and shading
    // variables.
    QList<QString>              varNameList;

    // ********************* Structure of front properties *********************
    // |-> actor properties
    //      |-- update
    //      |-- auto update
    //      |-> input variables
    //          |-- detection var
    //          |-- eastward wind (u)
    //          |-- northward wind (v)
    //          |-- geopotential height (only 3D fronts)
    //      |-> display options
    //          |-- 3D fronts
    //          |-- 3D normal curves
    //          |-- 2D fronts
    //          |-- 2D normal curves
    //      |-> front filter options
    //          |-- FPMA Distance of next axis
    //          |-- TFP transfer function
    //          |-> generic filter
    //              |-- add filter
    //              |-> filter #1
    //                  |-- variable
    //                  |-- transfer function
    //                  |-- type
    //                  |-- remove
    //          |-> optional filter
    //              |-- ABZ transfer function
    //              |-- breadth of frontal zone transfer function
    //              |-- slope of frontal surface (3D) transfer function
    //              |-- normal curve over-tracing
    //              |-- minimum front length (2D)
    //              |-- iso value front location equation
    //              |-- TFP threshold (2D)
    //      |-> front appearance
    //          |-> shading
    //              |-- mode
    //              |-- variable
    //              |-- type
    //              |-- transfer function
    //          |-> appearance 2D fronts
    //              |-- elevation
    //              |-- tube radius
    //          |-- bounding box
    //          |-- show warm/cold fronts
    //          |-- show cold side of fronts
    //          |-> light and shadow
    //              |-- shadow elevation
    //              |-- shadow color
    //              |-- light mode

    //      |-> normal curve appearance
    // TODO: (AB, 2021-03-09) change options for normal curve appearance


    // ************************ MAIN PROPERTIES ********************************
    //      |-- update
    QtProperty* applySettingsClickProperty;

    //      |-- auto update
    QtProperty* autoComputeProperty;
    bool        autoCompute;

    QtProperty* computeOnGPUProperty;
    bool        computeOnGPU;

    // ************************ INPUT VAR PROPERTIES ***************************
    //      |-> input variables
    QtProperty* inputVarGroupProperty;

    //          |-- detection var
    QtProperty*                 detectionVariableIndexProperty;
    uint8_t                     detectionVariableIndex;
    MNWP3DVolumeActorVariable*  detectionVariable;

    //          |-- eastward wind (u)
    // at the moment only for 3D fronts
    QtProperty*                 windUVarIndexProp;
    uint8_t                     windUVariableIndex;
    MNWP3DVolumeActorVariable*  windUVar;

    //          |-- northward wind (v)
    // at the moment only for 3D fronts
    QtProperty*                 windVVarIndexProp;
    uint8_t                     windVVariableIndex;
    MNWP3DVolumeActorVariable*  windVVar;

    //          |-- geopotential height (only 3D fronts)
    QtProperty*                 zVarIndexProp;
    uint8_t                     zVariableIndex;
    MNWP3DVolumeActorVariable*  zVar;


    // ************************ DISPLAY OPTIONS ********************************
    //      |-> display options
    QtProperty* displayOptionsGroupProperty;

    //          |-- 3D fronts
    QtProperty* render3DFrontProperty;
    bool        render3DFront;

    //          |-- 3D normal curves
    QtProperty* render3DNCProperty;
    bool        render3DNC;

    //          |-- 2D fronts
    QtProperty* render2DFrontProperty;
    bool        render2DFront;

    //          |-- 2D normal curves
    QtProperty* render2DNCProperty;
    bool        render2DNC;


    // ************************ FILTER OPTIONS *********************************
    //      |-> front filter options
    QtProperty* filterGroupProperty;

    //          |-- FPMA Distance of next axis
    QtProperty* fpmaDistanceProperty;
    double      fpmaDistanceValue_km;

    //          |-- TFP transfer function
    QtProperty*                 transferFunctionTFPProperty;
    MTransferFunction1D*        transferFunctionTFP;
    int                         transferFunctionTFPTexUnit;

    //          |-> Frontal strength transfer function
    QtProperty*                 transferFunctionFSProperty;
    MTransferFunction1D*        transferFunctionFS;
    int                         transferFunctionFSTexUnit;

    //          |-> generic filter
    QtProperty* genericFilterGroupProperty;

    //              |-- add filter
    //              |-> filter #1
    //                  |-- variable
    //                  |-- transfer function
    //                  |-- type
    //                  |-- remove
    QtProperty* addNormalCurveFilterProperty;
    QVector<NormalCurvesFilterSettings*> normalCurvesFilterSetList;

    //          |-> optional filter
    QtProperty* optionalFilterProperty;

    //              |-- ABZ transfer function
    QtProperty*                 transferFunctionABZProperty;
    MTransferFunction1D*        transferFunctionABZ;
    int                         transferFunctionABZTexUnit;

    //              |-- breadth of frontal zone transfer function
    QtProperty*                 transferFunctionBreadthProperty;
    MTransferFunction1D*        transferFunctionBreadth;
    int                         transferFunctionBreadthTexUnit;

    //              |-- slope of frontal surface (3D) transfer function
    QtProperty*                 transferFunctionSlopeProperty;
    MTransferFunction1D*        transferFunctionSlope;
    int                         transferFunctionSlopeTexUnit;

    //              |-- normal curve over-tracing
    QtProperty     *ncOvertracingProperty;
    int             ncOvertracing;

    //              |-- iso value front location equation
    QtProperty* fleIsoProperty;
    double      fleIsoValue;

    // ************************ FRONT APPEARANCE *******************************
    //      |-> front appearance
    QtProperty* appearanceGroupProperty;

    //          |-> shading
    QtProperty* shadingGroupProperty;

    //              |-- mode
    QtProperty*  shadingModeProperty;
    int          shadingMode;

    //              |-- variable
    QtProperty*                 shadingVariableIndexProperty;
    uint8_t                     shadingVariableIndex;
    MNWP3DVolumeActorVariable*  shadingVariable;

    //              |-- type
    QtProperty*  shadingVarModeProperty;
    filterType   shadingVarMode;

    //              |-- transfer function
    QtProperty*                 transferFunctionShadingProperty;
    MTransferFunction1D*        transferFunctionShading;
    int                         transferFunctionShadingTexUnit;

    //          |-> appearance 2D fronts
    QtProperty* shading2dGroupProperty;

    //              |-- elevation
    QtProperty* frontElevationProperty;
    double      frontElevation2d_hPa;

    //              |-- tube radius
    QtProperty*  frontsTubeRadiusProp;
    float        frontsTubeRadius;

    //          |-- bounding box
    // Initialized in cpp-file

    //          |-- show warm/cold fronts
    QtProperty* showFrontTypesProperty;
    bool        showFrontTypes;

    //          |-- show cold side of fronts
    QtProperty* showColdSideFrontProperty;
    bool        showColdSideFront;

    //          |-> light and shadow
    QtProperty* lightShadowGroupProperty;

    //              |-- shadow elevation
    QtProperty*                 shadowHeightProperty;
    float                       shadowHeight;

    //              |-- shadow color
    QtProperty*                 shadowColorProperty;
    QColor                      shadowColor;

    //              |-- light mode
    QtProperty*                 lightingModeProperty;
    GLint                       lightingMode;

    // ************************ NORMAL CURVES PROPERTIES ***********************
    //      |-> normal curve appearance
    // TODO: (AB, 2021-03-09) change options for normal curve appearance

    // TODO (AB,01.2021): do we need this extra computation of normal curves?
    QtProperty*                 normalCurvesGroupProp;

    bool        trigger3DNormalCurveFilter;
    bool        trigger2DNormalCurveFilter;

    QtProperty*                 normalCurvesTubeRadiusProp;
    float                       normalCurvesTubeRadius;

    QtProperty*                 ncShadingVarIndexProp;
    uint8_t                     ncShadingVariableIndex;
    MNWP3DVolumeActorVariable*  ncShadingVar;

    QtProperty*                 seedPointSpacingProp;
    float                       seedPointSpacing;

    QtProperty*                 seedPointSpacingZProp;
    int                         seedPointSpacingZ;

    QtProperty*                 ncShadingModeProperty;
    int                         ncShadingMode;

    QtProperty*                 transferFunctionNCShadingProperty;
    MTransferFunction1D*        transferFunctionNCShading;
    int                         transferFunctionNCShadingTexUnit;

// ********************** DEPRECREATED PROPERTIES ******************************
    bool            useTFPTransferFunction;
    bool            useFSTransferFunction;
    bool            useABZTransferFunction;
    bool            useSlopeTransferFunction;
    bool            useBreadthTransferFunction;

    std::shared_ptr<MFrontLocationEquationSource>   fleSource;
    std::shared_ptr<MThermalFrontParameterSource>   tfpSource;
    std::shared_ptr<MAdjacentBaroclinicZoneSource>  abzSource;
    std::shared_ptr<MPartialDerivativeFilter>       partialDerivFilter;
    std::shared_ptr<MFrontDetection2DSource>        frontDetection2DSource;
    std::shared_ptr<MFrontDetection3DSource>        frontDetection3DSource;
    std::shared_ptr<MFrontDetection3DSourceGPU>     frontDetection3DSourceGPU;

    M3DFrontSelection*                              frontSurfaceSelection;
    M2DFrontSelection*                                 frontLineSelection;

    MTriangleMeshSelection*                         frontMesh3D;
    MLineSelection*                                 frontLines2D;
    MNormalCurvesSelection*                         normalCurves3D;
    MNormalCurvesSelection*                         normalCurves2D;

    // shaders to render fronts
    std::shared_ptr<GL::MShaderEffect>  fronts3DSurfaceShader;
    std::shared_ptr<GL::MShaderEffect>  frontTubeShaderOIT;
    std::shared_ptr<GL::MShaderEffect>  normalCurvesTubeShaderOIT;

    GL::MVertexBuffer*                  vbo3DFronts;
    GL::MIndexBuffer*                   ibo3DFronts;
    GL::MVertexBuffer*                  vbo3DFrontsShading;

    GL::MVertexBuffer*                  vbo3DNormalCurves;
    GL::MIndexBuffer*                   ibo3DNormalCurves;
    GL::MVertexBuffer*                  vbo3DNormalCurvesShading;

    GL::MVertexBuffer*                  vbo2DFronts;
    GL::MIndexBuffer*                   ibo2DFronts;
    GL::MVertexBuffer*                  vbo2DFrontsShading;

    GL::MVertexBuffer*                  vbo2DNormalCurves;
    GL::MIndexBuffer*                   ibo2DNormalCurves;
    GL::MVertexBuffer*                  vbo2DNormalCurvesShading;


    std::vector<QList<QString>>               sampleSubroutines;
    std::vector<QList<QString>>               normalCompSubroutines;

    QVector<QVector2D> alphaShading3DFronts;
    QVector<QVector2D> alphaShading3DNormalCurves;

    QVector<QVector2D> alphaShading2DFronts;
    QVector<QVector2D> alphaShading2DNormalCurves;


    bool prevComputedOnGPU;
    bool recomputeAlpha3DFronts;
    bool recomputeShading3DFronts;
    bool recomputeAlpha3DNormalCurves;
    bool recomputeShading3DNormalCurves;
    int num3DFrontVertices;

    bool recomputeAlpha2DFronts;
    bool recomputeShading2DFronts;
    bool recomputeAlpha2DNormalCurves;
    bool recomputeShading2DNormalCurves;
    int num2DFrontVertices;

    u_int32_t restartIndex;

};

class MFrontDetectionActorFactory : public MAbstractActorFactory
{
public:
    MFrontDetectionActorFactory() : MAbstractActorFactory() {}

protected:
    MActor* createInstance() override { return new MFrontDetectionActor(); }
};

} // namespace Met3D



#endif //MET_3D_FRONTDETECTIONACTOR_H
