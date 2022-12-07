/******************************************************************************
**
**  This file is part of Met.3D -- a research environment for the
**  three-dimensional visual exploration of numerical ensemble weather
**  prediction data.
**
**  Copyright 2015-2021 Marc Rautenhaus
**  Copyright 2022      Julian Münsterberg
**
**  Regional Computing Center, Visualization
**  Universitaet Hamburg, Hamburg, Germany
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

#include "tropopausedetectionactor.h"

#include <QFileDialog>

// related third party imports
#include <netcdf>

using namespace Met3D;
using namespace netCDF;
using namespace netCDF::exceptions;

#include <log4cplus/loggingmacros.h>

MTropopauseDetectionActor::MTropopauseDetectionActor():
          MNWPMultiVarActor(),
          MBoundingBoxInterface(this, MBoundingBoxConnectionType::VOLUME),
          suppressUpdates(false),
          isComputing(false),
          applySettingsClickProperty(nullptr),
          inputVarGroupProperty(nullptr),
          detectionVariableIndexProperty(nullptr),
          detectionVariableIndex(0),
          detectionVariable(nullptr),
          displayOptionsGroupProperty(nullptr),
          renderTropopauseProperty(nullptr),
          renderTropopause(false),
          transferFunctionFDProperty(nullptr),
          transferFunctionFD(nullptr),
          transferFunctionFDTexUnit(-1),
          tropopauseIsoProperty(nullptr),
          tropopauseIsoValue(0.0),
          tropopauseColourProperty(nullptr),
          tropopauseColour(QColor(0, 0, 0, 255)),
          lightShadowGroupProperty(nullptr),
          shadowHeightProperty(nullptr),
          shadowHeight(1045),
          shadowColorProperty(nullptr),
          shadowColor(QColor(70, 70, 70, 150)),
          lightingModeProperty(nullptr),
          lightingMode(0),
          useFDTransferFunction(false),
          tropopauseDetectionSource(nullptr),
          tropopauseMesh3D(nullptr),
          tropopauseSurfaceShader(nullptr),

          vboTropopause(nullptr),
          iboTropopause(nullptr),
          vboTropopauseShading(nullptr),

          alphaShadingTropopause(0),

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

    //Properties initialisieren
    beginInitialiseQtProperties();

    setActorType("Tropopause Detection Actor");
    setName(getActorType());

    /**************************************************************************
                            MAIN PROPERTIES
    **************************************************************************/
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
                "example: should be temperature");


    /**************************************************************************
                            Render OPTIONS
    **************************************************************************/

    renderTropopauseProperty = addProperty(BOOL_PROPERTY, "Tropopause",
                                        actorPropertiesSupGroup);
    properties->mBool()->setValue(renderTropopauseProperty, renderTropopause);


    tropopauseIsoProperty = addProperty(SCIENTIFICDOUBLE_PROPERTY, "iso value",
                                 actorPropertiesSupGroup);
    properties->setSciDouble(tropopauseIsoProperty, tropopauseIsoValue,
                             -1.00, 1.00, 2, 0.01);
    tropopauseIsoProperty->setToolTip("Iso value of Tropopause location. \n"
                                       "Default and recommendation value = 0.0");

    tropopauseColourProperty = addProperty(
                COLOR_PROPERTY, "tropopause colour", actorPropertiesSupGroup);
    properties->mColor()->setValue(tropopauseColourProperty, tropopauseColour);

    transferFunctionFDProperty = addProperty(ENUM_PROPERTY,
                                              "transfer function firstDeriv filter",
                                              actorPropertiesSupGroup);
    properties->mEnum()->setEnumNames(transferFunctionFDProperty, availableTFs);
    transferFunctionFDProperty->setToolTip(
                "Set alpha value of transfer function for fuzzy filtering.");

    actorPropertiesSupGroup->addSubProperty(bBoxConnection->getProperty());


    /**************************************************************************
                            TROPOPAUSE APPEARANCE PROPERTIES
    **************************************************************************/
    appearanceGroupProperty = addProperty(GROUP_PROPERTY,
                                                 "tropopause appearance",
                                                 actorPropertiesSupGroup);

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

    connect(glRM, SIGNAL(actorCreated(MActor*)),
            SLOT(onActorCreated(MActor*)));
    connect(glRM, SIGNAL(actorDeleted(MActor*)),
            SLOT(onActorDeleted(MActor*)));
    connect(glRM, SIGNAL(actorRenamed(MActor*,QString)),
            SLOT(onActorRenamed(MActor*,QString)));


    endInitialiseQtProperties();
}

MTropopauseDetectionActor::~MTropopauseDetectionActor()
{

}

/******************************************************************************
***                            PUBLIC METHODS                               ***
*******************************************************************************/

#define SHADER_VERTEX_ATTRIBUTE 0

void MTropopauseDetectionActor::reloadShaderEffects()
{
    LOG4CPLUS_DEBUG(mlog, "loading shader programs");

    beginCompileShaders(2);

    compileShadersFromFileWithProgressDialog(
                tropopauseSurfaceShader, "src/glsl/tropopause_geometry.fx.glsl");

    endCompileShaders();
}

void MTropopauseDetectionActor::saveConfiguration(QSettings *settings)
{
    LOG4CPLUS_DEBUG(mlog, "saveConfiguration");
    MNWPMultiVarActor::saveConfiguration(settings);

    settings->beginGroup(getSettingsID());

    MBoundingBoxInterface::saveConfiguration(settings);

    settings->setValue("renderTropopause", renderTropopause);
    settings->setValue("detectionVariableIndex", detectionVariableIndex);
    settings->setValue("tropopauseIsoValue", tropopauseIsoValue);
    settings->setValue("tropopauseColour", tropopauseColour);
    settings->setValue("useFDTransferFunction", useFDTransferFunction);


    settings->setValue("transferFunctionFDProperty",
                       properties->getEnumItem(
                               transferFunctionFDProperty));

    settings->setValue("shadowColor", shadowColor);
    settings->setValue("shadowHeight", shadowHeight);
    settings->setValue("lightingMode", lightingMode);
    settings->endGroup();
}


void MTropopauseDetectionActor::loadConfiguration(QSettings *settings)
{

    LOG4CPLUS_DEBUG(mlog, "loadConfiguration");
    MNWPMultiVarActor::loadConfiguration(settings);

    suppressUpdates = true;

    settings->beginGroup(getSettingsID());

    MBoundingBoxInterface::loadConfiguration(settings);

    renderTropopause = settings->value("renderTropopause").toBool();
    properties->mBool()->setValue(renderTropopauseProperty, renderTropopause);

    detectionVariableIndex = settings->value("detectionVariableIndex").toInt();
    properties->mInt()->setValue(detectionVariableIndexProperty,
                                 detectionVariableIndex);


    tropopauseIsoValue = settings->value("tropopauseIsoValue").toDouble();
    properties->mDouble()->setValue(tropopauseIsoProperty, tropopauseIsoValue);

    tropopauseColour = settings->value("tropopauseColour").value<QColor>();
    properties->mColor()->setValue(tropopauseColourProperty, tropopauseColour);

    useFDTransferFunction = settings->value("useFDTransferFunction").toBool();

    QString tfName = settings->value("transferFunctionFDProperty", "None").toString();
    if (!setTransferFunction(transferFunctionFDProperty, tfName))
    {
    }

    setTransferFunctionFromProperty(transferFunctionFDProperty,
                                    &transferFunctionFD);

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
    }

    settings->endGroup();

    suppressUpdates = false;
    emitActorChangedSignal();
}

bool MTropopauseDetectionActor::setTransferFunction(QtProperty* tfProp, const QString &tfName)
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

    return false;
}


void MTropopauseDetectionActor::setTransferFunctionFromProperty(
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

MNWPActorVariable* MTropopauseDetectionActor::createActorVariable(
        const MSelectableDataSource& dataSource)
{
    MNWP3DVolumeActorVariable *newVar = new MNWP3DVolumeActorVariable(this);

    newVar->dataSourceID = dataSource.dataSourceID;
    newVar->levelType = dataSource.levelType;
    newVar->variableName = dataSource.variableName;

    return newVar;
}

const QList<MVerticalLevelType> MTropopauseDetectionActor::supportedLevelTypes()
{
    return (QList<MVerticalLevelType>()
            << HYBRID_SIGMA_PRESSURE_3D << PRESSURE_LEVELS_3D << AUXILIARY_PRESSURE_3D);
}

void MTropopauseDetectionActor::onBoundingBoxChanged()
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

void MTropopauseDetectionActor::asynchronousTropopauseAvailable(MDataRequest request)
{
    if (tropopauseMesh3D)
    {
        tropopauseMesh3D->releaseIndexBuffer();
        tropopauseMesh3D->releaseVertexBuffer();

        tropopauseDetectionSource->releaseData(tropopauseMesh3D);
    }
    tropopauseMesh3D = tropopauseDetectionSource->getData(request);

    numtropopauseVertices = tropopauseMesh3D->getNumVertices();

    vboTropopause = tropopauseMesh3D->getVertexBuffer();
    iboTropopause = tropopauseMesh3D->getIndexBuffer();
    applySettingsClickProperty->setEnabled(true);
    actorPropertiesSupGroup->setEnabled(true);

    suppressUpdates = false;
    isComputing = false;

    emitActorChangedSignal();
}

void MTropopauseDetectionActor::onAddActorVariable(MNWPActorVariable *var)
{
    varNameList << var->variableName;

    // Temporarily save variable indices.
    int tmpVarIndex = detectionVariableIndex;

    // Update enum lists.
    properties->mEnum()->setEnumNames(detectionVariableIndexProperty, varNameList);
    properties->mEnum()->setValue(detectionVariableIndexProperty, tmpVarIndex);

    refreshEnumsProperties(nullptr);
}

void MTropopauseDetectionActor::onDeleteActorVariable(MNWPActorVariable *var)
{

    int i = variables.indexOf(var);

    // Update variableIndex and shadingVariableIndex if these point to
    // the removed variable or to one with a lower index.
    if (i <= detectionVariableIndex)
    {
        detectionVariableIndex = std::max(-1, detectionVariableIndex - 1);
    }

    // Temporarily save variable indices.
    int tmpVarIndex = detectionVariableIndex;

    // Remove the variable name from the enum lists.
    varNameList.removeAt(i);

    // Update enum lists.
    properties->mEnum()->setEnumNames(detectionVariableIndexProperty, varNameList);
    properties->mEnum()->setValue(detectionVariableIndexProperty, tmpVarIndex);

    refreshEnumsProperties(nullptr);
}

void MTropopauseDetectionActor::onActorCreated(MActor *actor)
{
    // If the new actor is a transfer function, add it to the list of
    // available transfer functions.
    if (MTransferFunction1D *tf = dynamic_cast<MTransferFunction1D *>(actor))
    {
        // Don't render while the properties are being updated.
        enableEmissionOfActorChangedSignal(false);

        int indexFD = properties->mEnum()->value(transferFunctionFDProperty);

        QStringList availableTFs = properties->mEnum()->enumNames(transferFunctionFDProperty);
        availableTFs << tf->transferFunctionName();

        properties->mEnum()->setEnumNames(transferFunctionFDProperty, availableTFs);

        properties->mEnum()->setValue(transferFunctionFDProperty, indexFD);

        enableEmissionOfActorChangedSignal(true);
    }
}

void MTropopauseDetectionActor::onActorDeleted(MActor *actor)
{
    // If the deleted actor is a transfer function, remove it from the list of
    // available transfer functions.
    if (MTransferFunction1D *tf = dynamic_cast<MTransferFunction1D *>(actor))
    {
        enableEmissionOfActorChangedSignal(false);

        int indexFD = properties->mEnum()->value(transferFunctionFDProperty);

        QStringList availableTFs =
                properties->mEnum()->enumNames(transferFunctionFDProperty);
        availableTFs << tf->transferFunctionName();

        // If the deleted transfer function is currently connected to this
        // variable, set current transfer function to "None" (index 0).
        if (availableTFs.at(indexFD) == tf->getName()) { indexFD = 0; }

        availableTFs.removeOne(tf->getName());

        properties->mEnum()->setEnumNames(transferFunctionFDProperty, availableTFs);

        properties->mEnum()->setEnumNames(transferFunctionFDProperty, availableTFs);

        enableEmissionOfActorChangedSignal(true);
    }
}

void MTropopauseDetectionActor::onActorRenamed(MActor *actor,
                                            QString oldName)
{
    // If the renamed actor is a transfer function, change its name in the list
    // of available transfer functions.
    if (MTransferFunction1D *tf = dynamic_cast<MTransferFunction1D *>(actor))
    {
        // Don't render while the properties are being updated.
        enableEmissionOfActorChangedSignal(false);

        int indexFD = properties->mEnum()->value(transferFunctionFDProperty);

        QStringList availableTFs = properties->mEnum()->enumNames(transferFunctionFDProperty);

        // Replace affected entry.
        availableTFs[availableTFs.indexOf(oldName)] = tf->getName();

        properties->mEnum()->setEnumNames(transferFunctionFDProperty, availableTFs);

        properties->mEnum()->setValue(transferFunctionFDProperty, indexFD);

        enableEmissionOfActorChangedSignal(true);
    }
}

void MTropopauseDetectionActor::updateShadow()
{
    emitActorChangedSignal();
}

/******************************************************************************
***                             PRIVATE METHODS                             ***
*******************************************************************************/
void MTropopauseDetectionActor::initializeActorResources()
{
    // Parent initialisation.
    MNWPMultiVarActor::initializeActorResources();

    MGLResourcesManager *glRM = MGLResourcesManager::getInstance();

    //Getting the variables
    if (!variables.empty())
    {
        detectionVariable = static_cast<MNWP3DVolumeActorVariable*>(
                variables.at(detectionVariableIndex));
    }

    varNameList.clear();
    for (int vi = 0; vi < variables.size(); vi++)
    {
        MNWPActorVariable* var = variables.at(vi);
        varNameList << var->variableName;
    }
    properties->mEnum()->setEnumNames(detectionVariableIndexProperty, varNameList);
    properties->mEnum()->setValue(detectionVariableIndexProperty, detectionVariableIndex);

    bool loadShaders = false;
    loadShaders |= glRM->generateEffectProgram("tropopause_geometry",
                                               tropopauseSurfaceShader);

    if (loadShaders) reloadShaderEffects();

    // Create filter / source pipeline
    MSystemManagerAndControl *sysMC = MSystemManagerAndControl::getInstance();
    MAbstractScheduler *scheduler = sysMC->getScheduler("MultiThread");
    MAbstractMemoryManager *memoryManager = sysMC->getMemoryManager("NWP");

    tropopauseDetectionSource = std::make_shared<MTropopauseDetectionSource>();
    tropopauseDetectionSource->setScheduler(scheduler);
    tropopauseDetectionSource->setMemoryManager(memoryManager);


    transferFunctionFDTexUnit = assignTextureUnit();

    sampleSubroutines.resize(MVerticalLevelType::SIZE_LEVELTYPES);


    bool connected = connect(tropopauseDetectionSource.get(),
            SIGNAL(dataRequestCompleted(MDataRequest)),
            this, SLOT(asynchronousTropopauseAvailable(MDataRequest)));
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
}

void MTropopauseDetectionActor::onQtPropertyChanged(QtProperty *property)
{
    if (suppressUpdates) { return; }

    MNWPMultiVarActor::onQtPropertyChanged(property);

    if (property == renderTropopauseProperty)
    {
        renderTropopause = properties->mBool()->value(renderTropopauseProperty);
        emitActorChangedSignal();
    }

    if (property == applySettingsClickProperty)
    {
        if (renderTropopause) triggerAsynchronousTropopauseRequest();
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
    else if (property == tropopauseIsoProperty)
    {
        tropopauseIsoValue = properties->mSciDouble()->value(tropopauseIsoProperty);
        emitActorChangedSignal();
    }
    else if (property == tropopauseColourProperty)
    {
        tropopauseColour = properties->mColor()->value(tropopauseColourProperty);
    }
    else if (property == transferFunctionFDProperty)
    {
        setTransferFunctionFromProperty(transferFunctionFDProperty,
                                        &transferFunctionFD);

        useFDTransferFunction = (transferFunctionFD != nullptr);
    }
}


// unused, since use OIT render methods
void MTropopauseDetectionActor::renderToCurrentContext(
        MSceneViewGLWidget *sceneView)
{
    Q_UNUSED(sceneView);
}

// for translucent objects, this method is called twice
void MTropopauseDetectionActor::renderTransparencyToCurrentContext(
        MSceneViewGLWidget *sceneView)
{
    // calling render methods with OIT
    renderTropopauseSurfaceOIT(sceneView);
}

// for translucent objects, this method is called twice
void MTropopauseDetectionActor::renderTropopauseSurfaceOIT(
        MSceneViewGLWidget *sceneView)
{
    if ( !tropopauseMesh3D
        || !renderTropopause
        || numtropopauseVertices < 1
            )
    {
        return;
    }
    if (alphaShadingTropopause.size() < 1
            || alphaShadingTropopause.size() != numtropopauseVertices)
    {
        alphaShadingTropopause = QVector<QVector2D>(numtropopauseVertices, QVector2D(1.0, 1.0));
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
            tropopauseSurfaceShader->bindProgram("TriangleFilteringShadowMap");
        } else
        {
            tropopauseSurfaceShader->bindProgram((i == 0) ? "TriangleFilteringShadow"
                                               : "TriangleFilteringOIT");
        }

        tropopauseSurfaceShader->setUniformValue("mvpMatrix",
                                      *(sceneView->getModelViewProjectionMatrix()));
        tropopauseSurfaceShader->setUniformValue("pToWorldZParams",
                                      sceneView->pressureToWorldZParameters());
        tropopauseSurfaceShader->setUniformValue("lightDirection",
                                      sceneView->getLightDirection());
        tropopauseSurfaceShader->setUniformValue("cameraPosition",
                                      sceneView->getCamera()->getOrigin());

        tropopauseSurfaceShader->setUniformValue("bboxMax", maxBBox);
        tropopauseSurfaceShader->setUniformValue("bboxMin", minBBox);

        tropopauseSurfaceShader->setUniformValue("colour", tropopauseColour);

        tropopauseSurfaceShader->setUniformValue("useFDFilter", useFDTransferFunction);

        float shadowHeightWorldZ = sceneView->worldZfromPressure(shadowHeight);
        tropopauseSurfaceShader->setUniformValue("shadowColor", shadowColor);
        tropopauseSurfaceShader->setUniformValue("shadowHeight",
                                               shadowHeightWorldZ);

        tropopauseSurfaceShader->setUniformValue("lightingMode", lightingMode);

        if (transferFunctionFD != nullptr)
        {
            transferFunctionFD->getTexture()->bindToTextureUnit(transferFunctionFDTexUnit);
            tropopauseSurfaceShader->setUniformValue("tfFD", transferFunctionFDTexUnit); CHECK_GL_ERROR;

            tropopauseSurfaceShader->setUniformValue("tfFDMinMax",
                                                   QVector2D(transferFunctionFD->getMinimumValue(),
                                                             transferFunctionFD->getMaximumValue()));
        }

        tropopauseSurfaceShader->setUniformValue(
                    "inShadowMappingMode", sceneView->shadowMappingInProgress());
        tropopauseSurfaceShader->setUniformValue(
                    "lightMatrix", *sceneView->getLightOrthoProjectionMatrix());

        sceneView->getShadowMap()->bindToTextureUnit(
                    static_cast<GLuint>(sceneView->getShadowMapTexUnit()));
        tropopauseSurfaceShader->setUniformValue(
                    "shadowMap", sceneView->getShadowMapTexUnit());
        CHECK_GL_ERROR;

        sceneView->setOITUniforms(tropopauseSurfaceShader);

        const QString vboIDa = QString("tropopausesurface_vbo_shading#%1").arg(myID);
        uploadVec2ToVertexBuffer(alphaShadingTropopause, vboIDa, &vboTropopauseShading,
                                 sceneView);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,
                     iboTropopause->getIndexBufferObject()); CHECK_GL_ERROR;
        glBindBuffer(GL_ARRAY_BUFFER,
                     vboTropopause->getVertexBufferObject()); CHECK_GL_ERROR;

        int tropopauseMeshVertexSize = sizeof(Geometry::TropopauseMeshVertex);
        long bytePosition = 0;

        // vertex
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, tropopauseMeshVertexSize,
                              nullptr); CHECK_GL_ERROR;
        // normals
        bytePosition += 3 * sizeof(float);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, tropopauseMeshVertexSize,
                              (const GLvoid*) bytePosition); CHECK_GL_ERROR;

        // firstDeriv
        bytePosition += 3 * sizeof(float);
        glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, tropopauseMeshVertexSize,
                              (const GLvoid*) bytePosition); CHECK_GL_ERROR;

        glEnableVertexAttribArray(0); CHECK_GL_ERROR;
        glEnableVertexAttribArray(1); CHECK_GL_ERROR;
        glEnableVertexAttribArray(2); CHECK_GL_ERROR;

        // Create linked list for fragments
        glPolygonMode(GL_FRONT_AND_BACK, (renderAsWireFrame) ? GL_LINE : GL_FILL); CHECK_GL_ERROR;
        glDrawElements(GL_TRIANGLES, tropopauseMesh3D->getNumTriangleIndices(),
                       GL_UNSIGNED_INT, nullptr); CHECK_GL_ERROR;
    }

    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

}

void MTropopauseDetectionActor::triggerAsynchronousTropopauseRequest()
{
    if (variables.size() < 1 || getViews().empty() || isComputing
        || !detectionVariable->hasData()
            ) { return; }

    //Getting Variables
    detectionVariable = static_cast<MNWP3DVolumeActorVariable*>(variables.at(detectionVariableIndex));

    if (detectionVariable->grid == nullptr) { return; }

    if (!detectionVariable->pendingRequests.empty() )
    { return; }

    suppressUpdates = true;
    isComputing = true;

    actorPropertiesSupGroup->setEnabled(false);
    applySettingsClickProperty->setEnabled(false);
    if (tropopauseMesh3D)
    {
        MGLResourcesManager *glRM = MGLResourcesManager::getInstance();
        glRM->releaseGPUItem(vboTropopause);
        glRM->releaseGPUItem(iboTropopause);
    }
    MDataRequestHelper rh = detectionVariable->constructAsynchronousDataRequest();
    rh.insert("TROPOPAUSE_ISOVALUE", QString::number(tropopauseIsoValue));

    tropopauseDetectionSource->setDetectionVariableSource(detectionVariable->dataSource);

    tropopauseDetectionSource->requestData(rh.request());
}

void MTropopauseDetectionActor::refreshEnumsProperties(MNWPActorVariable *var)
{
    Q_UNUSED(var);
}
