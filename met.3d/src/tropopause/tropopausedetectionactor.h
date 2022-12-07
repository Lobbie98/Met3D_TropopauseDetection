/******************************************************************************
**
**  This file is part of Met.3D -- a research environment for the
**  three-dimensional visual exploration of numerical ensemble weather
**  prediction data.
**
**  Copyright 2015-2020 Marc Rautenhaus
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

#ifndef MET_3D_TROPOPAUSEDETECTIONACTOR_H
#define MET_3D_TROPOPAUSEDETECTIONACTOR_H


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

#include "tropopausedetectionsource.h"


namespace Met3D
{

class MTropopauseDetectionActor: public MNWPMultiVarActor,
                                 public MBoundingBoxInterface
{
Q_OBJECT

public:
    explicit MTropopauseDetectionActor();
    ~MTropopauseDetectionActor() override;

    QString getSettingsID() override
    { return "TropopauseDetectionActor"; }

    void reloadShaderEffects() override;

    void saveConfiguration(QSettings *settings) override;
    void loadConfiguration(QSettings *settings) override;

    MNWPActorVariable* createActorVariable(
            const MSelectableDataSource& dataSource) override;

    const QList<MVerticalLevelType> supportedLevelTypes() override;

    void onBoundingBoxChanged() override;

public slots:
    void asynchronousTropopauseAvailable(MDataRequest request);

    void onAddActorVariable(MNWPActorVariable *var) override;

    void onDeleteActorVariable(MNWPActorVariable *var) override;

    void onActorCreated(MActor *actor);

    void onActorDeleted(MActor *actor);

    void onActorRenamed(MActor *actor, QString oldName);

    void updateShadow();

private:
    void initializeActorResources() override;

    void onQtPropertyChanged(QtProperty *property) override;

    void renderToCurrentContext(MSceneViewGLWidget *sceneView) override;

    void renderTropopauseSurfaceOIT(MSceneViewGLWidget *sceneView);//TODO: sehr viel Funktion umbauen?

    void triggerAsynchronousTropopauseRequest();

    bool setTransferFunction(QtProperty* tfProp, const QString &tfName);//TODO:brauch ich ne transferfunction?

    void setTransferFunctionFromProperty(
            QtProperty* tfProp, MTransferFunction1D** transferFunc);//TODO: same, wie darüber

    typedef enum {
        ABSOLUTE_CHANGE = 0,
        PER_100KM = 1,
        VALUE_AT_FRONT = 2,
    } filterType;//TODO:brauche wahrscheinlich etwas ähnliches, aber nicht das gleiche

    /**
    * Refresh all variables select boxes if a variable was added or removed.
     */
    void refreshEnumsProperties(MNWPActorVariable *var);

    static QString filterTypeToString(filterType type);

    static filterType stringToFilterType(QString type);

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
    //          |-- detection var TODO gucken, welche Variablen ich brauche
    //          |-- eastward wind (u)
    //          |-- northward wind (v)
    //          |-- geopotential height (only 3D fronts)
    //      |-> display options TODO: kann komplett weg?
    //          |-- 3D fronts
    //          |-- 3D normal curves
    //          |-- 2D fronts
    //          |-- 2D normal curves
    //      |-> front filter options TODO: Welche Filter müssen hier rein?
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
    //      |-> front appearance
    //          |-> shading
    //              |-- mode
    //              |-- variable
    //              |-- type
    //              |-- transfer function
    //          |-- bounding box
    //          |-- show warm/cold fronts TODO gucken, ob ich sowas auch brauche
    //          |-- show cold side of fronts
    //          |-> light and shadow
    //              |-- shadow elevation
    //              |-- shadow color
    //              |-- light mode

    // ************************ MAIN PROPERTIES ********************************
    //      |-- update
    QtProperty* applySettingsClickProperty;

    // ************************ INPUT VAR PROPERTIES ***************************
    //      |-> input variables
    QtProperty* inputVarGroupProperty;

    //          |-- detection var
    QtProperty*                 detectionVariableIndexProperty;
    uint8_t                     detectionVariableIndex;
    MNWP3DVolumeActorVariable*  detectionVariable;

    // ************************ DISPLAY OPTIONS ********************************
    //      |-> display options
    QtProperty* displayOptionsGroupProperty;

    //          |-- 3D tropopause
    QtProperty* renderTropopauseProperty;
    bool        renderTropopause;

    //          |-- TropopauseIsoValue
    QtProperty* tropopauseIsoProperty;
    double tropopauseIsoValue;

    //          |-- tropopauseColour
    QtProperty* tropopauseColourProperty;
    QColor tropopauseColour;



    // ************************ FILTER OPTIONS *********************************
    //      |-> front filter options
    QtProperty* filterGroupProperty;

    //          |-- firstDeriv transfer function
    QtProperty*                 transferFunctionFDProperty;
    MTransferFunction1D*        transferFunctionFD;
    int                         transferFunctionFDTexUnit;

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
    /*QtProperty*                 transferFunctionShadingProperty;
    MTransferFunction1D*        transferFunctionShading;
    int                         transferFunctionShadingTexUnit;*/

    //          |-- bounding box
    // Initialized in cpp-file

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

// ********************** DEPRECREATED PROPERTIES ******************************
    bool                                            useFDTransferFunction;

    std::shared_ptr<MTropopauseDetectionSource>       tropopauseDetectionSource;

    MTropopauseTriangleMeshSelection*                  tropopauseMesh3D;//TODO: für 3D fläche ausreichend?

    // shaders to render fronts
    std::shared_ptr<GL::MShaderEffect>  tropopauseSurfaceShader;//TODO: nur ein shadet? sollte reichen

    GL::MVertexBuffer*                  vboTropopause;//TODO: die 3 reichen?
    GL::MIndexBuffer*                   iboTropopause;
    GL::MVertexBuffer*                  vboTropopauseShading;

    std::vector<QList<QString>>               sampleSubroutines;//TODO: brauch ich sowas?

    QVector<QVector2D> alphaShadingTropopause;

    int numtropopauseVertices;

    u_int32_t restartIndex;
};



class MTropopauseDetectionActorFactory : public MAbstractActorFactory
{
public:
    MTropopauseDetectionActorFactory() : MAbstractActorFactory() {}

protected:
    MActor* createInstance() override { return new MTropopauseDetectionActor(); }
};

}

#endif // MET_3D_TROPOPAUSEDETECTIONACTOR_H
