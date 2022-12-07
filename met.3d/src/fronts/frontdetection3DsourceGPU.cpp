/******************************************************************************
**
**  This file is part of Met.3D -- a research environment for the
**  three-dimensional visual exploration of numerical ensemble weather
**  prediction data.
**
**  Copyright 2015-2017 Marc Rautenhaus
**  Copyright 2017      Michael Kern
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

#include "frontdetection3DsourceGPU.h"
#include "util/mmarchingcubes.h"
#include <unistd.h>


// standard library imports

// related third party imports
#include <omp.h>

// local application imports
#include <log4cplus/loggingmacros.h>

#define COMPUTE_PARALLEL
#define MEASURE_CPU_TIME
#define qNaN (std::numeric_limits<float>::quiet_NaN())

using namespace Met3D;

MFrontDetection3DSourceGPU::MFrontDetection3DSourceGPU()
        : M3DFrontFilterSource()
{
    pbot    = 1050.; // hPa
    logpbot = log(pbot);
    zbot    = 0.;
    ztop    = 36.;
    ptop    = 20.;
    slopePtoZ = (ztop - zbot) / (log(ptop) - log(pbot));
}


void MFrontDetection3DSourceGPU::setFLESource(MFrontLocationEquationSource *s)
{
    fleSource = s;
    registerInputSource(fleSource);
    enablePassThrough(fleSource);
}


void MFrontDetection3DSourceGPU::setTFPSource(MThermalFrontParameterSource *s)
{
    tfpSource = s;
    registerInputSource(tfpSource);
    enablePassThrough(tfpSource);
}


void MFrontDetection3DSourceGPU::setABZSource(MAdjacentBaroclinicZoneSource *s)
{
    abzSource = s;
    registerInputSource(abzSource);
    enablePassThrough(abzSource);
}


void MFrontDetection3DSourceGPU::setDetectionVariableSource(MWeatherPredictionDataSource *s)
{
    detectionVariableSource = s;
    registerInputSource(detectionVariableSource);
    enablePassThrough(detectionVariableSource);
}


void MFrontDetection3DSourceGPU::setDetectionVariablePartialDeriveSource(
        MPartialDerivativeFilter *s)
{
    detectionVarPartialDerivativeSource = s;
    registerInputSource(detectionVarPartialDerivativeSource);
    enablePassThrough(detectionVarPartialDerivativeSource);
}


void MFrontDetection3DSourceGPU::setWindUSource(MWeatherPredictionDataSource *s)
{
    windUSource = s;
    registerInputSource(windUSource);
    enablePassThrough(windUSource);
}


void MFrontDetection3DSourceGPU::setWindVSource(MWeatherPredictionDataSource *s)
{
    windVSource = s;
    registerInputSource(windVSource);
    enablePassThrough(windVSource);
}


void MFrontDetection3DSourceGPU::setZSource(MWeatherPredictionDataSource *s)
{
    zSource = s;
    registerInputSource(zSource);
    enablePassThrough(zSource);
}


M3DFrontSelection* MFrontDetection3DSourceGPU::produceData(MDataRequest request)
{
    assert(fleSource != nullptr);
    assert(tfpSource != nullptr);
    assert(abzSource != nullptr);
    assert(detectionVariableSource != nullptr);
    assert(detectionVarPartialDerivativeSource != nullptr);
    assert(windUSource != nullptr);
    assert(windVSource != nullptr);
    assert(zSource != nullptr);

    MDataRequestHelper rh(request);

    const double isovalue = rh.value("FRONTS_ISOVALUE").toFloat();
    const QString windUVar = rh.value("FRONTS_WINDU_VAR");
    const QString windVVar = rh.value("FRONTS_WINDV_VAR");
    const QString zVar = rh.value("FRONTS_Z_VAR");

    rh.removeAll(locallyRequiredKeys());

    MStructuredGrid* fleGrid = fleSource->getData(rh.request());
    MStructuredGrid* tfpGrid = tfpSource->getData(rh.request());
    MStructuredGrid* abzGrid = abzSource->getData(rh.request());

    MStructuredGrid* detectionVarGrid = detectionVariableSource->getData(rh.request());
    const int DLON = MGradientProperties::DLON;
    const int DLAT = MGradientProperties::DLAT;
    rh.insert("GRADIENT", DLON);
    MStructuredGrid* dDetectionVarDXGrid =
            detectionVarPartialDerivativeSource->getData(rh.request());
    rh.insert("GRADIENT", DLAT);
    MStructuredGrid* dDetectionVarDYGrid =
            detectionVarPartialDerivativeSource->getData(rh.request());
    rh.remove("GRADIENT");

    rh.insert("VARIABLE", windUVar);
    MStructuredGrid* windUGrid = windUSource->getData(rh.request());
    rh.insert("VARIABLE", windVVar);
    MStructuredGrid* windVGrid = windVSource->getData(rh.request());
    rh.insert("VARIABLE", zVar);
    MStructuredGrid* zGrid = zSource->getData(rh.request());

    // 1) Create voxel cells
    const int nx = fleGrid->getNumLons() - 1;
    const int ny = fleGrid->getNumLats() - 1;
    const int nz = fleGrid->getNumLevels() - 1;

    //const uint32_t numVoxelCells = nx * ny * nz;

    LOG4CPLUS_DEBUG(mlog, QString("[0] Number of voxel cells: %1 x %2 x %3").arg(nx).arg(ny).arg(nz).toStdString());

    // QVector max size is 2GB (2000000000 bytes)
    // Here we check if we can compute mc or if our data is too large
    u_int32_t numEdges = (3 * fleGrid->getNumValues()
                    - fleGrid->getNumLons() * fleGrid->getNumLevels()
                    - fleGrid->getNumLons() * fleGrid->getNumLats()
                    - fleGrid->getNumLats() * fleGrid->getNumLevels());
    u_int32_t sizefloatQVector3d = numEdges * 3 * 4;
    u_int32_t maxSizeQVector = 2000000000;
    if (sizefloatQVector3d > maxSizeQVector)
    {
        LOG4CPLUS_DEBUG(mlog, "Your data set is too large to compute 3D fronts");
        return new M3DFrontSelection();
    }

    LOG4CPLUS_DEBUG(mlog, "[1] Compute voxels and isosurface triangle geometry...");

    MMarchingCubes mc(fleGrid, zGrid);
    mc.computeMeshOnCPU(isovalue);

    LOG4CPLUS_DEBUG(mlog, "[1] \t->done.");

    QVector<Geometry::MTriangle>* triangles = mc.getFlattenTriangles();
    QVector<QVector3D>* positions = mc.getFlattenInterPoints();
    QVector<QVector3D>* normals   = mc.getFlattenInterNormals();
    QVector<QVector3D>* normalsZ  = mc.getFlattenInterNormalsZ();

    int p = positions->size();
    LOG4CPLUS_DEBUG(
                mlog,
                QString("[2] Compute integration values and create mesh of %1 values...").arg(p).toStdString());
#ifdef MEASURE_CPU_TIME
    auto start2 = std::chrono::system_clock::now();
#endif

    int numNormalCurves = positions->size();
    QVector<QVector4D> posTFP(numNormalCurves, QVector4D(0, 0, 0, 0));
    for (int i = 0; i < numNormalCurves; i++)
    {
        posTFP[i].setX(positions->at(i).x());
        posTFP[i].setY(positions->at(i).y());
        posTFP[i].setZ(positions->at(i).z());
        posTFP[i].setW(tfpGrid->interpolateValue(positions->at(i)));
    }

    double westBoundary = fleGrid->getNorthWestTopDataVolumeCorner_lonlatp().x();
    double northBoundary = fleGrid->getNorthWestTopDataVolumeCorner_lonlatp().y();
    double eastBoundary = fleGrid->getSouthEastBottomDataVolumeCorner_lonlatp().x();
    double southBoundary = fleGrid->getSouthEastBottomDataVolumeCorner_lonlatp().y();
    QVector2D minBBox = QVector2D(westBoundary, southBoundary);
    QVector2D maxBBox = QVector2D(eastBoundary, northBoundary);

    QVector<QVector4D> endPos(numNormalCurves, QVector4D(0, 0, 0, 0));

    std::vector<QVector4D> endPosGPU(numNormalCurves);

    const float deltaLatKM = 111.2; //km
    const float intStepSize = detectionVarGrid->getDeltaLon();
    const int maxNumIterations = 3000 / (deltaLatKM * intStepSize);
    const float minValue = fleGrid->min();

    if (numNormalCurves == 0)
    {
        LOG4CPLUS_ERROR(mlog, "Warning: could not find any normal curve init "
                              "points");

        fleSource->releaseData(fleGrid);
        tfpSource->releaseData(tfpGrid);
        abzSource->releaseData(abzGrid);
        detectionVariableSource->releaseData(detectionVarGrid);
        detectionVarPartialDerivativeSource->releaseData(dDetectionVarDXGrid);
        detectionVarPartialDerivativeSource->releaseData(dDetectionVarDYGrid);
        windUSource->releaseData(windUGrid);
        windVSource->releaseData(windVGrid);
        zSource->releaseData(zGrid);

        auto rawFrontSurfaces = new MTriangleMeshSelection(0, 0);
        auto rawNormalCurves = new MNormalCurvesSelection(0);
        M3DFrontSelection *raw3DFronts = new M3DFrontSelection;
        raw3DFronts->setTriangleMeshSelection(rawFrontSurfaces);
        raw3DFronts->setNormalCurvesSelection(rawNormalCurves);
        return raw3DFronts;
    }

    else
    {
        std::vector<QList<QString>> normalCompSubroutines;

        normalCompSubroutines.resize(MVerticalLevelType::SIZE_LEVELTYPES);

        normalCompSubroutines[PRESSURE_LEVELS_3D]
                << "samplePressureLevel"
                << "pressureLevelGradient";

        normalCompSubroutines[HYBRID_SIGMA_PRESSURE_3D]
                << "sampleHybridLevel"
                << "hybridLevelGradient";

        normalCompSubroutines[AUXILIARY_PRESSURE_3D]
                << "sampleAuxiliaryPressure"
                << "auxiliaryPressureGradient";

        QMatrix4x4 mvpMatrix(0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0);

        QVector2D pToWorldZParams = QVector2D(logpbot, slopePtoZ);

        // upload start position
        const QString ssboNCCurvesID =
                QString("nc_integration_ssbo_%1").arg(getID());
        MGLResourcesManager *glRM = MGLResourcesManager::getInstance();
        glRM->makeCurrent();

        auto ssboNCStart =  dynamic_cast<GL::MShaderStorageBufferObject *>(
                    glRM->getGPUItem(ssboNCCurvesID));

        if(ssboNCStart)
        {
            ssboNCStart->updateSize(numNormalCurves);
            ssboNCStart->upload(posTFP.data(), GL_DYNAMIC_COPY);
        }
        else
        {
            ssboNCStart =  new GL::MShaderStorageBufferObject(
                        ssboNCCurvesID, sizeof(QVector4D), numNormalCurves);
            if (glRM->tryStoreGPUItem(ssboNCStart))
            {
                // Upload normal curve properties
                ssboNCStart->upload(posTFP.data(), GL_DYNAMIC_COPY);
            }
            else
            {
                LOG4CPLUS_WARN(mlog, "WARNING: cannot store buffer for normal curves"
                                     " in GPU memory, skipping normal curves computation.");

                delete ssboNCStart;

                fleSource->releaseData(fleGrid);
                tfpSource->releaseData(tfpGrid);
                abzSource->releaseData(abzGrid);
                detectionVariableSource->releaseData(detectionVarGrid);
                detectionVarPartialDerivativeSource->releaseData(dDetectionVarDXGrid);
                detectionVarPartialDerivativeSource->releaseData(dDetectionVarDYGrid);
                windUSource->releaseData(windUGrid);
                windVSource->releaseData(windVGrid);
                zSource->releaseData(zGrid);

                auto rawFrontSurfaces = new MTriangleMeshSelection(0, 0);
                auto rawNormalCurves = new MNormalCurvesSelection(0);
                M3DFrontSelection *raw3DFronts = new M3DFrontSelection;
                raw3DFronts->setTriangleMeshSelection(rawFrontSurfaces);
                raw3DFronts->setNormalCurvesSelection(rawNormalCurves);
                return raw3DFronts;
            }
        }


        // ssbo for result values
        const QString ssboNCEndID =
                QString("nc_end_ssbo_%1").arg(getID());

        auto ssboNCEnd =  dynamic_cast<GL::MShaderStorageBufferObject *>(
                    glRM->getGPUItem(ssboNCEndID));

        if(ssboNCEnd)
        {
            ssboNCEnd->updateSize(numNormalCurves);
            ssboNCEnd->upload(endPos.data(), GL_DYNAMIC_COPY);
        }
        else
        {
            ssboNCEnd =  new GL::MShaderStorageBufferObject(
                        ssboNCEndID, sizeof(QVector4D), numNormalCurves);
            if (glRM->tryStoreGPUItem(ssboNCEnd))
            {
                // Upload normal curve properties
                ssboNCEnd->upload(endPos.data(), GL_DYNAMIC_COPY);
            }
            else
            {
                LOG4CPLUS_WARN(mlog, "WARNING: cannot store buffer for normal curves"
                                     " in GPU memory, skipping normal curves computation.");

                delete ssboNCEnd;

                fleSource->releaseData(fleGrid);
                tfpSource->releaseData(tfpGrid);
                abzSource->releaseData(abzGrid);
                detectionVariableSource->releaseData(detectionVarGrid);
                detectionVarPartialDerivativeSource->releaseData(dDetectionVarDXGrid);
                detectionVarPartialDerivativeSource->releaseData(dDetectionVarDYGrid);
                windUSource->releaseData(windUGrid);
                windVSource->releaseData(windVGrid);
                zSource->releaseData(zGrid);

                auto rawFrontSurfaces = new MTriangleMeshSelection(0, 0);
                auto rawNormalCurves = new MNormalCurvesSelection(0);
                M3DFrontSelection *raw3DFronts = new M3DFrontSelection;
                raw3DFronts->setTriangleMeshSelection(rawFrontSurfaces);
                raw3DFronts->setNormalCurvesSelection(rawNormalCurves);
                return raw3DFronts;
            }
        }


        normalCurvesComputeShader->bindProgram("SingleIntegration");
        normalCurvesComputeShader->setUniformValue("mvpMatrix", mvpMatrix);
        normalCurvesComputeShader->setUniformValue("pToWorldZParams",
                                      pToWorldZParams);


        setVarSpecificShaderVars(normalCurvesComputeShader, fleGrid, "dataExtent",
                                 "dataVolume",
                                 "pressureTable", "surfacePressure",
                                 "hybridCoefficients", "lonLatLevAxes",
                                 "pressureTexCoordTable2D", "minMaxAccel3D",
                                 "flagsVolume", "auxPressureField3D_hPa");

        setVarSpecificShaderVars(normalCurvesComputeShader, detectionVarGrid, "dataExtentShV",
                                 "dataVolumeShV",
                                 "pressureTableShV", "surfacePressureShV",
                                 "hybridCoefficientsShV", "lonLatLevAxesShV",
                                 "pressureTexCoordTable2DShV", "minMaxAccel3DShV",
                                 "flagsVolumeShV", "auxPressureField3DShV_hPa");

        normalCurvesComputeShader->setUniformSubroutineByName(
                GL_COMPUTE_SHADER,
                normalCompSubroutines[detectionVarGrid->getLevelType()]);

        normalCurvesComputeShader->setUniformValue(
                "integrationStepSize", intStepSize); CHECK_GL_ERROR;
        normalCurvesComputeShader->setUniformValue(
                "maxNumIterations", maxNumIterations); CHECK_GL_ERROR;
        normalCurvesComputeShader->setUniformValue(
                "bisectionSteps", GLint(5)); CHECK_GL_ERROR;
        normalCurvesComputeShader->setUniformValue("integrationMode", -1); CHECK_GL_ERROR;
        normalCurvesComputeShader->setUniformValue("minValue", minValue); CHECK_GL_ERROR;
        normalCurvesComputeShader->setUniformValue("numNormalCurves", numNormalCurves); CHECK_GL_ERROR;
        normalCurvesComputeShader->setUniformValue("minBBox", minBBox); CHECK_GL_ERROR;
        normalCurvesComputeShader->setUniformValue("maxBBox", maxBBox); CHECK_GL_ERROR;


        ssboNCStart->bindToIndex(0);

        ssboNCEnd->bindToIndex(1);

        //int dispatchX = maxNumIterations / 32 + 1;

        glDispatchCompute(numNormalCurves, 1, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
        //sleep(10);//sleeps for 10 second

        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboNCEnd->getBufferObject()); CHECK_GL_ERROR;

        GLint bufMask = GL_MAP_READ_BIT;
        QVector4D* verticesGPU = (QVector4D*)
                glMapBufferRange(GL_SHADER_STORAGE_BUFFER, 0,
                numNormalCurves * sizeof(QVector4D), bufMask); CHECK_GL_ERROR;
        for (GLuint i = 0; i < GLuint(numNormalCurves); ++i)
        {
            endPosGPU[i] = verticesGPU[i];
        }

        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0); CHECK_GL_ERROR;
        normalCurvesComputeShader.reset();
    }


    MTriangleMeshSelection *rawFrontSurfaces = new MTriangleMeshSelection(
                 positions->size(),
                 triangles->size());

    MNormalCurvesSelection *rawNormalCurves = new MNormalCurvesSelection(0);

#ifdef COMPUTE_PARALLEL
#pragma omp parallel for
#endif
    // set triangles
    for (int i = 0; i < triangles->size(); i++)
    {
        rawFrontSurfaces->setTriangle(i, triangles->at(
                                          i));
    }

#ifdef COMPUTE_PARALLEL
#pragma omp parallel for
#endif
    for (auto k = 0; k < positions->size(); ++k)
    {
        // get current position
        QVector3D position = positions->at(k);

        // Check if position is valid
        if (position.x() == qNaN) { continue; }

        // initialize front mesh vertex
        Geometry::FrontMeshVertex frontVertex;

        // get all relevant values at the position
        const float tfp = tfpGrid->interpolateValue(position);
        const float abz = abzGrid->interpolateValue(position);

        // compute front type (warm or cold)
        QVector2D gradTheta(dDetectionVarDXGrid->interpolateValue(position),
                            dDetectionVarDYGrid->interpolateValue(position));
        gradTheta.normalize();
        QVector2D wind(windUGrid->interpolateValue(position),
                       windVGrid->interpolateValue(position));
        wind.normalize();
        const float type = float(QVector2D::dotProduct(-gradTheta, wind) >= 0);

        // compute frontal slope
        QVector3D normalZ = normalsZ->at(k);
        const float lengthXY = std::sqrt(normalZ.x()
                              * normalZ.x() + normalZ.y() * normalZ.y());

        const float lengthZ = normalZ.z();
        const float slope = lengthXY / lengthZ;


        // compute frontal strength:
        QVector3D ncEnd = QVector3D(endPosGPU[k].x(),
                                    endPosGPU[k].y(),
                                    endPosGPU[k].z());


        float strength = detectionVarGrid->interpolateValue(position) -
                detectionVarGrid->interpolateValue(ncEnd);

        // all needed values are computed, fill front vertex
        frontVertex.position = position;
        frontVertex.normal = normals->at(k);
        frontVertex.nCEnd = ncEnd;
        frontVertex.tfp = tfp;
        frontVertex.abz = abz;
        frontVertex.strength = strength;
        frontVertex.type = type;
        frontVertex.breadth = endPosGPU[k].w();
        frontVertex.slope = slope;
        rawFrontSurfaces->setVertex(k, frontVertex);

    }


    LOG4CPLUS_DEBUG(mlog, "[2] \t->done.");
#ifdef MEASURE_CPU_TIME
    auto end2 = std::chrono::system_clock::now();
    auto elapsed2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2);
    LOG4CPLUS_DEBUG(mlog, "[2] \t->done in " << elapsed2.count() << "ms");
#endif

    fleSource->releaseData(fleGrid);
    tfpSource->releaseData(tfpGrid);
    abzSource->releaseData(abzGrid);
    detectionVariableSource->releaseData(detectionVarGrid);
    detectionVarPartialDerivativeSource->releaseData(dDetectionVarDXGrid);
    detectionVarPartialDerivativeSource->releaseData(dDetectionVarDYGrid);
    windUSource->releaseData(windUGrid);
    windVSource->releaseData(windVGrid);
    zSource->releaseData(zGrid);

    normalCurvesComputeShader.reset();

    M3DFrontSelection *raw3DFronts = new M3DFrontSelection;
    raw3DFronts->setTriangleMeshSelection(rawFrontSurfaces);
    raw3DFronts->setNormalCurvesSelection(rawNormalCurves);
    return raw3DFronts;
}


MTask* MFrontDetection3DSourceGPU::createTaskGraph(MDataRequest request)
{
    assert(fleSource != nullptr);
    assert(tfpSource != nullptr);
    assert(abzSource != nullptr);
    assert(detectionVariableSource != nullptr);
    assert(detectionVarPartialDerivativeSource != nullptr);
    assert(windUSource != nullptr);
    assert(windVSource != nullptr);
    assert(zSource != nullptr);

    MTask* task =  new MTask(request, this);
    task->setGPUTask();
    MDataRequestHelper rh(request);

    const QString windUVar = rh.value("FRONTS_WINDU_VAR");
    const QString windVVar = rh.value("FRONTS_WINDV_VAR");
    const QString zVar = rh.value("FRONTS_Z_VAR");

    rh.removeAll(locallyRequiredKeys());

    task->addParent(fleSource->getTaskGraph(rh.request()));
    task->addParent(tfpSource->getTaskGraph(rh.request()));
    task->addParent(abzSource->getTaskGraph(rh.request()));

    task->addParent(detectionVariableSource->getTaskGraph(rh.request()));
    const int DLON = MGradientProperties::DLON;
    const int DLAT = MGradientProperties::DLAT;
    rh.insert("GRADIENT", DLON);
    task->addParent(detectionVarPartialDerivativeSource->getTaskGraph(rh.request()));
    rh.insert("GRADIENT", DLAT);
    task->addParent(detectionVarPartialDerivativeSource->getTaskGraph(rh.request()));
    rh.remove("GRADIENT");

    rh.insert("VARIABLE", windUVar);
    task->addParent(windUSource->getTaskGraph(rh.request()));
    rh.insert("VARIABLE", windVVar);
    task->addParent(windVSource->getTaskGraph(rh.request()));
    rh.insert("VARIABLE", zVar);
    task->addParent(zSource->getTaskGraph(rh.request()));

    if (!normalCurvesComputeShader)
    {
        // compile shader
        MGLResourcesManager *glRM = MGLResourcesManager::getInstance();
        glRM->makeCurrent();

        GLint maxWorkGroupInvocations; //, maxWorkGroupCount, maxWorkGroupSize;

        glGetIntegerv(GL_MAX_COMPUTE_WORK_GROUP_INVOCATIONS, &maxWorkGroupInvocations); CHECK_GL_ERROR;
        //glGetIntegerv(GL_MAX_COMPUTE_WORK_GROUP_COUNT, &maxWorkGroupCount); CHECK_GL_ERROR;
        //glGetIntegerv(GL_MAX_COMPUTE_WORK_GROUP_SIZE, &maxWorkGroupSize); CHECK_GL_ERROR;

        bool reloadShader = glRM->generateEffectProgram("normal_curves_compute_shader",
                                                        normalCurvesComputeShader);

        if (reloadShader)
        {
            if (!normalCurvesComputeShader->compileFromFile_Met3DHome("src/glsl/normal_curves_compute_shader.fx.glsl"))
            {
                throw std::exception();
            }
        }
    }
    return task;
}


const QStringList MFrontDetection3DSourceGPU::locallyRequiredKeys()
{
    return (QStringList() << "FRONTS_ISOVALUE"
                          << "FRONTS_WINDU_VAR"
                          << "FRONTS_WINDV_VAR"
                          << "FRONTS_Z_VAR");
}


void MFrontDetection3DSourceGPU::setVarSpecificShaderVars(
        std::shared_ptr<GL::MShaderEffect>& shader,
        MStructuredGrid* grid,
        const QString& structName,
        const QString& volumeName,
        const QString& pressureTableName,
        const QString& surfacePressureName,
        const QString& hybridCoeffName,
        const QString& lonLatLevAxesName,
        const QString& pressureTexCoordTable2DName,
        const QString& minMaxAccelStructure3DName,
        const QString& dataFlagsVolumeName,
        const QString& auxPressureField3DName)
{
    // Reset optional textures to avoid draw errors.
    // =============================================

    GL::MTexture *textureDummy1D = new GL::MTexture(GL_TEXTURE_1D, GL_ALPHA32F_ARB, 1);
    GL::MTexture *textureDummy2D = new GL::MTexture(GL_TEXTURE_2D, GL_ALPHA32F_ARB, 1, 1);
    GL::MTexture *textureDummy3D = new GL::MTexture(GL_TEXTURE_3D, GL_ALPHA32F_ARB, 1, 1, 1);

    GLuint textureUnitUnusedTextures;
    glGenTextures(1, &textureUnitUnusedTextures);
    // 1D textures...
    textureDummy1D->bindToTextureUnit(textureUnitUnusedTextures);
    shader->setUniformValue(pressureTableName, textureUnitUnusedTextures); CHECK_GL_ERROR;
    shader->setUniformValue(hybridCoeffName, textureUnitUnusedTextures); CHECK_GL_ERROR;

    // 2D textures...
    textureDummy2D->bindToTextureUnit(textureUnitUnusedTextures);
    shader->setUniformValue(surfacePressureName, textureUnitUnusedTextures); CHECK_GL_ERROR;
#ifdef ENABLE_HYBRID_PRESSURETEXCOORDTABLE
    shader->setUniformValue(pressureTexCoordTable2DName, textureUnitUnusedTextures); CHECK_GL_ERROR;
#endif

    // 3D textures...
    textureDummy3D->bindToTextureUnit(textureUnitUnusedTextures);
    shader->setUniformValue(dataFlagsVolumeName, textureUnitUnusedTextures); CHECK_GL_ERROR;
    shader->setUniformValue(auxPressureField3DName, textureUnitUnusedTextures); CHECK_GL_ERROR;

    // Bind textures and set uniforms.
    // ===============================

    // Bind volume data

    // bind detection var
    GL::MTexture* textureFleGrid = grid->getTexture();
    GLuint textureUnitFleGrid;
    glGenTextures(1, &textureUnitFleGrid);
    textureFleGrid->bindToTextureUnit(textureUnitFleGrid);
    shader->setUniformValue(volumeName, textureUnitFleGrid);

    shader->setUniformValue(structName + ".tfMinimum", 0.); CHECK_GL_ERROR;
    shader->setUniformValue(structName + ".tfMaximum", 0.); CHECK_GL_ERROR;

    // bind detection grid (lon, lat, lev)
    GL::MTexture* textureDataField = grid->getLonLatLevTexture();
    GLuint textureUnitDataField;
    glGenTextures(1, &textureUnitDataField);
    textureDataField->bindToTextureUnit(textureUnitDataField);
    shader->setUniformValue(lonLatLevAxesName, textureUnitDataField);

#ifdef ENABLE_RAYCASTER_ACCELERATION
    // Bind acceleration grid.
    GLuint textureUnitMinMaxAccelStructure;
    glGenTextures(1, &textureUnitMinMaxAccelStructure);
    grid->getMinMaxAccelTexture3D()->bindToTextureUnit(
                textureUnitMinMaxAccelStructure);
    shader->setUniformValue(minMaxAccelStructure3DName,
                            textureUnitMinMaxAccelStructure); CHECK_GL_ERROR;
#endif

    if (grid->getFlagsTexture() != nullptr)
    {
        // The data flags texture will only be valid if the grid contains
        // a flags field and this actor's render mode requests the flags
        // bitfield.
        GLuint textureUnitDataFlags;
        glGenTextures(1, &textureUnitDataFlags);
        grid->getFlagsTexture()->bindToTextureUnit(textureUnitDataFlags); CHECK_GL_ERROR;
        shader->setUniformValue("flagsVolume", textureUnitDataFlags);
    }

    // Set uniforms specific to data var level type.
    // =============================================

    QVector3D dataNWCrnr = grid->getNorthWestTopDataVolumeCorner_lonlatp();
    dataNWCrnr.setZ(worldZfromPressure(dataNWCrnr.z()));
    QVector3D dataSECrnr = grid->getSouthEastBottomDataVolumeCorner_lonlatp();
    dataSECrnr.setZ(worldZfromPressure(dataSECrnr.z()));

    if (grid->getLevelType() == PRESSURE_LEVELS_3D)
    {
        shader->setUniformValue(structName + ".levelType", GLint(0)); CHECK_GL_ERROR;

        MRegularLonLatStructuredPressureGrid *pgrid =
                            dynamic_cast<MRegularLonLatStructuredPressureGrid*>(grid);

        // Bind pressure to texture coordinate LUT.
        GL::MTexture* texturePressureTexCoordTable = pgrid->getPressureTexCoordTexture1D();
        GLuint textureUnitPressureTexCoordTable;
        glGenTextures(1, &textureUnitPressureTexCoordTable);
        texturePressureTexCoordTable->bindToTextureUnit(
                    textureUnitPressureTexCoordTable); CHECK_GL_ERROR;

        // Helper variables for texture coordinate LUT.
        const GLint nPTable = texturePressureTexCoordTable->getWidth();
        const GLfloat deltaZ_PTable = abs(dataSECrnr.z() - dataNWCrnr.z()) / (nPTable - 1);
        const GLfloat upperPTableBoundary = dataNWCrnr.z() + deltaZ_PTable / 2.0f;
        const GLfloat vertPTableExtent = abs(dataNWCrnr.z() - dataSECrnr.z()) + deltaZ_PTable;
        // shader->setUniformValue(structName + ".nPTable", nPTable); CHECK_GL_ERROR;
        // shader->setUniformValue(structName + ".deltaZ_PTable", deltaZ_PTable); CHECK_GL_ERROR;
        shader->setUniformValue(structName + ".upperPTableBoundary", upperPTableBoundary); CHECK_GL_ERROR;
        shader->setUniformValue(structName + ".vertPTableExtent", vertPTableExtent); CHECK_GL_ERROR;
    }

    else if (grid->getLevelType() == LOG_PRESSURE_LEVELS_3D)
    {
        shader->setUniformValue(structName + ".levelType", GLint(2)); CHECK_GL_ERROR;
    }

    else if (grid->getLevelType() == HYBRID_SIGMA_PRESSURE_3D)
    {
        shader->setUniformValue(structName + ".levelType", GLint(1)); CHECK_GL_ERROR;

        // bind detection var
        MLonLatHybridSigmaPressureGrid *hgrid =
                dynamic_cast<MLonLatHybridSigmaPressureGrid*>(grid);
                // Bind pressure to texture coordinate LUT.
        GL::MTexture* textureHybridCoefficients = hgrid->getHybridCoeffTexture();

        // Bind hybrid coefficients
        GLuint textureUnitHybridCoefficients;
        glGenTextures(1, &textureUnitHybridCoefficients);
        textureHybridCoefficients->bindToTextureUnit(textureUnitHybridCoefficients);
        shader->setUniformValue(hybridCoeffName,
                                textureUnitHybridCoefficients); CHECK_GL_ERROR;

        // Bind surface pressure
        GL::MTexture* textureSurfacePressure = hgrid->getSurfacePressureGrid()->getTexture();
        GLuint textureUnitSurfacePressure;
        glGenTextures(1, &textureUnitSurfacePressure);
        textureSurfacePressure->bindToTextureUnit(textureUnitSurfacePressure);
        shader->setUniformValue(
                surfacePressureName,
                textureUnitSurfacePressure); CHECK_GL_ERROR;

#ifdef ENABLE_HYBRID_PRESSURETEXCOORDTABLE

        // Bind pressure to texture coordinate LUT.
        GL::MTexture* texturePressureTexCoordTable = hgrid->getPressureTexCoordTexture2D();
        GLuint textureUnitPressureTexCoordTable;
        glGenTextures(1, &textureUnitPressureTexCoordTable);
        texturePressureTexCoordTable->bindToTextureUnit(
                textureUnitPressureTexCoordTable);
        shader->setUniformValue(
                pressureTexCoordTable2DName,
                textureUnitPressureTexCoordTable); CHECK_GL_ERROR;
#endif
    }

    else if (grid->getLevelType() == AUXILIARY_PRESSURE_3D)
    {
        shader->setUniformValue(structName + ".levelType", GLint(2)); CHECK_GL_ERROR;

        // Bind pressure field.
        MLonLatAuxiliaryPressureGrid *apgrid =
                dynamic_cast<MLonLatAuxiliaryPressureGrid*>(grid);
        GL::MTexture* textureAuxiliaryPressure = apgrid->getAuxiliaryPressureFieldGrid()->getTexture();
        GLuint textureUnitAuxiliaryPressure;
        glGenTextures(1, &textureUnitAuxiliaryPressure);
        textureAuxiliaryPressure->bindToTextureUnit(textureUnitAuxiliaryPressure);
        shader->setUniformValue(auxPressureField3DName,
                                textureUnitAuxiliaryPressure); CHECK_GL_ERROR;
    }

    // Precompute data extent variables and store in uniform struct.
    // =============================================================
    const GLfloat westernBoundary = dataNWCrnr.x() - grid->getDeltaLon() / 2.0f;
    const GLfloat eastWestExtent = dataSECrnr.x() - dataNWCrnr.x() + grid->getDeltaLon();
    const GLfloat northernBoundary = dataNWCrnr.y() + grid->getDeltaLat() / 2.0f;
    const GLfloat northSouthExtent = dataNWCrnr.y() - dataSECrnr.y() + grid->getDeltaLat();

    const GLint nLon = grid->getNumLons();
    const GLint nLat = grid->getNumLats();
    const GLint nLev = grid->getNumLevels();
    const GLfloat deltaLnP = std::abs(dataSECrnr.z() - dataNWCrnr.z()) / (nLev-1);
    const GLfloat upperBoundary = dataNWCrnr.z() + deltaLnP /2.0f;
    const GLfloat verticalExtent = abs(dataNWCrnr.z() - dataSECrnr.z()) + deltaLnP;

    // Assume that lat/lon spacing is the same.
    shader->setUniformValue(structName + ".deltaLat", grid->getDeltaLat()); CHECK_GL_ERROR;
    shader->setUniformValue(structName + ".deltaLon", grid->getDeltaLon()); CHECK_GL_ERROR;
    shader->setUniformValue(structName + ".dataSECrnr", dataSECrnr); CHECK_GL_ERROR;
    shader->setUniformValue(structName + ".dataNWCrnr", dataNWCrnr); CHECK_GL_ERROR;
    shader->setUniformValue(structName + ".westernBoundary", westernBoundary); CHECK_GL_ERROR;
    shader->setUniformValue(structName + ".eastWestExtent", eastWestExtent); CHECK_GL_ERROR;
    shader->setUniformValue(structName + ".northernBoundary", northernBoundary); CHECK_GL_ERROR;
    shader->setUniformValue(structName + ".northSouthExtent", northSouthExtent); CHECK_GL_ERROR;
    shader->setUniformValue(
            structName + ".gridIsCyclicInLongitude",
            grid->gridIsCyclicInLongitude()); CHECK_GL_ERROR;
    shader->setUniformValue(structName + ".nLon", nLon); CHECK_GL_ERROR;
    shader->setUniformValue(structName + ".nLat", nLat); CHECK_GL_ERROR;
    shader->setUniformValue(structName + ".nLev", nLev); CHECK_GL_ERROR;
    shader->setUniformValue(structName + ".deltaLnP", deltaLnP); CHECK_GL_ERROR;
    shader->setUniformValue(structName + ".upperBoundary", upperBoundary); CHECK_GL_ERROR;
    shader->setUniformValue(structName + ".verticalExtent", verticalExtent); CHECK_GL_ERROR;
}


double MFrontDetection3DSourceGPU::worldZfromPressure(double p_hPa)
{
    return worldZfromPressure(p_hPa, logpbot, slopePtoZ);
}


double MFrontDetection3DSourceGPU::worldZfromPressure(
        double p_hPa, double log_pBottom_hPa, double deltaZ_deltaLogP)
{
    return (log(p_hPa)-log_pBottom_hPa) * deltaZ_deltaLogP;
}


