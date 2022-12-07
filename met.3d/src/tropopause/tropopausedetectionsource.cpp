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

#include "tropopausedetectionsource.h"
#include "util/mmarchingcubes.h"

// standard library imports

// related third party imports
#include <omp.h>

// local application imports
#include <log4cplus/loggingmacros.h>

#define COMPUTE_PARALLEL
#define MEASURE_CPU_TIME
#define qNaN (std::numeric_limits<float>::quiet_NaN())


using namespace Met3D;
MDummyFilter::MDummyFilter()
        : MSingleInputProcessingWeatherPredictionDataSource()
{

}
/******************************************************************************
***                            PUBLIC METHODS                               ***
*******************************************************************************/

MStructuredGrid *MDummyFilter::produceData(Met3D::MDataRequest request)
{
    assert(inputSource != nullptr);

    MDataRequestHelper rh(request);
    rh.removeAll(locallyRequiredKeys());

    const int DP = MGradientProperties::DP;
    rh.insert("GRADIENT", DP);

    return inputSource->getData(rh.request());
}

void MDummyFilter::setInputSource(MWeatherPredictionDataSource* s)
{
    inputSource = s;
    registerInputSource(inputSource);
    //enablePassThrough(s);
}

MTask *MDummyFilter::createTaskGraph(MDataRequest request)
{
    assert(inputSource != nullptr);
    MTask* task = new MTask(request, this);

    MDataRequestHelper rh(request);

    const int DP = MGradientProperties::DP;
    rh.insert("GRADIENT", DP);
    task->addParent(inputSource->getTaskGraph(rh.request()));
    return task;
}



/******************************************************************************
***                          PROTECTED METHODS                              ***
*******************************************************************************/

const QStringList MDummyFilter::locallyRequiredKeys()
{
    return (QStringList());
}


////////////////////////////////////////////////////
/// \brief MTropopauseDetectionSource::MTropopauseDetectionSource
///////////////////////////////////////////////////////////////

MTropopauseDetectionSource::MTropopauseDetectionSource()
        : detectionVarPartialDerivative1Source(new MPartialDerivativeFilter()),
          detectionVarPartialDerivative2Source(new MPartialDerivativeFilter()),
          dummy1(new MDummyFilter())
{
}



void MTropopauseDetectionSource::setDetectionVariableSource(
        MWeatherPredictionDataSource *s)
{
    detectionVariableSource = s;
    registerInputSource(detectionVariableSource);
    initializeTDPipeline();
}

MTropopauseTriangleMeshSelection* MTropopauseDetectionSource::produceData(MDataRequest request)
{
    assert(detectionVarPartialDerivative1Source != nullptr);
    assert(detectionVarPartialDerivative2Source != nullptr);

    MDataRequestHelper rh(request);

    const double isovalue = rh.value("TROPOPAUSE_ISOVALUE").toFloat();

    rh.removeAll(locallyRequiredKeys());

    const int DP = MGradientProperties::DP;
    rh.insert("GRADIENT", DP);

    MStructuredGrid* firstDerivativeGrid = detectionVarPartialDerivative1Source->getData(rh.request());

    MStructuredGrid* secondDerivativeGrid = detectionVarPartialDerivative2Source->getData(rh.request());


    rh.remove("GRADIENT");

    // 1) Create voxel cells
    const int nx = secondDerivativeGrid->getNumLons() - 1;
    const int ny = secondDerivativeGrid->getNumLats() - 1;
    const int nz = secondDerivativeGrid->getNumLevels() - 1;

    // QVector max size is 2GB (2000000000 bytes)
    // Here we check if we can compute mc or if our data is too large
    u_int32_t numEdges = (3 * secondDerivativeGrid->getNumValues()
                    - secondDerivativeGrid->getNumLons() * secondDerivativeGrid->getNumLevels()
                    - secondDerivativeGrid->getNumLons() * secondDerivativeGrid->getNumLats()
                    - secondDerivativeGrid->getNumLats() * secondDerivativeGrid->getNumLevels());
    u_int32_t sizefloatQVector3d = numEdges * 3 * 4;
    u_int32_t maxSizeQVector = 2000000000;
    if (sizefloatQVector3d > maxSizeQVector)
    {
        LOG4CPLUS_DEBUG(mlog, "Your data set is too large to compute Tropopause");
        return new MTropopauseTriangleMeshSelection(0,0);
    }

    LOG4CPLUS_DEBUG(mlog, "[1] Compute voxels and isosurface triangle geometry...");
    //MarchingCubes
    MMarchingCubes mc(secondDerivativeGrid);
    mc.computeMeshOnCPU(isovalue);

    LOG4CPLUS_DEBUG(mlog, "[1] \t->done.");

    QVector<Geometry::MTriangle>* triangles = mc.getFlattenTriangles();
    QVector<QVector3D>* positions = mc.getFlattenInterPoints();
    QVector<QVector3D>* normals   = mc.getFlattenInterNormals();

    int p = positions->size();
    LOG4CPLUS_DEBUG(
                mlog,
                QString("[2] Compute integration values and create mesh of %1 values...").arg(p).toStdString());
#ifdef MEASURE_CPU_TIME
    auto start2 = std::chrono::system_clock::now();
#endif

    MTropopauseTriangleMeshSelection *rawTropopause = new MTropopauseTriangleMeshSelection(
                 positions->size(),
                 triangles->size());

    //const float minValue = thirdDerivativeGrid->min();

#ifdef COMPUTE_PARALLEL
#pragma omp parallel for
#endif
    for (auto k = 0; k < positions->size(); ++k)
    {
        // get current position
        QVector3D position = positions->at(k);

        // Check if position is valid
        if (position.x() == qNaN) { continue; }

        // initialize tropopause mesh vertex
        Geometry::TropopauseMeshVertex tropopauseVertex;

        //compute first Derivative
        const float firstDeriv = firstDerivativeGrid->interpolateValue(position);

        // all needed values are computed, fill tropopause vertex
        tropopauseVertex.position = position;
        tropopauseVertex.normal = normals->at(k);
        tropopauseVertex.firstDeriv = firstDeriv;


        // set vertices to raw frontal surfaces
        rawTropopause->setVertex(k, tropopauseVertex);
    }

#ifdef COMPUTE_PARALLEL
#pragma omp parallel for
#endif
    // set triangles
    for (int i = 0; i < triangles->size(); i++)
    {
        rawTropopause->setTriangle(i, triangles->at(
                                          i));
    }

    LOG4CPLUS_DEBUG(mlog, "[2] \t->done.");
#ifdef MEASURE_CPU_TIME
    auto end2 = std::chrono::system_clock::now();
    auto elapsed2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2);
    LOG4CPLUS_DEBUG(mlog, "[2] \t->done in " << elapsed2.count() << "ms");
#endif

    detectionVarPartialDerivative1Source->releaseData(firstDerivativeGrid);
    detectionVarPartialDerivative2Source->releaseData(secondDerivativeGrid);

    return rawTropopause;

}

MTropopauseTriangleMeshSelection* MTropopauseDetectionSource::getData(MDataRequest request)
{
    return dynamic_cast<MTropopauseTriangleMeshSelection*>(MScheduledDataSource::getData(request));
}


MTask* MTropopauseDetectionSource::createTaskGraph(MDataRequest request)
{
    assert(detectionVariableSource != nullptr);
    assert(detectionVarPartialDerivative1Source != nullptr);
    assert(detectionVarPartialDerivative2Source != nullptr);
    assert(dummy1 != nullptr);

    MTask* task =  new MTask(request, this);
    MDataRequestHelper rh(request);

    rh.removeAll(locallyRequiredKeys());

    const int DP = MGradientProperties::DP;
    rh.insert("GRADIENT", DP);
    task->addParent(detectionVarPartialDerivative1Source->getTaskGraph(rh.request()));
    task->addParent(dummy1->getTaskGraph(rh.request()));
    task->addParent(detectionVarPartialDerivative2Source->getTaskGraph(rh.request()));
    rh.remove("GRADIENT");
    return task;
}

void MTropopauseDetectionSource::initializeTDPipeline()
{
    if (!isInizialized)
    {
        MSystemManagerAndControl *sysMC =
                MSystemManagerAndControl::getInstance();
        MAbstractScheduler *scheduler = sysMC->getScheduler("MultiThread");
        MAbstractMemoryManager *memoryManager = sysMC->getMemoryManager("NWP");

        detectionVarPartialDerivative1Source->setScheduler(scheduler);
        detectionVarPartialDerivative1Source->setMemoryManager(memoryManager);

        detectionVarPartialDerivative2Source->setScheduler(scheduler);
        detectionVarPartialDerivative2Source->setMemoryManager(memoryManager);

        dummy1->setScheduler(scheduler);
        dummy1->setMemoryManager(memoryManager);

        isInizialized = true;
    }
    detectionVarPartialDerivative1Source->setInputSource(detectionVariableSource);
    dummy1->setInputSource(detectionVarPartialDerivative1Source);
    detectionVarPartialDerivative2Source->setInputSource(dummy1);
}


const QStringList MTropopauseDetectionSource::locallyRequiredKeys()
{
    return (QStringList() << "TROPOPAUSE_ISOVALUE");
}
