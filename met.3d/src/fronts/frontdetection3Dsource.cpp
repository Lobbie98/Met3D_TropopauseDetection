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

#include "frontdetection3Dsource.h"
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

MFrontDetection3DSource::MFrontDetection3DSource()
        : M3DFrontFilterSource()
{
}


void MFrontDetection3DSource::setFLESource(MFrontLocationEquationSource *s)
{
    fleSource = s;
    registerInputSource(fleSource);
    enablePassThrough(fleSource);
}


void MFrontDetection3DSource::setTFPSource(MThermalFrontParameterSource *s)
{
    tfpSource = s;
    registerInputSource(tfpSource);
    enablePassThrough(tfpSource);
}


void MFrontDetection3DSource::setABZSource(MAdjacentBaroclinicZoneSource *s)
{
    abzSource = s;
    registerInputSource(abzSource);
    enablePassThrough(abzSource);
}


void MFrontDetection3DSource::setDetectionVariableSource(MWeatherPredictionDataSource *s)
{
    detectionVariableSource = s;
    registerInputSource(detectionVariableSource);
    enablePassThrough(detectionVariableSource);
}


void MFrontDetection3DSource::setDetectionVariablePartialDeriveSource(
        MPartialDerivativeFilter *s)
{
    detectionVarPartialDerivativeSource = s;
    registerInputSource(detectionVarPartialDerivativeSource);
    enablePassThrough(detectionVarPartialDerivativeSource);
}


void MFrontDetection3DSource::setWindUSource(MWeatherPredictionDataSource *s)
{
    windUSource = s;
    registerInputSource(windUSource);
    enablePassThrough(windUSource);
}


void MFrontDetection3DSource::setWindVSource(MWeatherPredictionDataSource *s)
{
    windVSource = s;
    registerInputSource(windVSource);
    enablePassThrough(windVSource);
}


void MFrontDetection3DSource::setZSource(MWeatherPredictionDataSource *s)
{
    zSource = s;
    registerInputSource(zSource);
    enablePassThrough(zSource);
}


M3DFrontSelection* MFrontDetection3DSource::produceData(MDataRequest request)
{
    LOG4CPLUS_DEBUG(mlog, "MFrontDetection3DSource::produceData");
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
    const float overtracing = rh.value("NC_OVER_TRACING").toFloat();

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

    MTriangleMeshSelection *rawFrontSurfaces = new MTriangleMeshSelection(
                 positions->size(),
                 triangles->size());

    MNormalCurvesSelection *rawNormalCurves = new MNormalCurvesSelection(
                positions->size());

    const float minValue = fleGrid->min();

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

        // compute normal curves
        Geometry::NormalCurve nc = integrateAlongThermalGradient(
                    fleGrid, dDetectionVarDXGrid, dDetectionVarDYGrid,
                    position, tfp, minValue, overtracing);

        // compute frontal strength:
        QVector3D ncEnd = nc.positions.last();
        float strength = detectionVarGrid->interpolateValue(nc.positions.first()) -
                detectionVarGrid->interpolateValue(ncEnd);

        // all needed values are computed, fill front vertex
        frontVertex.position = position;
        frontVertex.normal = normals->at(k);
        frontVertex.nCEnd = nc.positions.last();
        frontVertex.tfp = tfp;
        frontVertex.abz = abz;
        frontVertex.strength = strength;
        frontVertex.type = type;
        frontVertex.breadth = nc.breadth;
        frontVertex.slope = slope;

        // all needed values are computed, fill normal curve vertices
        nc.tfp = tfp;
        nc.abz = abz;
        nc.strength = strength;
        nc.type = type;

        // set vertices to raw frontal surfaces and raw normal curves
        rawFrontSurfaces->setVertex(k, frontVertex);
        rawNormalCurves->setNormalCurve(k, nc);
    }

#ifdef COMPUTE_PARALLEL
#pragma omp parallel for
#endif
    // set triangles
    for (int i = 0; i < triangles->size(); i++)
    {
        rawFrontSurfaces->setTriangle(i, triangles->at(
                                          i));
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

    M3DFrontSelection *raw3DFronts = new M3DFrontSelection;
    raw3DFronts->setTriangleMeshSelection(rawFrontSurfaces);
    raw3DFronts->setNormalCurvesSelection(rawNormalCurves);
    return raw3DFronts;
}


MTask* MFrontDetection3DSource::createTaskGraph(MDataRequest request)
{
    LOG4CPLUS_DEBUG(mlog, "MFrontDetection3DSource::createTaskGraph");
    assert(fleSource != nullptr);
    assert(tfpSource != nullptr);
    assert(abzSource != nullptr);
    assert(detectionVariableSource != nullptr);
    assert(detectionVarPartialDerivativeSource != nullptr);
    assert(windUSource != nullptr);
    assert(windVSource != nullptr);
    assert(zSource != nullptr);

    MTask* task =  new MTask(request, this);
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

    return task;
}


const QStringList MFrontDetection3DSource::locallyRequiredKeys()
{
    return (QStringList() << "FRONTS_ISOVALUE"
                          << "FRONTS_WINDU_VAR"
                          << "FRONTS_WINDV_VAR"
                          << "FRONTS_Z_VAR"
                          << "FRONTS_NC_OVER_TRACING"
                          << "FRONTS_NC_INTEGRATION_GPU");
}


Geometry::NormalCurve MFrontDetection3DSource::integrateAlongThermalGradient(
        MStructuredGrid* fleGrid,
        MStructuredGrid* ddxGrid,
        MStructuredGrid* ddyGrid,
        QVector3D startPosition,
        const float tfp,
        const float minValue,
        const float overtracing)
{
    Q_UNUSED(overtracing);

    Geometry::NormalCurve normalCurve;

    const float deltaLatKM = 111.2; //km
    const float deltaLatM = 1.12E5; //m
    const float intStepSize = fleGrid->getDeltaLon();
    const int maxNumIterations = 3000 / (deltaLatKM * intStepSize);

    float normalCurveLength_km = 0; // KM

    QVector3D currentPos = startPosition;
    QVector3D gradient(0, 0, 0);
    QVector3D prevPos;

    int direction = (tfp < 0) ? -1 : 1;

    double westBoundary = ddxGrid->getNorthWestTopDataVolumeCorner_lonlatp().x();
    double northBoundary = ddxGrid->getNorthWestTopDataVolumeCorner_lonlatp().y();
    double eastBoundary = ddxGrid->getSouthEastBottomDataVolumeCorner_lonlatp().x();
    double southBoundary = ddxGrid->getSouthEastBottomDataVolumeCorner_lonlatp().y();

    normalCurve.positions.append(startPosition);

    for (auto cc = 0; cc < maxNumIterations; ++cc)
    {
        gradient.setX(ddxGrid->interpolateValue(currentPos) * deltaLatM);
        gradient.setY(ddyGrid->interpolateValue(currentPos) * deltaLatM);
        gradient.normalize();

        prevPos = currentPos;
        currentPos -= gradient * intStepSize * direction;

        // if next position out of boundary or in an area where
        // fle is not defined return normal curve vertices.
        if(currentPos.x() < westBoundary || currentPos.x() > eastBoundary
         ||currentPos.y() < southBoundary || currentPos.y() > northBoundary
         || (fleGrid->interpolateValue(currentPos) < minValue))
        {
            normalCurve.breadth = normalCurveLength_km;
            return normalCurve;
        }

        // check if next position is out of frontal zone, that is, if the
        // fle field is larger or equal to zero
        float locator = fleGrid->interpolateValue(currentPos);

        if (locator <= 0.)
        {
            // compute a bisection correction
            currentPos = bisectionCorrection(currentPos,
                                             prevPos,
                                             fleGrid,
                                             locator);

            normalCurve.positions.push_back(currentPos);

            // calculate length of normal curve segment and add to total length
            QVector3D diffPos = currentPos - prevPos;
            float deltaLonKM = deltaLatKM * std::cos(currentPos.y() / 180.0 * M_PI);
            diffPos.setX(diffPos.x() * deltaLonKM);
            diffPos.setY(diffPos.y() * deltaLatKM);
            diffPos.setZ(0);

            normalCurveLength_km += diffPos.length();
            normalCurve.breadth = normalCurveLength_km;

            return normalCurve;

        }

        // add next vertex to normal curve
        normalCurve.positions.push_back(currentPos);

        // calculate length of normal curve segment and add to total length
        QVector3D diffPos = currentPos - prevPos;
        float deltaLonKM = deltaLatKM * std::cos(currentPos.y() / 180.0 * M_PI);
        diffPos.setX(diffPos.x() * deltaLonKM);
        diffPos.setY(diffPos.y() * deltaLatKM);
        diffPos.setZ(0);
        normalCurveLength_km += diffPos.length();
    }

    // return normal curve if max number of iterations is reached
    normalCurve.breadth = normalCurveLength_km;
    return normalCurve;
}


// Correct the position of any detected iso-surface by using the
// bisection algorithm.
QVector3D MFrontDetection3DSource::bisectionCorrection(QVector3D position,
                              QVector3D prevPosition,
                              MStructuredGrid* fleGrid,
                              float locator)
{
    QVector3D centerPosition;
    const int numBisectionSteps = 5;
    float locatorCenter = locator;

    for (int i = 0; i < numBisectionSteps; ++i)
    {
        centerPosition = (position + prevPosition) / 2.0;
        locatorCenter = fleGrid->interpolateValue(centerPosition);

        if (locatorCenter <= 0.)
        {
            position = centerPosition;
        } else {
            prevPosition = centerPosition;
        }
    }
    return centerPosition;
}


