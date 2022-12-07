/******************************************************************************
**
**  This file is part of Met.3D -- a research environment for the
**  three-dimensional visual exploration of numerical ensemble weather
**  prediction data.
**
**  Copyright 2015-2020 Marc Rautenhaus
**  Copyright 2020      Andreas Beckert
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
// standard library imports
#include "assert.h"
#include <chrono>

// related third party imports
#include <log4cplus/loggingmacros.h>

// local application imports
#include "frontlocationequationsource.h"
#include "gxfw/nwpactorvariableproperties.h"
#include "actors/nwphorizontalsectionactor.h"
#include "util/mutil.h"
#include "util/mexception.h"
#include "util/metroutines.h"

#define MEASURE_CPU_TIME

using namespace std;

namespace Met3D
{

/******************************************************************************
***                     CONSTRUCTOR / DESTRUCTOR                            ***
*******************************************************************************/

MFrontLocationEquationSource::MFrontLocationEquationSource()
    : MSingleInputProcessingWeatherPredictionDataSource(),
      partialDerivativeFilter1(new MPartialDerivativeFilter()),
      vectorMagnitudeFilter(new MVectorMagnitudeFilter()),
      partialDerivativeFilter2(new MPartialDerivativeFilter())
{
}

MFrontLocationEquationSource::~MFrontLocationEquationSource()
{
    delete partialDerivativeFilter1;
    delete vectorMagnitudeFilter;
    delete partialDerivativeFilter2;
}
/******************************************************************************
***                            PUBLIC METHODS                               ***
*******************************************************************************/

MStructuredGrid *MFrontLocationEquationSource::produceData(MDataRequest request)
{
    assert(inputSource != nullptr);

#ifdef MEASURE_CPU_TIME
    auto start = std::chrono::system_clock::now();
#endif
    LOG4CPLUS_DEBUG(mlog, "Start to solve the front location equation");

    //request = baseRequest;
    MDataRequestHelper rh(request);
    QVector<MStructuredGrid*> grids(2);

    float distance_km = rh.value("FPMA_DISTANCE").toFloat();
    //float distance_km = 15.;

    rh.removeAll(locallyRequiredKeys());
    int DLON = MGradientProperties::DLON;
    int DLAT = MGradientProperties::DLAT;


    QString vecMagKeyValue = QString("0/GRADIENT:%1/1/GRADIENT:%2"
                                     ).arg(DLON).arg(DLAT);

    rh.insert("VECMAG_INP_REQ", vecMagKeyValue);
    rh.insert("GRADIENT", DLON);
    grids[0] = partialDerivativeFilter2->getData(rh.request());

    rh.insert("GRADIENT", DLAT);
    grids[1] = partialDerivativeFilter2->getData(rh.request());

    MStructuredGrid* resultGrid = createAndInitializeResultGrid(grids[0]);

    // solve the FLE
    computeFivePointMeanAxis(grids, resultGrid, distance_km);

    partialDerivativeFilter2->releaseData(grids[0]);
    partialDerivativeFilter2->releaseData(grids[1]);

#ifdef MEASURE_CPU_TIME
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    LOG4CPLUS_DEBUG(mlog, "Front location equation: solved in "
                    << elapsed.count() << "ms");
#endif

    return resultGrid;
}


void MFrontLocationEquationSource::setInputSource(MWeatherPredictionDataSource* s)
{
    inputSource = s;
    registerInputSource(inputSource);
    initializeFLEPipeline();
    //enablePassThrough(s);
}

MTask* MFrontLocationEquationSource::createTaskGraph(MDataRequest request)
{
    LOG4CPLUS_DEBUG(mlog, "MFrontLocationEquationSource::createTaskGraph");
    assert(inputSource != nullptr);
    assert(partialDerivativeFilter1 != nullptr);
    assert(vectorMagnitudeFilter != nullptr);
    assert(partialDerivativeFilter2 != nullptr);

    MTask *task = new MTask(request, this);
    MDataRequestHelper rh(request);

    LOG4CPLUS_DEBUG(mlog, request.toStdString());

    // remove distance key for five point mean axis
    rh.removeAll(locallyRequiredKeys());
    int DLON = MGradientProperties::DLON;
    int DLAT = MGradientProperties::DLAT;

    QString vecMagKeyValue = QString("0/GRADIENT:%1/1/GRADIENT:%2"
                                     ).arg(DLON).arg(DLAT);

    rh.insert("VECMAG_INP_REQ", vecMagKeyValue);
    rh.insert("GRADIENT", DLON);
    task->addParent(partialDerivativeFilter2->getTaskGraph(rh.request()));

    rh.insert("GRADIENT", DLAT);
    task->addParent(partialDerivativeFilter2->getTaskGraph(rh.request()));

    return task;
}


void MFrontLocationEquationSource::initializeFLEPipeline()
{
    if (!isInizialized)
    {
        partialDerivativeFilter1->setScheduler(scheduler);
        partialDerivativeFilter1->setMemoryManager(memoryManager);

        vectorMagnitudeFilter->setScheduler(scheduler);
        vectorMagnitudeFilter->setMemoryManager(memoryManager);

        partialDerivativeFilter2->setScheduler(scheduler);
        partialDerivativeFilter2->setMemoryManager(memoryManager);

        isInizialized = true;
    }

    partialDerivativeFilter1->setInputSource(inputSource);
    vectorMagnitudeFilter->setInputSource(0, partialDerivativeFilter1);
    vectorMagnitudeFilter->setInputSource(1, partialDerivativeFilter1);
    partialDerivativeFilter2->setInputSource(vectorMagnitudeFilter);
}


const QStringList MFrontLocationEquationSource::locallyRequiredKeys()
{
    return (QStringList()) << "FPMA_DISTANCE";
//    return (QStringList());
}


void MFrontLocationEquationSource::computeFivePointMeanAxis(
        QVector<MStructuredGrid *> grids, MStructuredGrid *resultGrid,
        float distance_km)
{
    LOG4CPLUS_DEBUG(mlog, "Start computing 5-point-mean-axis");

    const double deltaY = grids[0]->getDeltaLat_km();// * 1000;
    const float deltaLat = (distance_km * 180.) /
            (M_PI * MetConstants::EARTH_RADIUS_km);
    const float lonMin = std::min(grids[0]->getLon(0) ,
                                  grids[0]->getLon(grids[0]->getNumLons() - 1));
    const float lonMax = std::max(grids[0]->getLon(0) ,
                                  grids[0]->getLon(grids[0]->getNumLons() - 1));

    const float latMin = std::min(grids[0]->getLat(0),
                                  grids[0]->getLat(grids[0]->getNumLats() - 1));
    const float latMax = std::max(grids[0]->getLat(0),
                                  grids[0]->getLat(grids[0]->getNumLats() - 1));

    const float gridDeltaLat = grids[0]->getDeltaLat();
    const float gridDeltaLon = grids[0]->getDeltaLon();

#pragma omp parallel for
    for (unsigned int k = 0; k < grids[0]->getNumLevels(); ++k)
    {
        for (unsigned int j = 0; j < grids[0]->getNumLats(); ++j)
        {
            float lat = grids[0]->getLat(j);
            float deltaLon = deltaLat / fabs(cos(degreesToRadians(lat)));

            const float latP = std::max(lat - deltaLat, latMin);
            const float latN = std::min(lat + deltaLat, latMax);

            double deltaX = grids[0]->getDeltaLon_km(j);// * 1000;

            for (unsigned int i = 0; i < grids[0]->getNumLons(); ++i)
            {
                QVector2D s(grids[0]->getValue(k, j, i),
                        grids[1]->getValue(k, j, i));

                float lon = grids[0]->getLon(i);
                const float lonP = std::max(lon - deltaLon, lonMin);
                const float lonN = std::min(lon + deltaLon, lonMax);

                float pressure = grids[0]->getPressure(k, j, i);

                // 1) Obtain gradients of surrounding grid points.
                QVector2D A(grids[0]->interpolateValue(lonP, lat, pressure),
                        grids[1]->interpolateValue(lonP, lat, pressure));
                QVector2D B(grids[0]->interpolateValue(lonN, lat, pressure),
                        grids[1]->interpolateValue(lonN, lat, pressure));
                QVector2D C(grids[0]->interpolateValue(lon, latP, pressure),
                        grids[1]->interpolateValue(lon, latP, pressure));
                QVector2D D(grids[0]->interpolateValue(lon, latN, pressure),
                        grids[1]->interpolateValue(lon, latN, pressure));

                QVector<QVector2D> vecs = { s, A, B, C, D };

                // 2) Compute 5-point-mean axis
                float P = 0;
                float Q = 0;

                for (auto &vec : vecs)
                {
                    float len = vec.length();

                    if (len <= 1E-20) { len = 1.0; } // why 1 not 0?

                    float arctan = std::atan2(vec.y(), vec.x()); // measure angle between vector and x-axis

                    float beta = 2 * arctan; // multiply angle by 2

                    P += len * std::cos(beta);   // x-component of doubled vector
                    Q += len * std::sin(beta);   // y-component of doubled vector
                }

                float arctan = std::atan2(Q, P); // [-PI, PI]

                float betaMean = arctan * 0.5; // [-PI / 2, PI / 2]

                QVector2D unitMeanS(std::cos(betaMean)/5., std::sin(betaMean)/5.);

                // 3) Resolve vectors into unit axis s
                float As = QVector2D::dotProduct(A, unitMeanS)
                        * unitMeanS.x();
                float Bs = QVector2D::dotProduct(B, unitMeanS)
                        * unitMeanS.x();
                float Cs = QVector2D::dotProduct(C, unitMeanS)
                        * unitMeanS.y();
                float Ds = QVector2D::dotProduct(D, unitMeanS)
                        * unitMeanS.y();

                // 4) Compute central differences in x/y direction of resolved vectors
                float ddx = (Bs - As) / ((lonN - lonP)/(gridDeltaLon * deltaX));
                float ddy = (Ds - Cs) / ((latN - latP)/(gridDeltaLat * deltaY));

                // 5) Compute total divergence
                float divergence = (ddx + ddy) * -1;
                resultGrid->setValue(k, j, i, divergence);
            }
        }
    }
}

}  // namespace Met3D
