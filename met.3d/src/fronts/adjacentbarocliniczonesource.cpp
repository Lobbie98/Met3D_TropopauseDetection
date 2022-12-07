/******************************************************************************
**
**  This file is part of Met.3D -- a research environment for the
**  three-dimensional visual exploration of numerical ensemble weather
**  prediction data.
**
**  Copyright 2015-2018 Marc Rautenhaus
**  Copyright 2018      Michael Kern
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
// Criteria M2
#include "adjacentbarocliniczonesource.h"
#include "gxfw/nwpactorvariableproperties.h"
#include "util/mutil.h"
#include "util/mexception.h"
#include "util/metroutines.h"

#include <log4cplus/loggingmacros.h>
#include <omp.h>

using namespace Met3D;

MAdjacentBaroclinicZoneSource::MAdjacentBaroclinicZoneSource()
    : MSingleInputProcessingWeatherPredictionDataSource(),
      partialDerivativeFilter(new MPartialDerivativeFilter()),
      vectorMagnitudeFilter(new MVectorMagnitudeFilter())
{
}

MAdjacentBaroclinicZoneSource::~MAdjacentBaroclinicZoneSource()
{
    delete partialDerivativeFilter;
    delete vectorMagnitudeFilter;
}


MStructuredGrid* MAdjacentBaroclinicZoneSource::produceData(MDataRequest request)
{
    assert(inputSource != nullptr);

    MDataRequestHelper rh(request);

    rh.removeAll(locallyRequiredKeys());

    const int DLON = MGradientProperties::DLON;
    const int DLAT = MGradientProperties::DLAT;
    QString vecMagKeyValue = QString("0/GRADIENT:%1/1/GRADIENT:%2"
                                     ).arg(DLON).arg(DLAT);

    rh.insert("GRADIENT", DLON);
    MStructuredGrid* gridU = partialDerivativeFilter->getData(rh.request());

    rh.insert("GRADIENT", DLAT);
    MStructuredGrid* gridV = partialDerivativeFilter->getData(rh.request());

    rh.insert("VECMAG_INP_REQ", vecMagKeyValue);
    MStructuredGrid* gridMag  = vectorMagnitudeFilter->getData(rh.request());

    MStructuredGrid* result = createAndInitializeResultGrid(gridU);

    // Actual geometric distance between two grids in latitudinal direction
    // ~111km.
    const float deltaLat = 111.2E3;

    const float gridLength = (gridMag->getLons()[1] - gridMag->getLons()[0]) * deltaLat;

    for (unsigned int k = 0; k < result->getNumLevels(); ++k)
    {
// #pragma omp parallel for collapse(2)
        for (unsigned int j = 0; j < result->getNumLats(); ++j)
        {
            for (unsigned int i = 0; i < result->getNumLons(); ++i)
            {

                QVector2D gradMag(gridU->getValue(k, j, i),
                                  gridV->getValue(k, j, i));
                const float magGradMag = gradMag.length();
                const float magGrad = gridMag->getValue(k, j, i);

                const float gridDistance = 1.0f / std::sqrt(2) * gridLength;
                const float abz = magGrad + magGradMag * gridDistance;

                result->setValue(k, j, i, abz);
            }
        }
    }

    vectorMagnitudeFilter->releaseData(gridMag);
    partialDerivativeFilter->releaseData(gridU);
    partialDerivativeFilter->releaseData(gridV);

    return result;
}


void MAdjacentBaroclinicZoneSource::setInputSource(MWeatherPredictionDataSource* s)
{
    inputSource = s;
    registerInputSource(inputSource);
    initializeAFPPipeline();
    //enablePassThrough(s);
}


MTask* MAdjacentBaroclinicZoneSource::createTaskGraph(MDataRequest request)
{
    assert(inputSource != nullptr);
    assert(partialDerivativeFilter != nullptr);
    assert(vectorMagnitudeFilter != nullptr);

    MTask* task =  new MTask(request, this);
    MDataRequestHelper rh(request);

    LOG4CPLUS_DEBUG(mlog, request.toStdString());

    // No keys required, only for completeness
    rh.removeAll(locallyRequiredKeys());

    // Gradient to enum
    int DLON = MGradientProperties::DLON;
    int DLAT = MGradientProperties::DLAT;
    QString vecMagKeyValue = QString("0/GRADIENT:%1/1/GRADIENT:%2"
                                     ).arg(DLON).arg(DLAT);

    rh.insert("GRADIENT", DLON);
    task->addParent(partialDerivativeFilter->getTaskGraph(rh.request()));
    rh.insert("GRADIENT", DLAT);
    task->addParent(partialDerivativeFilter->getTaskGraph(rh.request()));

    rh.insert("VECMAG_INP_REQ", vecMagKeyValue);
    task->addParent(vectorMagnitudeFilter->getTaskGraph(rh.request()));

    return task;
}


void MAdjacentBaroclinicZoneSource::initializeAFPPipeline()
{
    if (!isInizialized)
    {
        MSystemManagerAndControl *sysMC =
                MSystemManagerAndControl::getInstance();
        MAbstractScheduler *scheduler = sysMC->getScheduler("MultiThread");
        //MAbstractScheduler *scheduler = sysMC->getScheduler("SingleThread");
        MAbstractMemoryManager *memoryManager = sysMC->getMemoryManager("NWP");

        partialDerivativeFilter->setScheduler(scheduler);
        partialDerivativeFilter->setMemoryManager(memoryManager);

        vectorMagnitudeFilter->setScheduler(scheduler);
        vectorMagnitudeFilter->setMemoryManager(memoryManager);
        isInizialized = true;
    }
    partialDerivativeFilter->setInputSource(inputSource);
    vectorMagnitudeFilter->setInputSource(0, partialDerivativeFilter);
    vectorMagnitudeFilter->setInputSource(1, partialDerivativeFilter);

}


const QStringList MAdjacentBaroclinicZoneSource::locallyRequiredKeys()
{
    return (QStringList());
}
