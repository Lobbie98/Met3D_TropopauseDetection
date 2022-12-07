/******************************************************************************
**
**  This file is part of Met.3D -- a research environment for the
**  three-dimensional visual exploration of numerical ensemble weather
**  prediction data.
**
**  Copyright 2017 Marc Rautenhaus
**  Copyright 2017 Michael Kern
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
// Criteria M1
#include "thermalfrontparametersource.h"
#include "gxfw/nwpactorvariableproperties.h"
#include "util/mutil.h"
#include "util/mexception.h"
#include "util/metroutines.h"

#include <log4cplus/loggingmacros.h>
#include <omp.h>

using namespace Met3D;

MThermalFrontParameterSource::MThermalFrontParameterSource()
        : MSingleInputProcessingWeatherPredictionDataSource(),
          partialDerivativeFilter1(new MPartialDerivativeFilter()),
          vectorMagnitudeFilter(new MVectorMagnitudeFilter()),
          partialDerivativeFilter2(new MPartialDerivativeFilter())
{
}

MThermalFrontParameterSource::~MThermalFrontParameterSource()
{
//    delete partialDerivativeFilter1;
//    delete vectorMagnitudeFilter;
//    delete partialDerivativeFilter2;
}


MStructuredGrid* MThermalFrontParameterSource::produceData(MDataRequest request)
{
    assert(inputSource != nullptr);

    //request = baseRequest;
    MDataRequestHelper rh(request);

    rh.removeAll(locallyRequiredKeys()); // no keys

    const int DLON = MGradientProperties::DLON;
    const int DLAT = MGradientProperties::DLAT;
    QString vecMagKeyValue = QString("0/GRADIENT:%1/1/GRADIENT:%2"
                                     ).arg(DLON).arg(DLAT);

    rh.insert("GRADIENT", DLON);
    MStructuredGrid* gridU = partialDerivativeFilter1->getData(rh.request());

    rh.insert("GRADIENT", DLAT);
    MStructuredGrid* gridV = partialDerivativeFilter1->getData(rh.request());

    rh.insert("VECMAG_INP_REQ", vecMagKeyValue);
    rh.insert("GRADIENT", DLON);
    MStructuredGrid* gridUU  = partialDerivativeFilter2->getData(rh.request());

    rh.insert("GRADIENT", DLAT);
    MStructuredGrid* gridVV = partialDerivativeFilter2->getData(rh.request());

    MStructuredGrid* result = createAndInitializeResultGrid(gridU);

    LOG4CPLUS_DEBUG(mlog, "Start computing tfp filter");

    for (unsigned int k = 0; k < result->getNumLevels(); ++k)
    {
// #pragma omp parallel for collapse(2)
        for (unsigned int j = 0; j < result->getNumLats(); ++j)
        {
            for (unsigned int i = 0; i < result->getNumLons(); ++i)
            {
                QVector2D gradMag(gridUU->getValue(k, j, i),
                                  gridVV->getValue(k, j, i));
                QVector2D gradient(gridU->getValue(k, j, i),
                                   gridV->getValue(k, j, i));
                gradient.normalize();

                const float tfp = -10000 * QVector2D::dotProduct(gradMag, gradient);

                result->setValue(k, j, i, tfp);
            }
        }
    }

    partialDerivativeFilter1->releaseData(gridU);
    partialDerivativeFilter1->releaseData(gridV);
    partialDerivativeFilter2->releaseData(gridUU);
    partialDerivativeFilter2->releaseData(gridVV);

    return result;
}

void MThermalFrontParameterSource::setInputSource(MWeatherPredictionDataSource* s)
{
    inputSource = s;
    registerInputSource(inputSource);
    initializeTFPPipeline();
    //enablePassThrough(s);
}


MTask* MThermalFrontParameterSource::createTaskGraph(MDataRequest request)
{
    assert(inputSource != nullptr);
    assert(partialDerivativeFilter1 != nullptr);
    assert(vectorMagnitudeFilter != nullptr);
    assert(partialDerivativeFilter2 != nullptr);

    MTask *task = new MTask(request, this);
    MDataRequestHelper rh(request);

    LOG4CPLUS_DEBUG(mlog, request.toStdString());

    // No keys required, only for completeness
    rh.removeAll(locallyRequiredKeys());

    // Gradient to enum
    int DLON = MGradientProperties::DLON;
    int DLAT = MGradientProperties::DLAT;

    rh.insert("GRADIENT", DLON);
    task->addParent(partialDerivativeFilter1->getTaskGraph(rh.request()));

    rh.insert("GRADIENT", DLAT);
    task->addParent(partialDerivativeFilter1->getTaskGraph(rh.request()));

    QString vecMagKeyValue = QString("0/GRADIENT:%1/1/GRADIENT:%2"
                                     ).arg(DLON).arg(DLAT);

    rh.insert("VECMAG_INP_REQ", vecMagKeyValue);
    rh.insert("GRADIENT", DLON);
    task->addParent(partialDerivativeFilter2->getTaskGraph(rh.request()));
    rh.insert("GRADIENT", DLAT);
    task->addParent(partialDerivativeFilter2->getTaskGraph(rh.request()));

    return task;
}


void MThermalFrontParameterSource::initializeTFPPipeline()
{
    if (!isInizialized)
    {
        MSystemManagerAndControl *sysMC =
                MSystemManagerAndControl::getInstance();
        MAbstractScheduler *scheduler = sysMC->getScheduler("MultiThread");
        //MAbstractScheduler *scheduler = sysMC->getScheduler("SingleThread");
        MAbstractMemoryManager *memoryManager = sysMC->getMemoryManager("NWP");

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


const QStringList MThermalFrontParameterSource::locallyRequiredKeys()
{
    return (QStringList());
}
