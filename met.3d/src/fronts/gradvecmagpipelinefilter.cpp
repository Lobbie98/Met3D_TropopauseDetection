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

// related third party imports
#include <log4cplus/loggingmacros.h>

// local application imports
#include "gradvecmagpipelinefilter.h"
#include "gxfw/nwpactorvariableproperties.h"
#include "util/mutil.h"
#include "util/mexception.h"
#include "util/metroutines.h"



//using namespace std;

namespace Met3D
{

/******************************************************************************
***                     CONSTRUCTOR / DESTRUCTOR                            ***
*******************************************************************************/

MGradVecMagPipelineFilter::MGradVecMagPipelineFilter()
    : MSingleInputProcessingWeatherPredictionDataSource(),
      partialDerivativeFilter(new MPartialDerivativeFilter()),
      vectorMagnitudeFilter(new MVectorMagnitudeFilter())
{
}

MGradVecMagPipelineFilter::~MGradVecMagPipelineFilter()
{
    delete partialDerivativeFilter;
    delete vectorMagnitudeFilter;
}
/******************************************************************************
***                            PUBLIC METHODS                               ***
*******************************************************************************/

MStructuredGrid *MGradVecMagPipelineFilter::produceData(MDataRequest request)
{
    assert(inputSource != nullptr);

    //request = baseRequest;
    MDataRequestHelper rh(request);
    MStructuredGrid *inputGrid;

    // Gradient to enum
    int dlon;
    int dlat;
    dlon = MGradientProperties::DLON;
    dlat = MGradientProperties::DLAT;
    //dlon = MGradientProperties::DLON_SOBEL;
    //dlat = MGradientProperties::DLAT_SOBEL;



    QString vecMagKeyValue = QString("0/GRADIENT:%1/1/GRADIENT:%2"
                                     ).arg(dlon).arg(dlat);

    rh.insert("VECMAG_INP_REQ", vecMagKeyValue);
    rh.removeAll(locallyRequiredKeys());

    inputGrid = vectorMagnitudeFilter->getData(rh.request());

    return inputGrid;
}


void MGradVecMagPipelineFilter::setInputSource(MWeatherPredictionDataSource* s)
{
    inputSource = s;
    registerInputSource(inputSource);
    initializeVecMagPipeline();
    enablePassThrough(s);
}

MTask* MGradVecMagPipelineFilter::createTaskGraph(MDataRequest request)
{
    assert(inputSource != nullptr);
    assert(partialDerivativeFilter != nullptr);
    assert(vectorMagnitudeFilter != nullptr);

    MTask *task = new MTask(request, this);
    MDataRequestHelper rh(request);

    // Gradient to enum
    int dlon;
    int dlat;
    dlon = MGradientProperties::DLON;
    dlat = MGradientProperties::DLAT;
    //dlon = MGradientProperties::DLON_SOBEL;
    //dlat = MGradientProperties::DLAT_SOBEL;

    QString vecMagKeyValue = QString("0/GRADIENT:%1/1/GRADIENT:%2"
                                     ).arg(dlon).arg(dlat);

    rh.insert("VECMAG_INP_REQ", vecMagKeyValue);
    rh.removeAll(locallyRequiredKeys());

    task->addParent(vectorMagnitudeFilter->getTaskGraph(rh.request()));

    return task;
}


void MGradVecMagPipelineFilter::initializeVecMagPipeline()
{

    MSystemManagerAndControl *sysMC =
            MSystemManagerAndControl::getInstance();
    MAbstractScheduler *scheduler = sysMC->getScheduler("MultiThread");
    //MAbstractScheduler *scheduler = sysMC->getScheduler("SingleThread");
    MAbstractMemoryManager *memoryManager = sysMC->getMemoryManager("NWP");

    partialDerivativeFilter->setScheduler(scheduler);
    partialDerivativeFilter->setMemoryManager(memoryManager);
    partialDerivativeFilter->setInputSource(inputSource);

    vectorMagnitudeFilter->setScheduler(scheduler);
    vectorMagnitudeFilter->setMemoryManager(memoryManager);
    vectorMagnitudeFilter->setInputSource(0, partialDerivativeFilter);
    vectorMagnitudeFilter->setInputSource(1, partialDerivativeFilter);

}


const QStringList MGradVecMagPipelineFilter::locallyRequiredKeys()
{
    return QStringList();
}

}  // namespace Met3D
