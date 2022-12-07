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

#include "vectormagnitudefilter.h"

// standard library imports

// related third party imports
#include <log4cplus/loggingmacros.h>
#include <omp.h>

// local application imports

using namespace Met3D;

MVectorMagnitudeFilter::MVectorMagnitudeFilter()
    : MProcessingWeatherPredictionDataSource()
{
}


MVectorMagnitudeFilter::~MVectorMagnitudeFilter()
{
    //for (int i = 0; i < inputSource.size(); i++)
    //{
    //    delete inputSource[i];
    //}
}

/******************************************************************************
***                            PUBLIC METHODS                               ***
*******************************************************************************/

void MVectorMagnitudeFilter::setInputSource(int id, MWeatherPredictionDataSource* s)
{
    if (inputSource.size() <= id)
    {
        inputSource.insert(id, s);
    }
    else
    {
        inputSource[id] = s;
    }
    registerInputSource(inputSource[id]);
    //enablePassThrough(s);
}


MStructuredGrid* MVectorMagnitudeFilter::produceData(MDataRequest request)
{
    const int nDim = inputSource.size();
    for (int i = 0; i < nDim; i++)
    {
        assert(inputSource[i] != nullptr);
    }

    MDataRequestHelper rh(request);
    QVector<MStructuredGrid*> inputGrids;

    for (int i = 0; i < nDim; i++)
    {
        inputGrids.append(inputSource[i]->getData(
                          constructInputSourceRequestFromRequest(i, request)));
    } // inputGrids now contains nDim entries

    rh.removeAll(locallyRequiredKeys());

    for (int i = 0; i < nDim; i++)
    {
        if (inputGrids[i]->getDataType() == SINGLE)
        {
            inputGrids[i]->copyFloatDataToDouble();
        }
    }

    MStructuredGrid* result = createAndInitializeResultGrid(inputGrids[0]);
    result->initializeDoubleData();

    //result->setGeneratingRequest(request);

    LOG4CPLUS_DEBUG(mlog, "Vector magnitude filter: computing vector magnitude ");

#pragma omp parallel for
    for (unsigned  int k = 0; k < result->getNumLevels(); k++)
    {
        for (unsigned  int j = 0; j < result->getNumLats(); j++)
        {
            for (unsigned  int i = 0; i < result->getNumLons(); i++)
            {
                double magnitude = 0;

                // Compute the Euclidean distance of all vector values
                for (int cc = 0; cc < nDim; cc++)
                {
                    const double value = inputGrids[cc]->getValue_double(k, j, i);
                    magnitude += value * value;
                }

                magnitude = std::sqrt(magnitude);

                result->setValue_double(k, j, i, magnitude);
            }
        }
    }

    // TODO (ab, 23.07.2020) test if its possible to release the input grids here
    for (int i = 0; i < nDim; ++i)
    {
        inputSource[i]->releaseData(inputGrids[i]);
    }
    result->copyDoubleDataToFloat();
    return result;
}


MTask* MVectorMagnitudeFilter::createTaskGraph(MDataRequest request)
{
    const int nDim = inputSource.size();
    MTask *task = new MTask(request, this);
    MDataRequestHelper rh(request);

    for (int i = 0; i < nDim; i++)
    {
        assert(inputSource[i] != nullptr);
        task->addParent(inputSource[i]->getTaskGraph(
                    constructInputSourceRequestFromRequest(i, rh)));
    }

    rh.removeAll(locallyRequiredKeys());

    return task;
}


QList<MVerticalLevelType> MVectorMagnitudeFilter::availableLevelTypes()
{
//TODO (mr, 26Feb2018) -- see availableValidTimes()
    assert(inputSource[0] != nullptr);
    return inputSource[0]->availableLevelTypes();
}


QStringList MVectorMagnitudeFilter::availableVariables(
        MVerticalLevelType levelType)
{
//TODO (mr, 26Feb2018) -- see availableValidTimes()
    assert(inputSource[0] != nullptr);
    return inputSource[0]->availableVariables(levelType);
}


QSet<unsigned int> MVectorMagnitudeFilter::availableEnsembleMembers(
        MVerticalLevelType levelType, const QString& variableName)
{
//TODO (mr, 26Feb2018) -- see availableValidTimes()
    assert(inputSource[0] != nullptr);
    return inputSource[0]->availableEnsembleMembers(levelType, variableName);
}


QList<QDateTime> MVectorMagnitudeFilter::availableInitTimes(
        MVerticalLevelType levelType, const QString& variableName)
{
//TODO (mr, 26Feb2018) -- see availableValidTimes()
    assert(inputSource[0] != nullptr);
    return inputSource[0]->availableInitTimes(levelType, variableName);
}


QList<QDateTime> MVectorMagnitudeFilter::availableValidTimes(
        MVerticalLevelType levelType,
        const QString& variableName, const QDateTime& initTime)
{
//TODO (mr, 26Feb2018) -- needs to use values from both input sources, depending
// on further usage (i.e. mapping from requested to input times etc.).
    assert(inputSource[0] != nullptr);
    return inputSource[0]->availableValidTimes(levelType, variableName, initTime);
}


QString MVectorMagnitudeFilter::variableLongName(
        MVerticalLevelType levelType,
        const QString& variableName)
{
    assert(inputSource[0] != nullptr);
    return inputSource[0]->variableLongName(levelType, variableName);
}


QString MVectorMagnitudeFilter::variableStandardName(
        MVerticalLevelType levelType,
        const QString& variableName)
{
    assert(inputSource[0] != nullptr);
    return inputSource[0]->variableStandardName(levelType, variableName);
}


QString MVectorMagnitudeFilter::variableUnits(
        MVerticalLevelType levelType,
        const QString& variableName)
{
    assert(inputSource[0] != nullptr);
    return inputSource[0]->variableUnits(levelType, variableName);
}



/******************************************************************************
***                          PROTECTED METHODS                              ***
*******************************************************************************/


MDataRequest MVectorMagnitudeFilter::constructInputSourceRequestFromRequest(
        int id, MDataRequestHelper rh)
{
    // request from "downstream" pipeline
    //MDataRequestHelper rhInp = rh; // for input source "id"
    // in base request e.g.: VEGMAG:0="GRADIENT/DLON"
    QString locallyReqKey = locallyRequiredKeys()[0];
    const QStringList keyReq = rh.value(locallyReqKey).split("/");
    QString localReqString = keyReq[id * 2 + 1];
    const QStringList localReq = localReqString.split(":");
    rh.removeAll(locallyRequiredKeys());
    rh.insert(localReq[0], localReq[1]);
    //rhInp.insert(localReq[0], localReq[1]);
    //rh.removeAll(locallyRequiredKeys());
    //rhInp.removeAll(locallyRequiredKeys());
    return rh.request();
}

const QStringList MVectorMagnitudeFilter::locallyRequiredKeys()
{
    return (QStringList() << "VECMAG_INP_REQ");
}
