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
#ifndef MET_3D_TROPOPAUSEDETECTIONSOURCE_H
#define MET_3D_TROPOPAUSEDETECTIONSOURCE_H

// standard library imports

// related third party imports

// local application imports
#include "data/scheduleddatasource.h"
#include "data/structuredgrid.h"
#include "data/weatherpredictiondatasource.h"
#include "data/smoothfilter.h"
#include "data/partialderivativefilter.h"

#include "tropopausesurfacemesh.h"

#include "gxfw/gl/shadereffect.h"

namespace Met3D
{

class MDummyFilter
        : public MSingleInputProcessingWeatherPredictionDataSource
{
public:
    MDummyFilter();

    MStructuredGrid* produceData(MDataRequest request);

    void setInputSource(MWeatherPredictionDataSource* s);

    MTask* createTaskGraph(MDataRequest request);
protected:
    const QStringList locallyRequiredKeys();

};


class MTropopauseDetectionSource : public MScheduledDataSource
{
public:
    explicit MTropopauseDetectionSource();

    void setDetectionVariableSource(MWeatherPredictionDataSource* s);


    MTropopauseTriangleMeshSelection* produceData(MDataRequest request) override;

    MTropopauseTriangleMeshSelection* getData(MDataRequest request) override;


    MTask *createTaskGraph(MDataRequest request) override;

private:
    const QStringList locallyRequiredKeys() override;
    void initializeTDPipeline();

    MWeatherPredictionDataSource* detectionVariableSource;
    MPartialDerivativeFilter* detectionVarPartialDerivative1Source;
    MPartialDerivativeFilter* detectionVarPartialDerivative2Source;

    MDummyFilter* dummy1;
    MDummyFilter* dummy2;

    bool isInizialized = false;
};


}

#endif // MET_3D_TROPOPAUSEDETECTIONSOURCE_H
