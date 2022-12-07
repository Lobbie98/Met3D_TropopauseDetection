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

#ifndef THERMALFRONTPARAMETERSOURCE_H
#define THERMALFRONTPARAMETERSOURCE_H

#include "data/scheduleddatasource.h"
#include "data/structuredgrid.h"
#include "data/processingwpdatasource.h"
#include "data/partialderivativefilter.h"
#include "fronts/vectormagnitudefilter.h"

namespace Met3D
{
class MThermalFrontParameterSource : public
        MSingleInputProcessingWeatherPredictionDataSource
{
public:
    MThermalFrontParameterSource();
    ~MThermalFrontParameterSource();

    MStructuredGrid* produceData(MDataRequest request);

    MTask *createTaskGraph(MDataRequest request);

    void setInputSource(MWeatherPredictionDataSource* s);

private:
    const QStringList locallyRequiredKeys();
    void generateTaskGraph();

    void initializeTFPPipeline();

    MPartialDerivativeFilter* partialDerivativeFilter1;
    MVectorMagnitudeFilter* vectorMagnitudeFilter;
    MPartialDerivativeFilter* partialDerivativeFilter2;

    bool isInizialized = false;
};
}


#endif //THERMALFRONTPARAMETERSOURCE_H
