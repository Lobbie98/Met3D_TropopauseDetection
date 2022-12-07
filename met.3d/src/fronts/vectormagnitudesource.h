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
#ifndef DEBUGVECMAGFILTER_H
#define DEBUGVECMAGFILTER_H

// standard library imports
#include <math.h>

// related third party imports
#include <QtCore>

// local application imports
#include "data/datarequest.h"
#include "gxfw/nwpmultivaractor.h"
#include "gxfw/synccontrol.h"
#include "data/derivedvars/derivedmetvarsdatasource.h"
#include "data/processingwpdatasource.h"
#include "data/structuredgridensemblefilter.h"
#include "data/structuredgrid.h"
#include "data/datarequest.h"
#include "data/partialderivativefilter.h"
#include "vectormagnitudefilter.h"

namespace Met3D
{


class MVectorMagnitudeSource
        : public MSingleInputProcessingWeatherPredictionDataSource
{
public:
    MVectorMagnitudeSource();
    ~MVectorMagnitudeSource();

   MStructuredGrid* produceData(MDataRequest request);

   MTask* createTaskGraph(MDataRequest request);

   void setInputSource(MWeatherPredictionDataSource* s);

protected:
    const QStringList locallyRequiredKeys();


private:

    void generateTaskGraph();

    void initializeVecMagPipeline();

    MPartialDerivativeFilter* partialDerivativeFilter;
    MVectorMagnitudeFilter* vectorMagnitudeFilter;
};

} // namespace Met3D

#endif // DEBUGVECMAGFILTER_H
