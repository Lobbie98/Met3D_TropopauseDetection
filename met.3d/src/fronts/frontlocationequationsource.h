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
#ifndef FRONTLOCATIONEQUATIONFILTER_H
#define FRONTLOCATIONEQUATIONFILTER_H

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

/**
  @brief MFrontLocationEquationSource implements the computation the
front location equation (FLE) of T. Hewson. This class builds its own pipeline
to include all filters necessary to implement the FLE. The FLE itself is
computed together with the five point mean axis in the method
computeFivePointMeanAxis.
 */
class MFrontLocationEquationSource
        : public MSingleInputProcessingWeatherPredictionDataSource
{
public:
    MFrontLocationEquationSource();
    ~MFrontLocationEquationSource();

   MStructuredGrid * produceData(MDataRequest request);

   MTask* createTaskGraph(MDataRequest request);

   void setInputSource(MWeatherPredictionDataSource* s);

protected:
    const QStringList locallyRequiredKeys();


private:

    void generateTaskGraph();

    void computeFivePointMeanAxis(QVector<MStructuredGrid*> grids,
                                  MStructuredGrid *resultGrid,
                                  float distance_km);
    void computeDirectionalDerivative();

    MPartialDerivativeFilter    *partialDerivativeFilter1;
    MVectorMagnitudeFilter      *vectorMagnitudeFilter;
    MPartialDerivativeFilter    *partialDerivativeFilter2;

    void initializeFLEPipeline();
    bool isInizialized = false;
};

} // namespace Met3D

#endif // FRONTLOCATIONEQUATIONFILTER_H
