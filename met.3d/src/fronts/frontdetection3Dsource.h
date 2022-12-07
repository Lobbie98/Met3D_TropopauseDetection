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
#ifndef MET_3D_FRONTDETECTION3DSOURCE_H
#define MET_3D_FRONTDETECTION3DSOURCE_H

// standard library imports

// related third party imports

// local application imports
#include "data/scheduleddatasource.h"
#include "data/structuredgrid.h"
#include "data/weatherpredictiondatasource.h"
#include "data/smoothfilter.h"
#include "data/partialderivativefilter.h"

#include "frontsurfacemesh.h"
#include "frontlocationequationsource.h"
#include "thermalfrontparametersource.h"
#include "adjacentbarocliniczonesource.h"

#include "gxfw/gl/shadereffect.h"


namespace Met3D
{
class MFrontDetection3DSource : public M3DFrontFilterSource
{
public:
    explicit MFrontDetection3DSource();

    void setFLESource(MFrontLocationEquationSource *s);
    void setTFPSource(MThermalFrontParameterSource *s);
    void setABZSource(MAdjacentBaroclinicZoneSource *s);
    void setDetectionVariableSource(MWeatherPredictionDataSource* s);
    void setDetectionVariablePartialDeriveSource(MPartialDerivativeFilter *s);
    void setWindUSource(MWeatherPredictionDataSource* s);
    void setWindVSource(MWeatherPredictionDataSource* s);
    void setZSource(MWeatherPredictionDataSource* s);

    M3DFrontSelection* produceData(MDataRequest request) override;

    MTask *createTaskGraph(MDataRequest request) override;

    static Geometry::NormalCurve integrateAlongThermalGradient(
            MStructuredGrid* fleGrid,
            MStructuredGrid* ddxGrid,
            MStructuredGrid* ddyGrid,
            QVector3D startPosition,
            const float tfp,
            const float minValue,
            const float overtracing);

    static QVector3D bisectionCorrection(QVector3D position,
                                  QVector3D prevPosition,
                                  MStructuredGrid* fleGrid,
                                  float locator);

private:
    const QStringList locallyRequiredKeys() override;

    MFrontLocationEquationSource* fleSource;
    MThermalFrontParameterSource* tfpSource;
    MAdjacentBaroclinicZoneSource* abzSource;
    MWeatherPredictionDataSource* detectionVariableSource;
    MPartialDerivativeFilter* detectionVarPartialDerivativeSource;
    MWeatherPredictionDataSource* windUSource;
    MWeatherPredictionDataSource* windVSource;
    MWeatherPredictionDataSource* zSource;
};


}

#endif //MET_3D_FRONTDETECTION3DSOURCE_H
