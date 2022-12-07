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
#ifndef MET_2D_FRONTDETECTION2DSOURCE_H
#define MET_2D_FRONTDETECTION2DSOURCE_H

// standard library imports

// related third party imports

// local application imports
#include "data/scheduleddatasource.h"
#include "data/structuredgrid.h"
#include "data/weatherpredictiondatasource.h"
#include "data/smoothfilter.h"
#include "data/partialderivativefilter.h"

#include "frontline.h"
#include "frontlocationequationsource.h"
#include "thermalfrontparametersource.h"
#include "adjacentbarocliniczonesource.h"


namespace Met3D
{
class MFrontDetection2DSource : public M2DFrontFilterSource
{
public:
    explicit MFrontDetection2DSource();

    void setFLESource(MFrontLocationEquationSource *s);
    void setTFPSource(MThermalFrontParameterSource *s);
    void setABZSource(MAdjacentBaroclinicZoneSource *s);
    void setDetectionVariableSource(MWeatherPredictionDataSource* s);
    void setDetectionVariablePartialDeriveSource(MPartialDerivativeFilter *s);
    void setWindUSource(MWeatherPredictionDataSource* s);
    void setWindVSource(MWeatherPredictionDataSource* s);

    M2DFrontSelection* produceData(MDataRequest request) override;

    MTask *createTaskGraph(MDataRequest request) override;

private:
    const QStringList locallyRequiredKeys() override;

    QVector<QVector3D>* getLinesFromEnsembleMember(
    MStructuredGrid* grid, const float isovalue, const float pHa);

    QVector3D getLineVertex(const QVector3D& p0, const QVector3D& p1,
                            const float v0, const float v1,
                            const float isovalue);

    QVector<u_int32_t> sortFrontLineVertices(
            MLineSelection* frontLineVertices,
            int numVertices);

    QVector<QVector<u_int32_t>> tryToJoinLineStrip(
            MLineSelection* frontLineVertices,
            QVector<QVector<u_int32_t>> indexArray, int i, int j);

    void addLine(QVector<Geometry::FrontLineVertex>* frontLineVertices,
                 QVector<QVector<u_int32_t>>* indexArray,
                  u_int32_t indexP1, u_int32_t indexP2);


    MFrontLocationEquationSource* fleSource;
    MThermalFrontParameterSource* tfpSource;
    MAdjacentBaroclinicZoneSource* abzSource;
    MWeatherPredictionDataSource* detectionVariableSource;
    MPartialDerivativeFilter* detectionVarPartialDerivativeSource;
    MWeatherPredictionDataSource* windUSource;
    MWeatherPredictionDataSource* windVSource;
};


}

#endif //MET_2D_FRONTDETECTION2DSOURCE_H
