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
#ifndef MET_3D_FRONTDETECTION3DSOURCEGPU_H
#define MET_3D_FRONTDETECTION3DSOURCEGPU_H

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
class MFrontDetection3DSourceGPU : public M3DFrontFilterSource
{
public:
    explicit MFrontDetection3DSourceGPU();

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

private:
    const QStringList locallyRequiredKeys() override;

    void integrateAlongThermalGradient(Geometry::FrontMeshVertex* v,
                                       MStructuredGrid* fleGrid,
                                       MStructuredGrid* detectionVarGrid,
                                       MStructuredGrid* ddxGrid,
                                       MStructuredGrid* ddyGrid,
                                       const float overtracing,
                                       const float minValue);

    void setVarSpecificShaderVars(
            std::shared_ptr<GL::MShaderEffect>& shader,
            MStructuredGrid* grid,
            const QString& structName,
            const QString& volumeName,
            const QString& pressureTableName,
            const QString& surfacePressureName,
            const QString& hybridCoeffName,
            const QString& lonLatLevAxesName,
            const QString& pressureTexCoordTable2DName,
            const QString& minMaxAccelStructure3DName,
            const QString& dataFlagsVolumeName,
            const QString& auxPressureField3DName);

    double worldZfromPressure(double p_hPa);


    double worldZfromPressure(
            double p_hPa, double log_pBottom_hPa, double deltaZ_deltaLogP);

    MFrontLocationEquationSource* fleSource;
    MThermalFrontParameterSource* tfpSource;
    MAdjacentBaroclinicZoneSource* abzSource;
    MWeatherPredictionDataSource* detectionVariableSource;
    MPartialDerivativeFilter* detectionVarPartialDerivativeSource;
    MWeatherPredictionDataSource* windUSource;
    MWeatherPredictionDataSource* windVSource;
    MWeatherPredictionDataSource* zSource;

    std::shared_ptr<GL::MShaderEffect> normalCurvesComputeShader;

    double pbot; // hPa
    double logpbot;
    double zbot;
    double ztop;
    double ptop;
    double slopePtoZ;

};


}

#endif //MET_3D_FRONTDETECTION3DSOURCE_H
