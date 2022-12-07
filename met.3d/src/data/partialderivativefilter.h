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

#ifndef GRADIENTFILTER_H
#define GRADIENTFILTER_H

// standard library imports
#include <math.h>

// related third party imports
#include <QtCore>

// local application imports
#include "processingwpdatasource.h"
#include "structuredgridensemblefilter.h"
#include "structuredgrid.h"
#include "datarequest.h"

namespace Met3D
{

/**
  @brief MPartialDerivativeFilter implements gradient computation for gridded data.
 */
class MPartialDerivativeFilter
        : public MSingleInputProcessingWeatherPredictionDataSource
{
public:
    MPartialDerivativeFilter();

    MStructuredGrid* produceData(MDataRequest request);

    MTask* createTaskGraph(MDataRequest request);

    /**
     Set data source computing geopotential height. This method has to be
     called if the input field is derived wrt. z (height).
     @param s
     */
    void setGeoPotSource(MWeatherPredictionDataSource* s);

protected:
    const QStringList locallyRequiredKeys();

    MWeatherPredictionDataSource* geoPotSource;

private:
    /**
     * @brief computePartialDerivativeLongitude Computes the horizontal partial
     * derivative in longitudinal direction using central and at boundaries
     * forward or backward finite differences.
     * @param inputGrid pointer to the input grid
     * @param resultGrid pointer to the result grid
     */
    void computePartialDerivativeLongitude(
            MStructuredGrid *inputGrid, MStructuredGrid *resultGrid);

    /**
     * @brief computePartialDerivativeLatitude Computes the horizontal partial
     * derivative in latitudinal direction using central and at boundaries
     * forward or backward finite differences.
     * @param inputGrid pointer to the input grid
     * @param resultGrid pointer to the result grid
     */
    void computePartialDerivativeLatitude(
            MStructuredGrid *inputGrid, MStructuredGrid *resultGrid);

    /**
     * @brief computePartialDerivativeLongitude Computes the horizontal partial
     * derivative in longitudinal direction using central and at boundaries
     * forward or backward finite differences.
     * @param inputGrid pointer to the input grid
     * @param resultGrid pointer to the result grid
     */
    void computePartialDerivativeLongitudeSobel(
            MStructuredGrid *inputGrid, MStructuredGrid *resultGrid);

    /**
     * @brief computePartialDerivativeLatitude Computes the horizontal partial
     * derivative in latitudinal direction using central and at boundaries
     * forward or backward finite differences.
     * @param inputGrid pointer to the input grid
     * @param resultGrid pointer to the result grid
     */
    void computePartialDerivativeLatitudeSobel(
            MStructuredGrid *inputGrid, MStructuredGrid *resultGrid);

    /**
     * @brief computePartialDerivativeVertical Computes the vertical partial
     * derivative using central and at boundaries
     * forward or backward finite differences.
     * @param inputGrid pointer to the input grid
     * @param resultGrid pointer to the result grid
     */
    void computePartialDerivativeVertical(
            MStructuredGrid *inputGrid, MStructuredGrid *resultGrid);

    /**
     * @brief computePartialDerivativeVerticalGeometricHeight Computes the
     * vertical partial derivative on geometric height using central and at
     * boundaries forward or backward finite differences.
     * @param inputGrid pointer to the input grid
     * @param geoPot pointer to the geopotential grid
     * @param resultGrid pointer to the result grid
    */
    void computePartialDerivativeVerticalGeometricHeight(
            MStructuredGrid *inputGrid, MStructuredGrid *geoPot,
            MStructuredGrid *resultGrid);

    /**
     * @brief computePartialDerivativePressureLongitude Computes the horizontal
     * partial derivative of pressure in longitudinal direction using central
     * and at boundaries forward or backward finite differences.
     * @param inputGrid pointer to the input grid
     * @param resultGrid pointer to the result grid
     */
    void computePartialDerivativePressureLongitude(
            MStructuredGrid *inputGrid, MStructuredGrid *resultGrid);

    /**
     * @brief computePartialDerivativePressureLongitude Computes the horizontal
     * partial derivative of pressure in latitudinal direction using central
     * and at boundaries forward or backward finite differences.
     * @param inputGrid pointer to the input grid
     * @param resultGrid pointer to the result grid
     */
    void computePartialDerivativePressureLatitude(
            MStructuredGrid *inputGrid, MStructuredGrid *resultGrid);


    void computePressureCoordinateTransormationLongitude(
            MStructuredGrid *inputGrid, MStructuredGrid *resultGrid);

    void computePressureCoordinateTransormationLatitude(
            MStructuredGrid *inputGrid, MStructuredGrid *resultGrid);

    /**
     * @brief periodicBoundaryTreatment test if periodic boundaries are present.
     * @param inputGrid pointer to the input grid
     * @return QList of 3 boolean:
     *          index 0: longitudinal periodic boundary conditions,
     *          index 1: periodic boundary conditions at north pole,
     *          index 2: periodic boundary conditions at south pole
     */
    QList<bool> periodicBoundaryTreatment(MStructuredGrid *inputGrid);


    /**
     * @brief callLibcalvarGradientRoutine function interface to call the
     * ddh3 routine of the libcalvar LAGRANTO routine.
     * @param inputGrid pointer to the input grid
     * @param resultGrid pointer to the result grid
     * @param direction horizontal direction of the gradient ("lon" or "lat")
     */
    void callLibcalvarGradientRoutine(
            MStructuredGrid *inputGrid, MStructuredGrid *resultGrid,
            QString direction);


    /**
     * @brief computeSecondHorizontalDerivative Compute the second
     * horizontal derivative in direction lon or lat
     * @param inputGrid pointer to the input grid
     * @param resultGrid pointer to the result grid
     * @param direction horizontal direction of the gradient ("lon" or "lat")
     */
    void computeSecondHorizontalDerivative(
            MStructuredGrid *inputGrid, MStructuredGrid *resultGrid,
            QString direction);


    /**
     * @brief computeSecondVerticalDerivative Compute the second vertical
     * derivative
     * @param inputGrid pointer to the input grid
     * @param resultGrid pointer to the result grid
     */
    void computeSecondVerticalDerivative(
            MStructuredGrid *inputGrid, MStructuredGrid *resultGrid);


    /**
     * @brief computeSecondVerticalDerivativeGeometricHeight Compute the vertical
     * gradient using geometric height instead of pressure
     * @param inputGrid
     * @param geoPot
     * @param resultGrid
     */
    void computeSecondVerticalDerivativeGeometricHeight(
            MStructuredGrid *inputGrid, MStructuredGrid *geoPot,
            MStructuredGrid *resultGrid);


    /**
     * @brief computeHorizontalPressureLevelTransformationForSecondDerivativeLon
     * This function
     * computes a vertical coordination transformation of horizontal gradient.
     * It transforms the gradient from native vertical coordinate system
     * to pressure coordinate system.
     * @param inputGrid pointer to the input grid
     * @param resultGrid pointer to the result grid
     */
    void computeHorizontalPressureLevelTransformationForSecondDerivativeLon(
            MStructuredGrid *inputGrid, MStructuredGrid *resultGrid);


    /**
     * @brief computeHorizontalPressureLevelTransformationForSecondDerivativeLat
     * This function
     * computes a vertical coordination transformation of horizontal gradient.
     * It transforms the gradient from native vertical coordinate system
     * to pressure coordinate system.
     * @param inputGrid pointer to the input grid
     * @param resultGrid pointer to the result grid
     */
    void computeHorizontalPressureLevelTransformationForSecondDerivativeLat(
            MStructuredGrid *inputGrid, MStructuredGrid *resultGrid);


    inline double computeDf(const double xm1, const double xp1)
    {
        return (xp1 - xm1)/2.0;
    }

    inline double computeD2f(const double xm1, const double x0, const double xp1)
    {
        return xp1 - 2.0 * x0 + xm1;
    }

    inline double computeDf2(const double xm1, const double xp1)
    {
        return pow(computeDf(xm1, xp1), 2.0);
    }

    inline double computeDx2(const double dx)
    {
        return pow(dx, 2.0);
    }
};

} // namespace Met3D

#endif // GRADIENTFILTER_H
