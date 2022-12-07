/******************************************************************************
**
**  This file is part of Met.3D -- a research environment for the
**  three-dimensional visual exploration of numerical ensemble weather
**  prediction data.
**
**  Copyright 2015-2021 Marc Rautenhaus
**  Copyright 2020-2021 Andreas Beckert
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

#include "partialderivativefilter.h"

// standard library imports
#include "assert.h"

// related third party imports
#include <log4cplus/loggingmacros.h>

// local application imports
#include "gxfw/nwpactorvariableproperties.h"
#include "util/mutil.h"
#include "util/mexception.h"
#include "util/metroutines.h"

using namespace std;

namespace Met3D
{

/******************************************************************************
***                     CONSTRUCTOR / DESTRUCTOR                            ***
*******************************************************************************/

MPartialDerivativeFilter::MPartialDerivativeFilter()
        : MSingleInputProcessingWeatherPredictionDataSource()
{

}
/******************************************************************************
***                            PUBLIC METHODS                               ***
*******************************************************************************/

MStructuredGrid *MPartialDerivativeFilter::produceData(Met3D::MDataRequest request)
{
    assert(inputSource != nullptr);

    MDataRequestHelper rh(request);
    // Parse request.
    // Examples: DERIVATIVE=D/LON, =D2/LAT
    QStringList parameterList = rh.value("GRADIENT").split("/");
    rh.removeAll(locallyRequiredKeys());
    // The first parameter passes the filter type.
    MGradientProperties::GradientModeTypes filterType =
            static_cast<MGradientProperties::GradientModeTypes>(
                parameterList[0].toInt());

    MStructuredGrid* inputGrid = inputSource->getData(rh.request());
    if (inputGrid->getDataType() == SINGLE)
    {
        inputGrid->copyFloatDataToDouble();
    }
    MStructuredGrid *resultGrid = createAndInitializeResultGrid(inputGrid);
    resultGrid->initializeDoubleData();

    MStructuredGrid *geoPotGrid = nullptr;
    if (geoPotGrid != nullptr)
    {
        if (geoPotGrid->getDataType() == SINGLE)
        {
            geoPotGrid->initializeDoubleData();
            geoPotGrid->copyFloatDataToDouble();
        }
    }
    QString gradientModeName = MGradientProperties::gradientModeToString(
                filterType);
    LOG4CPLUS_DEBUG(mlog, "Gradient filter: computing "
                    << gradientModeName.toUtf8().constData());

    QString levelType = MStructuredGrid::verticalLevelTypeToString(
                inputGrid->getLevelType());
    LOG4CPLUS_DEBUG(mlog, "Vertical level type:"
                    << levelType.toUtf8().constData());

    switch (filterType)
    {
    case MGradientProperties::DLON:
    {
        computePartialDerivativeLongitude(inputGrid, resultGrid);
        computePressureCoordinateTransormationLongitude(inputGrid, resultGrid);
        break;
    }
    case MGradientProperties::DLON_LAGRANTO:
    {
        callLibcalvarGradientRoutine(inputGrid, resultGrid,
                                               "lon");
        break;
    }
    case MGradientProperties::DLON_SOBEL:
    {
        computePartialDerivativeLongitudeSobel(inputGrid, resultGrid);
        break;
    }
    case MGradientProperties::DLAT:
    {
        computePartialDerivativeLatitude(inputGrid, resultGrid);
        computePressureCoordinateTransormationLatitude(inputGrid, resultGrid);
        break;
    }
    case MGradientProperties::DLAT_LAGRANTO:
    {
        callLibcalvarGradientRoutine(inputGrid, resultGrid,
                                               "lat");
        break;
    }
    case MGradientProperties::DLAT_SOBEL:
    {
        computePartialDerivativeLatitudeSobel(inputGrid, resultGrid);
        break;
    }
    case MGradientProperties::DP:
    {
        computePartialDerivativeVertical(inputGrid, resultGrid);
        break;
    }
    case MGradientProperties::DZ:
    {
        computePartialDerivativeVerticalGeometricHeight(inputGrid, geoPotGrid,
                                               resultGrid);
        break;
    }
    case MGradientProperties::D2LON:
    {
        computeSecondHorizontalDerivative(inputGrid, resultGrid, "lon");
        computeHorizontalPressureLevelTransformationForSecondDerivativeLon(
                    inputGrid, resultGrid);
        break;
    }
    case MGradientProperties::D2LAT:
    {
        computeSecondHorizontalDerivative(inputGrid, resultGrid, "lat");
        computeHorizontalPressureLevelTransformationForSecondDerivativeLat(
                    inputGrid, resultGrid);
        break;
    }
    case MGradientProperties::D2P:
    {
        computeSecondVerticalDerivative(inputGrid, resultGrid);
        break;
    }
    case MGradientProperties::D2Z:
    {
        computeSecondVerticalDerivativeGeometricHeight(inputGrid, geoPotGrid,
                                                       resultGrid);
        break;
    }
    default:
        LOG4CPLUS_DEBUG(mlog, "This gradient filter does not exists."
                        << gradientModeName.toUtf8().constData());
    }
    inputSource->releaseData(inputGrid);
    resultGrid->copyDoubleDataToFloat();
    return resultGrid;
}



MTask *MPartialDerivativeFilter::createTaskGraph(MDataRequest request)
{
    assert(inputSource != nullptr);
    MTask* task = new MTask(request, this);
    // Simply request the variable that was requested from this data source
    //(we're requesting the unsmoothed field and pass on the smoothed
    //version).
    MDataRequestHelper rh(request);
    rh.removeAll(locallyRequiredKeys());
    task->addParent(inputSource->getTaskGraph(rh.request()));
    return task;
}

/******************************************************************************
***                          PROTECTED METHODS                              ***
*******************************************************************************/

const QStringList MPartialDerivativeFilter::locallyRequiredKeys()
{
    return (QStringList() << "GRADIENT");
}


void MPartialDerivativeFilter::computePartialDerivativeLongitude(
        MStructuredGrid *inputGrid, MStructuredGrid *resultGrid)
{
    const int nLon = inputGrid->getNumLons();
    const int nLat = inputGrid->getNumLats();
    const int nLev = inputGrid->getNumLevels();
    // compute flags for pole and periodic treatment
    const QList<bool> periodicBC = periodicBoundaryTreatment(inputGrid);
#pragma omp parallel for
    for (int k = 0; k < nLev; k++)
    {
        int iP, iN;
        double dx;
        double dfdx;
        for (int j = 0; j < nLat; j++)
        {
            dx = inputGrid->getDeltaLon_km(j);
            for (int i = 0; i < nLon; i++)
            {
                iP = std::max(i - 1, 0);
                iN = std::min(i + 1, nLon - 1);
                dfdx = (inputGrid->getValue_double(k, j, iN)
                      - inputGrid->getValue_double(k, j, iP))
                      / (dx * double(abs(iN - iP)));
                resultGrid->setValue_double(k, j, i, dfdx);
            }
            if (periodicBC[0])
            {
                dfdx = (inputGrid->getValue_double(k, j, 1)
                      - inputGrid->getValue_double(k, j, nLon - 1)) / (2 * dx);
                resultGrid->setValue_double(k, j, 0, dfdx);
                dfdx = (inputGrid->getValue_double(k, j, 0)
                      - inputGrid->getValue_double(k, j, nLon - 2))  / (2 * dx);
                resultGrid->setValue_double(k, j, nLon - 1, dfdx);
            }
        }
    }
}


void MPartialDerivativeFilter::computePartialDerivativeLatitude(
        MStructuredGrid *inputGrid, MStructuredGrid *resultGrid)
{
    const int nLon = inputGrid->getNumLons();
    const int nLat = inputGrid->getNumLats();
    const int nLev = inputGrid->getNumLevels();
    const double dy = inputGrid->getDeltaLat_km();
    // compute flags for pole and periodic treatment
    const QList<bool> periodicBC = periodicBoundaryTreatment(inputGrid);

#pragma omp parallel for
    for (int k = 0; k < nLev; k++)
    {
        int jP, jN, iOpp;
        double dfdy;
        for (int j = 0; j < nLat; j++)
        {
            jP = std::max(j - 1, 0);
            jN = std::min(j + 1, nLat - 1);
            for (int i = 0; i < nLon; i++)
            {
                dfdy = (inputGrid->getValue_double(k, jN, i)
                      - inputGrid->getValue_double(k, jP, i))
                      / (dy * double(abs(jN - jP)));
                resultGrid->setValue_double(k, j, i, -dfdy);
            }
        }
        if (periodicBC[1])
        {
            for (int i = 0; i < nLon; i++)
            {
                iOpp = ((i + (nLon / 2)) % (nLon));
                dfdy = (inputGrid->getValue_double(k, 1, i)
                    - inputGrid->getValue_double(k, 0, iOpp)) / (2 * dy);
                resultGrid->setValue_double(k, 0, i, -dfdy);
            }
        }
        if (periodicBC[2])
        {
            for (int i = 0; i < nLon; i++)
            {
                iOpp = ((i + (nLon / 2)) % (nLon));
                dfdy = (inputGrid->getValue_double(k, nLat - 1, iOpp)
                    - inputGrid->getValue_double(k, nLat -  2, i)) / (2 * dy);
                resultGrid->setValue_double(k, nLat - 1, i, -dfdy);
            }
        }
    }
}


void MPartialDerivativeFilter::computePartialDerivativeLongitudeSobel(
        MStructuredGrid *inputGrid, MStructuredGrid *resultGrid)
{
    const int nLon = inputGrid->getNumLons();
    const int nLat = inputGrid->getNumLats();
    const int nLev = inputGrid->getNumLevels();
    const double sobel_x[5][5] =
    { {-2., -1., 0., 1., 2},
      {-2., -1., 0., 1., 2},
      {-4., -2., 0., 2., 4},
      {-2., -1., 0., 1., 2},
      {-2., -1., 0., 1., 2} };
#pragma omp parallel for
    for (int k = 0; k < nLev; k++)
    {
        int iP1, iP2, iN1, iN2, jP1, jP2, jN1, jN2;
        double df, dx, dfdx;
        for (int j = 2; j < nLat - 2; j++)
        {
            dx = inputGrid->getDeltaLon_km(j);
            for (int i = 2; i < nLon - 2; i++)
            {
                iP2 = max(0, i - 2);
                iP1 = max(0, i - 1);
                iN2 = min(nLon - 1, i + 2);
                iN1 = min(nLon - 1, i + 1);
                jP2 = max(0, j - 2);
                jP1 = max(0, j - 1);
                jN2 = min(nLat - 1, j + 2);
                jN1 = min(nLat - 1, j + 1);
                df = inputGrid->getValue_double(k, jP2, iP2) * sobel_x[0][0]
                   + inputGrid->getValue_double(k, jP2, iP1) * sobel_x[0][1]
                   + inputGrid->getValue_double(k, jP2, i    ) * sobel_x[0][2]
                   + inputGrid->getValue_double(k, jP2, iN1) * sobel_x[0][3]
                   + inputGrid->getValue_double(k, jP2, iN2) * sobel_x[0][4]

                   + inputGrid->getValue_double(k, jP1, iP2) * sobel_x[1][0]
                   + inputGrid->getValue_double(k, jP1, iP1) * sobel_x[1][1]
                   + inputGrid->getValue_double(k, jP1, i    ) * sobel_x[1][2]
                   + inputGrid->getValue_double(k, jP1, iN1) * sobel_x[1][3]
                   + inputGrid->getValue_double(k, jP1, iN2) * sobel_x[1][4]

                   + inputGrid->getValue_double(k, j    , iP2) * sobel_x[2][0]
                   + inputGrid->getValue_double(k, j    , iP1) * sobel_x[2][1]
                   + inputGrid->getValue_double(k, j    , i    ) * sobel_x[2][2]
                   + inputGrid->getValue_double(k, j    , iN1) * sobel_x[2][3]
                   + inputGrid->getValue_double(k, j    , iN2) * sobel_x[2][4]

                   + inputGrid->getValue_double(k, jN1, iP2) * sobel_x[3][0]
                   + inputGrid->getValue_double(k, jN1, iP1) * sobel_x[3][1]
                   + inputGrid->getValue_double(k, jN1, i    ) * sobel_x[3][2]
                   + inputGrid->getValue_double(k, jN1, iN1) * sobel_x[3][3]
                   + inputGrid->getValue_double(k, jN1, iN2) * sobel_x[3][4]

                   + inputGrid->getValue_double(k, jN2, iP2) * sobel_x[4][0]
                   + inputGrid->getValue_double(k, jN2, iP1) * sobel_x[4][1]
                   + inputGrid->getValue_double(k, jN2, i    ) * sobel_x[4][2]
                   + inputGrid->getValue_double(k, jN2, iN1) * sobel_x[4][3]
                   + inputGrid->getValue_double(k, jN2, iN2) * sobel_x[4][4];
                df = df / 25.;
                dfdx = df / (dx * 2);
                resultGrid->setValue_double(k, j, i, dfdx);
            }
        }
    }
}


void MPartialDerivativeFilter::computePartialDerivativeLatitudeSobel(
        MStructuredGrid *inputGrid, MStructuredGrid *resultGrid)
{
    const int nLon = inputGrid->getNumLons();
    const int nLat = inputGrid->getNumLats();
    const int nLev = inputGrid->getNumLevels();
    const double sobel_x[5][5] =
    { { -2., -2., -4., -2., -2.},
      { -1., -1., -2., -1., -1.},
      {  0.,  0.,  0.,  0.,  0.},
      {  2.,  2.,  4.,  2.,  2.},
      {  1.,  1.,  2.,  1.,  1.},};

#pragma omp parallel for
    for (int k = 0; k < nLev; k++)
    {
        int iP1, iP2, iN1, iN2, jP1, jP2, jN1, jN2;
        double dy = inputGrid->getDeltaLat_km();
        double df, dfdy;
        for (int j = 2; j < nLat - 2; j++)
        {
            for (int i = 2; i < nLon - 2; i++)
            {
                iP2 = max(0, i - 2);
                iP1 = max(0, i - 1);
                iN2 = min(nLon - 1, i + 2);
                iN1 = min(nLon - 1, i + 1);
                jP2 = max(0, j - 2);
                jP1 = max(0, j - 1);
                jN2 = min(nLat - 1, j + 2);
                jN1 = min(nLat - 1, j + 1);
                df = inputGrid->getValue_double(k, jP2, iP2) * sobel_x[0][0]
                   + inputGrid->getValue_double(k, jP2, iP1) * sobel_x[0][1]
                   + inputGrid->getValue_double(k, jP2, i    ) * sobel_x[0][2]
                   + inputGrid->getValue_double(k, jP2, iN1) * sobel_x[0][3]
                   + inputGrid->getValue_double(k, jP2, iN2) * sobel_x[0][4]

                   + inputGrid->getValue_double(k, jP1, iP2) * sobel_x[1][0]
                   + inputGrid->getValue_double(k, jP1, iP1) * sobel_x[1][1]
                   + inputGrid->getValue_double(k, jP1, i    ) * sobel_x[1][2]
                   + inputGrid->getValue_double(k, jP1, iN1) * sobel_x[1][3]
                   + inputGrid->getValue_double(k, jP1, iN2) * sobel_x[1][4]

                   + inputGrid->getValue_double(k, j    , iP2) * sobel_x[2][0]
                   + inputGrid->getValue_double(k, j    , iP1) * sobel_x[2][1]
                   + inputGrid->getValue_double(k, j    , i    ) * sobel_x[2][2]
                   + inputGrid->getValue_double(k, j    , iN1) * sobel_x[2][3]
                   + inputGrid->getValue_double(k, j    , iN2) * sobel_x[2][4]

                   + inputGrid->getValue_double(k, jN1, iP2) * sobel_x[3][0]
                   + inputGrid->getValue_double(k, jN1, iP1) * sobel_x[3][1]
                   + inputGrid->getValue_double(k, jN1, i    ) * sobel_x[3][2]
                   + inputGrid->getValue_double(k, jN1, iN1) * sobel_x[3][3]
                   + inputGrid->getValue_double(k, jN1, iN2) * sobel_x[3][4]

                   + inputGrid->getValue_double(k, jN2, iP2) * sobel_x[4][0]
                   + inputGrid->getValue_double(k, jN2, iP1) * sobel_x[4][1]
                   + inputGrid->getValue_double(k, jN2, i    ) * sobel_x[4][2]
                   + inputGrid->getValue_double(k, jN2, iN1) * sobel_x[4][3]
                   + inputGrid->getValue_double(k, jN2, iN2) * sobel_x[4][4];
                df = df / 25.;
                dfdy = df / (dy * 2);
                resultGrid->setValue_double(k, j, i, -dfdy);
            }
        }
    }
}


void MPartialDerivativeFilter::computePartialDerivativeVertical(
        MStructuredGrid *inputGrid, MStructuredGrid *resultGrid)
{
    const int nLon = inputGrid->getNumLons();
    const int nLat = inputGrid->getNumLats();
    const int nLev = inputGrid->getNumLevels();

    #pragma omp parallel for
    for (int k = 0; k < nLev; k++)
    {
        int kP, kN;
        double df, dp;
        for (int j = 0; j < nLat; j++)
        {
            for (int i = 0; i < nLon; i++)
            {
                kP = std::max(k - 1, 0);
                kN = std::min(k + 1, nLev - 1);
                df = inputGrid->getValue_double(kN, j, i)
                        - inputGrid->getValue_double(kP, j, i);
                dp = inputGrid->getPressure(kN, j, i)
                        - inputGrid->getPressure(kP, j, i);
                //When df is zero, then the result should be zero, else it wouldn't be defined
                float result;
                if(df == 0)
                {
                    result = 0;
                }
                else
                {
                    result = df/dp;
                }

                resultGrid->setValue_double(k, j, i, df / dp);

            }
        }
    }
}


void MPartialDerivativeFilter::computePartialDerivativePressureLongitude(
            MStructuredGrid *inputGrid, MStructuredGrid *resultGrid)
{
    const int nLon = inputGrid->getNumLons();
    const int nLat = inputGrid->getNumLats();
    const int nLev = inputGrid->getNumLevels();
    // compute flags for pole and periodic treatment
    const QList<bool> periodicBC = periodicBoundaryTreatment(inputGrid);

#pragma omp parallel for
    for (int k = 0; k < nLev; k++)
    {
        int iP, iN;
        double dx;
        double dpdx;
        for (int j = 0; j < nLat; j++)
        {
            dx = inputGrid->getDeltaLon_km(j);
            for (int i = 0; i < nLon; i++)
            {
                iP = std::max(i - 1, 0);
                iN = std::min(i + 1, nLon - 1);
                dpdx = (inputGrid->getPressure(k, j, iN)
                      - inputGrid->getPressure(k, j, iP))
                      / (dx * double(abs(iN - iP)));
                resultGrid->setValue_double(k, j, i, dpdx);
            }
            if (periodicBC[0])
            {
                dpdx = (inputGrid->getPressure(k, j, 1)
                      - inputGrid->getPressure(k, j, nLon - 1)) / (2 * dx);
                resultGrid->setValue_double(k, j, 0, dpdx);
                dpdx = (inputGrid->getPressure(k, j, 0)
                      - inputGrid->getPressure(k, j, nLon - 2))  / (2 * dx);
                resultGrid->setValue_double(k, j, nLon - 1, dpdx);
            }
        }
    }
}


void MPartialDerivativeFilter::computePartialDerivativePressureLatitude(
        MStructuredGrid *inputGrid, MStructuredGrid *resultGrid)
{
    const int nLon = inputGrid->getNumLons();
    const int nLat = inputGrid->getNumLats();
    const int nLev = inputGrid->getNumLevels();
    const double dy = inputGrid->getDeltaLat_km();
    // compute flags for pole and periodic treatment
    const QList<bool> periodicBC = periodicBoundaryTreatment(inputGrid);

#pragma omp parallel for
    for (int k = 0; k < nLev; k++)
    {
        double dpdy;
        int jP, jN, iOpp;
        for (int j = 0; j < nLat; j++)
        {
            jP = std::max(j - 1, 0);
            jN = std::min(j + 1, nLat - 1);
            for (int i = 0; i < nLon; i++)
            {
                dpdy = (inputGrid->getPressure(k, jN, i)
                      - inputGrid->getPressure(k, jP, i))
                      / (dy * double(abs(jN - jP)));
                resultGrid->setValue_double(k, j, i, dpdy);
            }
        }
        if (periodicBC[1])
        {
            for (int i = 0; i < nLon; i++)
            {
                iOpp = ((i + (nLon / 2)) % (nLon));
                dpdy = (inputGrid->getPressure(k, 1, i)
                    - inputGrid->getPressure(k, 0, iOpp)) / (2 * dy);
                resultGrid->setValue_double(k, 0, i, dpdy);
            }
        }
        if (periodicBC[2])
        {
            for (int i = 0; i < nLon; i++)
            {
                iOpp = ((i + (nLon / 2)) % (nLon));
                dpdy = (inputGrid->getPressure(k, nLat - 1, iOpp)
                    - inputGrid->getPressure(k, nLat -  2, i)) / (2 * dy);
                resultGrid->setValue_double(k, nLat - 1, i, dpdy);
            }
        }
    }
}


void MPartialDerivativeFilter::computePressureCoordinateTransormationLongitude(
        MStructuredGrid *inputGrid, MStructuredGrid *resultGrid)
{
    MStructuredGrid *dfdp = createAndInitializeResultGrid(inputGrid);
    dfdp->copyFloatDataToDouble();
    MStructuredGrid *dpdx = createAndInitializeResultGrid(inputGrid);
    dpdx->copyFloatDataToDouble();


    const int nLon = inputGrid->getNumLons();
    const int nLat = inputGrid->getNumLats();
    const int nLev = inputGrid->getNumLevels();

    computePartialDerivativeVertical(inputGrid, dfdp);
    computePartialDerivativePressureLongitude(inputGrid, dpdx);

#pragma omp parallel for
    for (int k = 0; k < nLev; k++)
    {
        double dfdxPGrid;
        for (int j = 0; j < nLat; j++)
        {
            for (int i = 0; i < nLon; i++)
            {
                dfdxPGrid = resultGrid->getValue_double(k, j, i)
                        - dfdp->getValue_double(k, j, i) * dpdx->getValue_double(k, j, i);
                resultGrid->setValue_double(k, j, i, dfdxPGrid);
            }
        }
    }
    delete dfdp;
    delete dpdx;
}


void MPartialDerivativeFilter::computePressureCoordinateTransormationLatitude(
        MStructuredGrid *inputGrid, MStructuredGrid *resultGrid)
{
    MStructuredGrid *dfdp = createAndInitializeResultGrid(inputGrid);
    dfdp->copyFloatDataToDouble();
    MStructuredGrid *dpdy = createAndInitializeResultGrid(inputGrid);
    dpdy->copyFloatDataToDouble();

    const int nLon = inputGrid->getNumLons();
    const int nLat = inputGrid->getNumLats();
    const int nLev = inputGrid->getNumLevels();

    computePartialDerivativeVertical(inputGrid, dfdp);
    computePartialDerivativePressureLatitude(inputGrid, dpdy);

#pragma omp parallel for
    for (int k = 0; k < nLev; k++)
    {
        double dfdyPGrid;
        for (int j = 0; j < nLat; j++)
        {
            for (int i = 0; i < nLon; i++)
            {
                dfdyPGrid = resultGrid->getValue_double(k, j, i)
                        - dfdp->getValue_double(k, j, i) * dpdy->getValue_double(k, j, i);
                resultGrid->setValue_double(k, j, i, dfdyPGrid);
            }
        }
    }
    delete dfdp;
    delete dpdy;
}


void MPartialDerivativeFilter::computePartialDerivativeVerticalGeometricHeight(
        MStructuredGrid *inputGrid, MStructuredGrid *geoPot,
        MStructuredGrid *resultGrid)
{
    const int nLon = inputGrid->getNumLons();
    const int nLat = inputGrid->getNumLats();
    const int nLev = inputGrid->getNumLevels();

#pragma omp parallel for
    for (int k = 0; k < nLev; k++)
    {
        int kP, kN;
        double df, dz;
        for (int j = 0; j < nLat; j++)
        {
            for (int i = 0; i < nLon; i++)
            {
                kP = std::max(k - 1, 0);
                kN = std::min(k + 1, nLev - 1);
                df = inputGrid->getValue_double(kN, j, i)
                        - inputGrid->getValue_double(kP, j, i);
                dz = geoPot->getPressure(kN, j, i)
                        - geoPot->getPressure(kP, j, i);
                resultGrid->setValue_double(k, j, i, - df / dz);
            }
        }
    }
}


// NO PERIODIC BOUNDARY THREATMENT
void MPartialDerivativeFilter::callLibcalvarGradientRoutine(
        MStructuredGrid *inputGrid, MStructuredGrid *resultGrid,
        QString direction)
{
    LOG4CPLUS_DEBUG(mlog, "You are using the libcalvar LAGRANTO library. "
                       << "This library is implemented for testing purpose only. "
                       << "No periodic boundary treatment!");
    // Cast the input grid to hybrid sigma pressure grid to access ak/bk
    // coefficients; convert the float arrays in MLonLatHybridSigmaPressureGrid
    // to float arrays.
    MLonLatHybridSigmaPressureGrid *hybridInputGrid =
            dynamic_cast<MLonLatHybridSigmaPressureGrid*>(inputGrid);

    // This method uses the LAGRANTO.ECMWF libcalvar function "ddh3" to
    // compute gradients. To call the FORTRAN function "ddh3",
    // the data contained in the MStructuredGrid classes needs to be
    // restructured:
    // * libcalvar requires float arrays that contain the full 3D variable
    //   fields.
    // * libcalvar requires the lat dimension to be reversed (in increasing
    //   lat order) in all spatial fields.
    // * surface pressure needs to be passed in hPa.
    // * ak and bk coefficients need to be passed as float arrays.
    // * dps surface pressure gradient needs to be precomputed and passed in
    //   hPa
    // Also compare to libcalvar usage in ppecmwf.py in the met.dp repository.

    // Grid sizes.
    int nlev = inputGrid->getNumLevels();
    int nlat = inputGrid->getNumLats();
    int nlon = inputGrid->getNumLons();
    int nlatnlon = nlat*nlon;

    // Convert surface pressure from Pa to hPa; reverse lat dimension.
    float *psfc_hPa_revLat = new float[nlat*nlon];
    for (int j = 0; j < nlat; j++)
        for (int i = 0; i < nlon; i++)
        {
            psfc_hPa_revLat[INDEX2yx(j, i, nlon)] =
                    hybridInputGrid->getSurfacePressure(nlat-1-j, i) / 100.;
        }


    float *ak_hPa_float = new float[nlev];
    float *bk_float = new float[nlev];
    for (int k = 0; k < nlev; k++)
    {
        ak_hPa_float[k] = hybridInputGrid->getAkCoeff(k);
        bk_float[k] = hybridInputGrid->getBkCoeff(k);
    }


    // compute the inputGrid with reverse lat dimension as a for "ddh3"
    float *a_revLat = new float[nlev*nlat*nlon];
    for (int k = 0; k < nlev; k++)
        for (int j = 0; j < nlat; j++)
            for (int i = 0; i < nlon; i++)
            {
                a_revLat[INDEX3zyx_2(k, j, i, nlatnlon, nlon)] =
                        inputGrid->getValue(k, nlat-1-j, i);

            }

    // Compute 2D field of cos(lat). Reverse lat
    // dimension.
    float *coslat_revLat = new float[nlat*nlon];
    for (int j = 0; j < nlat; j++)
        for (int i = 0; i < nlon; i++)
        {
            coslat_revLat[INDEX2yx(j, i, nlon)] =
                    cos(inputGrid->getLats()[nlat-1-j] / 180. * M_PI);
        }


    // "ddh3" requires two 4-element vectors that contain the lon/lat range.
    float *varmin = new float[4];
    varmin[0] = inputGrid->getLons()[0]; // min lon
    varmin[1] = inputGrid->getLats()[nlat-1]; // min lat
    varmin[2] = 0.;
    varmin[3] = 0.;
    float *varmax = new float[4];
    varmax[0] = inputGrid->getLons()[nlon-1]; // max lon
    varmax[1] = inputGrid->getLats()[0]; // max lat
    varmax[2] = 0.;
    varmax[3] = 0.;

    int dir;
    float df;
    float dx, dy;
    int i, j;


    // compute "dps" surface pressure gradient according to the
    // chosen direction and convert pressure from Pa to hPa
    // convert direction to an integer, reverse lat and direction
    // of latitudinal gradient
    float *dps_revLat = new float[nlat*nlon];
    if (direction == "lon")
    {
        dir = 0;
        for (j = 0; j < nlat; j++)
        {
            //dx = inputGrid->getDeltaLon_km(nlat-1-j) * 1000;
            // do not convert in m
            dx = inputGrid->getDeltaLon_km(nlat-1-j);
            // central differences
            for (i = 1; i < nlon - 1; i++)
            {
                df = (hybridInputGrid->getSurfacePressure(nlat-1-j, i + 1)
                        - hybridInputGrid->getSurfacePressure(nlat-1-j, i - 1))
                        / (2. * 100.);
                dps_revLat[INDEX2yx(j, i, nlon)] = df / dx;
            };
            // forward differences for the first index
            i = 0;
            df = (hybridInputGrid->getSurfacePressure(nlat-1-j, i + 1)
                    - hybridInputGrid->getSurfacePressure(nlat-1-j, i))
                    / (1. * 100.);
            dps_revLat[INDEX2yx(j, i, nlon)] = df / dx;
            // backward differences for the last index
            i = nlon - 1;
            df = (hybridInputGrid->getSurfacePressure(nlat-1-j, i)
                    - hybridInputGrid->getSurfacePressure(nlat-1-j, i - 1))
                    / (1. * 100.);
            dps_revLat[INDEX2yx(j, i, nlon)] = df / dx;
        };
    }
    else if (direction == "lat")
    {
        dir = 1;
        //dy = inputGrid->getDeltaLat_km() * 1000;
        //do not convert in m
        dy = inputGrid->getDeltaLat_km();
        for (i = 0; i < nlon; i++)
        {
            // central differences
            for (j = 1; j < nlat - 1; j++)
            {
                df = (hybridInputGrid->getSurfacePressure(nlat-1-j + 1, i)
                        - hybridInputGrid->getSurfacePressure(nlat-1-j - 1, i))
                        / (2. * 100.);
                dps_revLat[INDEX2yx(j, i, nlon)] = -(df / dy);
            };
            // backward differences for the last index, revers lat
            j = 0;
            df = (hybridInputGrid->getSurfacePressure(nlat - 1, i)
                    - hybridInputGrid->getSurfacePressure(nlat - 2, i))
                    / (1. * 100.);
            dps_revLat[INDEX2yx(j, i, nlon)] = -(df / dy);
            // forward differences for the first index, revers lat
            j = nlat - 1;
            df = (hybridInputGrid->getSurfacePressure(1, i)
                    - hybridInputGrid->getSurfacePressure(0, i))
                    / (1. * 100.);
            dps_revLat[INDEX2yx(j, i, nlon)] = -(df / dy);
        };
    }
    else
    {
        LOG4CPLUS_DEBUG(mlog, "Upsi, something went wrong! This direction to "
                           << "calculate the gradients does not exists."
                           << "Please try with another direction again.");
        exit(1);
    }

    // Call the "ddh3" LAGRANTO function inside the FORTRAN libcalvar library.
    float *gradient_revLat = new float[nlev*nlat*nlon];
    horizontalGradient_calvar(
                a_revLat, gradient_revLat, psfc_hPa_revLat, dps_revLat,
                coslat_revLat, dir, nlon, nlat, nlev, varmin, varmax,
                ak_hPa_float, bk_float);

    // Reverse the lat dimension of the computed gradient field and store the
    // gradient field in the result grid
    // convert the pressure from hPa to Pa
    for (int k = 0; k < nlev; k++)
        for (int j = 0; j < nlat; j++)
            for (int i = 0; i < nlon; i++)
            {
                resultGrid->setValue(k, j, i,
                                      gradient_revLat[INDEX3zyx_2(
                            k, nlat-1-j, i, nlatnlon, nlon)]);
            }


    // Delete temporary memory.
    delete[] psfc_hPa_revLat;
    delete[] ak_hPa_float;
    delete[] bk_float;
    delete[] coslat_revLat;
    delete[] a_revLat;
    delete[] varmin;
    delete[] varmax;
    delete[] gradient_revLat;
    delete[] dps_revLat;
}


QList<bool> MPartialDerivativeFilter::periodicBoundaryTreatment(
        MStructuredGrid *inputGrid)
{
    // compute flags for pole and periodic treatment
    bool southpl = false;
    bool northpl = false;
    bool lonper  = false;
    QList<bool> periodic;
    double phiTotal = inputGrid->getLons()[inputGrid->getNumLons() - 1]
            - inputGrid->getLons()[0]
            + inputGrid->getDeltaLon();
    // periodic boundaries test (x-dir)
    // TODO (ab, Mar2020) test periodic boundary conditions
    if (phiTotal >= 360.)
    {
        lonper = true;
        if ((inputGrid->getLats()[inputGrid->getNumLats() - 1]
             - inputGrid->getDeltaLat()) <= -90.)
        {
            southpl = true;
        }
        if ((inputGrid->getLats()[0] + inputGrid->getDeltaLat()) >= 90)
        {
            northpl = true;
        }
    }
    periodic.append(lonper);
    periodic.append(northpl);
    periodic.append(southpl);
    return periodic;
}

void MPartialDerivativeFilter::computeSecondHorizontalDerivative(
        MStructuredGrid *inputGrid, MStructuredGrid *resultGrid,
        QString direction)
{
    int nLat = inputGrid->getNumLats();
    int nLon = inputGrid->getNumLons();
    int iOpp;
    double dx, dy;
    double d2f;
    double dx2, dy2;

    // compute flags for pole and periodic treatment
    QList<bool> periodicBC = periodicBoundaryTreatment(inputGrid);

    // compute horizontal gradient in longitudinal direction.
    if (direction == "lon")
    {
        for (unsigned int k = 0; k < inputGrid->getNumLevels(); k++)
        {
            for (unsigned int j = 0; j < inputGrid->getNumLats(); j++)
            {
                // convert km to m
                //dx = inputGrid->getDeltaLon_km(j) * 1000;
                // do not convert in m
                dx = inputGrid->getDeltaLon_km(j);
                dx2 = computeDx2(dx);
                for (int i = 1; i < nLon - 1; i++)
                {
                    d2f = computeD2f(inputGrid->getValue_double(k, j, i - 1),
                                     inputGrid->getValue_double(k, j, i),
                                     inputGrid->getValue_double(k, j, i + 1));
                    resultGrid->setValue_double(k, j, i, d2f / dx2);
                };
                // i=0
                d2f = computeD2f(inputGrid->getValue_double(k, j, 2),
                                 inputGrid->getValue_double(k, j, 1),
                                 inputGrid->getValue_double(k, j, 0));
                resultGrid->setValue_double(k, j, 0, d2f / dx2);
                // i=nlon
                d2f = computeD2f(inputGrid->getValue_double(k, j, nLon - 1),
                                 inputGrid->getValue_double(k, j, nLon - 2),
                                 inputGrid->getValue_double(k, j, nLon - 3));
                resultGrid->setValue_double(k, j, nLon - 1, d2f / dx2);

                if (periodicBC[0])
                {
                    d2f = computeD2f(inputGrid->getValue_double(k, j, nLon - 1),
                                     inputGrid->getValue_double(k, j, 0),
                                     inputGrid->getValue_double(k, j, 1));
                    resultGrid->setValue_double(k, j, 0, d2f / dx2);

                    d2f = computeD2f(inputGrid->getValue_double(k, j, nLon - 2),
                                     inputGrid->getValue_double(k, j, nLon - 1),
                                     inputGrid->getValue_double(k, j, 0));
                    resultGrid->setValue_double(k, j, nLon - 1, d2f / dx2);
                }
            }
        }
    }
    // compute horizontal gradient in latitudinal direction.
    // reverse the direction of the gradient
    if (direction == "lat")
    {
        // convert km to m
        //dy = inputGrid->getDeltaLat_km() * 1000;
        // do not convert in m
        dy = inputGrid->getDeltaLat_km();
        dy2 = computeDx2(dy);
        for (unsigned int k = 0; k < inputGrid->getNumLevels(); k++)
        {
            for (unsigned int i = 0; i < inputGrid->getNumLons(); i++)
            {
                for (int j = 1; j < nLat - 1; j++)
                {
                    d2f = computeD2f(inputGrid->getValue_double(k, j - 1, i),
                                     inputGrid->getValue_double(k, j, i),
                                     inputGrid->getValue_double(k, j + 1, i));
                    resultGrid->setValue_double(k, j, i, d2f / dy2);
                };

                // j=0
                d2f = computeD2f(inputGrid->getValue_double(k, 2, i),
                                 inputGrid->getValue_double(k, 1, i),
                                 inputGrid->getValue_double(k, 0, i));
                resultGrid->setValue_double(k, 0, i, d2f / dy2);
                if (periodicBC[1])
                {
                    iOpp = ((i + (nLon / 2)) % (nLon - 1));
                    d2f = computeD2f(inputGrid->getValue_double(k, 0, iOpp),
                                     inputGrid->getValue_double(k, 0, i),
                                     inputGrid->getValue_double(k, 1, i));
                    resultGrid->setValue_double(k, 0, i, d2f / dy2);
                }

                // j= nLat-1
                d2f = computeD2f(inputGrid->getValue_double(k, nLat - 1, i),
                                 inputGrid->getValue_double(k, nLat - 2, i),
                                 inputGrid->getValue_double(k, nLat - 3, i));
                resultGrid->setValue_double(k, nLat - 1, i, d2f / dy2);
                if (periodicBC[2])
                {
                    iOpp = ((i + (nLon / 2)) % (nLon - 1));
                    d2f = computeD2f(inputGrid->getValue_double(k, nLat - 1, iOpp),
                                     inputGrid->getValue_double(k, nLat - 1, i),
                                     inputGrid->getValue_double(k, nLat - 2, i));
                    resultGrid->setValue_double(k, nLat - 1, i, d2f / dy2);
                }
            };
        };
    };
}


void MPartialDerivativeFilter::computeSecondVerticalDerivative(
        MStructuredGrid *inputGrid, MStructuredGrid *resultGrid)
{
    double d2f, dp2;
    int nLevel = inputGrid->getNumLevels();
    for (unsigned int j = 0; j < inputGrid->getNumLats(); j++)
    {
        for (unsigned int i = 0; i < inputGrid->getNumLons(); i++)
        {
            // upper and lower boundary
            dp2 = computeDx2(inputGrid->getPressure(0, j, i)
                             - inputGrid->getPressure(1, j, i));
            d2f = computeD2f(inputGrid->getValue_double(2, j, i),
                             inputGrid->getValue_double(1, j, i),
                             inputGrid->getValue_double(0, j, i));
            resultGrid->setValue_double(0, j, i, (d2f / dp2));

            dp2 = computeDx2(inputGrid->getPressure(nLevel - 2, j, i)
                             - inputGrid->getPressure(nLevel - 1, j, i));
            d2f = computeD2f(inputGrid->getValue_double(nLevel - 1, j, i),
                             inputGrid->getValue_double(nLevel - 2, j, i),
                             inputGrid->getValue_double(nLevel - 3, j, i));
            resultGrid->setValue_double(nLevel - 1, j, i, (d2f / dp2));

            for (int k = 1; k < nLevel - 1; k++)
            {
                dp2 = computeDf2(inputGrid->getPressure(k - 1, j, i),
                                 inputGrid->getPressure(k + 1, j, i));
                d2f = computeD2f(inputGrid->getValue_double(k - 1, j, i),
                                 inputGrid->getValue_double(k, j, i),
                                 inputGrid->getValue_double(k + 1, j, i));
                resultGrid->setValue_double(k, j, i, (d2f / dp2));
            }
        }
    }
}

void MPartialDerivativeFilter::computeSecondVerticalDerivativeGeometricHeight(
        MStructuredGrid *inputGrid, MStructuredGrid *geoPot,
        MStructuredGrid *resultGrid)
{
    double d2f, dz2;
    int nLevel = inputGrid->getNumLevels();
    for (unsigned int j = 0; j < inputGrid->getNumLats(); j++)
    {
        for (unsigned int i = 0; i < inputGrid->getNumLons(); i++)
        {
            // upper and lower boundary
            dz2 = computeDx2(geoPot->getValue_double(0, j, i)
                             - geoPot->getValue_double(1, j, i));
            d2f = computeD2f(inputGrid->getValue_double(2, j, i),
                             inputGrid->getValue_double(1, j, i),
                             inputGrid->getValue_double(0, j, i));
            resultGrid->setValue_double(0, j, i, (d2f / dz2));

            dz2 = computeDx2(geoPot->getValue_double(nLevel - 2, j, i)
                             - geoPot->getValue_double(nLevel - 1, j, i));
            d2f = computeD2f(inputGrid->getValue_double(nLevel - 1, j, i),
                             inputGrid->getValue_double(nLevel - 2, j, i),
                             inputGrid->getValue_double(nLevel - 3, j, i));
            resultGrid->setValue_double(nLevel - 1, j, i, (d2f / dz2));

            for (int k = 1; k < nLevel - 1; k++)
            {
                dz2 = computeDf2(geoPot->getValue_double(k - 1, j, i),
                                 geoPot->getValue_double(k + 1, j, i));
                d2f = computeD2f(inputGrid->getValue_double(k - 1, j, i),
                                 inputGrid->getValue_double(k, j, i),
                                 inputGrid->getValue_double(k + 1, j, i));
                resultGrid->setValue_double(k, j, i, (d2f / dz2));
            }
        }
    }
}


void MPartialDerivativeFilter::computeHorizontalPressureLevelTransformationForSecondDerivativeLon(
        MStructuredGrid *inputGrid, MStructuredGrid *resultGrid)
{
    unsigned int nLat = inputGrid->getNumLats();
    unsigned int nLon = inputGrid->getNumLons();
    unsigned int nLev = inputGrid->getNumLevels();
    double df_p, dp_x, dp_p, dx;
    double d2f_p, d2p_x, dp2, dx2;
    double dfdp, dpdx;
    double d2fdp2, d2pdx2;
    double result;

    // compute flags for pole and periodic treatment
    QList<bool> periodicBC = periodicBoundaryTreatment(inputGrid);

    // compute horizontal gradient in longitudinal direction.
    for (unsigned int k = 1; k < nLev - 1; k++) // ToDo: Treat k=0 and k=nlev-1
    {
        for (unsigned int j = 0; j < nLat; j++)
        {
            // convert km to m
            //dx = inputGrid->getDeltaLon_km(j) * 1000;
            // do not convert in m
            dx = inputGrid->getDeltaLon_km(j);
            dx2 = computeDx2(dx);
            for (unsigned int i = 1; i < nLon - 1; i++)
            {
                // deviation of f in p dir
                df_p = computeDf(inputGrid->getValue_double(k-1, j, i),
                                 inputGrid->getValue_double(k+1, j, i));
                dp_p = computeDf(inputGrid->getPressure(k-1, j, i),
                                 inputGrid->getPressure(k+1, j, i));
                dfdp = df_p/dp_p;

                // deviation of p in x dir
                dp_x = computeDf(inputGrid->getPressure(k, j, i - 1),
                                 inputGrid->getPressure(k, j, i + 1));
                dpdx = dp_x/dx;

                // second deviaion of p in x dir
                d2p_x = computeD2f(inputGrid->getPressure(k, j, i - 1),
                                   inputGrid->getPressure(k, j, i),
                                   inputGrid->getPressure(k, j, i + 1));
                d2pdx2 = d2p_x/dx2;

                // second deviation of f in p dir
                d2f_p = computeD2f(inputGrid->getValue_double(k - 1, j, i),
                                 inputGrid->getValue_double(k, j, i),
                                 inputGrid->getValue_double(k + 1, j, i));
                dp2 = computeDx2(dp_p);
                d2fdp2 = d2f_p/dp2;

                result = resultGrid->getValue_double(k, j, i)
                        - d2fdp2 * (dpdx * dpdx)
                        - dfdp * d2pdx2
                        ;
                resultGrid->setValue_double(k, j, i, result);
            };
            //i = 0
            df_p = computeDf(inputGrid->getValue_double(k-1, j, 0),
                             inputGrid->getValue_double(k+1, j, 0));
            dp_p = computeDf(inputGrid->getPressure(k-1, j, 0),
                             inputGrid->getPressure(k+1, j, 0));
            dfdp = df_p/dp_p;

            // deviation of p in x dir
            dp_x = computeDf(inputGrid->getPressure(k, j, 0),
                             inputGrid->getPressure(k, j, 1)) * 2.0;
            dpdx = dp_x/dx;

            // second deviaion of p in x dir
            d2p_x = computeD2f(inputGrid->getPressure(k, j, 0),
                               inputGrid->getPressure(k, j, 1),
                               inputGrid->getPressure(k, j, 2));
            d2pdx2 = d2p_x/dx2;

            // second deviation of f in p dir
            d2f_p = computeD2f(inputGrid->getValue_double(k - 1, j, 0),
                             inputGrid->getValue_double(k, j, 0),
                             inputGrid->getValue_double(k + 1, j, 0));
            dp2 = computeDx2(dp_p);
            d2fdp2 = d2f_p/dp2;

            if (periodicBC[0])
            {
                // deviation of p in x dir
                dp_x = computeDf(inputGrid->getPressure(k, j, nLon -1),
                                 inputGrid->getPressure(k, j, 1));
                dpdx = dp_x/dx;

                // second deviaion of p in x dir
                d2p_x = computeD2f(inputGrid->getPressure(k, j, nLon - 1),
                                   inputGrid->getPressure(k, j, 0),
                                   inputGrid->getPressure(k, j, 1));
                d2pdx2 = d2p_x/dx2;
            }
            result = resultGrid->getValue_double(k, j, 0)
                    - d2fdp2 * (dpdx * dpdx)
                    - dfdp * d2pdx2
                    ;
            resultGrid->setValue_double(k, j, 0, result);

            //i = nLon - 1
            df_p = computeDf(inputGrid->getValue_double(k-1, j, nLon - 1),
                             inputGrid->getValue_double(k+1, j, nLon - 1));
            dp_p = computeDf(inputGrid->getPressure(k-1, j, nLon - 1),
                             inputGrid->getPressure(k+1, j, nLon - 1));
            dfdp = df_p/dp_p;

            // deviation of p in x dir
            dp_x = computeDf(inputGrid->getPressure(k, j, nLon - 2),
                             inputGrid->getPressure(k, j, nLon - 1)) * 2.0;
            dpdx = dp_x/dx;

            // second deviaion of p in x dir
            d2p_x = computeD2f(inputGrid->getPressure(k, j, nLon - 3),
                               inputGrid->getPressure(k, j, nLon - 2),
                               inputGrid->getPressure(k, j, nLon - 1));
            d2pdx2 = d2p_x/dx2;

            // second deviation of f in p dir
            d2f_p = computeD2f(inputGrid->getValue_double(k - 1, j, nLon - 1),
                             inputGrid->getValue_double(k, j, nLon - 1),
                             inputGrid->getValue_double(k + 1, j, nLon - 1));
            dp2 = computeDx2(dp_p);
            d2fdp2 = d2f_p/dp2;

            if (periodicBC[0])
            {
                // deviation of p in x dir
                dp_x = computeDf(inputGrid->getPressure(k, j, 0),
                                 inputGrid->getPressure(k, j, nLon - 2));
                dpdx = dp_x/dx;

                // second deviaion of p in x dir
                d2p_x = computeD2f(inputGrid->getPressure(k, j, 0),
                                   inputGrid->getPressure(k, j, nLon - 1),
                                   inputGrid->getPressure(k, j, nLon - 2));
                d2pdx2 = d2p_x/dx2;
            }
            result = resultGrid->getValue_double(k, j, nLon - 1)
                    - d2fdp2 * (dpdx * dpdx)
                    - dfdp * d2pdx2
                    ;
            resultGrid->setValue_double(k, j, nLon - 1, result);

        }
    }
    // k = 0
    for (unsigned int j = 0; j < nLat; j++)
    {
        // convert km to m
        //dx = inputGrid->getDeltaLon_km(j) * 1000;
        // do not convert in m
        dx = inputGrid->getDeltaLon_km(j);
        dx2 = computeDx2(dx);
        for (unsigned int i = 1; i < nLon - 1; i++)
        {
            // deviation of f in p dir
            df_p = computeDf(inputGrid->getValue_double(0, j, i),
                             inputGrid->getValue_double(1, j, i)) * 2;
            dp_p = computeDf(inputGrid->getPressure(0, j, i),
                             inputGrid->getPressure(1, j, i)) * 2;
            dfdp = df_p/dp_p;

            // deviation of p in x dir
            dp_x = computeDf(inputGrid->getPressure(0, j, i - 1),
                             inputGrid->getPressure(0, j, i + 1));
            dpdx = dp_x/dx;

            // second deviaion of p in x dir
            d2p_x = computeD2f(inputGrid->getPressure(0, j, i - 1),
                               inputGrid->getPressure(0, j, i),
                               inputGrid->getPressure(0, j, i + 1));
            d2pdx2 = d2p_x/dx2;

            // second deviation of f in p dir
            d2f_p = computeD2f(inputGrid->getValue_double(0, j, i),
                             inputGrid->getValue_double(1, j, i),
                             inputGrid->getValue_double(2, j, i));
            dp2 = computeDx2(dp_p);
            d2fdp2 = d2f_p/dp2;

            result = resultGrid->getValue_double(0, j, i)
                    - d2fdp2 * (dpdx * dpdx)
                    - dfdp * d2pdx2
                    ;
            resultGrid->setValue_double(0, j, i, result);
        };
        //i = 0
        df_p = computeDf(inputGrid->getValue_double(0, j, 0),
                         inputGrid->getValue_double(1, j, 0)) * 2.0;
        dp_p = computeDf(inputGrid->getPressure(0, j, 0),
                         inputGrid->getPressure(1, j, 0)) * 2.0;
        dfdp = df_p/dp_p;

        // deviation of p in x dir
        dp_x = computeDf(inputGrid->getPressure(0, j, 0),
                         inputGrid->getPressure(0, j, 1)) * 2.0;
        dpdx = dp_x/dx;

        // second deviaion of p in x dir
        d2p_x = computeD2f(inputGrid->getPressure(0, j, 0),
                           inputGrid->getPressure(0, j, 1),
                           inputGrid->getPressure(0, j, 2));
        d2pdx2 = d2p_x/dx2;

        // second deviation of f in p dir
        d2f_p = computeD2f(inputGrid->getValue_double(0, j, 0),
                         inputGrid->getValue_double(1, j, 0),
                         inputGrid->getValue_double(2, j, 0));
        dp2 = computeDx2(dp_p);
        d2fdp2 = d2f_p/dp2;

        if (periodicBC[0])
        {
            // deviation of p in x dir
            dp_x = computeDf(inputGrid->getPressure(0, j, nLon -1),
                             inputGrid->getPressure(0, j, 1));
            dpdx = dp_x/dx;

            // second deviaion of p in x dir
            d2p_x = computeD2f(inputGrid->getPressure(0, j, nLon - 1),
                               inputGrid->getPressure(0, j, 0),
                               inputGrid->getPressure(0, j, 1));
            d2pdx2 = d2p_x/dx2;
        }
        result = resultGrid->getValue_double(0, j, 0)
                - d2fdp2 * (dpdx * dpdx)
                - dfdp * d2pdx2
                ;
        resultGrid->setValue_double(0, j, 0, result);

        //i = nLon - 1
        df_p = computeDf(inputGrid->getValue_double(0, j, nLon - 1),
                         inputGrid->getValue_double(1, j, nLon - 1)) * 2.0;
        dp_p = computeDf(inputGrid->getPressure(0, j, nLon - 1),
                         inputGrid->getPressure(1, j, nLon - 1)) * 2.0;
        dfdp = df_p/dp_p;

        // deviation of p in x dir
        dp_x = computeDf(inputGrid->getPressure(0, j, nLon - 2),
                         inputGrid->getPressure(0, j, nLon - 1)) * 2.0;
        dpdx = dp_x/dx;

        // second deviaion of p in x dir
        d2p_x = computeD2f(inputGrid->getPressure(0, j, nLon - 3),
                           inputGrid->getPressure(0, j, nLon - 2),
                           inputGrid->getPressure(0, j, nLon - 1));
        d2pdx2 = d2p_x/dx2;

        // second deviation of f in p dir
        d2f_p = computeD2f(inputGrid->getValue_double(0, j, nLon - 1),
                         inputGrid->getValue_double(1, j, nLon - 1),
                         inputGrid->getValue_double(2, j, nLon - 1));
        dp2 = computeDx2(dp_p);
        d2fdp2 = d2f_p/dp2;

        if (periodicBC[0])
        {
            // deviation of p in x dir
            dp_x = computeDf(inputGrid->getPressure(0, j, 0),
                             inputGrid->getPressure(0, j, nLon - 2));
            dpdx = dp_x/dx;

            // second deviaion of p in x dir
            d2p_x = computeD2f(inputGrid->getPressure(0, j, 0),
                               inputGrid->getPressure(0, j, nLon - 1),
                               inputGrid->getPressure(0, j, nLon - 2));
            d2pdx2 = d2p_x/dx2;
        }
        result = resultGrid->getValue_double(0, j, nLon - 1)
                - d2fdp2 * (dpdx * dpdx)
                - dfdp * d2pdx2
                ;
        resultGrid->setValue_double(0, j, nLon - 1, result);
    }

    // k = nLev - 1
    for (unsigned int j = 0; j < nLat; j++)
    {
        // convert km to m
        //dx = inputGrid->getDeltaLon_km(j) * 1000;
        // do not convert in m
        dx = inputGrid->getDeltaLon_km(j);
        dx2 = computeDx2(dx);
        for (unsigned int i = 1; i < nLon - 1; i++)
        {
            // deviation of f in p dir
            df_p = computeDf(inputGrid->getValue_double(nLev - 2, j, i),
                             inputGrid->getValue_double(nLev - 1, j, i)) * 2;
            dp_p = computeDf(inputGrid->getPressure(nLev - 2, j, i),
                             inputGrid->getPressure(nLev - 1, j, i)) * 2;
            dfdp = df_p/dp_p;

            // deviation of p in x dir
            dp_x = computeDf(inputGrid->getPressure(nLev - 1, j, i - 1),
                             inputGrid->getPressure(nLev - 1, j, i + 1));
            dpdx = dp_x/dx;

            // second deviaion of p in x dir
            d2p_x = computeD2f(inputGrid->getPressure(nLev - 1, j, i - 1),
                               inputGrid->getPressure(nLev - 1, j, i),
                               inputGrid->getPressure(nLev - 1, j, i + 1));
            d2pdx2 = d2p_x/dx2;

            // second deviation of f in p dir
            d2f_p = computeD2f(inputGrid->getValue_double(nLev - 3, j, i),
                             inputGrid->getValue_double(nLev - 2, j, i),
                             inputGrid->getValue_double(nLev - 1, j, i));
            dp2 = computeDx2(dp_p);
            d2fdp2 = d2f_p/dp2;

            result = resultGrid->getValue_double(nLev - 1, j, i)
                    - d2fdp2 * (dpdx * dpdx)
                    - dfdp * d2pdx2
                    ;
            resultGrid->setValue_double(nLev - 1, j, i, result);
        };
        //i = 0
        df_p = computeDf(inputGrid->getValue_double(nLev - 2, j, 0),
                         inputGrid->getValue_double(nLev - 1, j, 0)) * 2.0;
        dp_p = computeDf(inputGrid->getPressure(nLev - 2, j, 0),
                         inputGrid->getPressure(nLev - 1, j, 0)) * 2.0;
        dfdp = df_p/dp_p;

        // deviation of p in x dir
        dp_x = computeDf(inputGrid->getPressure(nLev - 1, j, 0),
                         inputGrid->getPressure(nLev - 1, j, 1)) * 2.0;
        dpdx = dp_x/dx;

        // second deviaion of p in x dir
        d2p_x = computeD2f(inputGrid->getPressure(nLev - 1, j, 0),
                           inputGrid->getPressure(nLev - 1, j, 1),
                           inputGrid->getPressure(nLev - 1, j, 2));
        d2pdx2 = d2p_x/dx2;

        // second deviation of f in p dir
        d2f_p = computeD2f(inputGrid->getValue_double(nLev - 3, j, 0),
                         inputGrid->getValue_double(nLev - 2, j, 0),
                         inputGrid->getValue_double(nLev - 1, j, 0));
        dp2 = computeDx2(dp_p);
        d2fdp2 = d2f_p/dp2;

        if (periodicBC[0])
        {
            // deviation of p in x dir
            dp_x = computeDf(inputGrid->getPressure(nLev - 1, j, nLon -1),
                             inputGrid->getPressure(nLev - 1, j, 1));
            dpdx = dp_x/dx;

            // second deviaion of p in x dir
            d2p_x = computeD2f(inputGrid->getPressure(nLev - 1, j, nLon - 1),
                               inputGrid->getPressure(nLev - 1, j, 0),
                               inputGrid->getPressure(nLev - 1, j, 1));
            d2pdx2 = d2p_x/dx2;
        }
        result = resultGrid->getValue_double(nLev - 1, j, 0)
                - d2fdp2 * (dpdx * dpdx)
                - dfdp * d2pdx2;
        resultGrid->setValue_double(nLev - 1, j, 0, result);

        //i = nLon - 1
        df_p = computeDf(inputGrid->getValue_double(nLev - 2, j, nLon - 1),
                         inputGrid->getValue_double(nLev - 1, j, nLon - 1)) * 2.0;
        dp_p = computeDf(inputGrid->getPressure(nLev - 2, j, nLon - 1),
                         inputGrid->getPressure(nLev - 1, j, nLon - 1)) * 2.0;
        dfdp = df_p/dp_p;

        // deviation of p in x dir
        dp_x = computeDf(inputGrid->getPressure(nLev - 1, j, nLon - 2),
                         inputGrid->getPressure(nLev - 1, j, nLon - 1)) * 2.0;
        dpdx = dp_x/dx;

        // second deviaion of p in x dir
        d2p_x = computeD2f(inputGrid->getPressure(nLev - 1, j, nLon - 3),
                           inputGrid->getPressure(nLev - 1, j, nLon - 2),
                           inputGrid->getPressure(nLev - 1, j, nLon - 1));
        d2pdx2 = d2p_x/dx2;

        // second deviation of f in p dir
        d2f_p = computeD2f(inputGrid->getValue_double(nLev - 3, j, nLon - 1),
                         inputGrid->getValue_double(nLev - 2, j, nLon - 1),
                         inputGrid->getValue_double(nLev - 1, j, nLon - 1));
        dp2 = computeDx2(dp_p);
        d2fdp2 = d2f_p/dp2;

        if (periodicBC[0])
        {
            // deviation of p in x dir
            dp_x = computeDf(inputGrid->getPressure(nLev - 1, j, 0),
                             inputGrid->getPressure(nLev - 1, j, nLon - 2));
            dpdx = dp_x/dx;

            // second deviaion of p in x dir
            d2p_x = computeD2f(inputGrid->getPressure(nLev - 1, j, 0),
                               inputGrid->getPressure(nLev - 1, j, nLon - 1),
                               inputGrid->getPressure(nLev - 1, j, nLon - 2));
            d2pdx2 = d2p_x/dx2;
        }
        result = resultGrid->getValue_double(nLev - 1, j, nLon - 1)
                - d2fdp2 * (dpdx * dpdx)
                - dfdp * d2pdx2;
        resultGrid->setValue_double(nLev - 1, j, nLon - 1, result);
    }
}


void MPartialDerivativeFilter::computeHorizontalPressureLevelTransformationForSecondDerivativeLat(
        MStructuredGrid *inputGrid, MStructuredGrid *resultGrid)
{
    unsigned int nLat = inputGrid->getNumLats();
    unsigned int nLon = inputGrid->getNumLons();
    unsigned int nLev = inputGrid->getNumLevels();
    int iOpp;
    double df_p, dp_y, dp_p, dy;
    double d2f_p, d2p_y, dp2, dy2;
    double dfdp, dpdy;
    double d2fdp2, d2pdy2;
    double result;
    // compute flags for pole and periodic treatment
    QList<bool> periodicBC = periodicBoundaryTreatment(inputGrid);
    // compute horizontal gradient in latitudinal direction.
    // reverse the direction of the gradient
    // convert km to m
    //dy = inputGrid->getDeltaLat_km() * 1000;
    // do not convert in m
    dy = inputGrid->getDeltaLat_km();
    dy2 = computeDx2(dy);

    for (unsigned int k = 1; k < nLev - 1; k++)
    {
        for (unsigned int i = 0; i < nLon; i++)
        {
            for (unsigned int j = 0; j < nLat; j++)
            {
                // deviation of f in p dir
                df_p = computeDf(inputGrid->getValue_double(k-1, j, i),
                                 inputGrid->getValue_double(k+1, j, i));
                dp_p = computeDf(inputGrid->getPressure(k-1, j, i),
                                 inputGrid->getPressure(k+1, j, i));
                dfdp = df_p/dp_p;

                // deviation of p in y dir
                dp_y = computeDf(inputGrid->getPressure(k, j - 1, i),
                                 inputGrid->getPressure(k, j - 1, i));
                dpdy = dp_y/dy;

                // second deviaion of p in y dir
                d2p_y = computeD2f(inputGrid->getPressure(k, j - 1, i),
                                   inputGrid->getPressure(k, j, i),
                                   inputGrid->getPressure(k, j + 1, i));
                d2pdy2 = d2p_y/dy2;

                // second deviation of f in p dir
                d2f_p = computeD2f(inputGrid->getValue_double(k - 1, j, i),
                                 inputGrid->getValue_double(k, j, i),
                                 inputGrid->getValue_double(k + 1, j, i));
                dp2 = computeDx2(dp_p);
                d2fdp2 = d2f_p/dp2;

                result = resultGrid->getValue_double(k, j, i)
                        - d2fdp2 * (dpdy * dpdy)
                        - dfdp * d2pdy2
                        ;
                resultGrid->setValue_double(k, j, i, result);
            };
            //j = 0
            df_p = computeDf(inputGrid->getValue_double(k - 1, 0, i),
                             inputGrid->getValue_double(k + 1, 0, i));
            dp_p = computeDf(inputGrid->getPressure(k - 1, 0, i),
                             inputGrid->getPressure(k + 1, 0, i));
            dfdp = df_p/dp_p;

            // deviation of p in y dir
            dp_y = computeDf(inputGrid->getPressure(k, 0, i),
                             inputGrid->getPressure(k, 1, i)) * 2.0;
            dpdy = dp_y/dy;

            // second deviaion of p in y dir
            d2p_y = computeD2f(inputGrid->getPressure(k, 0, i),
                               inputGrid->getPressure(k, 1, i),
                               inputGrid->getPressure(k, 2, i));
            d2pdy2 = d2p_y/dy2;

            // second deviation of f in p dir
            d2f_p = computeD2f(inputGrid->getValue_double(k - 1, 0, i),
                             inputGrid->getValue_double(k, 0, i),
                             inputGrid->getValue_double(k + 1, 0, i));
            dp2 = computeDx2(dp_p);
            d2fdp2 = d2f_p/dp2;

            if (periodicBC[1])
            {
                iOpp = ((i + (nLon / 2)) % (nLon - 1));
                // deviation of p in y dir
                dp_y = computeDf(inputGrid->getPressure(k, 0, iOpp),
                                 inputGrid->getPressure(k, 1, i));
                dpdy = dp_y/dy;

                // second deviaion of p in y dir
                d2p_y = computeD2f(inputGrid->getPressure(k, 0, iOpp),
                                   inputGrid->getPressure(k, 0, i),
                                   inputGrid->getPressure(k, 1, i));
                d2pdy2 = d2p_y/dy2;
            }
            result = resultGrid->getValue_double(k, 0, i)
                    - d2fdp2 * (dpdy * dpdy)
                    - dfdp * d2pdy2
                    ;
            resultGrid->setValue_double(k, 0, i, result);

            //i = nLon - 1
            df_p = computeDf(inputGrid->getValue_double(k-1, nLat - 1, i),
                             inputGrid->getValue_double(k+1, nLat - 1, i));
            dp_p = computeDf(inputGrid->getPressure(k-1, nLat - 1, i),
                             inputGrid->getPressure(k+1, nLat - 1, i));
            dfdp = df_p/dp_p;

            // deviation of p in y dir
            dp_y = computeDf(inputGrid->getPressure(k, nLat - 2, i),
                             inputGrid->getPressure(k, nLat - 1, i)) * 2.0;
            dpdy = dp_y/dy;

            // second deviaion of p in y dir
            d2p_y = computeD2f(inputGrid->getPressure(k, nLat - 3, i),
                               inputGrid->getPressure(k, nLat - 2, i),
                               inputGrid->getPressure(k, nLat - 1, i));
            d2pdy2 = d2p_y/dy2;

            // second deviation of f in p dir
            d2f_p = computeD2f(inputGrid->getValue_double(k - 1, nLat - 1, i),
                               inputGrid->getValue_double(k, nLat - 1, i),
                               inputGrid->getValue_double(k + 1, nLat - 1, i));
            dp2 = computeDx2(dp_p);
            d2fdp2 = d2f_p/dp2;

            if (periodicBC[2])
            {
                iOpp = ((i + (nLon / 2)) % (nLon - 1));
                // deviation of p in y dir
                dp_y = computeDf(inputGrid->getPressure(k, nLat - 2, i),
                                 inputGrid->getPressure(k, nLat - 1, iOpp));
                dpdy = dp_y/dy;

                // second deviaion of p in y dir
                d2p_y = computeD2f(inputGrid->getPressure(k, nLat - 2, i),
                                   inputGrid->getPressure(k, nLat - 1, i),
                                   inputGrid->getPressure(k, nLat - 1, iOpp));
                d2pdy2 = d2p_y/dy2;
            }
            result = resultGrid->getValue_double(k, nLat - 1, i)
                    - d2fdp2 * (dpdy * dpdy)
                    - dfdp * d2pdy2
                    ;
            resultGrid->setValue_double(k, nLat - 1, i, result);

        }
    }
    // k = 0
    for (unsigned int i = 0; i < nLon; i++)
    {
        for (unsigned int j = 1; j < nLat - 1; j++)
        {
            // deviation of f in p dir
            df_p = computeDf(inputGrid->getValue_double(0, j, i),
                             inputGrid->getValue_double(1, j, i)) * 2;
            dp_p = computeDf(inputGrid->getPressure(0, j, i),
                             inputGrid->getPressure(1, j, i)) * 2;
            dfdp = df_p/dp_p;

            // deviation of p in y dir
            dp_y = computeDf(inputGrid->getPressure(0, j - 1, i),
                             inputGrid->getPressure(0, j + 1, i));
            dpdy = dp_y/dy;

            // second deviaion of p in y dir
            d2p_y = computeD2f(inputGrid->getPressure(0, j - 1, i),
                               inputGrid->getPressure(0, j, i),
                               inputGrid->getPressure(0, j + 1, i));
            d2pdy2 = d2p_y/dy2;

            // second deviation of f in p dir
            d2f_p = computeD2f(inputGrid->getValue_double(0, j, i),
                             inputGrid->getValue_double(1, j, i),
                             inputGrid->getValue_double(2, j, i));
            dp2 = computeDx2(dp_p);
            d2fdp2 = d2f_p/dp2;

            result = resultGrid->getValue_double(0, j, i)
                    - d2fdp2 * (dpdy * dpdy)
                    - dfdp * d2pdy2;
            resultGrid->setValue_double(0, j, i, result);
        };
        //j = 0
        df_p = computeDf(inputGrid->getValue_double(0, 0, i),
                         inputGrid->getValue_double(1, 0, i)) * 2.0;
        dp_p = computeDf(inputGrid->getPressure(0, 0, i),
                         inputGrid->getPressure(1, 0, i)) * 2.0;
        dfdp = df_p/dp_p;

        // deviation of p in y dir
        dp_y = computeDf(inputGrid->getPressure(0, 0, i),
                         inputGrid->getPressure(0, 1, i)) * 2.0;
        dpdy = dp_y/dy;

        // second deviaion of p in y dir
        d2p_y = computeD2f(inputGrid->getPressure(0, 0, i),
                           inputGrid->getPressure(0, 1, i),
                           inputGrid->getPressure(0, 2, i));
        d2pdy2 = d2p_y/dy2;

        // second deviation of f in p dir
        d2f_p = computeD2f(inputGrid->getValue_double(0, 0, i),
                         inputGrid->getValue_double(1, 0, i),
                         inputGrid->getValue_double(2, 0, i));
        dp2 = computeDx2(dp_p);
        d2fdp2 = d2f_p/dp2;

        if (periodicBC[1])
        {
            iOpp = ((i + (nLon / 2)) % (nLon - 1));
            // deviation of p in y dir
            dp_y = computeDf(inputGrid->getPressure(0, 0, iOpp),
                             inputGrid->getPressure(0, 1, i));
            dpdy = dp_y/dy;

            // second deviaion of p in y dir
            d2p_y = computeD2f(inputGrid->getPressure(0, 0, iOpp),
                               inputGrid->getPressure(0, 0, i),
                               inputGrid->getPressure(0, 1, i));
            d2pdy2 = d2p_y/dy2;
        }
        result = resultGrid->getValue_double(0, 0, i)
                - d2fdp2 * (dpdy * dpdy)
                - dfdp * d2pdy2
                ;
        resultGrid->setValue_double(0, 0, i, result);

        //j = nLat - 1
        df_p = computeDf(inputGrid->getValue_double(0, nLat - 1, i),
                         inputGrid->getValue_double(1, nLat - 1, i)) * 2.0;
        dp_p = computeDf(inputGrid->getPressure(0, nLat - 1, i),
                         inputGrid->getPressure(1, nLat - 1, i)) * 2.0;
        dfdp = df_p/dp_p;

        // deviation of p in y dir
        dp_y = computeDf(inputGrid->getPressure(0, nLat - 2, i),
                         inputGrid->getPressure(0, nLat - 1, i)) * 2.0;
        dpdy = dp_y/dy;

        // second deviaion of p in y dir
        d2p_y = computeD2f(inputGrid->getPressure(0, nLat - 3, i),
                           inputGrid->getPressure(0, nLat - 3, i),
                           inputGrid->getPressure(0, nLat - 3, i));
        d2pdy2 = d2p_y/dy2;

        // second deviation of f in p dir
        d2f_p = computeD2f(inputGrid->getValue_double(0, nLat - 1, i),
                         inputGrid->getValue_double(1, nLat - 1, i),
                         inputGrid->getValue_double(2, nLat - 1, i));
        dp2 = computeDx2(dp_p);
        d2fdp2 = d2f_p/dp2;

        if (periodicBC[2])
        {
            iOpp = ((i + (nLon / 2)) % (nLon - 1));
            // deviation of p in y dir
            dp_y = computeDf(inputGrid->getPressure(0, nLat - 2, i),
                             inputGrid->getPressure(0, nLat - 1, iOpp));
            dpdy = dp_y/dy;

            // second deviaion of p in y dir
            d2p_y = computeD2f(inputGrid->getPressure(0, nLat - 2, i),
                               inputGrid->getPressure(0, nLat - 1, i),
                               inputGrid->getPressure(0, nLat - 1, iOpp));
            d2pdy2 = d2p_y/dy2;
        }

        result = resultGrid->getValue_double(0, nLat - 1, i)
                - d2fdp2 * (dpdy * dpdy)
                - dfdp * d2pdy2
                ;
        resultGrid->setValue_double(0, nLat - 1, i, result);
    }

    // k = nLev - 1
    for (unsigned int i = 0; i < nLon; i++)
    {
        for (unsigned int j = 1; j < nLat - 1; j++)
        {
            // deviation of f in p dir
            df_p = computeDf(inputGrid->getValue_double(nLev - 2, j, i),
                             inputGrid->getValue_double(nLev - 1, j, i)) * 2;
            dp_p = computeDf(inputGrid->getPressure(nLev - 2, j, i),
                             inputGrid->getPressure(nLev - 1, j, i)) * 2;
            dfdp = df_p/dp_p;

            // deviation of p in y dir
            dp_y = computeDf(inputGrid->getPressure(nLev - 1, j - 1, i),
                             inputGrid->getPressure(nLev - 1, j + 1, i));
            dpdy = dp_y/dy;

            // second deviaion of p in y dir
            d2p_y = computeD2f(inputGrid->getPressure(nLev - 1, j - 1, i),
                               inputGrid->getPressure(nLev - 1, j, i),
                               inputGrid->getPressure(nLev - 1, j + 1, i));
            d2pdy2 = d2p_y/dy2;

            // second deviation of f in p dir
            d2f_p = computeD2f(inputGrid->getValue_double(nLev - 3, j, i),
                             inputGrid->getValue_double(nLev - 2, j, i),
                             inputGrid->getValue_double(nLev - 1, j, i));
            dp2 = computeDx2(dp_p);
            d2fdp2 = d2f_p/dp2;

            result = resultGrid->getValue_double(nLev - 1, j, i)
                    - d2fdp2 * (dpdy * dpdy)
                    - dfdp * d2pdy2
                    ;
            resultGrid->setValue_double(nLev - 1, j, i, result);
        };
        //j = 0
        df_p = computeDf(inputGrid->getValue_double(nLev - 2, 0, i),
                         inputGrid->getValue_double(nLev - 1, 0, i)) * 2.0;
        dp_p = computeDf(inputGrid->getPressure(nLev - 2, 0, i),
                         inputGrid->getPressure(nLev - 1, 0, i)) * 2.0;
        dfdp = df_p/dp_p;

        // deviation of p in y dir
        dp_y = computeDf(inputGrid->getPressure(nLev - 1, 0, i),
                         inputGrid->getPressure(nLev - 1, 1, i)) * 2.0;
        dpdy = dp_y/dy;

        // second deviaion of p in y dir
        d2p_y = computeD2f(inputGrid->getPressure(nLev - 1, 0, i),
                           inputGrid->getPressure(nLev - 1, 1, i),
                           inputGrid->getPressure(nLev - 1, 2, i));
        d2pdy2 = d2p_y/dy2;

        // second deviation of f in p dir
        d2f_p = computeD2f(inputGrid->getValue_double(nLev - 3, 0, i),
                         inputGrid->getValue_double(nLev - 2, 0, i),
                         inputGrid->getValue_double(nLev - 1, 0, i));
        dp2 = computeDx2(dp_p);
        d2fdp2 = d2f_p/dp2;

        if (periodicBC[1])
        {
            iOpp = ((i + (nLon / 2)) % (nLon - 1));
            // deviation of p in y dir
            dp_y = computeDf(inputGrid->getPressure(nLev - 1, 0, iOpp),
                             inputGrid->getPressure(nLev - 1, 1, i));
            dpdy = dp_y/dy;

            // second deviaion of p in y dir
            d2p_y = computeD2f(inputGrid->getPressure(nLev - 1, 0, iOpp),
                               inputGrid->getPressure(nLev - 1, 0, i),
                               inputGrid->getPressure(nLev - 1, 1, i));
            d2pdy2 = d2p_y/dy2;
        }
        result = resultGrid->getValue_double(nLev - 1, 0, i)
                - d2fdp2 * (dpdy * dpdy)
                - dfdp * d2pdy2
                ;
        resultGrid->setValue_double(nLev - 1, 0, i, result);

        //j = nLat - 1
        df_p = computeDf(inputGrid->getValue_double(nLev - 2, nLat -1, i),
                         inputGrid->getValue_double(nLev - 1, nLat - 1, i)) * 2.0;
        dp_p = computeDf(inputGrid->getPressure(nLev - 2, nLat - 1, i),
                         inputGrid->getPressure(nLev - 1, nLat - 1, i)) * 2.0;
        dfdp = df_p/dp_p;

        // deviation of p in y dir
        dp_y = computeDf(inputGrid->getPressure(nLev - 1, nLat - 2, i),
                         inputGrid->getPressure(nLev - 1, nLat - 1, i)) * 2.0;
        dpdy = dp_y/dy;

        // second deviaion of p in y dir
        d2p_y = computeD2f(inputGrid->getPressure(nLev - 1, nLat - 3, i),
                           inputGrid->getPressure(nLev - 1, nLat - 2, i),
                           inputGrid->getPressure(nLev - 1, nLat - 2, i));
        d2pdy2 = d2p_y/dy2;

        // second deviation of f in p dir
        d2f_p = computeD2f(inputGrid->getValue_double(nLev - 3, nLat - 1, i),
                           inputGrid->getValue_double(nLev - 2, nLat - 1, i),
                           inputGrid->getValue_double(nLev - 1, nLat - 1, i));
        dp2 = computeDx2(dp_p);
        d2fdp2 = d2f_p/dp2;

        if (periodicBC[2])
        {
            iOpp = ((i + (nLon / 2)) % (nLon - 1));
            // deviation of p in y dir
            dp_y = computeDf(inputGrid->getPressure(nLev - 1, nLat - 2, i),
                             inputGrid->getPressure(nLev - 1, nLat - 1, iOpp));
            dpdy = dp_y/dy;

            // second deviaion of p in y dir
            d2p_y = computeD2f(inputGrid->getPressure(nLev - 1, nLat - 2, i),
                               inputGrid->getPressure(nLev - 1, nLat - 1, i),
                               inputGrid->getPressure(nLev - 1, nLat - 1, iOpp));
            d2pdy2 = d2p_y/dy2;
        }

        result = resultGrid->getValue_double(nLev - 1, nLat - 1, i)
                - d2fdp2 * (dpdy * dpdy)
                - dfdp * d2pdy2
                ;
        resultGrid->setValue_double(nLev - 1, nLat - 1, i, result);
    }
}

}  // namespace Met3D
