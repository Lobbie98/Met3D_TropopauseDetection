/******************************************************************************
**
**  This file is part of Met.3D -- a research environment for the
**  three-dimensional visual exploration of numerical ensemble weather
**  prediction data.
**
**  Copyright 2020-2021 Andreas Beckert
**  Copyright 2020-2021 Marc Rautenhaus
**
**  Regional Computing Center, Visual Data Analysis Group
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
#include "smoothfilter.h"

// standard library imports
#include "assert.h"
#include<cmath>
#include <chrono>

// related third party imports
#include <log4cplus/loggingmacros.h>

// local application imports
#include "util/mutil.h"
#include "util/mexception.h"
#include "util/metroutines.h"

#define MEASURE_CPU_TIME

using namespace std;

namespace Met3D
{

/******************************************************************************
***                     CONSTRUCTOR / DESTRUCTOR                            ***
*******************************************************************************/

MSmoothFilter::MSmoothFilter()
    : MSingleInputProcessingWeatherPredictionDataSource()
{

}


/******************************************************************************
***                            PUBLIC METHODS                               ***
*******************************************************************************/

MStructuredGrid *MSmoothFilter::produceData(MDataRequest request)
{
    assert(inputSource != nullptr);

#ifdef MEASURE_CPU_TIME
    auto start = std::chrono::system_clock::now();
#endif

    MDataRequestHelper rh(request);

    // Parse request and initialize result grids.
    // ==========================================
    // "SMOOTH" = Filter type / std deviation (km) / std deviation (grid points)
    // filter type is uniform weights, std deviation is interpreted as radius.
    QString smoothParameter = rh.value("SMOOTH");
    QStringList parameterList = smoothParameter.split("/");
    rh.removeAll(locallyRequiredKeys());
    // The first parameter passes the filter type.
    MSmoothProperties::SmoothModeTypes filterType =
            static_cast<MSmoothProperties::SmoothModeTypes>(
                parameterList[0].toInt());

    MStructuredGrid *result = nullptr;
    //Julian M:tempresult, damit vertikal und horizontal gefiltert werden kann.
    MStructuredGrid *tempresult = nullptr;
    MStructuredGrid* inputGrid = inputSource->getData(rh.request());
    inputGrid->copyFloatDataToDouble();
    result = createAndInitializeResultGrid(inputGrid);
    result->initializeDoubleData();
    //TempResult ist für die zwischenlösung nach dem ExponentialGrid da. Dies sollte definitiv nicht so implementiert werden
    tempresult = createAndInitializeResultGrid(inputGrid);
    tempresult->initializeDoubleData();


    MStructuredGrid *smoothedSfcAuxGrid = nullptr;
    MStructuredGrid *sfcAuxInputGrid = nullptr;

    QString smoothModeName = MSmoothProperties::smoothModeToString(filterType);
    LOG4CPLUS_DEBUG(mlog, "Smooth filter: computing smoothed data fields using "
                    << "method " << smoothModeName.toStdString()
                    << "...");

    // For hybrid-sigma-pressure and aux-pressure grids, we also need to
    // smooth the surface pressure or aux pressure field. Since these can be
    // shared between multiple grid objects, check if the smoothed field
    // already is under memory mangement (i.e., has already been computed
    // previously). If not, compute and store.
    MDataRequest smoothedSfcAuxRequest;
    bool smoothedSfcAuxFieldNeedsToBeComputed = false;
    if ((result->getVerticalLevelType() == HYBRID_SIGMA_PRESSURE_3D)
            || (result->getVerticalLevelType() == AUXILIARY_PRESSURE_3D))
    {
        // Obtain pointer to input sfc/aux grid.
        switch (inputGrid->getVerticalLevelType())
        {
        case HYBRID_SIGMA_PRESSURE_3D:
            sfcAuxInputGrid =
                    dynamic_cast<MLonLatHybridSigmaPressureGrid*>(inputGrid)
                    ->getSurfacePressureGrid(); // does not increase ref count
            break;
        case AUXILIARY_PRESSURE_3D:
            sfcAuxInputGrid =
                    dynamic_cast<MLonLatAuxiliaryPressureGrid*>(inputGrid)
                    ->getAuxiliaryPressureFieldGrid(); // does not incr. ref c.
            break;
        default:
            break;
        }

        // Construct request for smoothed surface pressure field.
        MDataRequestHelper smoothedSfcAuxReqHelp(
                    sfcAuxInputGrid->getGeneratingRequest());
        smoothedSfcAuxReqHelp.insert("SMOOTH", smoothParameter);
        smoothedSfcAuxRequest = smoothedSfcAuxReqHelp.request();

        // Find out whether the sfc pressure field with the required request
        // has already been computed and thus is available in the memory manager,
        // or whether it needs to be computed.
        if (memoryManager->containsData(this, smoothedSfcAuxRequest))
        {
            LOG4CPLUS_DEBUG(mlog, "Smooth filter: required sfc/aux-p field "
                            << "is available in cache.");

            // containsData() increases the item's reference count, hence
            // get the data item (i.e. the sfc/aux pressure field) and
            // exchange it in the result grid.
            if (MLonLatHybridSigmaPressureGrid *hybridResult =
                    dynamic_cast<MLonLatHybridSigmaPressureGrid*>(result))
            {
                hybridResult->exchangeSurfacePressureGrid(
                            static_cast<MRegularLonLatGrid*>(
                                memoryManager->getData(
                                    this, smoothedSfcAuxRequest)));
            }
            else if (MLonLatAuxiliaryPressureGrid *auxPResult =
                     dynamic_cast<MLonLatAuxiliaryPressureGrid*>(result))
            {
                auxPResult->exchangeAuxiliaryPressureGrid(
                            static_cast<MLonLatAuxiliaryPressureGrid*>(
                                memoryManager->getData(
                                    this, smoothedSfcAuxRequest)));
            }
        }
        else
        {
            LOG4CPLUS_DEBUG(mlog, "Smooth filter: required sfc/aux-p field "
                            << "needs to be computed.");

            smoothedSfcAuxFieldNeedsToBeComputed = true;

            // Initialize new smoothed surface pressure or aux pressure field.
            switch (result->getVerticalLevelType())
            {
            case HYBRID_SIGMA_PRESSURE_3D:
                smoothedSfcAuxGrid = createAndInitializeResultGrid(
                            dynamic_cast<MLonLatHybridSigmaPressureGrid*>(result)
                            ->getSurfacePressureGrid());
                break;
            case AUXILIARY_PRESSURE_3D:
                smoothedSfcAuxGrid = createAndInitializeResultGrid(
                            dynamic_cast<MLonLatAuxiliaryPressureGrid*>(result)
                            ->getAuxiliaryPressureFieldGrid());
                // NOTE: As a special case, the aux-p grid references itself
                // as aux-p grid. The createAndInitializeResultGrid() method
                // copies the pointer to the *unsmoothed* aux-p grid. Hence,
                // at the end of this method, the smoothed grid in
                // "smoothedSfcAuxGrid" needs to be fixed to reference itself.
                break;
            default:
                break;
            }

            smoothedSfcAuxGrid->setGeneratingRequest(smoothedSfcAuxRequest);
        }
    }

    //Julian M:Das gleiche, wie darüber, aber für das TempResultgrid. Ansonsten sind die Daten komisch
    MDataRequest tempsmoothedSfcAuxRequest;
    bool tempsmoothedSfcAuxFieldNeedsToBeComputed = false;
    if ((tempresult->getVerticalLevelType() == HYBRID_SIGMA_PRESSURE_3D)
            || (tempresult->getVerticalLevelType() == AUXILIARY_PRESSURE_3D))
    {
        // Obtain pointer to input sfc/aux grid.
        switch (inputGrid->getVerticalLevelType())
        {
        case HYBRID_SIGMA_PRESSURE_3D:
            sfcAuxInputGrid =
                    dynamic_cast<MLonLatHybridSigmaPressureGrid*>(inputGrid)
                    ->getSurfacePressureGrid(); // does not increase ref count
            break;
        case AUXILIARY_PRESSURE_3D:
            sfcAuxInputGrid =
                    dynamic_cast<MLonLatAuxiliaryPressureGrid*>(inputGrid)
                    ->getAuxiliaryPressureFieldGrid(); // does not incr. ref c.
            break;
        default:
            break;
        }

        // Construct request for smoothed surface pressure field.
        MDataRequestHelper smoothedSfcAuxReqHelp(
                    sfcAuxInputGrid->getGeneratingRequest());
        smoothedSfcAuxReqHelp.insert("SMOOTH", smoothParameter);
        tempsmoothedSfcAuxRequest = smoothedSfcAuxReqHelp.request();

        // Find out whether the sfc pressure field with the required request
        // has already been computed and thus is available in the memory manager,
        // or whether it needs to be computed.
        if (memoryManager->containsData(this, tempsmoothedSfcAuxRequest))
        {
            LOG4CPLUS_DEBUG(mlog, "Smooth filter: required sfc/aux-p field "
                            << "is available in cache.");

            // containsData() increases the item's reference count, hence
            // get the data item (i.e. the sfc/aux pressure field) and
            // exchange it in the result grid.
            if (MLonLatHybridSigmaPressureGrid *hybridTempResult =
                    dynamic_cast<MLonLatHybridSigmaPressureGrid*>(tempresult))
            {
                hybridTempResult->exchangeSurfacePressureGrid(
                            static_cast<MRegularLonLatGrid*>(
                                memoryManager->getData(
                                    this, tempsmoothedSfcAuxRequest)));
            }
            else if (MLonLatAuxiliaryPressureGrid *auxPTempResult =
                     dynamic_cast<MLonLatAuxiliaryPressureGrid*>(tempresult))
            {
                auxPTempResult->exchangeAuxiliaryPressureGrid(
                            static_cast<MLonLatAuxiliaryPressureGrid*>(
                                memoryManager->getData(
                                    this, tempsmoothedSfcAuxRequest)));
            }
        }
        else
        {
            LOG4CPLUS_DEBUG(mlog, "Smooth filter: required sfc/aux-p field "
                            << "needs to be computed.");

            tempsmoothedSfcAuxFieldNeedsToBeComputed = true;

            // Initialize new smoothed surface pressure or aux pressure field.
            switch (tempresult->getVerticalLevelType())
            {
            case HYBRID_SIGMA_PRESSURE_3D:
                smoothedSfcAuxGrid = createAndInitializeResultGrid(
                            dynamic_cast<MLonLatHybridSigmaPressureGrid*>(tempresult)
                            ->getSurfacePressureGrid());
                break;
            case AUXILIARY_PRESSURE_3D:
                smoothedSfcAuxGrid = createAndInitializeResultGrid(
                            dynamic_cast<MLonLatAuxiliaryPressureGrid*>(tempresult)
                            ->getAuxiliaryPressureFieldGrid());
                // NOTE: As a special case, the aux-p grid references itself
                // as aux-p grid. The createAndInitializeResultGrid() method
                // copies the pointer to the *unsmoothed* aux-p grid. Hence,
                // at the end of this method, the smoothed grid in
                // "smoothedSfcAuxGrid" needs to be fixed to reference itself.
                break;
            default:
                break;
            }

            smoothedSfcAuxGrid->setGeneratingRequest(tempsmoothedSfcAuxRequest);
        }
    }


    // Compute smoothed data fields.
    // =============================
    switch (filterType)
    {
    // Original Gaussian blur filter with precomputed weights.
    case MSmoothProperties::GAUSS_DISTANCE:
    {
        double stdDev_km = parameterList[1].toFloat();
        computeHorizontalGaussianSmoothing_GCDistance(
                    inputGrid, result, stdDev_km);

        if (smoothedSfcAuxFieldNeedsToBeComputed)
        {
            sfcAuxInputGrid->copyFloatDataToDouble();
            smoothedSfcAuxGrid->initializeDoubleData();
            computeHorizontalGaussianSmoothing_GCDistance(
                        sfcAuxInputGrid, smoothedSfcAuxGrid, stdDev_km);
            smoothedSfcAuxGrid->copyDoubleDataToFloat();
            smoothedSfcAuxGrid->deleteDoubleData();
        }
        break;
    }
    // Box blur filter, where the box size is precalculated according to
    // the distance between grid points. Distance changes between longitudes
    // according to the latitude are considered.
    case MSmoothProperties::BOX_BLUR_DISTANCE_FAST:
    {
        //Julian M: Exponentialfilter anwenden
        computeExponential(inputGrid, tempresult);

        double stdDev_km = parameterList[1].toFloat();
        MSmoothProperties::BoundaryModeTypes boundaryType =
                static_cast<MSmoothProperties::BoundaryModeTypes>(
                    parameterList[3].toInt());

        //Julian M: Wenn die stdDev auf Null ist, dann soll nur vertikal geglättet werden
        if(stdDev_km != 0)
        {
            computeHorizontalBoxBlurSmoothing_GCDistanceFast(
                        tempresult, result, stdDev_km, boundaryType);

            if (smoothedSfcAuxFieldNeedsToBeComputed)
            {
                sfcAuxInputGrid->copyFloatDataToDouble();
                smoothedSfcAuxGrid->initializeDoubleData();
                computeHorizontalBoxBlurSmoothing_GCDistanceFast(
                            sfcAuxInputGrid, smoothedSfcAuxGrid, stdDev_km, boundaryType);
                smoothedSfcAuxGrid->copyDoubleDataToFloat();
                smoothedSfcAuxGrid->deleteDoubleData();
            }
        }
        else
        {
            result = tempresult;
        }

        break;
    }
    // Uniform weights of surrounding grid points.
    case MSmoothProperties::UNIFORM_WEIGHTED_GRIDPOINTS:
    {
        int radius_gp = parameterList[2].toInt();
        computeHorizontalUniformWeightedSmoothing_GCGridpoints_GPU(
                    inputGrid, result, radius_gp);

        if (smoothedSfcAuxFieldNeedsToBeComputed)
        {
            sfcAuxInputGrid->copyFloatDataToDouble();
            smoothedSfcAuxGrid->initializeDoubleData();
            computeHorizontalUniformWeightedSmoothing_GCGridpoints(
                        sfcAuxInputGrid, smoothedSfcAuxGrid, radius_gp);
            smoothedSfcAuxGrid->copyDoubleDataToFloat();
            smoothedSfcAuxGrid->deleteDoubleData();
        }

        break;
    }
    // Original Gaussian blur filter on grid points.
    case MSmoothProperties::GAUSS_GRIDPOINTS:
    {
        int stdDev_gp = parameterList[2].toInt();
        computeHorizontalGaussianSmoothing_GCGridpoints(
                    inputGrid, result, stdDev_gp);

        if (smoothedSfcAuxFieldNeedsToBeComputed)
        {
            sfcAuxInputGrid->copyFloatDataToDouble();
            smoothedSfcAuxGrid->initializeDoubleData();
            computeHorizontalGaussianSmoothing_GCGridpoints(
                        sfcAuxInputGrid, smoothedSfcAuxGrid, stdDev_gp);
            smoothedSfcAuxGrid->copyDoubleDataToFloat();
            smoothedSfcAuxGrid->deleteDoubleData();
        }

        break;
    }
    // Box blur filter. This is the original implementation of the box blur
    // filter but very slow.
    case MSmoothProperties::BOX_BLUR_GRIDPOINTS_SLOW:
    {
        int stdDev_gp = parameterList[2].toInt();
        computeHorizontalBoxBlurSmoothing_GCGridpointsSlow(
                    inputGrid, result, stdDev_gp);

        if (smoothedSfcAuxFieldNeedsToBeComputed)
        {
            sfcAuxInputGrid->copyFloatDataToDouble();
            smoothedSfcAuxGrid->initializeDoubleData();
            computeHorizontalBoxBlurSmoothing_GCGridpointsSlow(
                        sfcAuxInputGrid, smoothedSfcAuxGrid, stdDev_gp);
            smoothedSfcAuxGrid->copyDoubleDataToFloat();
            smoothedSfcAuxGrid->deleteDoubleData();
        }

        break;
    }
    // Fastest box blur filter, with the same result as the slow box blur filter.
    case MSmoothProperties::BOX_BLUR_GRIDPOINTS_FAST:
    {
        int stdDev_gp = parameterList[2].toInt();
        MSmoothProperties::BoundaryModeTypes boundaryType =
                static_cast<MSmoothProperties::BoundaryModeTypes>(
                    parameterList[3].toInt());

        computeHorizontalBoxBlurSmoothing_GCGridpointsFast(
                    inputGrid, result, stdDev_gp, boundaryType);

        if (smoothedSfcAuxFieldNeedsToBeComputed)
        {
            sfcAuxInputGrid->copyFloatDataToDouble();
            smoothedSfcAuxGrid->initializeDoubleData();
            computeHorizontalBoxBlurSmoothing_GCGridpointsFast(
                        sfcAuxInputGrid, smoothedSfcAuxGrid, stdDev_gp, boundaryType);
            smoothedSfcAuxGrid->copyDoubleDataToFloat();
            smoothedSfcAuxGrid->deleteDoubleData();
        }

        break;
    }
    default:
        LOG4CPLUS_ERROR(mlog, "ERROR: Requested smooth method '"
                        << smoothModeName.toStdString()
                        << "' does not exist. Returning nullptr data field.");
    }


    // For hybrid-sigma-pressure and aux-pressure grids, if new surface pressure
    // or aux pressure field has been computed:
    if (smoothedSfcAuxFieldNeedsToBeComputed)
    {
        // Special case (cf comments above where grids are initialized):
        // The 3D pressure field that acts as the aux-p grid references itself.
        // At this point, the smoothed field still references the unsmoothed
        // field. Fix this.
        if (MLonLatAuxiliaryPressureGrid *smoothedAuxP =
                dynamic_cast<MLonLatAuxiliaryPressureGrid*>(smoothedSfcAuxGrid))
        {
            smoothedAuxP->exchangeAuxiliaryPressureGrid(smoothedAuxP);
        }

        // Store sfc/aux grid in memory manager. The call to "storeData()" will
        // place an initial reference of "1" on the item, hence upon success
        // the corresponding fields in the result object can be replaced. In
        // case "storeData()" fails (e.g. in the unlikely event that another
        // thread has stored a field with the same request in the mean time),
        // "storeData()" also increases the reference count, hence the item
        // computed in this thread can be deleted and a reference to the
        // already stored field needs to be obtained.
        if ( !memoryManager->storeData(this, smoothedSfcAuxGrid) )
        {
            delete smoothedSfcAuxGrid;
        }

        // In both cases (storeData() fail and success) the item is now in the
        // memory manager. Hence exchange sfc pressure or aux pressure grid in
        // "result" field.
        if (MLonLatHybridSigmaPressureGrid *hybridResult =
                dynamic_cast<MLonLatHybridSigmaPressureGrid*>(result))
        {
            hybridResult->exchangeSurfacePressureGrid(
                        static_cast<MRegularLonLatGrid*>(
                            memoryManager->getData(
                                this, smoothedSfcAuxRequest)));
        }
        else if (MLonLatAuxiliaryPressureGrid *auxPResult =
                 dynamic_cast<MLonLatAuxiliaryPressureGrid*>(result))
        {
            auxPResult->exchangeAuxiliaryPressureGrid(
                        static_cast<MLonLatAuxiliaryPressureGrid*>(
                            memoryManager->getData(
                                this, smoothedSfcAuxRequest)));
        }
    }

    //LOG4CPLUS_DEBUG(mlog, "Smooth filter: computation finished.");

    result->copyDoubleDataToFloat();

#ifdef MEASURE_CPU_TIME
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    LOG4CPLUS_DEBUG(mlog, "Smooth filter: done in " << elapsed.count() << "ms");
#endif

    // That's it. Release input grid and return result.
    inputSource->releaseData(inputGrid);
    return result;
}


MTask *MSmoothFilter::createTaskGraph(MDataRequest request)
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

const QStringList MSmoothFilter::locallyRequiredKeys()
{
    return (QStringList() << "SMOOTH");
}


/******************************************************************************
***                          PRIVATE METHODS                                ***
*******************************************************************************/

//*************************** GAUSSIAN SMOOTHING *******************************

void MSmoothFilter::computeExponential(
        const MStructuredGrid *inputGrid, MStructuredGrid *resultGrid)
{
    unsigned int nlevs = inputGrid->getNumLevels();
    unsigned int nlats = inputGrid->getNumLats();
    unsigned int nlons = inputGrid->getNumLons();

    //Julian M:Alle Datenpunkte müssen neu berechnet werden
    for (unsigned int j = 0; j < nlats; ++j)
    {
        for (unsigned int i = 0; i < nlons; ++i)
        {
            for (unsigned int k = 0 ; k < nlevs; ++k)
            {
                //ExponentialFilter
                if(k == 0 ||k == 1 || k == nlevs-2 || k == nlevs-1)
                {
                    //Julian M: Die ersten beiden und letzten Werte aus dem Rohgrid kopieren
                    resultGrid->setValue_double( k, j, i, inputGrid->getValue_double(k,j,i) );
                }
                else
                {
                    //Die zwei Werte darüber, zwei Werte darunter und den aktuellen Wert holen
                    int kMinus = k-1;
                    int kPlus = k+1;
                    double vMinusMinus = inputGrid->getValue_double(kMinus-1,j,i);
                    double vMinus = inputGrid->getValue_double(kMinus,j,i);
                    double vakt = inputGrid->getValue_double(k, j, i);
                    double vPlus = inputGrid->getValue_double(kPlus, j, i);
                    double vPlusPlus = inputGrid->getValue_double(kPlus+1, j, i);

                    //Julian M:Gaußsche Verteilung
                    double weighted_val =  0.125f * (vMinusMinus+vPlusPlus) + 0.225f * (vMinus +vPlus) + 0.3f* vakt;
                    resultGrid->setValue_double(k,j,i, weighted_val);
                }
            }
        }
    }

}



// Kernel cannot be precomputed as distances between center point and
// surrounding points change depending on geographical position.
// Gauss is implemented as convolution of latitudinal and longitudinal Gauss
void MSmoothFilter::computeHorizontalGaussianSmoothing_GCDistance(
        const MStructuredGrid *inputGrid, MStructuredGrid *resultGrid,
        double stdDev_km)
{
    MStructuredGrid *resultGridTemp = nullptr;
    resultGridTemp = createAndInitializeResultGrid(resultGrid);
    resultGridTemp->initializeDoubleData();
    //double totalWeight;
    //double totalValue, currentValue;
    //int nWeights;
    const int nLons = inputGrid->getNumLons();
    const int nLats = inputGrid->getNumLats();
    const int nLev = inputGrid->getNumLevels();
    // int iMin, iMax, jMin, jMax;
    const QList<QList<double>> latDependentLonWeights =
            precomputeLatDependendDistanceWeightsOfLongitude(inputGrid,
                                                             stdDev_km);
    const QList<double> weightsLat = precomputeDistanceWeightsOfLatitude(
                inputGrid, stdDev_km);
#pragma omp parallel for
    for (int k = 0; k < nLev; k++)
    {

        // Longitudinal Gauss smoothing.
        for (int j = 0; j < nLats; j++)
        {
            int nWeights = latDependentLonWeights[j].size();
            for (int i = 0; i < nLons; i++)
            {
                double totalValue = 0;
                double totalWeight = 0;
                double currentValue = inputGrid->getValue_double(k, j, i);
                if (IS_MISSING(currentValue))
                {
                    resultGridTemp->setValue_double(k, j, i, M_MISSING_VALUE);
                }
                else
                {
                    int iMin = std::max(i - nWeights + 1, 0);
                    int iMax = std::min(i + nWeights, nLons);
                    for (int m = iMin; m < iMax; m++)
                    {
                        currentValue = inputGrid->getValue_double(k, j, m);
                        if (!IS_MISSING(currentValue))
                        {
                            totalValue += currentValue
                                    * latDependentLonWeights[j][abs(i - m)];
                            totalWeight += latDependentLonWeights[j][abs(i - m)];
                        }
                    }
                    resultGridTemp->setValue_double(k, j, i, totalValue / totalWeight);
                }
            }
        }
        // Latitudinal Gauss smoothing.
        int nWeights = weightsLat.size();
        for (int i = 0; i < nLons; i++)
        {
            for (int j = 0; j < nLats; j++)
            {
                double totalValue = 0;
                double totalWeight = 0;
                double currentValue = resultGridTemp->getValue_double(k, j, i);
                if (IS_MISSING(currentValue))
                {
                    resultGrid->setValue_double(k, j, i, M_MISSING_VALUE);
                }
                else
                {
                    int jMin = std::max(j - nWeights + 1, 0);
                    int jMax = std::min(j + nWeights, nLats);
                    for (int m = jMin; m < jMax; m++)
                    {
                        currentValue = resultGridTemp->getValue_double(k, m, i);
                        if (!IS_MISSING(currentValue))
                        {
                            totalValue += currentValue
                                    * weightsLat[abs(j - m)];
                            totalWeight += weightsLat[abs(j - m)];
                        }
                    }
                    resultGrid->setValue_double(k, j, i, totalValue / totalWeight);
                }
            }
        }
    }

    delete resultGridTemp;
}


void MSmoothFilter::computeHorizontalGaussianSmoothing_GCGridpoints(
        const MStructuredGrid *inputGrid, MStructuredGrid *resultGrid,
        int stdDev_gp)
{
    const unsigned int nlon = inputGrid->getNumLons();
    const unsigned int nlat = inputGrid->getNumLats();
    const double radius = 2. * pow(double(stdDev_gp), 2.);
    //Significant radius, all grid points within the 99% quantile
    //of a Gaussian distribution are considered. The 99% quantile is the
    //result of the standard deviation multiplied by 2.576. For more
    //information see: https://en.wikipedia.org/wiki/Normal_distribution
    const int sigRadius = ceil(double(stdDev_gp) * 2.576);
    int squaredDistance;
    double weight, currentValue, totalValue, addValue, weightSum;
//#pragma omp parallel for
    for (unsigned int k = 0; k < inputGrid->getNumLevels(); ++k)
        {
        for (unsigned int j = 0; j < nlat; ++j)
        {
            for (unsigned int i = 0; i < nlon; ++i)
            {
                currentValue = inputGrid->getValue_double(k, j, i);
                if (IS_MISSING(currentValue))
                {
                    resultGrid->setValue_double(k, j, i, M_MISSING_VALUE);
                }
                else
                {
                    totalValue = 0;
                    weightSum = 0;
                    int south = max(0, static_cast<int>(i) - sigRadius);
                    int north = min(static_cast<int>(nlon), static_cast<int>(i)
                                    + sigRadius);
                    int west = max(0, static_cast<int>(j) - sigRadius);
                    int east = min(static_cast<int>(nlat), static_cast<int>(j)
                                    + sigRadius);
                    for(int n = south; n < north; n++)
                    {
                        for(int m = west; m < east + 1; m++)
                        {
                            addValue = inputGrid->getValue_double(k, m, n);
                            if (!IS_MISSING(addValue))
                            {
                                squaredDistance = (m - j) * (m - j) + (n - i) * (n - i);
                                weight = exp(-squaredDistance / radius)
                                        / (M_PI * radius);
                                totalValue += addValue * weight;
                                weightSum += weight;
                            }
                        }
                    }
                    resultGrid->setValue_double(k, j, i, (totalValue / weightSum));
                }
            }
        }
    }
}


QList<QList<double>> MSmoothFilter::precomputeLatDependendDistanceWeightsOfLongitude(
        const MStructuredGrid *inputGrid, double stdDev_km)
{
    //Significant radius, all grid points within the 99% quantile
    //of a Gaussian distribution are considered. The 99% quantile is the
    //result of the standard deviation multiplied by 2.576. For more
    //information see: https://en.wikipedia.org/wiki/Normal_distribution
    double sigRadius = stdDev_km * 2.576;
    int sigRadius_gridpoints;
    double deltaGridpoint_km;
    double distance_km;
    QList<double> weights;
    QList<QList<double>> latDependendWeights;
    for (unsigned int i = 0; i<inputGrid->getNumLats(); i++)
    {
        //deltaGridpoint_km = inputGrid->getGeometricDeltaLat_km(i);
        deltaGridpoint_km = inputGrid->getDeltaLon_km(i);
        sigRadius_gridpoints = round(sigRadius/deltaGridpoint_km);
        for (int j = 0; j <= sigRadius_gridpoints; j++)
        {
            distance_km = j *  deltaGridpoint_km;
            weights.append(computeGaussWeight(stdDev_km, distance_km));
        }
        latDependendWeights.append((QList<double>) weights);
        weights.clear();
    }
    return latDependendWeights;
}


QList<double> MSmoothFilter::precomputeDistanceWeightsOfLatitude(
        const MStructuredGrid *inputGrid, double stdDev_km)
{
    //Significant radius, all grid points within the 99% quantile
    //of a Gaussian distribution are considered. The 99% quantile is the
    //result of the standard deviation multiplied by 2.576. For more
    //information see: https://en.wikipedia.org/wiki/Normal_distribution
    double significantRadius = stdDev_km * 2.576;
    int significantRadius_gridpoints;
    //double deltaGridpoints_km = inputGrid->getGeometricDeltaLat_km();
    double deltaGridpoints_km = inputGrid->getDeltaLat_km();
    double distance_km;
    QList<double> weights;
    significantRadius_gridpoints =
            round(significantRadius / deltaGridpoints_km);

    for (int j = 0; j <= significantRadius_gridpoints; j++)
    {
        distance_km = j * deltaGridpoints_km;
        weights.append(computeGaussWeight(stdDev_km, distance_km));
    }
    return weights;
}


double MSmoothFilter::computeGaussWeight(double stdDev_km, double distance_km)
{
    double weight = exp(-(pow(distance_km, 2.) / (2. * pow(stdDev_km, 2.))))
            / sqrt(1. / (2. * M_PI * pow(stdDev_km, 2.)));
    return weight;
}


//*************************** BOX BLUR SMOOTHING *******************************

// Longitudinal Box Blur smoothing.
void MSmoothFilter::computeHorizontalBoxBlurSmoothing_GCDistanceFast(
        const MStructuredGrid *inputGrid, MStructuredGrid *resultGrid,
        double stdDev_km, MSmoothProperties::BoundaryModeTypes boundaryType)
{

    QList<QList<int>> latDependentBoxRadii;
    int n = 3;
    latDependentBoxRadii = computeLatDependentBoxRadii(inputGrid, stdDev_km, n);
    // Distance between latitudes.
    const double deltaGp_km = inputGrid->getDeltaLat_km();
    int distanceInGridpoints = static_cast<int>(round(stdDev_km / deltaGp_km));
    QVector<int> lonBoxRadii = computeBoxRadii(distanceInGridpoints, n);
    MStructuredGrid *resultGridTemp = nullptr;
    resultGridTemp = createAndInitializeResultGrid(resultGrid);
    resultGridTemp->initializeDoubleData();
    boxBlurTotalFast(inputGrid, resultGrid, lonBoxRadii[0],
            latDependentBoxRadii[0], boundaryType);
    boxBlurTotalFast(resultGrid, resultGridTemp, lonBoxRadii[1],
            latDependentBoxRadii[1], boundaryType);
    boxBlurTotalFast(resultGridTemp, resultGrid, lonBoxRadii[2],
            latDependentBoxRadii[2], boundaryType);
    delete resultGridTemp;
}


void MSmoothFilter::computeHorizontalBoxBlurSmoothing_GCGridpointsFast(
        const MStructuredGrid *inputGrid, MStructuredGrid *resultGrid,
        int stdDev_gp, MSmoothProperties::BoundaryModeTypes boundaryType)
{
    QVector<int> boxRadii = computeBoxRadii(stdDev_gp, 3);
    MStructuredGrid *resultGridTemp = nullptr;
    resultGridTemp = createAndInitializeResultGrid(resultGrid);
    resultGridTemp->initializeDoubleData();
    boxBlurTotalFast(inputGrid, resultGrid, boxRadii[0], boundaryType);
    boxBlurTotalFast(resultGrid, resultGridTemp, boxRadii[1], boundaryType);
    boxBlurTotalFast(resultGridTemp, resultGrid, boxRadii[2], boundaryType);
    delete resultGridTemp;
}


void MSmoothFilter::computeHorizontalBoxBlurSmoothing_GCGridpointsSlow(
        const MStructuredGrid *inputGrid, MStructuredGrid *resultGrid,
        int stdDev_gp)
{
    QVector<int> boxRadii = computeBoxRadii(stdDev_gp, 3);
    MStructuredGrid *resultGridTemp = nullptr;
    resultGridTemp = createAndInitializeResultGrid(resultGrid);
    resultGridTemp->initializeDoubleData();
    boxBlurTotalSlow(inputGrid, resultGrid, boxRadii[0]);
    boxBlurTotalSlow(resultGrid, resultGridTemp, boxRadii[1]);
    boxBlurTotalSlow(resultGridTemp, resultGrid, boxRadii[2]);
    delete resultGridTemp;
}


QList<QList<int>> MSmoothFilter::computeLatDependentBoxRadii(
        const MStructuredGrid *inputGrid, double stdDev_km,
        int n)
{
    int stdDev_gp;
    QList<QList<int>> boxes;
    QList<int> t;
    for (int i = 0; i < n; i++)
    {
        boxes.append((QList<int>) t);
    }
    for (unsigned int i = 0; i < inputGrid->getNumLats(); i++)
    {
        stdDev_gp = numGridpointsSpannedByDistance(inputGrid, i, stdDev_km);
        // Ideal averaging filter width.
        double widthIdeal = sqrt((12. * pow(stdDev_gp, 2) / n) + 1.);
        // Filter width rounded to nearest odd integer less than widthIdeal.
        int widthIdealLess = floor(widthIdeal);
        // Check if widthIdealLess is odd.
        if ((widthIdealLess % 2) == 0)
        {
            widthIdealLess--;
        }
        // Nearest odd integer higher than ideal width.
        int widthIdealUp = widthIdealLess + 2;
        // mIdeal determines at which pass the nearest odd integer width higher
        // than the width ideal should be used. This should compensate the
        // rounding in the first place to the nearest integer less than
        // the ideal width.
        double mIdeal = (12. * pow(stdDev_gp, 2.) - n *
                        pow(widthIdealLess, 2.)
                        - 4. * double(n) * double(widthIdealLess) - 3. * n)
                /(-4. * double(widthIdealLess) - 4.);
        // m has to be rounded to nearest integer.
        int m = round(mIdeal);
        for (int j = 0; j < n; j++)
        {
            if (j < m)
            {
                boxes[j].append((widthIdealLess - 1) / 2);
            }
            else
            {
                boxes[j].append((widthIdealUp - 1) / 2);
            }
        }
    }
    return boxes;
}


int MSmoothFilter::numGridpointsSpannedByDistance(
        const MStructuredGrid *inputGrid, int iLat, double distance_km)
{
    double distanceInDeg;
    int distanceInGridpoints;
    double phi = abs(inputGrid->getLats()[iLat]) * M_PI / 180.;
    double latitudeCircleInKm = cos(phi) * 2. * M_PI
            * MetConstants::EARTH_RADIUS_km;
    // This should prevent an division by zero when phi is 90°.
    // TODO(ab, 13Feb2020) --: Rethink this error query, probably there are
    // smarter ways.
    // distanceInDegree should not be larger than 360°.
    if (latitudeCircleInKm > 0)
    {
        distanceInDeg = distance_km/latitudeCircleInKm * 360.;
        if (distanceInDeg > 360.) {distanceInDeg = 360.;}
    }
    else
    {
        distanceInDeg = 360.;
    }
    distanceInGridpoints = round(distanceInDeg / inputGrid->getDeltaLon());
    return distanceInGridpoints;
}


QVector<int> MSmoothFilter::computeBoxRadii(int stdDev_gp, int n)
{
    QVector<int> boxRadii(n);
    // Ideal averaging filter width.
    double widthIdeal = sqrt((12. * pow(stdDev_gp, 2.) / double(n)) + 1.);
    // Filter width rounded to nearest odd integer less than widthIdeal.
    int widthIdealLess = floor(widthIdeal);
    if ((widthIdealLess % 2) == 0)
    {
        widthIdealLess--;
    }
    // Check if widthIdealLess is odd.
    int widthIdealUp = widthIdealLess + 2;
    // mIdeal determines at which pass the nearest odd integer width is higher
    // than the ideal width should be used. This should compensate the
    // rounding in the first place to the nearest integer less than
    // the ideal width.
    double mIdeal = (12. * pow(stdDev_gp, 2.) - n
                    * pow(double(widthIdealLess), 2.) - 4.
                    * double(n) * widthIdealLess - 3. * double(n))
            / (-4. * double(widthIdealLess) - 4.);
    int m = round(mIdeal);
    for (int i = 0; i < n; i++)
    {
        if (i < m)
        {
            boxRadii[i] = (widthIdealLess - 1) / 2;
        }
        else
        {
            boxRadii[i] = (widthIdealUp - 1) / 2;
        }
    }
    return boxRadii;
}


void MSmoothFilter::boxBlurTotalFast(
        const MStructuredGrid *inputGrid, MStructuredGrid *resultGrid,
        int lonBoxRadius, QList<int> latDependentBoxRadii,
        MSmoothProperties::BoundaryModeTypes boundaryType)
{
    MStructuredGrid *resultGridTemp = nullptr;
    resultGridTemp = createAndInitializeResultGrid(resultGrid);
    resultGridTemp->initializeDoubleData();
    boxBlurLongitudinalFast(inputGrid, resultGridTemp, latDependentBoxRadii,
                            boundaryType);
    boxBlurLatitudinalFast(resultGridTemp, resultGrid, lonBoxRadius,
                           boundaryType);
    delete resultGridTemp;
}


void MSmoothFilter::boxBlurTotalFast(
        const MStructuredGrid *inputGrid, MStructuredGrid *resultGrid,
        int boxRadius, MSmoothProperties::BoundaryModeTypes boundaryType)
{
    MStructuredGrid *resultGridTemp = nullptr;
    resultGridTemp = createAndInitializeResultGrid(resultGrid);
    resultGridTemp->initializeDoubleData();
    boxBlurLongitudinalFast(inputGrid, resultGridTemp, boxRadius, boundaryType);
    boxBlurLatitudinalFast(resultGridTemp, resultGrid, boxRadius, boundaryType);
    delete resultGridTemp;
}


void MSmoothFilter::boxBlurLongitudinalFast(
        const MStructuredGrid *inputGrid, MStructuredGrid *resultGrid,
        const QList<int> latDependentBoxRadii,
        const MSmoothProperties::BoundaryModeTypes boundaryType)
{

    const int nLons = inputGrid->getNumLons();

    if (boundaryType == MSmoothProperties::BoundaryModeTypes::CONSTANT)
    {
#pragma omp parallel for
        for (unsigned int k = 0; k < inputGrid->getNumLevels(); k++)
        {
            double value, plusValue, minusValue;
            int boxRadius, nGridPoints, iMinus, iPlus;
            for (unsigned int j = 0; j < inputGrid->getNumLats(); j++)
            {
                boxRadius = latDependentBoxRadii[j];
                if IS_MISSING(inputGrid->getValue_double(k, j, 0))
                {
                    nGridPoints = 0;
                    value = 0.;
                }
                else
                {
                    value = inputGrid->getValue_double(k, j, 0) * double(boxRadius + 1);
                    nGridPoints = boxRadius + 1;
                }
                // add values until box Radius is reached
                for(int i = 1; i < boxRadius + 1; i++)
                {
                    plusValue = inputGrid->getValue_double(k, j, i);
                    if IS_MISSING(plusValue)
                    {
                        plusValue = 0.;
                    }
                    else
                    {
                        nGridPoints ++;
                    }
                    value += plusValue;
                }
                // set the first value
                if IS_MISSING(inputGrid->getValue_double(k, j, 0))
                {
                    resultGrid->setValue_double(k, j, 0, M_MISSING_VALUE);
                }
                else
                {
                    resultGrid->setValue_double(k, j, 0, value/double(nGridPoints));
                }

                // Get and set all other values starting from index i = 1.
                for (int i = 1; i < nLons; i++)
                {
                    iMinus = std::max(0, i - boxRadius - 1);
                    minusValue = inputGrid->getValue_double(k, j, iMinus);
                    if IS_MISSING(minusValue)
                    {
                        minusValue = 0.;
                        nGridPoints ++;
                    }
                    value -= minusValue;

                    iPlus = std::min(nLons - 1, i + boxRadius);
                    plusValue = inputGrid->getValue_double(k, j, iPlus);
                    if IS_MISSING(plusValue)
                    {
                        nGridPoints --;
                        plusValue = 0.;
                    }
                    value += plusValue;
                    if IS_MISSING(inputGrid->getValue_double(k, j, i))
                    {
                        resultGrid->setValue_double(k, j, i, M_MISSING_VALUE);
                    }
                    else
                    {
                        resultGrid->setValue_double(k, j, i, value/double(nGridPoints));
                    }
                }
            }
        }
    }
    else if (boundaryType == MSmoothProperties::BoundaryModeTypes::SYMMETRIC)
    {
#pragma omp parallel for
        for (unsigned int k = 0; k < inputGrid->getNumLevels(); k++)
        {
            double value, plusValue, minusValue;
            int boxRadius, nGridPoints, iMinus, iPlus;
            for (unsigned int j = 0; j < inputGrid->getNumLats(); j++)
            {
                boxRadius = latDependentBoxRadii[j];
                const QVector<int> indexList = createIndexList(inputGrid,
                                                               boxRadius,
                                                         QString("LON"));
                value = 0;
                nGridPoints = 2 * boxRadius + 1;
                // add values until box 2 * Radius + 1 is reached
                for (int i = 0; i < (2 * boxRadius + 1); i++)
                {
                    plusValue = inputGrid->getValue_double(k, j, indexList[i]);
                    if IS_MISSING(plusValue)
                    {
                        plusValue = 0.;
                        nGridPoints --;
                    }
                    value += plusValue;
                }
                if IS_MISSING(inputGrid->getValue_double(k, j, 0))
                {
                     resultGrid->setValue_double(k, j, 0, M_MISSING_VALUE);
                }
                else
                {
                    resultGrid->setValue_double(k, j, 0, value/double(nGridPoints));
                }
                for (int i = 1; i < nLons; i++)
                {
                    iMinus = i - 1;
                    minusValue = inputGrid->getValue_double(k, j, indexList[iMinus]);
                    if IS_MISSING(minusValue)
                    {
                        minusValue = 0.;
                        nGridPoints ++;
                    }
                    value -= minusValue;

                    iPlus = i + 2 * boxRadius;
                    plusValue = inputGrid->getValue_double(k, j, indexList[iPlus]);
                    if IS_MISSING(plusValue)
                    {
                        plusValue = 0.;
                        nGridPoints --;
                    }
                    value += plusValue;

                    if IS_MISSING(inputGrid->getValue_double(k, j, i))
                    {
                         resultGrid->setValue_double(k, j, i, M_MISSING_VALUE);
                    }
                    else
                    {
                        resultGrid->setValue_double(k, j, i, value/double(nGridPoints));
                    }
                }
            }
        }
    }
    else if (boundaryType == MSmoothProperties::BoundaryModeTypes::NANPADDING)
    {
#pragma omp parallel for
        for (unsigned int k = 0; k < inputGrid->getNumLevels(); k++)
        {
            double value, plusValue, minusValue;
            int boxRadius, nGridPoints;
            double iarr, currentValue;
            for (unsigned int j = 0; j < inputGrid->getNumLats(); j++)
            {
                boxRadius = latDependentBoxRadii[j];
                value = 0.;
                nGridPoints = 0;
                // Adds the values until radius is reached.
                // This is then the first value to be set on the result grid.
                for (int i = 0; i < boxRadius; i++)
                {
                    plusValue = inputGrid->getValue_double(k, j, i);
                    if (IS_MISSING(plusValue))
                    {
                        plusValue = 0.;
                        nGridPoints--;
                    }
                    value += plusValue;
                    nGridPoints++;
                }
                for(int i = 0; i < boxRadius + 1; i++)
                {
                    plusValue = inputGrid->getValue_double(k, j, i + boxRadius);
                    currentValue = inputGrid->getValue_double(k, j, i);
                    if (IS_MISSING(plusValue))
                    {
                        plusValue = 0.;
                        nGridPoints--;
                    }
                    value += plusValue;
                    nGridPoints++;
                    if (IS_MISSING(currentValue))
                    {
                        resultGrid->setValue_double(k, j, i, M_MISSING_VALUE);
                    }
                    else
                    {
                        iarr = 1. / double(nGridPoints);
                        resultGrid->setValue_double(k, j, i, value * iarr);
                    }
                }

                for(int i = boxRadius + 1; i < (nLons - boxRadius); i++)
                {
                    plusValue = inputGrid->getValue_double(k, j, i + boxRadius);
                    minusValue = inputGrid->getValue_double(k, j, i - boxRadius - 1);
                    currentValue = inputGrid->getValue_double(k, j, i);
                    if (IS_MISSING(plusValue))
                    {
                        plusValue = 0.;
                        nGridPoints--;
                    }
                    if (IS_MISSING(minusValue))
                    {
                        minusValue = 0.;
                        nGridPoints++;
                    }
                    value += plusValue - minusValue;
                    if (IS_MISSING(currentValue))
                    {
                        resultGrid->setValue_double(k, j, i, M_MISSING_VALUE);
                    }
                    else
                    {
                        iarr = 1. / double(nGridPoints);
                        resultGrid->setValue_double(k, j, i, value * iarr);
                    }
                }

                for(int i = (nLons - boxRadius); i < nLons; i++)
                {
                    minusValue = inputGrid->getValue_double(k, j, i - boxRadius - 1);
                    currentValue = inputGrid->getValue_double(k, j, i);
                    if (IS_MISSING(minusValue))
                    {
                        minusValue = 0.;
                        nGridPoints ++;
                    }
                    value -= minusValue;
                    nGridPoints--;
                    if (IS_MISSING(currentValue))
                    {
                        resultGrid->setValue_double(k, j, i, M_MISSING_VALUE);
                    }
                    else
                    {
                        iarr = 1. / double(nGridPoints);
                        resultGrid->setValue_double(k, j, i, value * iarr);
                    }
                }
            }
        }
    }
}


void MSmoothFilter::boxBlurLongitudinalFast(
        const MStructuredGrid *inputGrid, MStructuredGrid *resultGrid,
        int boxRadius, MSmoothProperties::BoundaryModeTypes boundaryType)
{
    const int nLons = inputGrid->getNumLons();

    if (boundaryType == MSmoothProperties::BoundaryModeTypes::CONSTANT)
    {
#pragma omp parallel for
        for (unsigned int k = 0; k < inputGrid->getNumLevels(); k++)
        {
            double value, plusValue, minusValue;
            int nGridPoints, iMinus, iPlus;
            for (unsigned int j = 0; j < inputGrid->getNumLats(); j++)
            {
                if IS_MISSING(inputGrid->getValue_double(k, j, 0))
                {
                    nGridPoints = 0;
                    value = 0.;
                }
                else
                {
                    value = inputGrid->getValue_double(k, j, 0) * double(boxRadius + 1);
                    nGridPoints = boxRadius + 1;
                }
                // add values until box Radius is reached
                for(int i = 1; i < boxRadius + 1; i++)
                {
                    plusValue = inputGrid->getValue_double(k, j, i);
                    if IS_MISSING(plusValue)
                    {
                        plusValue = 0.;
                    }
                    else
                    {
                        nGridPoints ++;
                    }
                    value += plusValue;
                }
                // set the first value
                if IS_MISSING(inputGrid->getValue_double(k, j, 0))
                {
                    resultGrid->setValue_double(k, j, 0, M_MISSING_VALUE);
                }
                else
                {
                    resultGrid->setValue_double(k, j, 0, value/double(nGridPoints));
                }

                // Get and set all other values starting from index i = 1.
                for (int i = 1; i < nLons; i++)
                {
                    iMinus = std::max(0, i - boxRadius - 1);
                    minusValue = inputGrid->getValue_double(k, j, iMinus);
                    if IS_MISSING(minusValue)
                    {
                        minusValue = 0.;
                        nGridPoints ++;
                    }
                    value -= minusValue;

                    iPlus = std::min(nLons - 1, i + boxRadius);
                    plusValue = inputGrid->getValue_double(k, j, iPlus);
                    if IS_MISSING(plusValue)
                    {
                        nGridPoints --;
                        plusValue = 0.;
                    }
                    value += plusValue;
                    if IS_MISSING(inputGrid->getValue_double(k, j, i))
                    {
                        resultGrid->setValue_double(k, j, i, M_MISSING_VALUE);
                    }
                    else
                    {
                        resultGrid->setValue_double(k, j, i, value/double(nGridPoints));
                    }
                }
            }
        }
    }
    else if (boundaryType == MSmoothProperties::BoundaryModeTypes::SYMMETRIC)
    {
#pragma omp parallel for
        for (unsigned int k = 0; k < inputGrid->getNumLevels(); k++)
        {
            double value, plusValue, minusValue;
            int nGridPoints, iMinus, iPlus;
            for (unsigned int j = 0; j < inputGrid->getNumLats(); j++)
            {
                const QVector<int> indexList = createIndexList(inputGrid, boxRadius,
                                                           QString("LON"));
                value = 0;
                nGridPoints = 2 * boxRadius + 1;
                // add values until box 2 * Radius + 1 is reached
                for (int i = 0; i < (2 * boxRadius + 1); i++)
                {
                    plusValue = inputGrid->getValue_double(k, j, indexList[i]);
                    if IS_MISSING(plusValue)
                    {
                        plusValue = 0.;
                        nGridPoints --;
                    }
                    value += plusValue;
                }
                if IS_MISSING(inputGrid->getValue_double(k, j, 0))
                {
                     resultGrid->setValue_double(k, j, 0, M_MISSING_VALUE);
                }
                else
                {
                    resultGrid->setValue_double(k, j, 0, value/double(nGridPoints));
                }
                for (int i = 1; i < nLons; i++)
                {
                    iMinus = i - 1;
                    minusValue = inputGrid->getValue_double(k, j, indexList[iMinus]);
                    if IS_MISSING(minusValue)
                    {
                        minusValue = 0.;
                        nGridPoints ++;
                    }
                    value -= minusValue;

                    iPlus = i + 2 * boxRadius;
                    plusValue = inputGrid->getValue_double(k, j, indexList[iPlus]);
                    if IS_MISSING(plusValue)
                    {
                        plusValue = 0.;
                        nGridPoints --;
                    }
                    value += plusValue;

                    if IS_MISSING(inputGrid->getValue_double(k, j, i))
                    {
                         resultGrid->setValue_double(k, j, i, M_MISSING_VALUE);
                    }
                    else
                    {
                        resultGrid->setValue_double(k, j, i, value/double(nGridPoints));
                    }
                }
            }
        }
    }
    else if (boundaryType == MSmoothProperties::BoundaryModeTypes::NANPADDING)
    {
#pragma omp parallel for
        for (unsigned int k = 0; k < inputGrid->getNumLevels(); k++)
        {
            double iarr, currentValue;
            double value, plusValue, minusValue;
            int nGridPoints;
            for (unsigned int j = 0; j < inputGrid->getNumLats(); j++)
            {
                value = 0.;
                nGridPoints = 0;
                // Adds the values until radius is reached.
                // This is then the first value to be set on the result grid.
                for (int i = 0; i < boxRadius; i++)
                {
                    plusValue = inputGrid->getValue_double(k, j, i);
                    if (IS_MISSING(plusValue))
                    {
                        plusValue = 0.;
                        nGridPoints--;
                    }
                    value += plusValue;
                    nGridPoints++;
                }
                for(int i = 0; i < boxRadius + 1; i++)
                {
                    plusValue = inputGrid->getValue_double(k, j, i + boxRadius);
                    currentValue = inputGrid->getValue_double(k, j, i);
                    if (IS_MISSING(plusValue))
                    {
                        plusValue = 0.;
                        nGridPoints--;
                    }
                    value += plusValue;
                    nGridPoints++;
                    if (IS_MISSING(currentValue))
                    {
                        resultGrid->setValue_double(k, j, i, M_MISSING_VALUE);
                    }
                    else
                    {
                        iarr = 1. / double(nGridPoints);
                        resultGrid->setValue_double(k, j, i, value * iarr);
                    }
                }

                for(int i = boxRadius + 1; i < (nLons - boxRadius); i++)
                {
                    plusValue = inputGrid->getValue_double(k, j, i + boxRadius);
                    minusValue = inputGrid->getValue_double(k, j, i - boxRadius - 1);
                    currentValue = inputGrid->getValue_double(k, j, i);
                    if (IS_MISSING(plusValue))
                    {
                        plusValue = 0.;
                        nGridPoints--;
                    }
                    if (IS_MISSING(minusValue))
                    {
                        minusValue = 0.;
                        nGridPoints++;
                    }
                    value += plusValue - minusValue;
                    if (IS_MISSING(currentValue))
                    {
                        resultGrid->setValue_double(k, j, i, M_MISSING_VALUE);
                    }
                    else
                    {
                        iarr = 1. / double(nGridPoints);
                        resultGrid->setValue_double(k, j, i, value * iarr);
                    }
                }

                for(int i = (nLons - boxRadius); i < nLons; i++)
                {
                    minusValue = inputGrid->getValue_double(k, j, i - boxRadius - 1);
                    currentValue = inputGrid->getValue_double(k, j, i);
                    if (IS_MISSING(minusValue))
                    {
                        minusValue = 0.;
                        nGridPoints ++;
                    }
                    value -= minusValue;
                    nGridPoints--;
                    if (IS_MISSING(currentValue))
                    {
                        resultGrid->setValue_double(k, j, i, M_MISSING_VALUE);
                    }
                    else
                    {
                        iarr = 1. / double(nGridPoints);
                        resultGrid->setValue_double(k, j, i, value * iarr);
                    }
                }
            }
        }
    }
}


void MSmoothFilter::boxBlurLatitudinalFast(
        const MStructuredGrid *inputGrid, MStructuredGrid *resultGrid,
        int boxRadius, MSmoothProperties::BoundaryModeTypes boundaryType)
{
    const int nLats = inputGrid->getNumLats();

    if (boundaryType == MSmoothProperties::BoundaryModeTypes::CONSTANT)
    {
#pragma omp parallel for
        for (unsigned int k = 0; k < inputGrid->getNumLevels(); k++)
        {
            double value, plusValue, minusValue;
            int nGridPoints, jMinus, jPlus;
            for (unsigned int i = 0; i < inputGrid->getNumLons(); i++)
            {
                if IS_MISSING(inputGrid->getValue_double(k, 0, i))
                {
                    nGridPoints = 0;
                    value = 0.;
                }
                else
                {
                    value = inputGrid->getValue_double(k, 0, i) * double(boxRadius + 1);
                    nGridPoints = boxRadius + 1;
                }
                // add values until box Radius is reached
                for(int j = 1; j < boxRadius + 1; j++)
                {
                    plusValue = inputGrid->getValue_double(k, j, i);
                    if IS_MISSING(plusValue)
                    {
                        plusValue = 0.;
                    }
                    else
                    {
                        nGridPoints ++;
                    }
                    value += plusValue;
                }
                // set the first value
                if IS_MISSING(inputGrid->getValue_double(k, 0, i))
                {
                    resultGrid->setValue_double(k, 0, i, M_MISSING_VALUE);
                }
                else
                {
                    resultGrid->setValue_double(k, 0, i, value/double(nGridPoints));
                }

                // Get and set all other values starting from index j = 1.
                for (int j = 1; j < nLats; j++)
                {
                    jMinus = std::max(0, j - boxRadius - 1);
                    minusValue = inputGrid->getValue_double(k, jMinus, i);
                    if IS_MISSING(minusValue)
                    {
                        minusValue = 0.;
                        nGridPoints ++;
                    }
                    value -= minusValue;

                    jPlus = std::min(nLats - 1, j + boxRadius);
                    plusValue = inputGrid->getValue_double(k, jPlus, i);
                    if IS_MISSING(plusValue)
                    {
                        nGridPoints --;
                        plusValue = 0.;
                    }
                    value += plusValue;
                    if IS_MISSING(inputGrid->getValue_double(k, j, i))
                    {
                        resultGrid->setValue_double(k, j, i, M_MISSING_VALUE);
                    }
                    else
                    {
                        resultGrid->setValue_double(k, j, i, value/double(nGridPoints));
                    }
                }
            }
        }
    }
    else if (boundaryType == MSmoothProperties::BoundaryModeTypes::SYMMETRIC)
    {
#pragma omp parallel for
        for (unsigned int k = 0; k < inputGrid->getNumLevels(); k++)
        {
            double value, plusValue, minusValue;
            int nGridPoints, jMinus, jPlus;
            for (unsigned int i = 0; i < inputGrid->getNumLons(); i++)
            {
                const QVector<int> indexList = createIndexList(inputGrid, boxRadius,
                                                           QString("LAT"));
                value = 0;
                nGridPoints = 2 * boxRadius + 1;
                // add values until box Radius + 1 is reached
                for (int j = 0; j < (2 * boxRadius + 1); j++)
                {
                    plusValue = inputGrid->getValue_double(k, indexList[j], i);
                    if IS_MISSING(plusValue)
                    {
                        plusValue = 0.;
                        nGridPoints --;
                    }
                    value += plusValue;
                }
                if IS_MISSING(inputGrid->getValue_double(k, 0, i))
                {
                     resultGrid->setValue_double(k, 0, i, M_MISSING_VALUE);
                }
                else
                {
                    resultGrid->setValue_double(k, 0, i, value/double(nGridPoints));
                }
                for (int j = 1; j < nLats; j++)
                {
                    jMinus = j - 1;
                    minusValue = inputGrid->getValue_double(k, indexList[jMinus], i);
                    if IS_MISSING(minusValue)
                    {
                        minusValue = 0.;
                        nGridPoints ++;
                    }
                    value -= minusValue;

                    jPlus = j + 2 * boxRadius;
                    plusValue = inputGrid->getValue_double(k, indexList[jPlus], i);
                    if IS_MISSING(plusValue)
                    {
                        plusValue = 0.;
                        nGridPoints --;
                    }
                    value += plusValue;

                    if IS_MISSING(inputGrid->getValue_double(k, j, i))
                    {
                         resultGrid->setValue_double(k, j, i, M_MISSING_VALUE);
                    }
                    else
                    {
                        resultGrid->setValue_double(k, j, i, value/double(nGridPoints));
                    }
                }
            }
        }
    }
    else if (boundaryType == MSmoothProperties::BoundaryModeTypes::NANPADDING)
    {
#pragma omp parallel for
        for (unsigned int k = 0; k < inputGrid->getNumLevels(); k++)
        {
            double iarr, currentValue;
            double value, plusValue, minusValue;
            int nGridPoints;
            for (unsigned int i = 0; i < inputGrid->getNumLons(); i++)
            {
                value = 0.;
                nGridPoints = 0;
                // Adds the values until radius is reached.
                // This is then the first value to be set on the result grid.
                for (int j = 0; j < boxRadius; j++)
                {
                    plusValue = inputGrid->getValue_double(k, j, i);
                    if (IS_MISSING(plusValue))
                    {
                        plusValue = 0.;
                        nGridPoints--;
                    }
                    value += plusValue;
                    nGridPoints++;
                }
                for(int j = 0; j < boxRadius + 1; j++)
                {
                    plusValue = inputGrid->getValue_double(k, j + boxRadius, i);
                    currentValue = inputGrid->getValue_double(k, j, i);
                    if (IS_MISSING(plusValue))
                    {
                        plusValue = 0.;
                        nGridPoints--;
                    }
                    value += plusValue;
                    nGridPoints++;
                    if (IS_MISSING(currentValue))
                    {
                        resultGrid->setValue_double(k, j, i, M_MISSING_VALUE);
                    }
                    else
                    {
                        iarr = 1. / double(nGridPoints);
                        resultGrid->setValue_double(k, j, i, value * iarr);
                    }
                }

                for(int j = boxRadius + 1; j < (nLats - boxRadius); j++)
                {
                    plusValue = inputGrid->getValue_double(k, j + boxRadius, i);
                    minusValue = inputGrid->getValue_double(k, j - boxRadius - 1, i);
                    currentValue = inputGrid->getValue_double(k, j, i);
                    if (IS_MISSING(plusValue))
                    {
                        plusValue = 0.;
                        nGridPoints--;
                    }
                    if (IS_MISSING(minusValue))
                    {
                        minusValue = 0.;
                        nGridPoints++;
                    }
                    value += plusValue - minusValue;
                    if (IS_MISSING(currentValue))
                    {
                        resultGrid->setValue_double(k, j, i, M_MISSING_VALUE);
                    }
                    else
                    {
                        iarr = 1. / static_cast<double>(nGridPoints);
                        resultGrid->setValue_double(k, j, i, value * iarr);
                    }
                }

                for(int j = (nLats - boxRadius); j < nLats; j++)
                {
                    minusValue = inputGrid->getValue_double(k, j - boxRadius - 1, i);
                    currentValue = inputGrid->getValue_double(k, j, i);
                    if (IS_MISSING(minusValue))
                    {
                        minusValue = 0.;
                        nGridPoints ++;
                    }
                    value -= minusValue;
                    nGridPoints--;
                    if (IS_MISSING(currentValue))
                    {
                        resultGrid->setValue_double(k, j, i, M_MISSING_VALUE);
                    }
                    else
                    {
                        iarr = 1. / double(nGridPoints);
                        resultGrid->setValue_double(k, j, i, value * iarr);
                    }
                }
            }
        }
    }
}


QVector<int> MSmoothFilter::createIndexList(const MStructuredGrid *inputGrid,
                                            const int boxRadius, const QString dir)
{
    int size;
    if (dir == "LON")
    {
        size = inputGrid->getNumLons();
    }
    else if (dir == "LAT")
    {
        size = inputGrid->getNumLats();
    }
    else
    {
        size = 0;
    }
    int n = boxRadius - 1;
    QVector<int> indexListStart(boxRadius);
    generate(indexListStart.begin(), indexListStart.end(), [&n]{ return n--;});

    n = 0;
    QVector<int> indexListMiddle(size);
    generate(indexListMiddle.begin(), indexListMiddle.end(), [&n]{ return n++;});


    n = size - 1;
    QVector<int> indexListEnd(boxRadius);
    generate(indexListEnd.begin(), indexListEnd.end(), [&n]{ return n--;});

    QVector<int> indexList;
    indexList.append(indexListStart);
    indexList.append(indexListMiddle);
    indexList.append(indexListEnd);

    return indexList;
}


void MSmoothFilter::boxBlurTotalSlow(const MStructuredGrid *inputGrid,
                                     MStructuredGrid *resultGrid, int boxRadius)
{
    // No missing value handling!
    // Same result as boxBlurTotalFast with constant boundaries.
    int nLat = inputGrid->getNumLats();
    int nLon = inputGrid->getNumLons();
//#pragma omp parallel for
    for (unsigned int k = 0; k < inputGrid->getNumLevels(); k++)
    {
        for (int j = 0; j < nLat; j++)
        {
            for (int i = 0; i < nLon; i++)
            {
                double value = 0;
                for(int iy = j - boxRadius; iy <= j + boxRadius; iy++)
                {
                    for(int ix = i - boxRadius; ix <= i + boxRadius; ix++)
                    {
                        int x = min(nLon - 1, max(0, ix));
                        int y = min(nLat - 1, max(0, iy));
                        value += inputGrid->getValue_double(k, y, x);
                    }
                }
                resultGrid->setValue_double(k, j, i,
                                     (value / ((2 * double(boxRadius) + 1)
                                               * (2 * double(boxRadius) + 1))));
            }
        }
    }
}


//*************************** UNIFORM WEIGHTED SMOOTHING ***********************

void MSmoothFilter::computeHorizontalUniformWeightedSmoothing_GCGridpoints(
        const MStructuredGrid *inputGrid, MStructuredGrid *resultGrid,
        int radius_gp)
{
    int nLon = inputGrid->getNumLons();
    int nlat = inputGrid->getNumLats();
    int nlev = inputGrid->getNumLevels();
#pragma omp parallel for
    for (int k = 0; k < nlev; ++k)
    {
        for (int j = 0; j < nlat; ++j)
        {
            for (int i = 0; i < nLon; ++i)
            {
                // Implementation of the simplest smoothing filter.
                // Take the sdtDev_gp parameter as grid distance and smooth
                // without accounting any weights.

                int iStart = max(0, static_cast<int>(i - radius_gp));
                int iEnd = min(static_cast<int>(nLon),
                               static_cast<int>(i + radius_gp));
                int jStart = max(0, static_cast<int>(j - radius_gp));
                int jEnd = min(static_cast<int>(nlat),
                               static_cast<int>(j + radius_gp));

                double totalValue = 0;
                int nSmoothPoints = 0;
                double currentValue = inputGrid->getValue_double(k, j, i);
                if (IS_MISSING(currentValue))
                {
                    resultGrid->setValue_double(k, j, i,
                                         M_MISSING_VALUE);
                }
                else
                {
                    for (int js = jStart; js < jEnd; js++)
                    {
                        for (int is = iStart; is < iEnd; is++)
                        {
                            double addValue = inputGrid->getValue_double(k, js, is);
                            if (!IS_MISSING(addValue))
                            {
                                totalValue += addValue;
                                nSmoothPoints ++;
                            }
                        }
                    }
                    resultGrid->setValue_double(k, j, i,
                                         totalValue/double(nSmoothPoints));
                 }
            }
        }
    }
}


//*************************** UNIFORM WEIGHTED SMOOTHING GPU TRY****************
//

void MSmoothFilter::computeHorizontalUniformWeightedSmoothing_GCGridpoints_GPU(
        const MStructuredGrid *inputGrid, MStructuredGrid *resultGrid,
        int radius_gp)
{
    int nlon = inputGrid->getNumLons();
    int nlat = inputGrid->getNumLats();
    unsigned int nlev = inputGrid->getNumLevels();
    int nlatlons = nlat * nlon;

    const int numValues = inputGrid->getNumValues();
    const double *data_double = inputGrid->getData_double();
    double *data_double_res = new double[nlev * nlat * nlon];

// #pragma omp target data map(to:radius_gp) map(to:nlon) map(to:nlat) map(to:nlev) map(to:data_double[0:numValues])  map(tofrom:data_double_res[0:numValues])
// #pragma omp target teams distribute parallel for
// #pragma omp parallel for
    for (unsigned int k = 0; k < nlev; ++k)
    {
        for (int j = 0; j < nlat; ++j)
        {
            for (int i = 0; i < nlon; ++i)
            {
                int index = k * nlatlons + j * nlon + i;

                int iStart = ((i - radius_gp) > 0) ? (i - radius_gp) : 0;
                int iEnd = ((i + radius_gp) < nlon) ? (i + radius_gp) : nlon;
                int jStart = ((j - radius_gp) > 0) ? (j - radius_gp) : 0;
                int jEnd = ((j + radius_gp) < nlat) ? (j + radius_gp) : nlat;

                double totalValue = 0;
                int nSmoothPoints = 0;
                double currentValue = data_double[((k)*(nlatlons) + (j)*(nlon) + (i))];
                if (IS_MISSING(currentValue))
                {
                    data_double_res[index] = M_MISSING_VALUE;
                }
                else
                {
                    for (int js = jStart; js < jEnd; js++)
                    {
                        for (int is = iStart; is < iEnd; is++)
                        {
                            double addValue = data_double[((k * nlatlons) + (js * nlon) + is)];
                            if (!IS_MISSING(addValue))
                            {
                                totalValue += addValue;
                                nSmoothPoints ++;
                            }
                        }
                    }
                    data_double_res[index] =
                                         totalValue/double(nSmoothPoints);
                }
            }
        }
    }
    for (int n = 0; n < numValues; ++n)
    {
        resultGrid->setValue_double(n, data_double_res[n]);
    }
    delete [] data_double_res;
}

}  // namespace Met3D
