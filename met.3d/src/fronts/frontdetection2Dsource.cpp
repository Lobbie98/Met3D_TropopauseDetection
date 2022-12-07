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

#include "frontdetection2Dsource.h"
#include "fronts/frontdetection3Dsource.h"
#include "util/mmarchingcubes.h"

// standard library imports

// related third party imports
#include <omp.h>

// local application imports
#include <log4cplus/loggingmacros.h>

#define COMPUTE_PARALLEL
#define MEASURE_CPU_TIME
#define qNaN (std::numeric_limits<float>::quiet_NaN())

using namespace Met3D;

MFrontDetection2DSource::MFrontDetection2DSource()
        : M2DFrontFilterSource()
{

}


void MFrontDetection2DSource::setFLESource(MFrontLocationEquationSource *s)
{
    fleSource = s;
    registerInputSource(fleSource);
    enablePassThrough(fleSource);
}


void MFrontDetection2DSource::setTFPSource(MThermalFrontParameterSource *s)
{
    tfpSource = s;
    registerInputSource(tfpSource);
    enablePassThrough(tfpSource);
}


void MFrontDetection2DSource::setABZSource(MAdjacentBaroclinicZoneSource *s)
{
    abzSource = s;
    registerInputSource(abzSource);
    enablePassThrough(abzSource);
}


void MFrontDetection2DSource::setDetectionVariableSource(MWeatherPredictionDataSource *s)
{
    detectionVariableSource = s;
    registerInputSource(detectionVariableSource);
    enablePassThrough(detectionVariableSource);
}


void MFrontDetection2DSource::setDetectionVariablePartialDeriveSource(
        MPartialDerivativeFilter *s)
{
    detectionVarPartialDerivativeSource = s;
    registerInputSource(detectionVarPartialDerivativeSource);
    enablePassThrough(detectionVarPartialDerivativeSource);
}


void MFrontDetection2DSource::setWindUSource(MWeatherPredictionDataSource *s)
{
    windUSource = s;
    registerInputSource(windUSource);
    enablePassThrough(windUSource);
}


void MFrontDetection2DSource::setWindVSource(MWeatherPredictionDataSource *s)
{
    windVSource = s;
    registerInputSource(windVSource);
    enablePassThrough(windVSource);
}


M2DFrontSelection* MFrontDetection2DSource::produceData(MDataRequest request)
{
    assert(fleSource != nullptr);
    assert(tfpSource != nullptr);
    assert(abzSource != nullptr);
    assert(detectionVariableSource != nullptr);
    assert(detectionVarPartialDerivativeSource != nullptr);
    assert(windUSource != nullptr);
    assert(windVSource != nullptr);

    MDataRequestHelper rh(request);

    const double isovalue = rh.value("FRONTS_ISOVALUE").toFloat();
    const QString windUVar = rh.value("FRONTS_WINDU_VAR");
    const QString windVVar = rh.value("FRONTS_WINDV_VAR");
    const float overtracing = rh.value("NC_OVER_TRACING").toFloat();
    const float elevation_hPa = rh.value("FRONTS_ELEVATION").toFloat();

    rh.removeAll(locallyRequiredKeys());

    MStructuredGrid* fleGrid = fleSource->getData(rh.request());
    MStructuredGrid* tfpGrid = tfpSource->getData(rh.request());
    MStructuredGrid* abzGrid = abzSource->getData(rh.request());

    MStructuredGrid* detectionVarGrid = detectionVariableSource->getData(rh.request());
    const int DLON = MGradientProperties::DLON;
    const int DLAT = MGradientProperties::DLAT;
    rh.insert("GRADIENT", DLON);
    MStructuredGrid* dDetectionVarDXGrid =
            detectionVarPartialDerivativeSource->getData(rh.request());
    rh.insert("GRADIENT", DLAT);
    MStructuredGrid* dDetectionVarDYGrid =
            detectionVarPartialDerivativeSource->getData(rh.request());
    rh.remove("GRADIENT");

    rh.insert("VARIABLE", windUVar);
    MStructuredGrid* windUGrid = windUSource->getData(rh.request());
    rh.insert("VARIABLE", windVVar);
    MStructuredGrid* windVGrid = windVSource->getData(rh.request());

    LOG4CPLUS_DEBUG(mlog, "[1] Start marching squares to get front line candidates");
#ifdef MEASURE_CPU_TIME
    auto start1 = std::chrono::system_clock::now();
#endif

    QVector<QVector3D>* positions = getLinesFromEnsembleMember(
                fleGrid, isovalue, elevation_hPa);
    LOG4CPLUS_DEBUG(mlog, "[1] \t->done.");

#ifdef MEASURE_CPU_TIME
    auto end1 = std::chrono::system_clock::now();
    auto elapsed1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1);
    LOG4CPLUS_DEBUG(mlog, "[1] \t->done in " << elapsed1.count() << "ms");
#endif


    MLineSelection *rawFrontLines = new MLineSelection(
                 positions->size());

    MNormalCurvesSelection *rawNormalCurves = new MNormalCurvesSelection(
                positions->size());


    const float minValue = fleGrid->min();


    LOG4CPLUS_DEBUG(mlog, "[2] Start integration along normal curves");
#ifdef MEASURE_CPU_TIME
    auto start2 = std::chrono::system_clock::now();
#endif

//#ifdef COMPUTE_PARALLEL
//#pragma omp parallel for
//#endif
    for (int k = 0; k < positions->size(); ++k)
    {
        // get current position
        QVector3D position = positions->at(k);

        // Check if position is valid
        if (position.x() == qNaN) { continue; }

        // initialize new front line new vertex
        Geometry::FrontLineVertex frontVertex;

        // get all relevant values at the position
        const float tfp = tfpGrid->interpolateValue(position);
        const float abz = abzGrid->interpolateValue(position);

        // compute front type (warm or cold)
        QVector2D gradTheta(dDetectionVarDXGrid->interpolateValue(position),
                            dDetectionVarDYGrid->interpolateValue(position));
        gradTheta.normalize();
        QVector2D wind(windUGrid->interpolateValue(position),
                       windVGrid->interpolateValue(position));
        wind.normalize();
        const float type = float(QVector2D::dotProduct(-gradTheta, wind) >= 0);

        // compute normal curves
        Geometry::NormalCurve nc = MFrontDetection3DSource::integrateAlongThermalGradient(
                    fleGrid, dDetectionVarDXGrid, dDetectionVarDYGrid,
                    position, tfp, minValue, overtracing);

        // compute frontal strength:
        QVector3D ncEnd = nc.positions.last();
        float strength = detectionVarGrid->interpolateValue(nc.positions.first()) -
                detectionVarGrid->interpolateValue(ncEnd);

        // all needed values are computed, fill front vertex
        frontVertex.position = position;
        frontVertex.nCEnd = nc.positions.last();
        frontVertex.tfp = tfp;
        frontVertex.abz = abz;
        frontVertex.strength = strength;
        frontVertex.type = type;
        frontVertex.breadth = nc.breadth;

        // all needed values are computed, fill normal curve vertices
        nc.tfp = tfp;
        nc.abz = abz;
        nc.strength = strength;
        nc.type = type;

        // set vertices to 2D fronts values
        rawFrontLines->setVertex(k, frontVertex);
        rawNormalCurves->setNormalCurve(k, nc);

    }

    LOG4CPLUS_DEBUG(mlog, "[2] \t->done.");
#ifdef MEASURE_CPU_TIME
    auto end2 = std::chrono::system_clock::now();
    auto elapsed2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2);
    LOG4CPLUS_DEBUG(mlog, "[2] \t->done in " << elapsed2.count() << "ms");
#endif


    QVector<u_int32_t> indexArray = sortFrontLineVertices(rawFrontLines,
                                                          positions->size());
    // set
    rawFrontLines->setIndexArray(indexArray);

    fleSource->releaseData(fleGrid);
    tfpSource->releaseData(tfpGrid);
    abzSource->releaseData(abzGrid);
    detectionVariableSource->releaseData(detectionVarGrid);
    detectionVarPartialDerivativeSource->releaseData(dDetectionVarDXGrid);
    detectionVarPartialDerivativeSource->releaseData(dDetectionVarDYGrid);
    windUSource->releaseData(windUGrid);
    windVSource->releaseData(windVGrid);

    M2DFrontSelection *raw2DFronts = new M2DFrontSelection;
    raw2DFronts->setLineSelection(rawFrontLines);
    raw2DFronts->setNormalCurvesSelection(rawNormalCurves);
    return raw2DFronts;
}


MTask* MFrontDetection2DSource::createTaskGraph(MDataRequest request)
{
    assert(fleSource != nullptr);
    assert(tfpSource != nullptr);
    assert(abzSource != nullptr);
    assert(detectionVariableSource != nullptr);
    assert(detectionVarPartialDerivativeSource != nullptr);
    assert(windUSource != nullptr);
    assert(windVSource != nullptr);

    MTask* task =  new MTask(request, this);
    MDataRequestHelper rh(request);

    const QString windUVar = rh.value("FRONTS_WINDU_VAR");
    const QString windVVar = rh.value("FRONTS_WINDV_VAR");

    rh.removeAll(locallyRequiredKeys());

    task->addParent(fleSource->getTaskGraph(rh.request()));
    task->addParent(tfpSource->getTaskGraph(rh.request()));
    task->addParent(abzSource->getTaskGraph(rh.request()));

    task->addParent(detectionVariableSource->getTaskGraph(rh.request()));
    const int DLON = MGradientProperties::DLON;
    const int DLAT = MGradientProperties::DLAT;
    rh.insert("GRADIENT", DLON);
    task->addParent(detectionVarPartialDerivativeSource->getTaskGraph(rh.request()));
    rh.insert("GRADIENT", DLAT);
    task->addParent(detectionVarPartialDerivativeSource->getTaskGraph(rh.request()));
    rh.remove("GRADIENT");

    rh.insert("VARIABLE", windUVar);
    task->addParent(windUSource->getTaskGraph(rh.request()));
    rh.insert("VARIABLE", windVVar);
    task->addParent(windVSource->getTaskGraph(rh.request()));

    return task;
}


const QStringList MFrontDetection2DSource::locallyRequiredKeys()
{
    return (QStringList() << "FRONTS_ISOVALUE"
                          << "FRONTS_WINDU_VAR"
                          << "FRONTS_WINDV_VAR"
                          << "FRONTS_NC_OVER_TRACING"
                          << "FRONTS_ELEVATION"
                          << "FRONTS_NC_INTEGRATION_GPU");
}


QVector<QVector3D>* MFrontDetection2DSource::getLinesFromEnsembleMember(
        MStructuredGrid* grid, const float isovalue, const float elevation_hPa)
{
    QVector<QVector3D>* lines = new  QVector<QVector3D>;

    //for (auto k = 31; k < 32; ++k)
    const int maxK = (elevation_hPa >= 0) ? 1 : grid->getNumLevels();

    //for (auto k = 58; k < 59; ++k)
    for (int k = 0; k < maxK; ++k)
    {
        for (unsigned int j = 0; j < grid->getNumLats() - 1; ++j)
        {
            float lat0 = grid->getLats()[j];
            float lat1 = grid->getLats()[j + 1];

            for (unsigned int i = 0; i < grid->getNumLons() - 1; ++i)
            {
                float lon0 = grid->getLons()[i];
                float lon1 = grid->getLons()[i + 1];

                QVector3D worldPos00(lon0, lat0, (elevation_hPa >= 0)
                                     ? elevation_hPa : grid->getPressure(k, j, i));
                QVector3D worldPos10(lon1, lat0, (elevation_hPa >= 0)
                                     ? elevation_hPa : grid->getPressure(k, j, i + 1));
                QVector3D worldPos01(lon0, lat1, (elevation_hPa >= 0)
                                     ? elevation_hPa : grid->getPressure(k, j + 1, i));
                QVector3D worldPos11(lon1, lat1, (elevation_hPa >= 0)
                                     ? elevation_hPa : grid->getPressure(k, j + 1, i + 1));

                const float value00 = grid->interpolateValue(worldPos00);
                const float value10 = grid->interpolateValue(worldPos10);
                const float value01 = grid->interpolateValue(worldPos01);
                const float value11 = grid->interpolateValue(worldPos11);

                uint8_t bits = 0;
                bits |= int(value00 > isovalue) * 0x1;
                bits |= int(value10 > isovalue) * 0x2;
                bits |= int(value01 > isovalue) * 0x4;
                bits |= int(value11 > isovalue) * 0x8;

                if (bits == 0) { continue; }
                if (bits > 7) { bits = 15 - bits; }

                switch (bits)
                {
                    case 1:
                        lines->push_back(getLineVertex(worldPos00, worldPos10,
                                                      value00, value10, isovalue));
                        lines->push_back(getLineVertex(worldPos00, worldPos01,
                                                      value00, value01, isovalue));
                        break;
                    case 2:
                        lines->push_back(getLineVertex(worldPos00, worldPos10,
                                                      value00, value10, isovalue));
                        lines->push_back(getLineVertex(worldPos10, worldPos11,
                                                      value10, value11, isovalue));
                        break;

                    case 3:
                        lines->push_back(getLineVertex(worldPos00, worldPos01,
                                                      value00, value01, isovalue));
                        lines->push_back(getLineVertex(worldPos10, worldPos11,
                                                      value10, value11, isovalue));
                        break;
                    case 4:
                        lines->push_back(getLineVertex(worldPos00, worldPos01,
                                                      value00, value01, isovalue));
                        lines->push_back(getLineVertex(worldPos01, worldPos11,
                                                      value01, value11, isovalue));
                        break;
                    case 5:
                        lines->push_back(getLineVertex(worldPos00, worldPos10,
                                                      value00, value10, isovalue));
                        lines->push_back(getLineVertex(worldPos01, worldPos11,
                                                      value01, value11, isovalue));
                        break;
                    case 6:
                    {
                        // TODO: handle ambiguous cases
                        float midValue =
                                (value00 + value10 + value01 + value11) / 4.0f;
                        bool sign00 = value00 > isovalue;
                        bool sign11 = value11 > isovalue;
                        bool signMid = midValue > isovalue;


                        if (signMid == sign00 && signMid == sign11)
                        {
                            lines->push_back(
                                    getLineVertex(worldPos00, worldPos10,
                                                  value00, value10, isovalue));
                            lines->push_back(
                                    getLineVertex(worldPos10, worldPos11,
                                                  value10, value11, isovalue));
                            lines->push_back(
                                    getLineVertex(worldPos00, worldPos01,
                                                  value00, value01, isovalue));
                            lines->push_back(
                                    getLineVertex(worldPos01, worldPos11,
                                                  value01, value11, isovalue));
                        }
                        else
                        {
                            lines->push_back(
                                    getLineVertex(worldPos00, worldPos10,
                                                  value00, value10, isovalue));
                            lines->push_back(
                                    getLineVertex(worldPos00, worldPos01,
                                                  value00, value01, isovalue));
                            lines->push_back(
                                    getLineVertex(worldPos10, worldPos11,
                                                  value10, value11, isovalue));
                            lines->push_back(
                                    getLineVertex(worldPos01, worldPos11,
                                                  value01, value11, isovalue));
                        }

                        break;
                    }
                    case 7:
                        lines->push_back(getLineVertex(worldPos10, worldPos11,
                                                      value10, value11, isovalue));
                        lines->push_back(getLineVertex(worldPos01, worldPos11,
                                                      value01, value11, isovalue));
                        break;
                    default:
                        break;
                }
            }
        }
    }

    return lines;
}


QVector3D MFrontDetection2DSource::getLineVertex(
        const QVector3D& p0, const QVector3D& p1,
        const float v0, const float v1, const float isovalue)
{
    const float t = (isovalue - v0) / (v1 - v0);

    if (t < 0 || t > 1) {
        std::cout << "Error?" << std::endl; }

    return p1 * t + p0 * (1.0f - t);
}


QVector<u_int32_t> MFrontDetection2DSource::sortFrontLineVertices(MLineSelection* frontLineVertices,
                                                                  int numVertices)
{
    LOG4CPLUS_DEBUG(mlog, "[3] Start sorting front line candidates and normal curves");
#ifdef MEASURE_CPU_TIME
    auto start3 = std::chrono::system_clock::now();
#endif

    // each line segment is defined by two vertices. The two corresponding
    // vertices are always ordered next to each other in the frontLineVertices array.
    // We know that the 1st and 2nd vertices form one line segment. The
    // 3rd and 4th vertices the next line segment, and so on...
    // add first line segment to empty index array as starting point
    QVector<QVector<u_int32_t>> indexArray;
    QVector<u_int32_t> firstIndexLineStrip;

    firstIndexLineStrip.push_back(0);
    firstIndexLineStrip.push_back(1);
    indexArray.append(firstIndexLineStrip);
    // start with 2nd line segment (vertex 3 and 4), since the first line segment is already
    // added to the index array.
    for (int v = 2; v + 1 < numVertices; v += 2)
    {
        int indexP1 = v;
        int indexP2 = v + 1;
        QVector3D p1(frontLineVertices->getVertex(indexP1).position);
        QVector3D p2(frontLineVertices->getVertex(indexP2).position);
        bool addedLineToPolygon = false;

        QVector<int> addedToPolygonNr;

        // Loop over all polygons available for this contour level. Start with the
        // last polygon, as this one is the most likely match for the new line
        // segment.
        for (int i = indexArray.size() - 1; i >= 0; i--)
        {
            // Get first and last point of polygon.
            QVector3D first = frontLineVertices->getVertex(indexArray.at(i).first()).position;
            QVector3D last = frontLineVertices->getVertex(indexArray.at(i).last()).position;

            // Does the new line connect to one end of the polygon?
            // Yes: Append line to polygon.
            if (p1 == first) {
                if (!addedLineToPolygon)
                {
                    indexArray[i].push_front(indexP2);
                    addedLineToPolygon = true;
                }
                addedToPolygonNr.append(i);
            }
            else if (p2 == first) {
                if (!addedLineToPolygon){
                    indexArray[i].push_front(indexP1);
                    addedLineToPolygon = true;
                }
                addedToPolygonNr.append(i);
            }
            else if (p1 == last) {
                if (!addedLineToPolygon){
                    indexArray[i].push_back(indexP2);
                    addedLineToPolygon = true;
                }
                addedToPolygonNr.append(i);
            }
            else if (p2 == last) {
                if (!addedLineToPolygon){
                    indexArray[i].push_back(indexP1);
                    addedLineToPolygon = true;
                }
                addedToPolygonNr.append(i);
            }
        }

        // If line couldn't be attached to any polygon create new polygon.
        if (!addedLineToPolygon)
        {
            QVector<u_int32_t> newIndexLineStrip;

            newIndexLineStrip.push_back(indexP1);
            newIndexLineStrip.push_back(indexP2);
            indexArray.append(newIndexLineStrip);
        }
        if(addedToPolygonNr.size() > 1)
        {
            if (addedToPolygonNr.size() > 2)
            {
                // print error message saying that there are more that two line strips
                // to join, will take the first and second one
            }
            if (addedToPolygonNr[0] != addedToPolygonNr[1])
            {
                indexArray = tryToJoinLineStrip(frontLineVertices,
                                                indexArray,
                                                addedToPolygonNr[0],
                                                addedToPolygonNr[1]);
            }
        }
        //addLine(frontLineVertices, indexArray, v, v+1);
    }

    // change sorted indices to a list of indices with a restart index
    // check that all lines have at least 10 elements. This is the minimum
    // length for a clean rendering of lines
    u_int32_t restartIndex = std::numeric_limits<u_int32_t>::max();
    QVector<u_int32_t> sortedIndeces;
    for (int i = 0; i < indexArray.size(); i++)
    {
        if (indexArray[i].size() > 9)
        {
            for (int j = 0; j < indexArray[i].size(); j++)
            {
                sortedIndeces.append(indexArray[i][j]);
            }
            sortedIndeces.append(restartIndex);

        }
    }

    LOG4CPLUS_DEBUG(mlog, "[3] \t->done.");
#ifdef MEASURE_CPU_TIME
    auto end3 = std::chrono::system_clock::now();
    auto elapsed3 = std::chrono::duration_cast<std::chrono::milliseconds>(end3 - start3);
    LOG4CPLUS_DEBUG(mlog, "[3] \t->done in " << elapsed3.count() << "ms");
#endif

    return sortedIndeces;
}


QVector<QVector<u_int32_t>> MFrontDetection2DSource::tryToJoinLineStrip(
        MLineSelection* frontLineVertices,
        QVector<QVector<u_int32_t>> indexArray, int i, int j)
{
    // Test if polygons #i and #k of contour level k can be merged (i.e. if one
    // of their end points is identical). If yes, merge them.
    QVector3D line1First =  frontLineVertices->getVertex(indexArray.at(i).first()).position;
    QVector3D line1Last =  frontLineVertices->getVertex(indexArray.at(i).last()).position;
    QVector3D line2First =  frontLineVertices->getVertex(indexArray.at(j).first()).position;
    QVector3D line2Last =  frontLineVertices->getVertex(indexArray.at(j).last()).position;
    if (line1Last == line2First)
    {
        // Last point of polygon i equals first point of polygon j. Append all
        // but the first point of polygon j to polygon i. Delete polygon j.
        indexArray[j].removeFirst();
        indexArray[i].append(indexArray[j]);
        indexArray.remove(j);
        return indexArray;
    }
    else if (line1First == line2Last) {
        // Last point of polygon j equals first point of polygon i. Append all
        // but the first point of polygon i to polygon j. Delete polygon i.
        indexArray[i].removeFirst();
        indexArray[j].append(indexArray[i]);
        indexArray.remove(i);
        return indexArray;
    }
    else if (line1First == line2First)
    {
        // reverse indexArray[i]
        QVector<u_int32_t> indexReversed(indexArray[i].size());
        for (int r = 0; r < indexArray[i].size(); r++)
        {
            indexReversed[r] = indexArray[i][indexArray[i].size() - 1 - r];
        }
        // First points of both polygons equal. Insert elements of polygon j
        // in the front of polygon i. Delete polygon j.
        indexArray[i] = indexReversed;
        indexArray[j].removeFirst();
        indexArray[i].append(indexArray[j]);
        indexArray.remove(j);
        return indexArray;
    }
    else if (line1Last == line2Last)
    {
        // reverse indexArray[j]
        QVector<u_int32_t> indexReversed(indexArray[j].size());
        for (int r = 0; r < indexArray[j].size(); r++)
        {
            indexReversed[r] = indexArray[j][indexArray[j].size() - 1 - r];
        }
        // Last points of both polygons are equal. Append polygon j in
        // reversed order (except for the last (after reverse first) element) to polygon i.
        // Delete polygon j.
        indexReversed.removeFirst();
        indexArray[i].append(indexReversed);
        indexArray.remove(j);
        return indexArray;
    }
    return indexArray;
}
