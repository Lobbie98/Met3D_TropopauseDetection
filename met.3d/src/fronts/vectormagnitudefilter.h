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

#ifndef MET_3D_VECTORMAGNITUDEFILTER_H
#define MET_3D_VECTORMAGNITUDEFILTER_H

#include "data/processingwpdatasource.h"
#include "data/structuredgrid.h"

namespace Met3D
{
class MVectorMagnitudeFilter
        : public MProcessingWeatherPredictionDataSource
{
public:
    MVectorMagnitudeFilter();
    ~MVectorMagnitudeFilter();

    void setInputSource(int id, MWeatherPredictionDataSource* s);

    MStructuredGrid* produceData(MDataRequest request) override;

    MTask* createTaskGraph(MDataRequest request) override;

    QList<MVerticalLevelType> availableLevelTypes() override;

    QStringList availableVariables(MVerticalLevelType levelType) override;

    QSet<unsigned int> availableEnsembleMembers(MVerticalLevelType levelType,
                                                const QString& variableName) override;

    QList<QDateTime> availableInitTimes(MVerticalLevelType levelType,
                                        const QString& variableName) override;

    QList<QDateTime> availableValidTimes(MVerticalLevelType levelType,
                                         const QString& variableName,
                                         const QDateTime& initTime) override;

    QString variableLongName(MVerticalLevelType levelType,
                             const QString&     variableName) override;

    QString variableStandardName(MVerticalLevelType levelType,
                             const QString&     variableName) override;

    QString variableUnits(MVerticalLevelType levelType,
                             const QString&     variableName) override;

protected:
    const QStringList locallyRequiredKeys() override;

    MDataRequest constructInputSourceRequestFromRequest(int id,
                                                        MDataRequestHelper rh);

    QVector<MWeatherPredictionDataSource*> inputSource;
};

} // namespace Met3D

#endif //MET_3D_VECTORMAGNITUDEFILTER_H
