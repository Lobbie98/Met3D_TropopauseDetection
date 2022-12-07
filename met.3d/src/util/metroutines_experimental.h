/******************************************************************************
**
**  This file is part of Met.3D -- a research environment for the
**  three-dimensional visual exploration of numerical ensemble weather
**  prediction data.
**
**  Copyright 2018 Marc Rautenhaus
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
#ifndef METROUTINES_EXPERIMENTAL_H
#define METROUTINES_EXPERIMENTAL_H

// standard library imports

// related third party imports

// local application imports
#include "util/mutil.h"

/**
  NOTE: Use this file to implement experimental new meteorological computations
  that need testing before going into "metroutines.h/.cpp".
 */

namespace Met3D
{

/******************************************************************************
***                               CONSTANTS                                 ***
*******************************************************************************/

namespace MetConstants
{
}


namespace MetConstantsExperimental
{
    const double TEMP_ZERO_KELVIN = 273.15; // K
    const double PRESSURE_ZERO_HPA = 1000.0; // hPa
    const double GAS_CONSTANT_DRY_AIR = 287.053; // R_d in J K^-1 kg^-1 [S]
    const double GAS_CONSTANT_RATIO = 0.62197; // eps in g g^-1 / kg kg^-1 [S]
}


/******************************************************************************
***                                METHODS                                  ***
*******************************************************************************/

namespace MetRoutinesExperimental
{
    /**
     * Compute the mixing ratio according to eq. 4.4 in Stull
     * @param P
     * @param e
     * @return mixing ratio in g g^-1 / kg kg^-1
     */
    double mixingRatio(const double P, const double e);

    // Code obtained from Tim Hewson (Metview macro)
    double vaporPressureFromTK(const double T);
    double latentHeat(const double T);
    double slopeOfVaporPressure(const double T, const double eS);
    double pseudoAdiabaticCurve(const double T, const double P, const double r);
    double computeWTIt(const double T, const double P, const double e);
    double computeWPTIt(const double T, const double r, const double aCurve);
    double funW(double T, double P);
    double computeWPTItSecant(const double Fw, const double P);

    double computeWPT(const double T, const double P, const double q);
}


} // namespace Met3D

#endif // METROUTINES_EXPERIMENTAL_H
