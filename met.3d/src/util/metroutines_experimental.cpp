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
#include "metroutines_experimental.h"

// standard library imports
#include <cmath>

// related third party imports
#include <log4cplus/loggingmacros.h>

// local application imports
#include "util/mutil.h"

using namespace std;


namespace Met3D
{

namespace MetRoutinesExperimental
{
    double mixingRatio(const double P, const double e)
    {
        using namespace MetConstantsExperimental;

        const double r = (GAS_CONSTANT_RATIO * e) / (P - e);

        return r;
    }


    //###############################################################
    //# YTVAP: computes water vapour pressure EP in hPa
    //#       TC>0 a polynomial of degree 6
    //#       TC<0 empirical Magnus formula
    //# input: TK = temperature in Kelvin
    //# output: EP = water vapour pressure in hPa
    //###############################################################
    double vaporPressureFromTK(const double T)
    {
        //A0=6.107799961
        //A1=4.436518521E-01
        //A2=1.428945805E-02
        //A3=2.650648471E-04
        //A4=3.031240396E-06
        //A5=2.034080948E-08
        //A6=6.136820929E-11
        const double A0 = 6.107799961;
        const double A1 = 4.436518521E-01;
        const double A2 = 1.428945805E-02;
        const double A3 = 2.650648471E-04;
        const double A4 = 3.031240396E-06;
        const double A5 = 2.034080948E-08;
        const double A6 = 6.136820929E-11;

        // Temperature in celsius
        //Z=TK-273.15
        double TC = T - 273.15;
        double e = 0;

        if (TC <= 0.0)
        {
            // A=7.4475*Z/(234.07+Z)
            // EP=6.108*(10.0**A)
            const double A = 7.4475 * TC / (234.07 + TC);
            e = 6.108 * std::pow(10.0, A);
        }
        else
        {
            //EP= A0+Z*(A1+Z*(A2+Z*(A3+Z*(A4+Z*(A5+A6*Z)))))
            e = A0 + TC * (A1 + TC * (A2 + TC * (A3 + TC * (A4 + TC *(A5 + A6 * TC)))));
        }

        return e;
    }

    //###############################################################
    //# BALC: computes the latent heat of condensation of water vapour
    //#       using the expression DL/DT=CPV-CW
    //# input: TK = temperature in Kelvin
    //# output: ALT = latent heat in cal*deg-1
    //###############################################################
    double latentHeat(const double T)
    {
        using namespace MetConstantsExperimental;

        //TC=TK-273.15
        const double TC = T - TEMP_ZERO_KELVIN;
        //ALT=597.3-(0.56*TC)
        double Lw = 597.3 - (0.56 * TC);

        return Lw;
    }

    //###############################################################
    //# BDES function
    //# BDES computes the slope of the saturation vapour pressure
    //#      curve using the Clausius-Clapeyron relation
    //# input: TK = temperature in K
    //#         EP = saturated water vapour pressure in hPa
    //# output: DEP=DE/DT in hPa*K-1
    //###############################################################
    double slopeOfVaporPressure(const double T, const double e)
    {
        using namespace MetConstantsExperimental;

        //Z=TK-273.15
        const double TC = T - TEMP_ZERO_KELVIN;
        //A=597.3-(0.56*Z)
        const double A = 597.3 - (0.56 * TC);

        //DEP=(A*EP)/(0.1101*TK*TK)
        double derivES = (A * e) / (0.1101 * T * T);

        return derivES;
    }

    //###############################################################
    //# BFF function
    //# BFF determines the value on a pseudoadiabatic curve passing
    //#     through a point of a given temperature TK and pressure P
    //# inputs: TK = temperature in K
    //#         P = pressure in hPa
    //#         AC =C*= CP+Q*CW EN CAL*GR-1*K-1
    //# output: AK =C*LN(TK)-R*LN(P-E)+L*Q/TK
    //###############################################################
    double pseudoAdiabaticCurve(const double Tw, const double P, const double r)
    {
        //ESAT=YTVAP(TK)
        double eS = vaporPressureFromTK(Tw);
        //QSAT=BQQ(P,ESAT)
        double rS = mixingRatio(P, eS);
        //ALAT=BALC(TK)
        double Lw = latentHeat(Tw);
        //A=AC*log(TK)
        double A = r * std::log(Tw);
        //B=0.0685*log(P-ESAT)
        double B = 0.0685 * std::log(P - eS);
        //D=ALAT*QSAT/TK
        double D = Lw * rS / Tw;
        //AK=A-B+D
        double AK = A - B + D;

        return AK;
    }

    //###############################################################
    //# BFTW: iteration procedure to compute the wet-bulb temperature
    //#       WTW=(R) temp. of the wet-bulb thermometer in K
    //# input: TK=(R) temperature in K
    //#         P=(R) pressure in hPa
    //#         EV=(R) water vapour pressure in hPa
    //# output: WTW=(R) wet-bulb temperature in K
    //###############################################################
    double computeWTIt(const double T, const double P, const double e)
    {
        //TWI=TK
        //TWF=TK
        //ERR=1.0
        double TwI = T;
        double TwF = T;
        double err = 1.0;

        //while (ERR >= 0.1) do
        while (err >= 0.1)
        {
            //EW=YTVAP(TWI)
            double eW = vaporPressureFromTK(TwI);
            //DEW=BDES(TWI,EW)
            double derivEW = slopeOfVaporPressure(TwI, eW);
            //AL=BALC(TWI)
            double Lw = latentHeat(TwI);
            //E1=(EW-EV)/(P-EW)
            double E1 = (eW - e) / (P - eW);
            //Y=0.2405*(TK-TWI)-0.62197*AL*E1
            double Y = 0.2405 * (T - TwI) - 0.62197 * Lw * E1;
            //Z=-0.2405-0.34833*E1-0.62197*DEW*AL*(P-EV)/(P-EW)**2
            double Z = -0.2405 - 0.34833 * E1 - 0.62197 * derivEW * Lw * (P - e)
                                                / pow((P - eW), 2);
            //TWF=TWI-Y/Z
            TwF = TwI - Y / Z;
            //ERR=maxvalue(abs(TWF-TWI))
            err = std::abs(TwF - TwI);

            //if (ERR>=0.1) then
            //        TWI=TWF
            if (err >= 0.1) { TwI = TwF; }
        }

        return TwF;
    }

    //###############################################################
    //# BCATPW: iteration procedure to compute
    //#         the wet-bulb potential temperature
    //# THW=(R) wet-bulb potential temperature in K
    //###############################################################
    double computeWPTIt(const double Tw, const double r, const double aCurve)
    {
        //TI=TKW
        //P=1000.0
        //ERR=1.0
        double Ti = Tw;
        double P = MetConstantsExperimental::PRESSURE_ZERO_HPA;
        double err = 1.0;

        //while (ERR >= 0.1) do
        while (err >= 0.1)
        {
            //AK=BFF(TI,P,AC)
            double aC = pseudoAdiabaticCurve(Ti, P, r);
            //AK=AK-ACT
            aC -= aCurve;
            //ESAT=YTVAP(TI)
            double eS = vaporPressureFromTK(Ti);
            //DE=BDES(TI,ESAT)
            double derivEs = slopeOfVaporPressure(Ti, eS);
            //ALAT=BALC(TI)
            double Lw = latentHeat(Ti);
            //QSAT=BQQ(P,ESAT)
            double rS = mixingRatio(P, eS);
            //DQ=621.97*DE/(P-ESAT)**2
            double derivRs = 621.97 * derivEs / std::pow((P - eS), 2);
            //A1=AC/TI+0.0685*DE/(P-ESAT)
            double A1 = r / Ti + 0.0685 * derivEs / (P - eS);
            //A2=((-0.56*QSAT+ALAT*DQ)*TI-ALAT*QSAT)/TI**2
            double A2 = ((-0.56 * rS + Lw * derivRs) * Ti - Lw * rS) / (Ti * Ti);
            //AD=A1+A2
            double AD = A1 + A2;
            //TF=TI-AK/AD
            double Tf = Ti - aC / AD;
            //ERR=maxvalue(abs(TF-TI))
            err = std::abs(Tf - Ti);

            //TI=TF
            Ti = Tf;
        }

        //THW=TI
        double thetaW = Ti;
        return thetaW;
    }

    //###############################################################
    //# FUNW function
    //###############################################################
    double funW(double T, double P)
    {
        //EVS=YTVAP(TEK)
        double eS = vaporPressureFromTK(T);
        //QS=BQQ(PR,EVS)
        double rS = mixingRatio(P, eS);
        //AL=BALC(TEK)
        double Lw = latentHeat(T);
        //A=0.2405*log(TEK)-0.0685*log(PR-EVS)+(AL*QS)/TEK
        double A = 0.2405 * std::log(T) - 0.0685 * std::log(P - eS) + (Lw * rS) / T;
        //FW=A
        double fW = A;

        return fW;
    }

    //##############################################################
    //# ITERW function
    //# ITERW: iteration procedure based on the secant method
    //#      : used for TPH calculations at upper levels
    //##############################################################
    double computeWPTItSecant(const double Fw, const double P)
    {
        using namespace MetConstantsExperimental;

        double XI1 = 0;
        double XI2 = 0;
        double XI = 0;

        //if (PR<50.0) then
        //  XI1=45.0+273.15
        //  XI2=55.0+273.15
        if (P < 50.0)
        {
            XI1 = 45.0 + TEMP_ZERO_KELVIN;
            XI2 = 55.0 + TEMP_ZERO_KELVIN;
        }
            //else if (PR<100.0) then
            //   XI1=35.0+273.15
            //   XI2=45.0+273.15
        else if (P < 100.0)
        {
            XI1 = 35.0 + TEMP_ZERO_KELVIN;
            XI2 = 45.0 + TEMP_ZERO_KELVIN;
        }
            //else
            //  XI1=20.0+273.15
            //  XI2=30.0+273.15
        else
        {
            XI1 = 20.0 + TEMP_ZERO_KELVIN;
            XI2 = 30.0 + TEMP_ZERO_KELVIN;
        }

        //ERR=1.0
        double err = 1.0;
        //while (ERR >= 0.1) do
        while (err >= 0.1)
        {
            //FI=FUNW(XI1,1000.0)
            double FwI = funW(XI1, 1000.0);
            //FI1=FI-FWI
            double FI1 = FwI - Fw;
            //FI=FUNW(XI2,1000.0)
            FwI = funW(XI2, 1000.0);
            //FI2=FI-FWI
            double FI2 = FwI - Fw;
            //XI=XI2-FI2*(XI2-XI1)/(FI2-FI1)
            XI = XI2 - FI2 * (XI2 - XI1) / (FI2 - FI1);
            //XI1=XI2
            XI1 = XI2;
            //XI2=XI
            XI2 = XI;
            //ERR=maxvalue(abs(XI1-XI2))
            err = std::abs(XI1 - XI2);
        }

        //FWF=XI
        double thetaW = XI;
        return thetaW;
    };

    //###############################################################
    //# TPH: computes the wet-bulb temperature and wet-bulb potential
    //#      temperature by an iterative process
    //# input:
    //#        PR = (R) pressure in hPa
    //#        TEK = (R) temperature in K
    //#        HR = (R) relative humidity in %
    //# output:
    //#        (XTW = (R) wet-bulb temperature in K)
    //#        XTHW = (R) wet-bulb POT. temperature in K
    //###############################################################
    double computeWPT(const double T, const double P, const double q)
    {
        using namespace MetConstantsExperimental;

        // g * g^-1
        //double r = mixingRatio(q);
        double thetaW = 0;

        // if (PR > 200.0) then
        if (P > 200)
        {
            // EVS=YTVAP(TEK)
            //double eS = vaporPressureFromTK(T);

            // Stull Clausius-Clapeyron eq. 4.1a on page 88
    //        double eS2 = VAPOR_PRESSURE * std::exp(LATENT_HEAT_VAPORIZATION
    //                      / GAS_CONSTANT_WATER_VAPOR
    //                      * (1.0 / TEMP_ZERO_KELVIN - 1.0 / T));

            // Stull Teten's formula eq 4.2 on page 89
    //        double eS = VAPOR_PRESSURE * std::exp(17.2694 * (T - TEMP_ZERO_KELVIN)
    //                      / (T - 35.86));
    //
    //        //QS=BQQ(PR,EVS)
    //        //double rS = mixingRatio(P, eS);
    //        //double r = mixingRatio(q);
    //        double qS = (GAS_CONSTANT_RATIO * eS) / (P - eS * (1 - GAS_CONSTANT_RATIO));
    //
    //        double RH = q / qS;
    //        //double RH = r / rS;
    //
    //        // EP=0.01*HR*EVS
    //        double e2 = RH * eS;

            // Stull used specific humidity to vapor pressure relationship, eq 4.7 on page 91
            double e = (q * P) / (GAS_CONSTANT_RATIO + q * (1 - GAS_CONSTANT_RATIO));
            //double e = vaporPressure(P, r);

            // XTW=BFTW(TEK,PR,EP)
            double Tw = computeWTIt(T, P, e);
            // AC=BQQ(PR,EP)
            // AC=0.2405+AC
            double rW = mixingRatio(P, e);
            rW += 0.2405;

            // ACT=BFF(XTW,PR,AC)
            double aCurve = pseudoAdiabaticCurve(Tw, P, rW);
            // XTHW=BCATPW(XTW,AC,ACT)
            thetaW = computeWPTIt(Tw, rW, aCurve);
        }
        else
        {
            // XTW=TEK
            double Tw = T;
            // FWI=FUNW(XTW,PR)
            double Fw = funW(Tw, P);
            // XTHW=ITERW(FWI,PR)
            thetaW = computeWPTItSecant(Fw, P);
        }

        return thetaW;
    }
} // namespace MetRoutinesExperimental

} // namespace Met3D
