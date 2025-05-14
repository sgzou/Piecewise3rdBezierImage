#include <iostream>
#include <fstream>
#include <math.h>
#include <filesystem>
#include "../OptimizeOperation.h"

class SystemCalculator : public OptimizeOperation
{
    public:

        /*  ==============================================
            The residual of Bezier approximate solution for the equation
            ==============================================  */    
        double EquationResidual(BezierImage img, std::vector<int> samples) override
        {
            double energy = 0.0;
            double stepGap = (img.DomainInterval[0][1]-img.DomainInterval[0][0])/samples[0];
            std::array<double, 4> t1 = {{0.0, 0.0, 0.0, 0.0}};
            std::array<double, 4> t2 = {{0.0, 0.0, 0.0, 0.0}};            
            coordinate v;
            double tangentLength[samples[0]+1];
            double arcLengthFrom0[samples[0]+1];

            v.coord.reserve(img.TargetDimension);
            for (int k = 0; k < img.TargetDimension; k++)
            {
                v.coord.push_back(0.0);
            }

            t1[0] = img.DomainInterval[0][0];
            for (int i = 0; i < samples[0]+1; i++)
            {
                v = img.DiffBezierMap(t1,0);
                tangentLength[i] = sqrt(VectorLenSquare(v.coord));  // |γ'(t1)|
                if (i==0) // record arc length from 0 to t1
                {
                    arcLengthFrom0[i] = tangentLength[i]*stepGap;
                }
                else
                {
                    arcLengthFrom0[i] = arcLengthFrom0[i-1] + tangentLength[i]*stepGap;
                }
                t1[0] = t1[0] + stepGap;
            }

            t1[0] = img.DomainInterval[0][0];
            t2[0] = img.DomainInterval[0][0];
            for (int i = 0; i < samples[0]+1; i++)
            {
                for (int j = 0; j < samples[0]+1; j++)
                {
                    if (i != j & !( i==0 & j == samples[0]) & !( i==samples[0] & j == 0))
                    {
                        // O'Hara energy integral
                        energy = energy + (( 1.0/img.DistanceSquare(t1,t2) ) - ( 1.0/pow(std::min(abs(arcLengthFrom0[i]-arcLengthFrom0[j]),arcLengthFrom0[samples[0]]-abs(arcLengthFrom0[i]-arcLengthFrom0[j])), 2) )) * tangentLength[i] * tangentLength[j];                        
                    }

                    if ((i==0 & j==samples[0]) || (j ==0 && i==samples[0]))
                    {
                        energy = energy + 1.0E+14*img.DistanceSquare(t1,t2);
                    }

                    t2[0] = t2[0] + stepGap;
                }
                t1[0] = t1[0] + stepGap;
                t2[0] = img.DomainInterval[0][0];
            }
            return energy * stepGap * stepGap;
        }


        /*  ==============================================
            The gradient of the residual of Bezier approximate solution for the equation
            ==============================================  */
        coordinate GradientEquationResidual(BezierImage img, int cpt, std::vector<int> samples) override
        {
            coordinate result;
            double energy = 0.0;
            double stepGap = (img.DomainInterval[0][1]-img.DomainInterval[0][0])/samples[0];
            std::array<double, 4> t1 = {{0.0, 0.0, 0.0, 0.0}};
            std::array<double, 4> t2 = {{0.0, 0.0, 0.0, 0.0}};            
            coordinate v;
            double vLength[samples[0]+1];
            double arcLengthFrom0[samples[0]+1];
            coordinate arcLengthFrom0Derivative[samples[0]+1];

            coordinate bmapt1;
            coordinate bmapt2;
            coordinate bmapDifft1;
            coordinate bmapDifft2;
            double bmapCptDiffFactort1;
            double bmapCptDiffFactort2;
            double disCptFactor;

            double distanceSquare = 0.0;
            double arcLength = 0.0;
            coordinate arcLengthDerivative;
            double OHaraEnergyElement = 0.0;
            
            v.coord.reserve(img.TargetDimension);
            result.coord.reserve(img.TargetDimension);
            arcLengthDerivative.coord.reserve(img.TargetDimension);
            for (int k = 0; k < img.TargetDimension; k++)
            {
                v.coord.push_back(0.0);
                result.coord.push_back(0.0);
                arcLengthDerivative.coord.push_back(0.0);
            }

            if (cpt == 0 || cpt == img.ControlPointsTotalNumber-1)
            {
                return result;
            }

            t1[0] = img.DomainInterval[0][0];
            for (int i = 0; i < samples[0]+1; i++)
            {
                v = img.DiffBezierMap(t1,0);
                vLength[i] = sqrt(VectorLenSquare(v.coord));  // |γ'(t1)|
                if (i==0) // record arc length from 0 to t1
                {
                    arcLengthFrom0[i] = vLength[i]*stepGap;
                    for (int j = 0; j < img.TargetDimension; j++)
                    {
                        arcLengthFrom0Derivative[i].coord.push_back( v.coord[j]*img.CPtFactorOfDiffBezierMap(cpt,t1,0)*stepGap / vLength[i] );
                    }
                }
                else
                {
                    arcLengthFrom0[i] = arcLengthFrom0[i-1] + vLength[i]*stepGap;
                    for (int j = 0; j < img.TargetDimension; j++)
                    {
                        arcLengthFrom0Derivative[i].coord.push_back( arcLengthFrom0Derivative[i-1].coord[j] + v.coord[j]*img.CPtFactorOfDiffBezierMap(cpt,t1,0)*stepGap / vLength[i] );
                    }
                }
                t1[0] = t1[0] + stepGap;
            }

            t1[0] = img.DomainInterval[0][0];
            t2[0] = img.DomainInterval[0][0];
            for (int i = 0; i < samples[0]+1; i++)
            {
                for (int j = 0; j < samples[0]+1; j++)
                {
                    if (i != j & !( i==0 & j == samples[0]) & !( i==samples[0] & j == 0))
                    {   
                        bmapt1 = img.BezierMap(t1);
                        bmapt2 = img.BezierMap(t2);
                        bmapDifft1 = img.DiffBezierMap(t1,0);
                        bmapDifft2 = img.DiffBezierMap(t2,0);
                        bmapCptDiffFactort1 = img.CPtFactorOfDiffBezierMap(cpt,t1,0);
                        bmapCptDiffFactort2 = img.CPtFactorOfDiffBezierMap(cpt,t2,0);
                        disCptFactor = img.CPtFactorOfBezierMap(cpt,t1) - img.CPtFactorOfBezierMap(cpt,t2);

                        distanceSquare = img.DistanceSquare(t1,t2);
                        if ( abs(arcLengthFrom0[i]-arcLengthFrom0[j]) < arcLengthFrom0[samples[0]]-abs(arcLengthFrom0[i]-arcLengthFrom0[j]) )
                        {
                            arcLength = abs(arcLengthFrom0[i]-arcLengthFrom0[j]);
                            for (int c = 0; c < img.TargetDimension; c++)
                            {
                                arcLengthDerivative.coord[c] = copysign(1.0, arcLengthFrom0[i]-arcLengthFrom0[j]) * (arcLengthFrom0Derivative[i].coord[c] - arcLengthFrom0Derivative[j].coord[c]);
                            }
                        }
                        else
                        {
                            arcLength = arcLengthFrom0[samples[0]]-abs(arcLengthFrom0[i]-arcLengthFrom0[j]);
                            for (int c = 0; c < img.TargetDimension; c++)
                            {
                                arcLengthDerivative.coord[c] = arcLengthFrom0Derivative[samples[0]].coord[c] - copysign(1.0, arcLengthFrom0[i]-arcLengthFrom0[j])*(arcLengthFrom0Derivative[i].coord[c] - arcLengthFrom0Derivative[j].coord[c]);
                            }
                        }

                        
                        OHaraEnergyElement =  (1.0/distanceSquare) - (1.0/pow(arcLength,2));
                        
                        result.coord[0] = result.coord[0] - ( ( 2.0/pow(distanceSquare,2) ) * ( bmapt1.coord[0] - bmapt2.coord[0] ) * disCptFactor - (2.0/pow(arcLength,3))* arcLengthDerivative.coord[0] ) * vLength[i] * vLength[j] 
                        + ( OHaraEnergyElement / (vLength[i]+1.0E-14) * vLength[j]) * ( bmapDifft1.coord[0] * bmapCptDiffFactort1 ) 
                        + ( OHaraEnergyElement * vLength[i] / (vLength[j]+1.0E-14)) * ( bmapDifft2.coord[0] * bmapCptDiffFactort2 );

                        result.coord[1] = result.coord[1] - ( ( 2.0/pow(distanceSquare,2) ) * ( bmapt1.coord[1] - bmapt2.coord[1] ) * disCptFactor - (2.0/pow(arcLength,3))* arcLengthDerivative.coord[1] ) * vLength[i] * vLength[j] 
                        + ( OHaraEnergyElement / (vLength[i]+1.0E-14) * vLength[j]) * ( bmapDifft1.coord[1] * bmapCptDiffFactort1 ) 
                        + ( OHaraEnergyElement * vLength[i] / (vLength[j]+1.0E-14)) * ( bmapDifft2.coord[1] * bmapCptDiffFactort2 );

                        result.coord[2] = result.coord[2] - ( ( 2.0/pow(distanceSquare,2) ) * ( bmapt1.coord[2] - bmapt2.coord[2] ) * disCptFactor - (2.0/pow(arcLength,3))* arcLengthDerivative.coord[2] ) * vLength[i] * vLength[j] 
                        + ( OHaraEnergyElement / (vLength[i]+1.0E-14) * vLength[j]) * ( bmapDifft1.coord[2] * bmapCptDiffFactort1 ) 
                        + ( OHaraEnergyElement * vLength[i] / (vLength[j]+1.0E-14)) * ( bmapDifft2.coord[2] * bmapCptDiffFactort2 );

                    }

                    if ((i==0 & j==samples[0]) || (j ==0 && i==samples[0]))
                    {
                        bmapt1 = img.BezierMap(t1);
                        bmapt2 = img.BezierMap(t2);
                        bmapCptDiffFactort1 = (img.CPtFactorOfBezierMap(cpt,t1) - img.CPtFactorOfBezierMap(cpt,t2));

                        result.coord[0] = result.coord[0] + 2.0E+14*(bmapt1.coord[0] - bmapt2.coord[0])*bmapCptDiffFactort1;
                        result.coord[1] = result.coord[1] + 2.0E+14*(bmapt1.coord[1] - bmapt2.coord[1])*bmapCptDiffFactort1;
                        result.coord[2] = result.coord[2] + 2.0E+14*(bmapt1.coord[2] - bmapt2.coord[2])*bmapCptDiffFactort1;
                    }

                    t2[0] = t2[0] + stepGap;
                }
                t1[0] = t1[0] + stepGap;
                t2[0] = img.DomainInterval[0][0];
            }

            for (int k = 0; k < img.TargetDimension; k++)
            {
                result.coord[k] = result.coord[k] * stepGap * stepGap;
            }
            return result;
        }

};


/*  ==============================================
    The main calculation flow
    ==============================================  */
int main(int argc, char * argv[])  
{
    std::vector<int> cPtsN = {21};

    std::vector<int> samples = {120};

    coordinate inipt;      // the initial point of the curve
    inipt.coord.push_back(0.2);
    inipt.coord.push_back(0.2);
    inipt.coord.push_back(0.2);

    std::map<std::array<int,4>, coordinate> fixPts;
    fixPts[{{0,0,0,0}}] = inipt;
    fixPts[{{21,0,0,0}}] = inipt;


    std::array<double,2> intvl = {{0.0, 1.0}};
    std::vector<std::array<double,2>> interval;
    interval.push_back(intvl);
    
    BezierImage initImg;
    //initImg.ReadFile("data/1.81877_param.txt");
    initImg.RandGenerator(cPtsN, fixPts, interval);

    SystemCalculator calculator;

    BezierImage optimizedResult = calculator.ConjugateGradientOptimization(initImg, samples, 40.0, 1.0E-5, true, 1.0E-6, true, true ,1.0);
    double residual = calculator.EquationResidual(optimizedResult, samples);

    std::filesystem::create_directory("data");

    std::ostringstream tempStream;
    tempStream << "data/" << residual << "_coord.txt";

    std::string fileName;
    fileName = tempStream.str();

    optimizedResult.WriteConfigFile(fileName, samples);

    tempStream.clear();
    tempStream.str("");
    tempStream << "data/" << residual << "_param.txt";
    fileName = tempStream.str();

    optimizedResult.WriteParamFile(fileName);

    return 0;
}