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
            double stepGap = 1.0/samples[0];
            std::array<double, 4> t1 = {{0.0, 0.0, 0.0, 0.0}};
            std::array<double, 4> t2 = {{0.0, 0.0, 0.0, 0.0}};            
            coordinate v;
            double vLength[samples[0]+1];

            double temp = 0.0;

            v.coord.reserve(img.TargetDimension);
            for (int k = 0; k < img.TargetDimension; k++)
            {
                v.coord.push_back(0.0);
            }

            t1[0] = img.DomainInterval[0][0];
            for (int i = 0; i < samples[0]+1; i++)
            {
                v = img.DiffBezierMap(t1,0);
                vLength[i] = sqrt(VectorLenSquare(v.coord));  // |γ'(t1)|
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
                        temp = img.DistanceSquare(t1,t2);  // |γ(t1)-γ(t2)|^2
                        if (isnan(temp) | temp < 1.0E-80)
                        {
                            std::cout << " distance square too small: " << temp << std::endl;
                            energy = 1.9E+200;
                            break;
                        }
                        
                        // O'Hara energy integral
                        energy = energy + (( 1.0/temp )  * vLength[i] * vLength[j] - ( 1.0/pow(std::min(abs(t1[0]-t2[0]),1.0-abs(t1[0]-t2[0])), 2) ));
                    }

                    if ((i==0 & j==samples[0]) || (j ==0 && i==samples[0]))
                    {
                        energy = energy + 1.0E+20*img.DistanceSquare(t1,t2);
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
            double stepGap = 1.0/samples[0];
            std::array<double, 4> t1 = {{0.0, 0.0, 0.0, 0.0}};
            std::array<double, 4> t2 = {{0.0, 0.0, 0.0, 0.0}};            
            coordinate v;
            double vLength[samples[0]+1];

            coordinate bmapt1;
            coordinate bmapt2;
            coordinate bmapDifft1;
            coordinate bmapDifft2;
            double bmapCptDiffFactort1;
            double bmapCptDiffFactort2;
            double disCptFactor;

            double distanceSquare = 0.0;
            double arcDistance = 0.0;
            double OHaraEnergySimple = 0.0;
            
            v.coord.reserve(img.TargetDimension);
            result.coord.reserve(img.TargetDimension);
            for (int k = 0; k < img.TargetDimension; k++)
            {
                v.coord.push_back(0.0);
                result.coord.push_back(0.0);
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
                        OHaraEnergySimple =  (1.0/distanceSquare); // the second term integral int int 1/arcLength^2 |γ'(t1)||γ'(t2)| dt1 dt2 = const. + int 1/arcLength |γ'(t1)| dt1 = const. + c. ln ( arcLength(t1,t2) ) is independent on γ
                        
                        result.coord[0] = result.coord[0] - ( ( 2.0/pow(distanceSquare,2) ) * ( bmapt1.coord[0] - bmapt2.coord[0] ) * disCptFactor ) * vLength[i] * vLength[j] 
                        + ( OHaraEnergySimple / (vLength[i]+1.0E-14) * vLength[j]) * ( bmapDifft1.coord[0] * bmapCptDiffFactort1 ) 
                        + ( OHaraEnergySimple * vLength[i] / (vLength[j]+1.0E-14)) * ( bmapDifft2.coord[0] * bmapCptDiffFactort2 );

                        result.coord[1] = result.coord[1] - ( ( 2.0/pow(distanceSquare,2) ) * ( bmapt1.coord[1] - bmapt2.coord[1] ) * disCptFactor ) * vLength[i] * vLength[j] 
                        + ( OHaraEnergySimple / (vLength[i]+1.0E-14) * vLength[j]) * ( bmapDifft1.coord[1] * bmapCptDiffFactort1 ) 
                        + ( OHaraEnergySimple * vLength[i] / (vLength[j]+1.0E-14)) * ( bmapDifft2.coord[1] * bmapCptDiffFactort2 );

                        result.coord[2] = result.coord[2] - ( ( 2.0/pow(distanceSquare,2) ) * ( bmapt1.coord[2] - bmapt2.coord[2] ) * disCptFactor ) * vLength[i] * vLength[j] 
                        + ( OHaraEnergySimple / (vLength[i]+1.0E-14) * vLength[j]) * ( bmapDifft1.coord[2] * bmapCptDiffFactort1 ) 
                        + ( OHaraEnergySimple * vLength[i] / (vLength[j]+1.0E-14)) * ( bmapDifft2.coord[2] * bmapCptDiffFactort2 );

                    }

                    if ((i==0 & j==samples[0]) || (j ==0 && i==samples[0]))
                    {
                        bmapt1 = img.BezierMap(t1);
                        bmapt2 = img.BezierMap(t2);
                        bmapCptDiffFactort1 = (img.CPtFactorOfBezierMap(cpt,t1) - img.CPtFactorOfBezierMap(cpt,t2));

                        result.coord[0] = result.coord[0] + 2.0E+20*(bmapt1.coord[0] - bmapt2.coord[0])*bmapCptDiffFactort1;
                        result.coord[1] = result.coord[1] + 2.0E+20*(bmapt1.coord[1] - bmapt2.coord[1])*bmapCptDiffFactort1;
                        result.coord[2] = result.coord[2] + 2.0E+20*(bmapt1.coord[2] - bmapt2.coord[2])*bmapCptDiffFactort1;
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
    std::vector<int> cPtsN = {15};

    std::vector<int> samples = {90};

    coordinate inipt;      // 曲線初始位置
    inipt.coord.push_back(0.2);
    inipt.coord.push_back(0.2);
    inipt.coord.push_back(0.2);

    std::map<std::array<int,4>, coordinate> fixPts;
    fixPts[{{0,0,0,0}}] = inipt;
    fixPts[{{15,0,0,0}}] = inipt;


    std::array<double,2> intvl = {{0.0, 1.0}};
    std::vector<std::array<double,2>> interval;
    interval.push_back(intvl);
    
    BezierImage initImg;
    //initImg.ReadFile("data/9.11058_param.txt");
    initImg.RandGenerator(cPtsN, fixPts, interval);

    SystemCalculator calculator;

    BezierImage optimizedResult = calculator.ConjugateGradientOptimization(initImg, samples, 5.0, 1.E-5, true, 1.0E-6, true, true, 1.0);
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