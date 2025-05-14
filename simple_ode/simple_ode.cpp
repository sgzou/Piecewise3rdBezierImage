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
            double residual = 0.0;
            double timeGap = (img.DomainInterval[0][1]-img.DomainInterval[0][0])/(double)samples[0];
            std::array<double, 4> ts = {{0.0, 0.0, 0.0, 0.0}};

            coordinate v;
            coordinate v0;
            coordinate v00;
            
            ts[0] = img.DomainInterval[0][0];

            for (int i = 0; i < samples[0]+1; i++)
            {
                v00 = img.DoubleDiffBezierMap(ts,0,0);

                residual = residual + pow(v00.coord[0] - 2.0*pow(ts[0],5),2);

                // initial (boundary) data
                if (i==0)
                {
                    v = img.BezierMap(ts);
                    v0 = img.DiffBezierMap(ts,0);

                    residual = residual + 1.0E5*pow(v.coord[0] - 2.0,2) + 1.0E5*pow(v0.coord[0] - 3.0,2);
                }
                
                ts[0] = ts[0] + timeGap;
            }

            return residual * timeGap;
        }


        /*  ==============================================
            The gradient of the residual of Bezier approximate solution for the equation
            ==============================================  */
        coordinate GradientEquationResidual(BezierImage img, int cpt, std::vector<int> samples) override
        {
            coordinate result;
            double timeGap = (img.DomainInterval[0][1]-img.DomainInterval[0][0])/(double)samples[0];
            double innerProduct = 0.0;
            double cptFactorOfDiffBezier;
            double cptFactorOfDiffBezier0;

            std::array<double, 4> ts = {{0.0, 0.0, 0.0, 0.0}};

            coordinate v;
            coordinate v0;
            coordinate v00;

            ts[0] = img.DomainInterval[0][0];

            for (int i = 0; i < img.TargetDimension; i++)
            {
                result.coord.push_back(0.0);
            }

            for (int i = 0; i < samples[0]+1; i++)
            {
                v00 = img.DoubleDiffBezierMap(ts,0,0);

                result.coord[0] = result.coord[0] + 2.0*(v00.coord[0] - 2.0*pow(ts[0],5))*img.CPtFactorOfDoubleDiffBezierMap(cpt, ts, 0 ,0);

                if (i==0)
                {
                    v = img.BezierMap(ts);
                    v0 = img.DiffBezierMap(ts,0);

                    result.coord[0] = result.coord[0] + 2.0E5*(v.coord[0] - 2.0)*img.CPtFactorOfBezierMap(cpt, ts) + 2.0E5*(v0.coord[0] - 3.0)*img.CPtFactorOfDiffBezierMap(cpt, ts, 0);
                }

                ts[0] = ts[0] + timeGap;
            }

            for (int k = 0; k < img.TargetDimension; k++)
            {
                result.coord[k] = result.coord[k]*timeGap;
            }
            return result;
        }

};


/*  ==============================================
    The main calculation flow
    ==============================================  */
int main(int argc, char * argv[])  
{    
    std::vector<int> samples = {60};
    
    BezierImage initImg;
    initImg.ReadFile("data/simple_param.txt");

    SystemCalculator calculator;

    BezierImage optimizedResult = calculator.ConjugateGradientOptimization(initImg, samples, 0.5, 1.0E-6);
    double residual = calculator.EquationResidual(optimizedResult, samples);

    std::filesystem::create_directory("data");

    std::ostringstream tempStream;
    tempStream << "data/" << residual << "_coord.txt";

    std::string fileName;
    fileName = tempStream.str();

    optimizedResult.WriteConfigWithVariableFile(fileName, samples);

    tempStream.clear();
    tempStream.str("");
    tempStream << "data/" << residual << "_param.txt";
    fileName = tempStream.str();

    optimizedResult.WriteParamFile(fileName);

    return 0;
}