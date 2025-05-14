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
            double loss = 0.0;
            double timeInterval = (img.DomainInterval[0][1]-img.DomainInterval[0][0]);
            double timeGap = timeInterval/(double)samples[0];
            std::array<double, 4> ts = {{0.0, 0.0, 0.0, 0.0}};

            coordinate v;
            coordinate v0;
            
            ts[0] = img.DomainInterval[0][0];
            for (int i = 0; i < samples[0]+1; i++)
            {
                v = img.BezierMap(ts);
                v0 = img.DiffBezierMap(ts,0);
    
                loss = loss + exp(-ts[0]/timeInterval*2.0)*( pow(v0.coord[0] - 10.0*(v.coord[1]-v.coord[0]),2) + pow(v0.coord[1] - v.coord[0]*(28.0-v.coord[2])+v.coord[1],2)  + pow(v0.coord[2] - v.coord[0]*v.coord[1] + 2.66666667*v.coord[2],2) );
    
                if (i==0)
                {
                    loss = loss + 1.0E6*pow(v.coord[0]-10.0,2) + 1.0E6*pow(v.coord[1]-10.0,2) + 1.0E6*pow(v.coord[2]-10.0,2);
                }                    
                ts[0] = ts[0] + timeGap;
            }
            
            return loss*timeGap;
        }


        /*  ==============================================
            The gradient of the residual of Bezier approximate solution for the equation
            ==============================================  */
        coordinate GradientEquationResidual(BezierImage img, int cpt, std::vector<int> samples) override
        {
            coordinate result;
            double timeInterval = (img.DomainInterval[0][1]-img.DomainInterval[0][0]);
            double timeGap = timeInterval/(double)samples[0];
            double cptFactorOfBezier;
            double cptFactorOfDiffBezier0;

            double temp0;
            double temp1;
            double temp2;
            double tempfactor;

            std::array<double, 4> ts = {{0.0, 0.0, 0.0, 0.0}};

            coordinate v;
            coordinate v0;

            result.coord.reserve(img.TargetDimension);
            for (int i = 0; i < img.TargetDimension; i++)
            {
                result.coord.push_back(0.0);
            }

            ts[0] = img.DomainInterval[0][0];
            for (int i = 0; i < samples[0]+1; i++)
            {
                v = img.BezierMap(ts);
                v0 = img.DiffBezierMap(ts,0);
                cptFactorOfBezier = img.CPtFactorOfBezierMap(cpt, ts);
                cptFactorOfDiffBezier0 = img.CPtFactorOfDiffBezierMap(cpt, ts, 0);
                tempfactor = 2.0*exp(-ts[0]/timeInterval*2.0);

                temp0 = tempfactor*(v0.coord[0] - 10.0*(v.coord[1]-v.coord[0]));
                temp1 = tempfactor*(v0.coord[1] - v.coord[0]*(28.0-v.coord[2])+v.coord[1]);
                temp2 = tempfactor*(v0.coord[2] - v.coord[0]*v.coord[1] + 2.66666667*v.coord[2]);

                result.coord[0] = result.coord[0] + temp0*(10.0*cptFactorOfBezier + cptFactorOfDiffBezier0) - temp1*(28.0-v.coord[2])*cptFactorOfBezier - temp2*v.coord[1]*cptFactorOfBezier;

                result.coord[1] = result.coord[1] - temp0*(10.0*cptFactorOfBezier) + temp1*(cptFactorOfDiffBezier0 + cptFactorOfBezier) - temp2*v0.coord[0]*cptFactorOfBezier;

                result.coord[2] = result.coord[2] + temp1*v.coord[0]*cptFactorOfBezier + temp2*(cptFactorOfDiffBezier0 + 2.66666667*cptFactorOfBezier);

                if (i==0)
                {
                    result.coord[0] = result.coord[0] + 2.0E6*(v.coord[0]-10.0)*cptFactorOfBezier;
                    result.coord[1] = result.coord[1] + 2.0E6*(v.coord[1]-10.0)*cptFactorOfBezier;
                    result.coord[2] = result.coord[2] + 2.0E6*(v.coord[2]-10.0)*cptFactorOfBezier;
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
    std::vector<int> samples = {1200};
    
    BezierImage initImg;
    initImg.ReadFile("data/18.5126_param.txt");

    SystemCalculator calculator;

    BezierImage optimizedResult = calculator.ConjugateGradientOptimization(initImg, samples, 0.1);
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