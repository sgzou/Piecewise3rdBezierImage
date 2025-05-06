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
        double EquationResidual(BezierImage obj, std::vector<int> samples) override
        {
            double loss = 0.0;
            double timeGap = (obj.DomainInterval[0][1]-obj.DomainInterval[0][0])/(double)samples[0];
            std::array<double, 4> ts = {{0.0, 0.0, 0.0, 0.0}};

            coordinate v;
            coordinate v0;
            coordinate v00;
            
            ts[0] = obj.DomainInterval[0][0];

            for (int i = 0; i < samples[0]+1; i++)
            {
                v = obj.BezierMap(ts);
                v0 = obj.DiffBezierMap(ts,0);
                v00 = obj.DoubleDiffBezierMap(ts,0,0);

                loss = loss + pow(v00.coord[0] - 36.0*sin(6.0*ts[0]),2);

                if (i==0)
                {
                    loss = loss + 100.0*pow(v.coord[0] - 2.0,2) + 100.0*pow(v0.coord[0] + 3.0,2);
                }
                
                ts[0] = ts[0] + timeGap;
            }

            return loss * timeGap;
        }


        /*  ==============================================
            The gradient of the residual of Bezier approximate solution for the equation
            ==============================================  */
        coordinate GradientEquationResidual(BezierImage obj, int cpt, std::vector<int> samples) override
        {
            coordinate result;
            double timeGap = (obj.DomainInterval[0][1]-obj.DomainInterval[0][0])/(double)samples[0];
            double innerProduct = 0.0;
            double cptFactorOfDiffBezier;
            double cptFactorOfDiffBezier0;

            std::array<double, 4> ts = {{0.0, 0.0, 0.0, 0.0}};

            coordinate v;
            coordinate v0;
            coordinate v00;

            ts[0] = obj.DomainInterval[0][0];

            for (int i = 0; i < obj.TargetDimeion; i++)
            {
                result.coord.push_back(0.0);
            }

            for (int i = 0; i < samples[0]+1; i++)
            {
                v = obj.BezierMap(ts);
                v0 = obj.DiffBezierMap(ts,0);
                v00 = obj.DoubleDiffBezierMap(ts,0,0);

                result.coord[0] = result.coord[0] + 2.0*(v00.coord[0] - 36.0*sin(6.0*ts[0]))*obj.CPtFactorOfDoubleDiffBezierMap(cpt, ts, 0 ,0);

                if (i==0)
                {
                    result.coord[0] = result.coord[0] + 200.0*(v.coord[0] - 2.0)*obj.CPtFactorOfBezierMap(cpt, ts) + 200.0*(v0.coord[0] + 3.0)*obj.CPtFactorOfDiffBezierMap(cpt, ts, 0);
                }

                ts[0] = ts[0] + timeGap;
            }

            for (int k = 0; k < obj.TargetDimeion; k++)
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
    
    BezierImage obj;
    double correctionFactor = 1.0;
    obj.ReadFile("data/0_1.9947_param.txt");

    SystemCalculator calculator;

    BezierImage optimizedResult = calculator.ConjugateGradientOptimization(obj, samples, 0.1);
    double residual = calculator.EquationResidual(optimizedResult, samples);

    std::filesystem::create_directory("data");

    std::ostringstream tempStream;
    tempStream << "data/" << residual << "_coord.txt";

    std::string fileName;
    fileName = tempStream.str();

    obj.WriteConfigWithVariableFile(fileName, samples);

    tempStream.clear();
    tempStream.str("");
    tempStream << "data/" << residual << "_param.txt";
    fileName = tempStream.str();

    obj.WriteParamFile(fileName);

    return 0;
}