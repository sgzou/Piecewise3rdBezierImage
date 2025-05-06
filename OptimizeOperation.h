#include <iostream>
#include <math.h>
#include "Piecewise3rdBezierImage.h"

class OptimizeOperation
{
    public:

        /*  ==============================================
            The residual of the equation for the Bezier approximate solution
            ==============================================  */
        virtual double EquationResidual(BezierImage img, std::vector<int> samples){
            return 0.0;
        }


        /*  ==============================================
            The gradient of the residual of the equation for the Bezier approximate solution
            ==============================================  */
        virtual coordinate GradientEquationResidual(BezierImage img, int cpt, std::vector<int> samples){
            coordinate cc;
            return cc;
        }


        /*  ==============================================
            The (control point) vector length square (in target space).
            ==============================================  */
        double VectorLenSquare(std::vector<double> vect){
            double result = 0.0;
            for (int i = 0; i < vect.size(); i++)
            {
                result = result + pow(vect[i], 2);
            }
            return result;
        };


        /*  ==============================================
            All control points with their target space vector form a large vector.
            This function calculate the total length square.
            ==============================================  */
        double CtrlPtsTotalLenSquare(std::vector<coordinate> vects){
            double totalLen = 0;
            
            for (int i = 0; i < vects.size(); i++)
            {
                totalLen = totalLen + VectorLenSquare(vects[i].coord);
            }
            return totalLen;
        };


        /*  ==============================================
            The main optimization process
            ==============================================  */
        BezierImage ConjugateGradientOptimization(BezierImage bImage, std::vector<int> samples, double goalResidual, bool autoIncreateSamples = true, bool autoPerturbation = true, double minimalValidDeform = 1.0E-4)
        {

            double residual = EquationResidual(bImage, samples);
            double newResidual;
        
            double startDeform = 1.0E-15;
            double realDeformScale;
            double changeStepDeformScale = minimalValidDeform*0.01;
            int testStep = 0;
        
            std::vector<coordinate> gradient;
            std::vector<coordinate> previousGradient;
            double grandientNorm;
            double previousGrandientNorm = 1.0E+300;
            double conjugateGrandientFactor;
        
            coordinate zeroGradient;
            coordinate gradientFlow;
        
            for (int i = 0; i < bImage.TargetDimeion; i++)
            {
                zeroGradient.coord.push_back( 0.0 );
                gradientFlow.coord.push_back( 0.0 );
            }
        
            for (int cPt = 0; cPt < bImage.CtrlPtsTotalNum; cPt++)
            {
                previousGradient.push_back( zeroGradient );
            }
        
            while (residual > goalResidual)
            {
        
                testStep = testStep + 1;
                
                BezierImage deformedObj = BezierImage(bImage);
                
                realDeformScale = startDeform;
                
                std::cout << "testStep = " << testStep << ", #sample = " << samples[0]  << std::endl;
        
                bool isCalculate;
                std::array<int, 4> index;
        
                gradient.clear();
                for (int cPt = 0; cPt < bImage.CtrlPtsTotalNum; cPt++)
                {
        
                    isCalculate = false;
                    index = bImage.CtrlPtSerialNumtoIndex[cPt];
                        
                    for (int i = 0; i < bImage.DomainDimeion; i++)
                    {
                        if (index[i] < 3 || index[i]%3 == 0 || index[i]%3 == 2)
                        {
                            isCalculate = true;
                        }
                    }
                    if (isCalculate)
                    {
                        gradient.push_back( GradientEquationResidual(bImage, cPt, samples) );
                    }
                    else
                    {
                        gradient.push_back( zeroGradient );
                    }
        
                    printf("\r calculate gradient for control points [%d%%] ",cPt*100/(bImage.CtrlPtsTotalNum-1));
                }
        
                grandientNorm = CtrlPtsTotalLenSquare(gradient);
        
                conjugateGrandientFactor = grandientNorm / previousGrandientNorm;
        
                for (int cPt = 0; cPt < bImage.CtrlPtsTotalNum; cPt++)
                {
                    for (int i = 0; i < bImage.TargetDimeion; i++)
                    {
                        gradient[cPt].coord[i] = gradient[cPt].coord[i] + conjugateGrandientFactor*previousGradient[cPt].coord[i];
                    }
                }
        
                while (realDeformScale < 1.0)
                {
                    for (int cPt = 0; cPt < bImage.CtrlPtsTotalNum; cPt++)
                    {
                        for (int i = 0; i < bImage.TargetDimeion; i++)
                        {
                            gradientFlow.coord[i] = - realDeformScale*gradient[cPt].coord[i];
                        }
                        deformedObj.SmoothDeform(cPt, gradientFlow);
                    }
        
                    newResidual = EquationResidual(deformedObj, samples);
                    
                    if (newResidual < residual + 1.0E-8)  // The numerical error margin may allow to continue the calculation.
                    {
                        bImage = BezierImage(deformedObj);
                        residual = newResidual;
        
                        if (realDeformScale < changeStepDeformScale)
                        {
                            realDeformScale = realDeformScale * 10.0;
                        }
                        else
                        {
                            realDeformScale = realDeformScale * 1.5;
                        }
                    }
                    else
                    {
                        break;
                    }
                }
        
                previousGradient = gradient;
                previousGrandientNorm = grandientNorm;
        
                /* If the deformation is not far enough, we perturb the control points and restart the Conjugate Gradient algorithm. */
                if (realDeformScale < minimalValidDeform || realDeformScale > 1.0)
                {
                    if (autoPerturbation == true)
                    {
                        coordinate tempDeform;
                        for (int i = 0; i < bImage.TargetDimeion; i++)
                        {
                            tempDeform.coord.push_back(0.0);
                        }
            
                        for (int cPt = 0; cPt < bImage.CtrlPtsTotalNum; cPt++)
                        {
                            bool isBoundary = false;
                            for (int d = 0; d < bImage.DomainDimeion; d++)
                            {
                                if (bImage.CtrlPtSerialNumtoIndex[cPt][d]==0 || bImage.CtrlPtSerialNumtoIndex[cPt][d] == bImage.CtrlPtsNumInDAxis[d]-1)
                                {
                                    isBoundary = true;
                                }
                            }

                            if (isBoundary == false)
                            {
                                for (int i = 0; i < bImage.TargetDimeion; i++)
                                {
                                    tempDeform.coord[i] = 1.0E-6*(static_cast <double> (rand()) / static_cast <double> (RAND_MAX));
                                    previousGradient[cPt].coord[i] = 0.0;
                                }
                                bImage.SmoothDeform(cPt, tempDeform);
                            }
                        }
                    }

                    if (autoIncreateSamples == true)
                    {
                        for (int i = 0; i < samples.size(); i++)
                        {
                            samples[i] = samples[i] + (int) 2*(static_cast <double> (rand()) / static_cast <double> (RAND_MAX)) + 1;
                        }
                    }
        
                    previousGrandientNorm = 1.0E+300;
                    residual = EquationResidual(bImage, samples);
                }
        
                std::cout << ", realDeformScale = " << realDeformScale << ", residual = " << residual << std::endl;
        
            }

            return bImage;
        }


        /*  ==============================================
            The main optimization process with rescale to < maxsize
            ==============================================  */
        BezierImage ConjugateGradientOptimizationRescaled(BezierImage bImage, std::vector<int> samples, double goalResidual, bool autoIncreateSamples = true, bool autoPerturbation = true, double minimalValidDeform = 1.0E-4, double maxsize = 1.0)
        {
            double maxlength = 0.0;
            double residual = EquationResidual(bImage, samples);
            double newResidual;
        
            double startDeform = 1.0E-15;
            double realDeformScale;
            double changeStepDeformScale = minimalValidDeform*0.01;
            int testStep = 0;
        
            std::vector<coordinate> gradient;
            std::vector<coordinate> previousGradient;
            double grandientNorm;
            double previousGrandientNorm = 1.0E+300;
            double conjugateGrandientFactor;
        
            coordinate zeroGradient;
            coordinate gradientFlow;
        
            for (int i = 0; i < bImage.TargetDimeion; i++)
            {
                zeroGradient.coord.push_back( 0.0 );
                gradientFlow.coord.push_back( 0.0 );
            }
        
            for (int cPt = 0; cPt < bImage.CtrlPtsTotalNum; cPt++)
            {
                previousGradient.push_back( zeroGradient );
            }
        
            while (residual > goalResidual)
            {
                testStep = testStep + 1;
                
                BezierImage deformedObj = BezierImage(bImage);
                
                realDeformScale = startDeform;
                
                std::cout << "testStep = " << testStep << ", #sample = " << samples[0]  << std::endl;
        
                bool isCalculate;
                std::array<int, 4> index;
        
                gradient.clear();
                for (int cPt = 0; cPt < bImage.CtrlPtsTotalNum; cPt++)
                {
        
                    isCalculate = false;
                    index = bImage.CtrlPtSerialNumtoIndex[cPt];
                        
                    for (int i = 0; i < bImage.DomainDimeion; i++)
                    {
                        if (index[i] < 3 || index[i]%3 == 0 || index[i]%3 == 2)
                        {
                            isCalculate = true;
                        }
                    }
                    if (isCalculate)
                    {
                        gradient.push_back( GradientEquationResidual(bImage, cPt, samples) );
                    }
                    else
                    {
                        gradient.push_back( zeroGradient );
                    }
        
                    printf("\r calculate gradient for control points [%d%%] ",cPt*100/(bImage.CtrlPtsTotalNum-1));
                }
        
                grandientNorm = CtrlPtsTotalLenSquare(gradient);
        
                conjugateGrandientFactor = grandientNorm / previousGrandientNorm;
        
                for (int cPt = 0; cPt < bImage.CtrlPtsTotalNum; cPt++)
                {
                    for (int i = 0; i < bImage.TargetDimeion; i++)
                    {
                        gradient[cPt].coord[i] = gradient[cPt].coord[i] + conjugateGrandientFactor*previousGradient[cPt].coord[i];
                    }
                }
        
                while (realDeformScale < 1.0)
                {
                    for (int cPt = 0; cPt < bImage.CtrlPtsTotalNum; cPt++)
                    {
                        for (int i = 0; i < bImage.TargetDimeion; i++)
                        {
                            gradientFlow.coord[i] = - realDeformScale*gradient[cPt].coord[i];
                        }
                        deformedObj.SmoothDeform(cPt, gradientFlow);
                    }
        
                    newResidual = EquationResidual(deformedObj, samples);
                    
                    if (newResidual < residual + 1.0E-8)  // The numerical error margin may allow to continue the calculation.
                    {
                        bImage = BezierImage(deformedObj);
                        residual = newResidual;
        
                        if (realDeformScale < changeStepDeformScale)
                        {
                            realDeformScale = realDeformScale * 10.0;
                        }
                        else
                        {
                            realDeformScale = realDeformScale * 1.5;
                        }
                    }
                    else
                    {
                        break;
                    }
                }
                
                previousGradient = gradient;
                previousGrandientNorm = grandientNorm;

                maxlength = 0.0;
                for (int cPt = 0; cPt < bImage.CtrlPtsTotalNum; cPt++)
                {
                    for (int i = 0; i < bImage.TargetDimeion; i++)
                    {
                        if (maxlength < bImage.CtrlPoint[cPt].coord[i])
                        {
                            maxlength = bImage.CtrlPoint[cPt].coord[i];
                        }                        
                    }
                }

                if (maxlength > maxsize)
                {
                    double ratio = maxsize / maxlength;
                    previousGrandientNorm = previousGrandientNorm * ratio;

                    for (int cPt = 0; cPt < bImage.CtrlPtsTotalNum; cPt++)
                    {
                        for (int i = 0; i < bImage.TargetDimeion; i++)
                        {
                            previousGradient[cPt].coord[i] = previousGradient[cPt].coord[i] * ratio;
                            bImage.CtrlPoint[cPt].coord[i] = bImage.CtrlPoint[cPt].coord[i] * ratio;
                        }
                    }
                }

                /* If the deformation is not far enough, we perturb the control points and restart the Conjugate Gradient algorithm. */
                if (realDeformScale < minimalValidDeform || realDeformScale > 1.0)
                {
                    if (autoPerturbation == true)
                    {
                        coordinate tempDeform;
                        for (int i = 0; i < bImage.TargetDimeion; i++)
                        {
                            tempDeform.coord.push_back(0.0);
                        }
            
                        for (int cPt = 0; cPt < bImage.CtrlPtsTotalNum; cPt++)
                        {
                            bool isBoundary = false;
                            for (int d = 0; d < bImage.DomainDimeion; d++)
                            {
                                if (bImage.CtrlPtSerialNumtoIndex[cPt][d]==0 || bImage.CtrlPtSerialNumtoIndex[cPt][d] == bImage.CtrlPtsNumInDAxis[d]-1)
                                {
                                    isBoundary = true;
                                }
                            }

                            if (isBoundary == false)
                            {
                                for (int i = 0; i < bImage.TargetDimeion; i++)
                                {
                                    tempDeform.coord[i] = 1.0E-6*(static_cast <double> (rand()) / static_cast <double> (RAND_MAX));
                                    previousGradient[cPt].coord[i] = 0.0;
                                }
                                bImage.SmoothDeform(cPt, tempDeform);
                            }
                        }
                    }

                    if (autoIncreateSamples == true)
                    {
                        for (int i = 0; i < samples.size(); i++)
                        {
                            samples[i] = samples[i] + (int) 2*(static_cast <double> (rand()) / static_cast <double> (RAND_MAX)) + 1;
                        }
                    }
        
                    previousGrandientNorm = 1.0E+300;
                    residual = EquationResidual(bImage, samples);
                }
        
                std::cout << ", realDeformScale = " << realDeformScale << ", residual = " << residual << std::endl;
        
            }

            return bImage;
        }

};