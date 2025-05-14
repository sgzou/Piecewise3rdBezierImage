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
        BezierImage ConjugateGradientOptimization(BezierImage bImage, std::vector<int> samples, double goalResidual, double extraSetAtSmallDeform = 1.0E-5, bool autoPerturbation = true, double pertubationSize = 1.0E-6, bool autoIncreateSamples = true, bool autoRescale = false, double rescaleMaxSize = 1.0)
        {
            const int MAX_ITERATIONS = 1000; 
            double maxlength = 0.0;
            double residual = EquationResidual(bImage, samples);
            double newResidual;
        
            double startDeform = 1.0E-15;
            double realDeformScale;
            double changeStepDeformScale = extraSetAtSmallDeform*0.01;
            int iterationCount = 0;
        
            std::vector<coordinate> gradient;
            std::vector<coordinate> previousGradient;
            double gradientNorm;
            double previousGradientNorm = 1.0E+300;
            double conjugateGradientFactor;
        
            coordinate zeroGradient;
            coordinate gradientFlow;
        
            for (int i = 0; i < bImage.TargetDimension; i++)
            {
                zeroGradient.coord.push_back( 0.0 );
                gradientFlow.coord.push_back( 0.0 );
            }
        
            for (int cPt = 0; cPt < bImage.ControlPointsTotalNumber; cPt++)
            {
                previousGradient.push_back( zeroGradient );
            }
        
            while (residual > goalResidual && iterationCount < MAX_ITERATIONS)
            {
                iterationCount = iterationCount + 1;
                
                BezierImage deformedObj = BezierImage(bImage);
                
                realDeformScale = startDeform;
                
                std::cout << "iterationCount = " << iterationCount << ", #sample = " << samples[0]  << std::endl;
        
                bool isCalculate;
                std::array<int, 4> index;
        
                gradient.clear();
                gradient.reserve(bImage.ControlPointsTotalNumber);
                for (int cPt = 0; cPt < bImage.ControlPointsTotalNumber; cPt++)
                {
        
                    isCalculate = false;
                    index = bImage.ControlPointSerialNumToIndex[cPt];
                        
                    for (int i = 0; i < bImage.DomainDimension; i++)
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
        
                    printf("\r calculate gradient for control points [%d%%] ",cPt*100/(bImage.ControlPointsTotalNumber-1));
                }
        
                gradientNorm = CtrlPtsTotalLenSquare(gradient);
        
                conjugateGradientFactor = gradientNorm / previousGradientNorm;
        
                for (int cPt = 0; cPt < bImage.ControlPointsTotalNumber; cPt++)
                {
                    for (int i = 0; i < bImage.TargetDimension; i++)
                    {
                        gradient[cPt].coord[i] = gradient[cPt].coord[i] + conjugateGradientFactor*previousGradient[cPt].coord[i];
                    }
                }
        
                while (realDeformScale < 1.0)
                {
                    for (int cPt = 0; cPt < bImage.ControlPointsTotalNumber; cPt++)
                    {
                        for (int i = 0; i < bImage.TargetDimension; i++)
                        {
                            gradientFlow.coord[i] = - realDeformScale*gradient[cPt].coord[i];
                        }
                        deformedObj.SmoothDeform(cPt, gradientFlow);
                    }
        
                    newResidual = EquationResidual(deformedObj, samples);
                    
                    if (newResidual < residual + 1.0E-10)  // The numerical error margin may allow to continue the calculation.
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
                previousGradientNorm = gradientNorm;

                if (autoRescale == true)
                {
                    maxlength = 0.0;
                    for (int cPt = 0; cPt < bImage.ControlPointsTotalNumber; cPt++)
                    {
                        for (int i = 0; i < bImage.TargetDimension; i++)
                        {
                            if (maxlength < bImage.ControlPoints[cPt].coord[i])
                            {
                                maxlength = bImage.ControlPoints[cPt].coord[i];
                            }                        
                        }
                    }

                    if (maxlength > rescaleMaxSize)
                    {
                        double ratio = rescaleMaxSize / maxlength;
                        previousGradientNorm = previousGradientNorm * ratio;

                        for (int cPt = 0; cPt < bImage.ControlPointsTotalNumber; cPt++)
                        {
                            for (int i = 0; i < bImage.TargetDimension; i++)
                            {
                                previousGradient[cPt].coord[i] = previousGradient[cPt].coord[i] * ratio;
                                bImage.ControlPoints[cPt].coord[i] = bImage.ControlPoints[cPt].coord[i] * ratio;
                            }
                        }
                    }
                }

                /* If the deformation of Conjugate Gradient method can not change enough, we perturb the control points and restart the Conjugate Gradient algorithm. */
                if (realDeformScale < extraSetAtSmallDeform || realDeformScale > 1.0)
                {
                    if (autoPerturbation == true)
                    {
                        coordinate tempDeform;
                        tempDeform.coord.reserve(bImage.TargetDimension);
                        for (int i = 0; i < bImage.TargetDimension; i++)
                        {
                            tempDeform.coord.push_back(0.0);
                        }
            
                        for (int cPt = 0; cPt < bImage.ControlPointsTotalNumber; cPt++)
                        {
                            bool isBoundary = false;
                            for (int d = 0; d < bImage.DomainDimension; d++)
                            {
                                if (bImage.ControlPointSerialNumToIndex[cPt][d]==0 || bImage.ControlPointSerialNumToIndex[cPt][d] == bImage.ControlPointsNumberInDomainAxis[d]-1)
                                {
                                    isBoundary = true;
                                }
                            }

                            if (isBoundary == false)
                            {
                                for (int i = 0; i < bImage.TargetDimension; i++)
                                {
                                    tempDeform.coord[i] = pertubationSize*(static_cast <double> (rand()) / static_cast <double> (RAND_MAX));
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
                            samples[i] = samples[i] + 1;
                        }
                    }
        
                    previousGradientNorm = 1.0E+300;
                    residual = EquationResidual(bImage, samples);
                }
        
                std::cout << ", realDeformScale = " << realDeformScale << ", residual = " << residual << std::endl;

                if (iterationCount >= MAX_ITERATIONS) {  
                    std::cout << "Warning: Reached maximum iterations without convergence." << std::endl;
                }
            }

            return bImage;
        }

};