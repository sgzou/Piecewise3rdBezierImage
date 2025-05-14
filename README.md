# Piecewise3rdBezierImage

## Approximate solutions of integro-differential equations by piecewise 3rd Bezier functions and conjugate gradient method

This is an algorithm for finding the numerical solution of the general (2nd order) integro-differential equation by the piecewise polynomial (3rd Bezier) functions.
The main improvements here are the consideration of the piecewise 3rd Bezier functions with $C^1$ continuity and the use of the conjugate gradient method for Bezier control points to find the approximate solution of integro-differential equations.

See the [document](/doc/Piecewise_3rd_Bezier_Image_approximation_method.pdf) for detail mathematical idea.

## Usage

1. Include "OptimizeOperation.h" and inherit the OptimizeOperation class, then override two methods EquationResidual() and GradientEquationResidual()
```c++
#include "../OptimizeOperation.h"

class SystemCalculator : public OptimizeOperation
{
    public:
      double EquationResidual(BezierImage img, std::vector<int> samples) override { ... }
      coordinate GradientEquationResidual(BezierImage img, int cpt, std::vector<int> samples) override { ... }
```
Note, the `BezierImage` class function `BezierMap(t)` always takes a parameter `t` of type `std::array<double, 4>`.
If the parameter is one-dimensional, we will just use the first index of this array, and the same applies to higher dimensions.
Therefore, the maximum dimension of the domain space is 4.

2. Create a inintial Bezer function
```c++
BezierImage initImg;
initImg.ReadFile("parameter.txt"); 
//initImg.RandGenerator(cPtsN, fixPts, interval); // see the example of knot_energy
```
where the structure of the parameter file "parameter.txt" is
```
1 3
0 0.2
0 10 10 10 
30 10 10 10 
```
they are
```
DomainSpaceDimension TargetSpaceDimension
DomainSpaceStartPointOfFirstDimension DomainSpaceEndPointOfFirstDimension
...
ZerothControlPointIndex .. ZerothControlPointPositionInTargetSpace ..
FirstControlPointIndex .. FirstControlPointPositionInTargetSpace ..
SecondControlPointIndex .. SecondControlPointPositionInTargetSpace ..
...
```
If some lost control points are lost, the function will automatically fill them in using information from previous points.

3. Implement and call ConjugateGradientOptimization
```c++
SystemCalculator calculator;
BezierImage optimizedResult = calculator.ConjugateGradientOptimization(initImg, samples, 1.0);
```

## Lorenz system example
![example fig](/lorenz_system/data/lorenz_system_ex.png)

## O'Hara knot energy example
![example fig](/knot_energy/data/knot_energy_ex.png)

## Wave equation example
![example fig](/wave_equation/data/wave_equation_ex.png)
