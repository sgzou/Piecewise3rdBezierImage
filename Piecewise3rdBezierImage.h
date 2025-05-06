#include <functional>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <map>
#include <array>
#include <chrono>
#include <thread>

class coordinate
{
    public:
        std::vector<double> coord;
};

/* Bernstein Polynomial */
static double Bernstein(int i, double start, double end, double t){
    switch (i)
    {
    case 0:
        return pow((end-t), 3);
    case 1:
        return 3.0*pow((end-t), 2)*(t-start);
    case 2:
        return 3.0*(end-t)*pow((t-start), 2);
    case 3:
      return 1.0*pow((t-start), 3);
    default:
        return 0.0;
    }
}

/* differential Bernstein Polynomial */
static double DiffBernstein(int i, double start, double end, double t){
    switch (i)
    {
    case 0:
        return -3.0*pow((end-t), 2);
    case 1:
        return -6.0*(end-t)*(t-start) + 3.0*pow((end-t), 2);
    case 2:
        return -3.0*pow((t-start), 2) + 6.0*(end-t)*(t-start);
    case 3:
      return 3.0*pow((t-start), 2);
    default:
        return 0.0;
    }
}

/* double differential Bernstein Polynomial */
static double DoubleDiffBernstein(int i, double start, double end, double t){
    switch (i)
    {
    case 0:
        return 6.0*(end-t);
    case 1:
        return 6.0*(t-start) - 12.0*(end-t);
    case 2:
        return -12.0*(t-start) + 6.0*(end-t);
    case 3:
      return 6.0*(t-start);
    default:
        return 0.0;
    }
}

/*
    Bezier geometric object : f(u,v,w) = P_{ijk} B_i(u) B_j(v) B_k(w) , where P_{ijk} are control points and B_i(u) is Bernstein polynomials.
    differentiable condition : P_{(3i+1)} = 2P_{(3i)} - P_{(3i-1)}
*/
class BezierImage
{
    public:
        /*  ==============================================
            Basic information of spaces and control points 
            ==============================================  */
        int TargetDimeion;                     // The dimension of target space
        int DomainDimeion;                     // The dimension of domain space
        int CtrlPtsTotalNum;                   // The total number of control points
        std::array<int,4> CtrlPtsNumInDAxis;   // The number of control points in domain axis
        std::vector<coordinate> CtrlPoint;     // The control point set (the order corresponding to the serial number), and their coordinate in target space

        /*  ==============================================
            Detail information for control point index and domain mesh
            ==============================================  */
        std::vector<std::array<int,4>> CtrlPtSerialNumtoIndex;      // Record the index (ijk) of the n-th control point P_{ijk}
        std::map<std::array<int,4>, int> CtrlPtIndextoSerialNum;    // The dictionary for index (ijk) of control point P_{ijk} map to its serial number

        std::array<std::array<double, 2>, 4> DomainInterval;        // Domain interval in individual dimension
        std::array<double, 4> DomainMeshSpacing;                    // The space size of domain in individual dimension
        std::array<std::vector<double>, 4> DomainMeshZoneStart;     // Record the start point in domain for the mesh zones in individual dimension and zone ; DomainMeshZoneStart[0][3] is the start point value in the first dimension and the third zone


        /*  ==============================================
            Generate the informations for Bezier image via few information
            ==============================================  */
        void Generator( std::vector<std::array<double,2>> domainIntvl, std::map<std::array<int,4>, coordinate> cPoint ){

            TargetDimeion = cPoint.begin()->second.coord.size();
            DomainDimeion = domainIntvl.size();
            for (int i = 0; i < 4; i++)
            {
                DomainInterval[i] = {0.0, 0.0};
            }
            for (int i = 0; i < domainIntvl.size(); i++)
            {
                DomainInterval[i][0] = domainIntvl[i][0];
                DomainInterval[i][1] = domainIntvl[i][1];
            }
            
            CtrlPtsTotalNum = cPoint.size();
            CtrlPtsNumInDAxis = {{1,1,1,1}};
            int index=0;
            for(std::map<std::array<int,4>, coordinate>::iterator it = cPoint.begin(); it != cPoint.end(); ++it) {
                std::array<int,4> tempIndex = it->first;
                for (int i = 0; i < DomainDimeion; i++)
                {
                    if (tempIndex[i]+1 > CtrlPtsNumInDAxis[i])
                    {
                        CtrlPtsNumInDAxis[i] = tempIndex[i]+1;
                    }
                }
                CtrlPtSerialNumtoIndex.push_back(tempIndex);
                CtrlPtIndextoSerialNum[tempIndex]=index;
                ++index;
                CtrlPoint.push_back(it->second);
            }
            for (int i = 0; i < 4; i++)
            {
                // Check the number of control points in each dimension is 3N+1. If not, change it to 3N+1, later will be filled the missing points.
                int remnant = CtrlPtsNumInDAxis[i]%3;
                if (remnant != 1)
                {
                    CtrlPtsNumInDAxis[i] = CtrlPtsNumInDAxis[i] + (2*remnant+1)%3;
                }

                // Set the domain mesh information
                int latticeNum = (int)CtrlPtsNumInDAxis[i]/3;
                std::vector<double> temp;
                DomainMeshSpacing[i] = abs(domainIntvl[i][1]-domainIntvl[i][0])/latticeNum;
                for (int j = 0; j < latticeNum; j++)
                {
                    temp.push_back(domainIntvl[i][0] + (double)j*DomainMeshSpacing[i]);
                }
                DomainMeshZoneStart[i] = temp;
            }

            // fill the missing control points
            for (int d0 = 0; d0 < CtrlPtsNumInDAxis[0]; d0++)
            {
                for (int d1 = 0; d1 < CtrlPtsNumInDAxis[1]; d1++)
                {
                    for (int d2 = 0; d2 < CtrlPtsNumInDAxis[2]; d2++)
                    {
                        for (int d3 = 0; d3 < CtrlPtsNumInDAxis[3]; d3++)
                        {
                            // key is present, it returns 1, otherwise, it returns 0
                            if (CtrlPtIndextoSerialNum.count({d0,d1,d2,d3}) == 0)
                            {
                                CtrlPtIndextoSerialNum[{d0,d1,d2,d3}] = CtrlPtsTotalNum;
                                CtrlPtSerialNumtoIndex.push_back({d0,d1,d2,d3});
                                // use the previous point target coordinate to fill the missing point
                                coordinate temp( CtrlPoint[CtrlPtIndextoSerialNum[{d0-1,d1,d2,d3}]] );
                                for (int i = 0; i < TargetDimeion; i++)
                                {
                                    temp.coord[i] = temp.coord[i];
                                }                                
                                CtrlPoint.push_back(temp);
                                CtrlPtsTotalNum = CtrlPtSerialNumtoIndex.size();
                            }
                        }
                    }
                }
            }

            Smoothize();
        };


        /*  ==============================================
            Smoothize by the differentiable condition P[3i+1] = 2P[3i] - P[3i-1]
            ==============================================  */
        void Smoothize(){
            std::array<int,4> cptIndex;
            std::array<int,4> temp1;
            std::array<int,4> temp2;
            std::array<int,4> temp3;
            std::array<int,4> temp4;

            for (int pt = 0; pt < CtrlPtsTotalNum; pt++)
            {
                for (int cr = 0; cr < TargetDimeion; cr++)
                {
                    cptIndex = CtrlPtSerialNumtoIndex[pt];
                    for (int i = 0; i < DomainDimeion; i++)
                    {
                        if (cptIndex[i]%3==1 && cptIndex[i]!=1)
                        {
                            std::copy(cptIndex.begin(), cptIndex.end(), temp1.begin());
                            temp1[i] = temp1[i]-1;
                            std::copy(cptIndex.begin(), cptIndex.end(), temp2.begin());
                            temp2[i] = temp2[i]-2;

                            CtrlPoint[pt].coord[cr] = 2.0*CtrlPoint[ CtrlPtIndextoSerialNum[temp1] ].coord[cr] - 1.0*CtrlPoint[ CtrlPtIndextoSerialNum[temp2] ].coord[cr];
                        }
                    }
                }
            }
        }


        /*  ==============================================
            deform the control point cPt in the specific vector, and impose the differentiable condition P[3i+1] = 2P[3i] - P[3i-1]
            ==============================================  */
        void SmoothDeform( int cPt, coordinate deformVect ){
            std::array<int,4> temp1p;
            std::array<int,4> temp2p;
            std::array<int,4> temp1m;
            std::array<int,4> temp2m;

            std::array<int,4> cptIndex = CtrlPtSerialNumtoIndex[cPt];

            for (int j = 0; j < deformVect.coord.size(); j++)
            {
                CtrlPoint[cPt].coord[j] = CtrlPoint[cPt].coord[j] + deformVect.coord[j];
            }

            for (int i = 0; i < DomainDimeion; i++)
            {
                if (cptIndex[i]%3==1 && cptIndex[i]!=1 && cptIndex[i] < CtrlPtsNumInDAxis[i]-2)
                {
                    std::copy(cptIndex.begin(), cptIndex.end(), temp1m.begin());
                    temp1m[i] = temp1m[i]-1;
                    std::copy(cptIndex.begin(), cptIndex.end(), temp2m.begin());
                    temp2m[i] = temp2m[i]-2;

                    for (int j = 0; j < deformVect.coord.size(); j++)
                    {
                        CtrlPoint[ CtrlPtIndextoSerialNum[temp1m] ].coord[j] = 0.5*( CtrlPoint[ CtrlPtIndextoSerialNum[temp2m] ].coord[j] + CtrlPoint[cPt].coord[j] );
                    }
                }

                if (cptIndex[i]%3==2 && cptIndex[i] < CtrlPtsNumInDAxis[i]-3)
                {
                    std::copy(cptIndex.begin(), cptIndex.end(), temp1p.begin());
                    temp1p[i] = temp1p[i]+1;
                    std::copy(cptIndex.begin(), cptIndex.end(), temp2p.begin());
                    temp2p[i] = temp2p[i]+2;

                    for (int j = 0; j < deformVect.coord.size(); j++)
                    {
                        // Changing the nearest point is better
                        CtrlPoint[ CtrlPtIndextoSerialNum[temp1p] ].coord[j] = 0.5*( CtrlPoint[ CtrlPtIndextoSerialNum[temp2p] ].coord[j] + CtrlPoint[cPt].coord[j] );
                    }
                }

                if (cptIndex[i]%3==0 && cptIndex[i]!=0 && cptIndex[i] < CtrlPtsNumInDAxis[i]-2)
                {
                    std::copy(cptIndex.begin(), cptIndex.end(), temp1m.begin());
                    temp1m[i] = temp1m[i]-1;
                    std::copy(cptIndex.begin(), cptIndex.end(), temp1p.begin());
                    temp1p[i] = temp1p[i]+1;

                    for (int j = 0; j < deformVect.coord.size(); j++)
                    {
                        CtrlPoint[ CtrlPtIndextoSerialNum[temp1p] ].coord[j] = 2.0*CtrlPoint[cPt].coord[j] - CtrlPoint[ CtrlPtIndextoSerialNum[temp1m] ].coord[j];
                    }
                }
            }
        }


        /*  ==============================================
            Read file
            ==============================================  */
        void ReadFile( std::string fileName ){
            std::ifstream myfile(fileName);
            std::string myline;
            std::string temp;
            
            int targetDim;
            int domainDim;
            std::map<std::array<int,4>, coordinate> pts;
            std::vector<std::array<double,2>> vregs;

            if ( myfile.is_open() ) {
                std::vector<int> dim;
                std::getline(myfile, myline);
                std::stringstream firstline(myline);
                while ( std::getline(firstline, temp, ' ') ) {
                    dim.push_back(std::stod(temp));
                }
                if (dim.size() <2)
                {
                    std::cout << "this file formation is not correct." << std::endl;
                }
                
                domainDim = dim[0];
                targetDim = dim[1];

                for (int i = 0; i < domainDim; i++)
                {
                    std::vector<double> regs;
                    std::getline(myfile, myline);
                    std::stringstream regiondata(myline);
                    while ( std::getline(regiondata, temp, ' ') ) {
                        regs.push_back(std::stod(temp));
                    }
                    vregs.push_back({regs[0], regs[1]});
                }                

                while ( myfile ) {
                    std::vector<double> data;

                    std::getline(myfile, myline);
                    std::stringstream ss(myline);
                    while ( std::getline(ss, temp, ' ') ) {
                        data.push_back(std::stod(temp));
                    }
                    if ( data.size()<2 ) { break; }

                    std::vector<double> val(targetDim);
                    std::copy(data.begin()+domainDim, data.end(), val.begin());

                    std::array<int,4> key={{0,0,0,0}};
                    for (int i = 0; i < domainDim; i++)
                    {
                        key[i] = (int) data[i];
                    }

                    pts[key] = coordinate{ val };
                }

                Generator(vregs, pts);
            }
            else
            {
                std::cout << "read error; please check the file name and path !" << std::endl;
            }
        };


        /*  ==============================================
            Random generate the control points with some specific fixed control point
            ==============================================  */
        void RandGenerator(std::vector<int> cPtN, std::map<std::array<int,4>, coordinate> fixedPt, std::vector<std::array<double,2>> domainIntvl){
            
            std::map<std::array<int,4>, coordinate> points(fixedPt);
            int totalPtN = 1;
            int dim = domainIntvl.size();
            std::array<int,4> cPtsNofDim = {{1,1,1,1}};

            // Let the number of control points in each dimension be 3N+1
            for (int i = 0; i < dim; i++)
            {
                cPtsNofDim[i] = cPtN[i] + (2*(cPtN[i]%3)+1)%3;
                totalPtN = totalPtN * cPtsNofDim[i];
            }

            int targetDim = fixedPt.begin()->second.coord.size();
            srand(time(0));
            // Generate the remaining points with random numbers
            for (int d0 = 0; d0 < cPtsNofDim[0]; d0++)
            {
                for (int d1 = 0; d1 < cPtsNofDim[1]; d1++)
                {
                    for (int d2 = 0; d2 < cPtsNofDim[2]; d2++)
                    {
                        for (int d3 = 0; d3 < cPtsNofDim[3]; d3++)
                        {
                            // key is present, it returns 1, otherwise, it returns 0
                            if (points.count({d0,d1,d2,d3}) == 0)
                            {
                                coordinate pt;
                                for (int i = 0; i < targetDim; i++)
                                {                                    
                                    pt.coord.push_back( static_cast <double> (rand()) / static_cast <double> (RAND_MAX) );                                    
                                }
                                points[{d0,d1,d2,d3}] = pt;
                            }
                        }
                    }
                }
            }
            Generator(domainIntvl, points);
        }
        

        /*  ==============================================
            Bezier map function
            ==============================================  */
        coordinate BezierMap( std::array<double,4> t ){
            coordinate result;
            std::array<int,4> startLatticeIndex = {{0,0,0,0}};
            std::array<double,4> startLatticePt = {{-1.0, -1.0, -1.0, -1.0}};
            std::array<double,4> endLatticePt = {{1.0, 1.0, 1.0, 1.0}};
            std::array<int,4> sCtrlPt = {{0,0,0,0}};
            std::array<int,4> lRange = {{1,1,1,1}};
            double latticeSpacing = 1.0;
            int pt;
            double productTemp;

            for (int i = 0; i < TargetDimeion; i++)
            {
                result.coord.push_back(0.0);
            }
            
            for (int i = 0; i < DomainDimeion; i++)
            {
                startLatticeIndex[i] = (int)((t[i]-DomainInterval[i][0])/(DomainMeshSpacing[i] + 1.0E-14)); // The double float point error is 10^-16
                startLatticePt[i] = DomainMeshZoneStart[i][startLatticeIndex[i]];
                endLatticePt[i] = startLatticePt[i] + DomainMeshSpacing[i];
                latticeSpacing = latticeSpacing * DomainMeshSpacing[i];
                sCtrlPt[i] = 3*startLatticeIndex[i];    // the first control point in the mesh zone
                lRange[i] = 4;
            }

            for (int d0 = 0; d0 < lRange[0]; d0++)
            {
                for (int d1 = 0; d1 < lRange[1]; d1++)
                {
                    for (int d2 = 0; d2 < lRange[2]; d2++)
                    {
                        for (int d3 = 0; d3 < lRange[3]; d3++)
                        {
                            pt = CtrlPtIndextoSerialNum[{d0 + sCtrlPt[0], d1 + sCtrlPt[1], d2 + sCtrlPt[2], d3 + sCtrlPt[3]}];
                            productTemp = Bernstein(d0, startLatticePt[0], endLatticePt[0], t[0])
                                    *Bernstein(d1, startLatticePt[1], endLatticePt[1], t[1])
                                    *Bernstein(d2, startLatticePt[2], endLatticePt[2], t[2])
                                    *Bernstein(d3, startLatticePt[3], endLatticePt[3], t[3]);
                            
                            for (int j = 0; j < TargetDimeion; j++)
                            {
                                result.coord[j] = result.coord[j] + CtrlPoint[pt].coord[j]*productTemp;
                            }
                        }
                    }
                }
            }

            productTemp = pow(latticeSpacing, 3);

            for (int j = 0; j < TargetDimeion; j++)
            {
                result.coord[j] = result.coord[j] / productTemp;
            }

            return result;
        };


        /*  ==============================================
            differential Bezier map function
            ==============================================  */
        coordinate DiffBezierMap( std::array<double,4> t, int k ){
            coordinate result;
            std::array<int,4> startLatticeIndex = {{0,0,0,0}};
            std::array<double,4> startLatticePt = {{-1.0, -1.0, -1.0, -1.0}};
            std::array<double,4> endLatticePt = {{1.0, 1.0, 1.0, 1.0}};
            std::array<int,4> sCtrlPt = {{0,0,0,0}};
            std::array<int,4> lRange = {{1,1,1,1}};
            double latticeSpacing = 1.0;
            int pt;
            double productTemp;

            for (int i = 0; i < TargetDimeion; i++)
            {
                result.coord.push_back(0.0);
            }

            for (int i = 0; i < DomainDimeion; i++)
            {
                startLatticeIndex[i] = (int)((t[i]-DomainInterval[i][0])/(DomainMeshSpacing[i] + 1.0E-14));
                startLatticePt[i] = DomainMeshZoneStart[i][startLatticeIndex[i]];
                endLatticePt[i] = startLatticePt[i] + DomainMeshSpacing[i];
                latticeSpacing = latticeSpacing * DomainMeshSpacing[i];
                sCtrlPt[i] = 3*startLatticeIndex[i];    // the first control point in the mesh zone
                lRange[i] = 4;
            }

            std::function<double(int, double, double, double)> func0 = Bernstein;
            std::function<double(int, double, double, double)> func1 = Bernstein;
            std::function<double(int, double, double, double)> func2 = Bernstein;
            std::function<double(int, double, double, double)> func3 = Bernstein;

            switch ( k )
            {
                case 0:
                    func0 = DiffBernstein;
                    break;
                case 1:
                    func1 = DiffBernstein;
                    break;
                case 2:
                    func2 = DiffBernstein;
                    break;
                case 3:
                    func3 = DiffBernstein;
                    break;
                default:
                    break;
            }


            for (int d0 = 0; d0 < lRange[0]; d0++)
            {
                for (int d1 = 0; d1 < lRange[1]; d1++)
                {
                    for (int d2 = 0; d2 < lRange[2]; d2++)
                    {
                        for (int d3 = 0; d3 < lRange[3]; d3++)
                        {
                            pt = CtrlPtIndextoSerialNum[{d0 + sCtrlPt[0], d1 + sCtrlPt[1], d2 + sCtrlPt[2], d3 + sCtrlPt[3]}];
                            productTemp = func0(d0, startLatticePt[0], endLatticePt[0], t[0])
                                    *func1(d1, startLatticePt[1], endLatticePt[1], t[1])
                                    *func2(d2, startLatticePt[2], endLatticePt[2], t[2])
                                    *func3(d3, startLatticePt[3], endLatticePt[3], t[3]);

                            for (int j = 0; j < TargetDimeion; j++)
                            {
                                result.coord[j] = result.coord[j] + CtrlPoint[pt].coord[j]*productTemp;
                            }
                        }
                    }
                }
            }

            productTemp = pow(latticeSpacing, 3);

            for (int j = 0; j < TargetDimeion; j++)
            {
                result.coord[j] = result.coord[j] / productTemp;
            }

            return result;
        };


        /*  ==============================================
            Double differential Bezier map function
            ==============================================  */
        coordinate DoubleDiffBezierMap( std::array<double,4> t, int k, int l){
            coordinate result;
            std::array<int,4> startLatticeIndex = {{0,0,0,0}};
            std::array<double,4> startLatticePt = {{-1.0, -1.0, -1.0, -1.0}};
            std::array<double,4> endLatticePt = {{1.0, 1.0, 1.0, 1.0}};
            std::array<int,4> sCtrlPt = {{0,0,0,0}};
            std::array<int,4> lRange = {{1,1,1,1}};
            double latticeSpacing = 1.0;
            int pt;
            double productTemp;

            for (int i = 0; i < TargetDimeion; i++)
            {
                result.coord.push_back(0.0);
            }

            for (int i = 0; i < DomainDimeion; i++)
            {
                startLatticeIndex[i] = (int)((t[i]-DomainInterval[i][0])/(DomainMeshSpacing[i] + 1.0E-14));
                startLatticePt[i] = DomainMeshZoneStart[i][startLatticeIndex[i]];
                endLatticePt[i] = startLatticePt[i] + DomainMeshSpacing[i];
                latticeSpacing = latticeSpacing * DomainMeshSpacing[i];
                sCtrlPt[i] = 3*startLatticeIndex[i];    // the first control point in the mesh zone
                lRange[i] = 4;
            }

            std::function<double(int, double, double, double)> func0 = Bernstein;
            std::function<double(int, double, double, double)> func1 = Bernstein;
            std::function<double(int, double, double, double)> func2 = Bernstein;
            std::function<double(int, double, double, double)> func3 = Bernstein;

            switch( k )
            {
                case 0:
                    switch( l )
                    {
                        case 0:
                            func0 = DoubleDiffBernstein;
                            break;
                        case 1:
                            func0 = DiffBernstein;
                            func1 = DiffBernstein;
                            break;
                        case 2:
                            func0 = DiffBernstein;
                            func2 = DiffBernstein;
                            break;
                        case 3:
                            func0 = DiffBernstein;
                            func3 = DiffBernstein;
                            break;
                        default:
                            break;
                    }
                    break;
                case 1:
                    switch( l )
                    {
                        case 0:
                            func1 = DiffBernstein;
                            func0 = DiffBernstein;
                            break;
                        case 1:
                            func1 = DoubleDiffBernstein;
                            break;
                        case 2:
                            func1 = DiffBernstein;
                            func2 = DiffBernstein;
                            break;
                        case 3:
                            func1 = DiffBernstein;
                            func3 = DiffBernstein;
                            break;
                        default:
                            break;
                    }
                    break;
                case 2:
                    switch( l )
                    {
                        case 0:
                            func2 = DiffBernstein;
                            func0 = DiffBernstein;
                            break;
                        case 1:
                            func2 = DiffBernstein;
                            func1 = DiffBernstein;
                            break;
                        case 2:
                            func2 = DoubleDiffBernstein;
                            break;
                        case 3:
                            func2 = DiffBernstein;
                            func3 = DiffBernstein;
                            break;
                        default:
                            break;
                    }
                    break;
                case 3:
                    switch( l )
                    {
                        case 0:
                            func3 = DiffBernstein;
                            func0 = DiffBernstein;
                            break;
                        case 1:
                            func3 = DiffBernstein;
                            func1 = DiffBernstein;
                            break;
                        case 2:
                            func3 = DiffBernstein;
                            func2 = DiffBernstein;
                            break;
                        case 3:
                            func3 = DoubleDiffBernstein;
                            break;
                        default:
                            break;
                    }
                    break;
                default:
                    break;
            }

            for (int d0 = 0; d0 < lRange[0]; d0++)
            {
                for (int d1 = 0; d1 < lRange[1]; d1++)
                {
                    for (int d2 = 0; d2 < lRange[2]; d2++)
                    {
                        for (int d3 = 0; d3 < lRange[3]; d3++)
                        {
                            pt = CtrlPtIndextoSerialNum[{d0 + sCtrlPt[0], d1 + sCtrlPt[1], d2 + sCtrlPt[2], d3 + sCtrlPt[3]}];
                            productTemp = func0(d0, startLatticePt[0], endLatticePt[0], t[0])
                                    *func1(d1, startLatticePt[1], endLatticePt[1], t[1])
                                    *func2(d2, startLatticePt[2], endLatticePt[2], t[2])
                                    *func3(d3, startLatticePt[3], endLatticePt[3], t[3]);

                            for (int j = 0; j < TargetDimeion; j++)
                            {
                                result.coord[j] = result.coord[j] + CtrlPoint[pt].coord[j]*productTemp;
                            }
                        }
                    }
                }
            }

            productTemp = pow(latticeSpacing, 3);

            for (int j = 0; j < TargetDimeion; j++)
            {
                result.coord[j] = result.coord[j] / productTemp;
            }

            return result;
        };


        /*  ==============================================
            the partial derivative of Bezier function w.r.t. control point cpt, i.e. (d f) / (d cpt)
            ==============================================  */
        double CPtFactorOfBezierMap(int cpt, std::array<double,4> t, double correctionFactor = 1.0){
            double result = 0.0;
            std::array<int,4> startLatticeIndex = {{0,0,0,0}};
            std::array<double,4> startLatticePt = {{-1.0, -1.0, -1.0, -1.0}};
            std::array<double,4> endLatticePt = {{1.0, 1.0, 1.0, 1.0}};
            std::array<int,4> sCtrlPt = {{0,0,0,0}};
            double latticeSpacing = 1.0;
            int pt;
            std::array<int, 4> index = CtrlPtSerialNumtoIndex[cpt];
            std::array<int, 4> tempIndex;
            int d0;
            int d1;
            int d2;
            int d3;

            for (int i = 0; i < DomainDimeion; i++)
            {
                startLatticeIndex[i] = (int)((t[i]-DomainInterval[i][0])/(DomainMeshSpacing[i] + 1.0E-14));
                startLatticePt[i] = DomainMeshZoneStart[i][startLatticeIndex[i]];
                endLatticePt[i] = startLatticePt[i] + DomainMeshSpacing[i];
                latticeSpacing = latticeSpacing * DomainMeshSpacing[i];
                sCtrlPt[i] = 3*startLatticeIndex[i];    // the first control point in the mesh zone
            }

            if (index[0] >= sCtrlPt[0] && index[0] <= sCtrlPt[0] + 3 && 
            index[1] >= sCtrlPt[1] && index[1] <= sCtrlPt[1] + 3 && 
            index[2] >= sCtrlPt[2] && index[2] <= sCtrlPt[2] + 3 && 
            index[3] >= sCtrlPt[3] && index[3] <= sCtrlPt[3] + 3)
            {
                d0 = index[0] - sCtrlPt[0];
                d1 = index[1] - sCtrlPt[1];
                d2 = index[2] - sCtrlPt[2];
                d3 = index[3] - sCtrlPt[3];

                result = Bernstein(d0, startLatticePt[0], endLatticePt[0], t[0])
                        *Bernstein(d1, startLatticePt[1], endLatticePt[1], t[1])
                        *Bernstein(d2, startLatticePt[2], endLatticePt[2], t[2])
                        *Bernstein(d3, startLatticePt[3], endLatticePt[3], t[3]);
            }

            /* The aim of the following steps is trying to correct the partial derivative w.r.t. control point by the differential condition: P[3i+1] = 2P[3i] - P[3i-1] , i.e. the control points are not all independent.  */
            for (int i = 0; i < DomainDimeion; i++)
            {
                std::copy(index.begin(), index.end(), tempIndex.begin());

                /* (d P[3k+1]) / (d P[3k]) = 2, this means that when we take partial derivative P[3k] we also need to count the contribution of  (d f) / (d P[3k+1]) x2 */
                if (index[i]%3 == 0 && index[i] != 0 && index[i] < CtrlPtsNumInDAxis[i]-3)
                {
                    tempIndex[i] = index[i] + 1;

                    if (tempIndex[0] >= sCtrlPt[0] && tempIndex[0] <= sCtrlPt[0] + 3 && 
                        tempIndex[1] >= sCtrlPt[1] && tempIndex[1] <= sCtrlPt[1] + 3 && 
                        tempIndex[2] >= sCtrlPt[2] && tempIndex[2] <= sCtrlPt[2] + 3 && 
                        tempIndex[3] >= sCtrlPt[3] && tempIndex[3] <= sCtrlPt[3] + 3)
                    {
                        d0 = tempIndex[0] - sCtrlPt[0];
                        d1 = tempIndex[1] - sCtrlPt[1];
                        d2 = tempIndex[2] - sCtrlPt[2];
                        d3 = tempIndex[3] - sCtrlPt[3];

                        result = result + 2.0*correctionFactor*( 
                                Bernstein(d0, startLatticePt[0], endLatticePt[0], t[0])
                                *Bernstein(d1, startLatticePt[1], endLatticePt[1], t[1])
                                *Bernstein(d2, startLatticePt[2], endLatticePt[2], t[2])
                                *Bernstein(d3, startLatticePt[3], endLatticePt[3], t[3]) );
                    }
                }

                /* because the differential condition will reduce one degree of freedom, so we don't consider the third control point P[3k+1] and its (d ...) / (d P[3k+1]) */
                // if (index[i]%3 == 1 && index[i] != 1)
                // {
                //     tempIndex[i] = index[i] - 1;

                //     if (tempIndex[0] >= sCtrlPt[0] && tempIndex[0] <= sCtrlPt[0] + 3 && 
                //         tempIndex[1] >= sCtrlPt[1] && tempIndex[1] <= sCtrlPt[1] + 3 && 
                //         tempIndex[2] >= sCtrlPt[2] && tempIndex[2] <= sCtrlPt[2] + 3 && 
                //         tempIndex[3] >= sCtrlPt[3] && tempIndex[3] <= sCtrlPt[3] + 3)
                //     {
                //         d0 = tempIndex[0] - sCtrlPt[0];
                //         d1 = tempIndex[1] - sCtrlPt[1];
                //         d2 = tempIndex[2] - sCtrlPt[2];
                //         d3 = tempIndex[3] - sCtrlPt[3];

                //         result = result + 0.5*correctionFactor*( 
                //                 Bernstein(d0, startLatticePt[0], endLatticePt[0], t[0])
                //                 *Bernstein(d1, startLatticePt[1], endLatticePt[1], t[1])
                //                 *Bernstein(d2, startLatticePt[2], endLatticePt[2], t[2])
                //                 *Bernstein(d3, startLatticePt[3], endLatticePt[3], t[3]) );
                //     }
                // }

                /* (d P[3k]) / (d P[3k-1]) = 0.5, this means that when we take partial derivative P[3k-1] we also need to count the contribution of  (d f) / (d P[3k]) x0.5. 
                Here we use the transmission property to recover the differentiable continuity; see the document for details. */
                if (index[i]%3 == 2 && index[i] < CtrlPtsNumInDAxis[i]-3)
                {
                    tempIndex[i] = index[i] + 1;

                    if (tempIndex[0] >= sCtrlPt[0] && tempIndex[0] <= sCtrlPt[0] + 3 && 
                        tempIndex[1] >= sCtrlPt[1] && tempIndex[1] <= sCtrlPt[1] + 3 && 
                        tempIndex[2] >= sCtrlPt[2] && tempIndex[2] <= sCtrlPt[2] + 3 && 
                        tempIndex[3] >= sCtrlPt[3] && tempIndex[3] <= sCtrlPt[3] + 3)
                    {
                        d0 = tempIndex[0] - sCtrlPt[0];
                        d1 = tempIndex[1] - sCtrlPt[1];
                        d2 = tempIndex[2] - sCtrlPt[2];
                        d3 = tempIndex[3] - sCtrlPt[3];

                        result = result + 0.5*correctionFactor*( 
                                Bernstein(d0, startLatticePt[0], endLatticePt[0], t[0])
                                *Bernstein(d1, startLatticePt[1], endLatticePt[1], t[1])
                                *Bernstein(d2, startLatticePt[2], endLatticePt[2], t[2])
                                *Bernstein(d3, startLatticePt[3], endLatticePt[3], t[3]) );
                    }
                }

            }

            return result / pow(latticeSpacing, 3);
        }


        /*  ==============================================
            the partial derivative of differential Bezier function w.r.t. control point cpt, i.e. (d d_uf ) / (d cpt).
            the similar idea can look the comment in the function CPtFactorOfBezierMap()
            ==============================================  */
        double CPtFactorOfDiffBezierMap(int cpt, std::array<double,4> t, int k, double correctionFactor = 1.0){
            double result = 0.0;
            std::array<int,4> startLatticeIndex = {{0,0,0,0}};
            std::array<double,4> startLatticePt = {{-1.0, -1.0, -1.0, -1.0}};
            std::array<double,4> endLatticePt = {{1.0, 1.0, 1.0, 1.0}};
            std::array<int,4> sCtrlPt = {{0,0,0,0}};
            double latticeSpacing = 1.0;
            int pt;
            std::array<int, 4> index = CtrlPtSerialNumtoIndex[cpt];
            std::array<int, 4> tempIndex;
            int d0;
            int d1;
            int d2;
            int d3;

            for (int i = 0; i < DomainDimeion; i++)
            {
                startLatticeIndex[i] = (int)((t[i]-DomainInterval[i][0])/(DomainMeshSpacing[i] + 1.0E-14));
                startLatticePt[i] = DomainMeshZoneStart[i][startLatticeIndex[i]];
                endLatticePt[i] = startLatticePt[i] + DomainMeshSpacing[i];
                latticeSpacing = latticeSpacing * DomainMeshSpacing[i];
                sCtrlPt[i] = 3*startLatticeIndex[i];
            }

            std::function<double(int, double, double, double)> func0 = Bernstein;
            std::function<double(int, double, double, double)> func1 = Bernstein;
            std::function<double(int, double, double, double)> func2 = Bernstein;
            std::function<double(int, double, double, double)> func3 = Bernstein;

            switch ( k )
            {
                case 0:
                    func0 = DiffBernstein;
                    break;
                case 1:
                    func1 = DiffBernstein;
                    break;
                case 2:
                    func2 = DiffBernstein;
                    break;
                case 3:
                    func3 = DiffBernstein;
                    break;
                default:
                    break;
            }

            if (index[0] >= sCtrlPt[0] && index[0] <= sCtrlPt[0] + 3 && 
            index[1] >= sCtrlPt[1] && index[1] <= sCtrlPt[1] + 3 && 
            index[2] >= sCtrlPt[2] && index[2] <= sCtrlPt[2] + 3 && 
            index[3] >= sCtrlPt[3] && index[3] <= sCtrlPt[3] + 3)
            {
                d0 = index[0] - sCtrlPt[0];
                d1 = index[1] - sCtrlPt[1];
                d2 = index[2] - sCtrlPt[2];
                d3 = index[3] - sCtrlPt[3];

                result = func0(d0, startLatticePt[0], endLatticePt[0], t[0])
                        *func1(d1, startLatticePt[1], endLatticePt[1], t[1])
                        *func2(d2, startLatticePt[2], endLatticePt[2], t[2])
                        *func3(d3, startLatticePt[3], endLatticePt[3], t[3]);
            }

            /* The aim of the following steps is trying to correct the partial derivative w.r.t. control point by the differential condition: P[3i+1] = 2P[3i] - P[3i-1] , i.e. the control points are not all independent.  */
            for (int i = 0; i < DomainDimeion; i++)
            {
                std::copy(index.begin(), index.end(), tempIndex.begin());

                if (index[i]%3 == 0 && index[i] != 0 && index[i] < CtrlPtsNumInDAxis[i]-3)
                {
                    tempIndex[i] = index[i] + 1;

                    if (tempIndex[0] >= sCtrlPt[0] && tempIndex[0] <= sCtrlPt[0] + 3 && 
                        tempIndex[1] >= sCtrlPt[1] && tempIndex[1] <= sCtrlPt[1] + 3 && 
                        tempIndex[2] >= sCtrlPt[2] && tempIndex[2] <= sCtrlPt[2] + 3 && 
                        tempIndex[3] >= sCtrlPt[3] && tempIndex[3] <= sCtrlPt[3] + 3)
                    {
                        d0 = tempIndex[0] - sCtrlPt[0];
                        d1 = tempIndex[1] - sCtrlPt[1];
                        d2 = tempIndex[2] - sCtrlPt[2];
                        d3 = tempIndex[3] - sCtrlPt[3];

                        result = result + 2.0*correctionFactor*( 
                                func0(d0, startLatticePt[0], endLatticePt[0], t[0])
                                *func1(d1, startLatticePt[1], endLatticePt[1], t[1])
                                *func2(d2, startLatticePt[2], endLatticePt[2], t[2])
                                *func3(d3, startLatticePt[3], endLatticePt[3], t[3]) );
                    }
                }

                // if (index[i]%3 == 1 && index[i] != 1)
                // {
                //     tempIndex[i] = index[i] - 1;

                //     if (tempIndex[0] >= sCtrlPt[0] && tempIndex[0] <= sCtrlPt[0] + 3 && 
                //         tempIndex[1] >= sCtrlPt[1] && tempIndex[1] <= sCtrlPt[1] + 3 && 
                //         tempIndex[2] >= sCtrlPt[2] && tempIndex[2] <= sCtrlPt[2] + 3 && 
                //         tempIndex[3] >= sCtrlPt[3] && tempIndex[3] <= sCtrlPt[3] + 3)
                //     {
                //         d0 = tempIndex[0] - sCtrlPt[0];
                //         d1 = tempIndex[1] - sCtrlPt[1];
                //         d2 = tempIndex[2] - sCtrlPt[2];
                //         d3 = tempIndex[3] - sCtrlPt[3];

                //         result = result + 0.5*correctionFactor*( 
                //                 Bernstein(d0, startLatticePt[0], endLatticePt[0], t[0])
                //                 *Bernstein(d1, startLatticePt[1], endLatticePt[1], t[1])
                //                 *Bernstein(d2, startLatticePt[2], endLatticePt[2], t[2])
                //                 *Bernstein(d3, startLatticePt[3], endLatticePt[3], t[3]) );
                //     }
                // }

                if (index[i]%3 == 2 && index[i] < CtrlPtsNumInDAxis[i]-3)
                {
                    tempIndex[i] = index[i] + 1;

                    if (tempIndex[0] >= sCtrlPt[0] && tempIndex[0] <= sCtrlPt[0] + 3 && 
                        tempIndex[1] >= sCtrlPt[1] && tempIndex[1] <= sCtrlPt[1] + 3 && 
                        tempIndex[2] >= sCtrlPt[2] && tempIndex[2] <= sCtrlPt[2] + 3 && 
                        tempIndex[3] >= sCtrlPt[3] && tempIndex[3] <= sCtrlPt[3] + 3)
                    {
                        d0 = tempIndex[0] - sCtrlPt[0];
                        d1 = tempIndex[1] - sCtrlPt[1];
                        d2 = tempIndex[2] - sCtrlPt[2];
                        d3 = tempIndex[3] - sCtrlPt[3];

                        result = result + 0.5*correctionFactor*( 
                                func0(d0, startLatticePt[0], endLatticePt[0], t[0])
                                *func1(d1, startLatticePt[1], endLatticePt[1], t[1])
                                *func2(d2, startLatticePt[2], endLatticePt[2], t[2])
                                *func3(d3, startLatticePt[3], endLatticePt[3], t[3]) );
                    }
                }
            }

            return result / pow(latticeSpacing, 3);
        }


        /*  ==============================================
            the partial derivative of double differential Bezier function w.r.t. control point cpt, i.e. (d d_uf ) / (d cpt).
            the similar idea can look the comment in the function CPtFactorOfBezierMap()
            ==============================================  */
        double CPtFactorOfDoubleDiffBezierMap(int cpt, std::array<double,4> t, int k, int l, double correctionFactor = 1.0){
            double result = 0.0;
            std::array<int,4> startLatticeIndex = {{0,0,0,0}};
            std::array<double,4> startLatticePt = {{-1.0, -1.0, -1.0, -1.0}};
            std::array<double,4> endLatticePt = {{1.0, 1.0, 1.0, 1.0}};
            std::array<int,4> sCtrlPt = {{0,0,0,0}};
            double latticeSpacing = 1.0;
            int pt;
            std::array<int, 4> index = CtrlPtSerialNumtoIndex[cpt];
            std::array<int, 4> tempIndex;
            int d0;
            int d1;
            int d2;
            int d3;

            for (int i = 0; i < DomainDimeion; i++)
            {
                startLatticeIndex[i] = (int)((t[i]-DomainInterval[i][0])/(DomainMeshSpacing[i] + 1.0E-14));
                startLatticePt[i] = DomainMeshZoneStart[i][startLatticeIndex[i]];
                endLatticePt[i] = startLatticePt[i] + DomainMeshSpacing[i];
                latticeSpacing = latticeSpacing * DomainMeshSpacing[i];
                sCtrlPt[i] = 3*startLatticeIndex[i];
            }

            std::function<double(int, double, double, double)> func0 = Bernstein;
            std::function<double(int, double, double, double)> func1 = Bernstein;
            std::function<double(int, double, double, double)> func2 = Bernstein;
            std::function<double(int, double, double, double)> func3 = Bernstein;
            
            switch( k )
            {
                case 0:
                    switch( l )
                    {
                        case 0:
                            func0 = DoubleDiffBernstein;
                            break;
                        case 1:
                            func0 = DiffBernstein;
                            func1 = DiffBernstein;
                            break;
                        case 2:
                            func0 = DiffBernstein;
                            func2 = DiffBernstein;
                            break;
                        case 3:
                            func0 = DiffBernstein;
                            func3 = DiffBernstein;
                            break;
                        default:
                            break;
                    }
                    break;
                case 1:
                    switch( l )
                    {
                        case 0:
                            func1 = DiffBernstein;
                            func0 = DiffBernstein;
                            break;
                        case 1:
                            func1 = DoubleDiffBernstein;
                            break;
                        case 2:
                            func1 = DiffBernstein;
                            func2 = DiffBernstein;
                            break;
                        case 3:
                            func1 = DiffBernstein;
                            func3 = DiffBernstein;
                            break;
                        default:
                            break;
                    }
                    break;
                case 2:
                    switch( l )
                    {
                        case 0:
                            func2 = DiffBernstein;
                            func0 = DiffBernstein;
                            break;
                        case 1:
                            func2 = DiffBernstein;
                            func1 = DiffBernstein;
                            break;
                        case 2:
                            func2 = DoubleDiffBernstein;
                            break;
                        case 3:
                            func2 = DiffBernstein;
                            func3 = DiffBernstein;
                            break;
                        default:
                            break;
                    }
                    break;
                case 3:
                    switch( l )
                    {
                        case 0:
                            func3 = DiffBernstein;
                            func0 = DiffBernstein;
                            break;
                        case 1:
                            func3 = DiffBernstein;
                            func1 = DiffBernstein;
                            break;
                        case 2:
                            func3 = DiffBernstein;
                            func2 = DiffBernstein;
                            break;
                        case 3:
                            func3 = DoubleDiffBernstein;
                            break;
                        default:
                            break;
                    }
                    break;
                default:
                    break;
            }
            
            if (index[0] >= sCtrlPt[0] && index[0] <= sCtrlPt[0] + 3 && 
            index[1] >= sCtrlPt[1] && index[1] <= sCtrlPt[1] + 3 && 
            index[2] >= sCtrlPt[2] && index[2] <= sCtrlPt[2] + 3 && 
            index[3] >= sCtrlPt[3] && index[3] <= sCtrlPt[3] + 3)
            {
                d0 = index[0] - sCtrlPt[0];
                d1 = index[1] - sCtrlPt[1];
                d2 = index[2] - sCtrlPt[2];
                d3 = index[3] - sCtrlPt[3];

                result = func0(d0, startLatticePt[0], endLatticePt[0], t[0])
                        *func1(d1, startLatticePt[1], endLatticePt[1], t[1])
                        *func2(d2, startLatticePt[2], endLatticePt[2], t[2])
                        *func3(d3, startLatticePt[3], endLatticePt[3], t[3]);
            }

            /* The aim of the following steps is trying to correct the partial derivative w.r.t. control point by the differential condition: P[3i+1] = 2P[3i] - P[3i-1] , i.e. the control points are not all independent.  */
            for (int i = 0; i < DomainDimeion; i++)
            {
                std::copy(index.begin(), index.end(), tempIndex.begin());

                if (index[i]%3 == 0 && index[i] != 0 && index[i] < CtrlPtsNumInDAxis[i]-3)
                {
                    tempIndex[i] = index[i] + 1;

                    if (tempIndex[0] >= sCtrlPt[0] && tempIndex[0] <= sCtrlPt[0] + 3 && 
                        tempIndex[1] >= sCtrlPt[1] && tempIndex[1] <= sCtrlPt[1] + 3 && 
                        tempIndex[2] >= sCtrlPt[2] && tempIndex[2] <= sCtrlPt[2] + 3 && 
                        tempIndex[3] >= sCtrlPt[3] && tempIndex[3] <= sCtrlPt[3] + 3)
                    {
                        d0 = tempIndex[0] - sCtrlPt[0];
                        d1 = tempIndex[1] - sCtrlPt[1];
                        d2 = tempIndex[2] - sCtrlPt[2];
                        d3 = tempIndex[3] - sCtrlPt[3];

                        result = result + 2.0*correctionFactor*( 
                                func0(d0, startLatticePt[0], endLatticePt[0], t[0])
                                *func1(d1, startLatticePt[1], endLatticePt[1], t[1])
                                *func2(d2, startLatticePt[2], endLatticePt[2], t[2])
                                *func3(d3, startLatticePt[3], endLatticePt[3], t[3]) );
                    }
                }

                // if (index[i]%3 == 1 && index[i] != 1)
                // {
                //     tempIndex[i] = index[i] - 1;

                //     if (tempIndex[0] >= sCtrlPt[0] && tempIndex[0] <= sCtrlPt[0] + 3 && 
                //         tempIndex[1] >= sCtrlPt[1] && tempIndex[1] <= sCtrlPt[1] + 3 && 
                //         tempIndex[2] >= sCtrlPt[2] && tempIndex[2] <= sCtrlPt[2] + 3 && 
                //         tempIndex[3] >= sCtrlPt[3] && tempIndex[3] <= sCtrlPt[3] + 3)
                //     {
                //         d0 = tempIndex[0] - sCtrlPt[0];
                //         d1 = tempIndex[1] - sCtrlPt[1];
                //         d2 = tempIndex[2] - sCtrlPt[2];
                //         d3 = tempIndex[3] - sCtrlPt[3];

                //         result = result + 0.5*correctionFactor*( 
                //                 Bernstein(d0, startLatticePt[0], endLatticePt[0], t[0])
                //                 *Bernstein(d1, startLatticePt[1], endLatticePt[1], t[1])
                //                 *Bernstein(d2, startLatticePt[2], endLatticePt[2], t[2])
                //                 *Bernstein(d3, startLatticePt[3], endLatticePt[3], t[3]) );
                //     }
                // }
                                
                if (index[i]%3 == 2 && index[i] < CtrlPtsNumInDAxis[i]-3)
                {
                    tempIndex[i] = index[i] + 1;

                    if (tempIndex[0] >= sCtrlPt[0] && tempIndex[0] <= sCtrlPt[0] + 3 && 
                        tempIndex[1] >= sCtrlPt[1] && tempIndex[1] <= sCtrlPt[1] + 3 && 
                        tempIndex[2] >= sCtrlPt[2] && tempIndex[2] <= sCtrlPt[2] + 3 && 
                        tempIndex[3] >= sCtrlPt[3] && tempIndex[3] <= sCtrlPt[3] + 3)
                    {
                        d0 = tempIndex[0] - sCtrlPt[0];
                        d1 = tempIndex[1] - sCtrlPt[1];
                        d2 = tempIndex[2] - sCtrlPt[2];
                        d3 = tempIndex[3] - sCtrlPt[3];

                        result = result + 0.5*correctionFactor*( 
                                func0(d0, startLatticePt[0], endLatticePt[0], t[0])
                                *func1(d1, startLatticePt[1], endLatticePt[1], t[1])
                                *func2(d2, startLatticePt[2], endLatticePt[2], t[2])
                                *func3(d3, startLatticePt[3], endLatticePt[3], t[3]) );
                    }
                }
            }

            return result / pow(latticeSpacing, 3);
        }


        /*  ==============================================
            The distance square of two point of Bezier image in target space.
            ============================================== */
        double DistanceSquare( std::array<double,4> t1, std::array<double,4> t2 ){
            double result = 0.0;
            coordinate pt1 = BezierMap(t1);
            coordinate pt2 = BezierMap(t2);

            for (int i = 0; i < TargetDimeion; i++)
            {
                result = result + pow( pt1.coord[i] - pt2.coord[i] , 2);
            }            
            return result;
        }


        /*  ==============================================
            Write the configuration data in the target space of the Bezier image into file.
            ============================================== */
        void WriteConfigFile(std::string fileName, std::vector<int> samples)
        {
            std::ofstream myfile;
            myfile.open (fileName);
            std::array<double, 4> t = {{0.0, 0.0, 0.0, 0.0}};
            std::vector<double> stepGap;
            std::array<int, 4> fSamples = {{1,1,1,1}};
            coordinate pt;

            for (int i = 0; i < samples.size(); i++)
            {
                stepGap.push_back( (DomainInterval[i][1]-DomainInterval[i][0])/(double)samples[i] );
                fSamples[i] = samples[i] + 1;
            }
            
            t[0] = DomainInterval[0][0];
            for (int d0 = 0; d0 < fSamples[0]; d0++)
            {
                t[1] = DomainInterval[1][0];
                for (int d1 = 0; d1 < fSamples[1]; d1++)
                {
                    t[2] = DomainInterval[2][0];
                    for (int d2 = 0; d2 < fSamples[2]; d2++)
                    {
                        t[3] = DomainInterval[3][0];
                        for (int d3 = 0; d3 < fSamples[3]; d3++)
                        {
                            pt = BezierMap(t);
                            for (int j = 0; j < TargetDimeion; j++)
                            {
                                myfile << std::setprecision(5) << pt.coord[j] << " ";
                            }
                            
                            if (DomainDimeion >= 4 && d3 != fSamples[3]-1)
                            {
                                myfile << std::endl;
                                t[3] = t[3] + stepGap[3];
                            }
                            
                        }

                        if (DomainDimeion >= 3 && d2 != fSamples[2]-1)
                        {
                            myfile << std::endl;
                            t[2] = t[2] + stepGap[2];
                        }

                    }

                    if (DomainDimeion >= 2 && d1 != fSamples[1]-1)
                    {
                        myfile << std::endl;
                        t[1] = t[1] + stepGap[1];
                    }

                }

                if (DomainDimeion >= 1)
                {
                    myfile << std::endl;
                    t[0] = t[0] + stepGap[0];
                }

            }

            myfile.close();
        }


        /*  ==============================================
            Write the control point data of the Bezier image into file.
            ============================================== */
        void WriteParamFile(std::string fileName)
        {
            std::ofstream myfile;
            myfile.open (fileName);
            
            myfile << DomainDimeion << " " << TargetDimeion << std::endl;

            for (int i = 0; i < DomainDimeion; i++)
            {
                myfile << DomainInterval[i][0] << " " << DomainInterval[i][1] << std::endl;
            }

            for (int i = 0; i < CtrlPtsTotalNum; i++)
            {
                for (int k = 0; k < DomainDimeion; k++)
                {
                    myfile << CtrlPtSerialNumtoIndex[i][k] << " ";
                }

                for (int j = 0; j < TargetDimeion; j++)
                {
                    myfile << std::setprecision(3) << CtrlPoint[i].coord[j] << " ";
                }
                myfile << std::endl;
            }
            myfile.close();
        }


        /*  ==============================================
            Write the configuration data in the target space of the Bezier image into file.
            And separate them by the first domain variable t[0].
            ============================================== */
        void WriteConfigFileEachTime(std::string fileName, std::vector<int> samples)
        {
            std::array<double, 4> t = {{0.0, 0.0, 0.0, 0.0}};
            std::vector<double> stepGap;
            std::array<int, 4> fSamples = {{1,1,1,1}};
            coordinate pt;

            for (int i = 0; i < samples.size(); i++)
            {
                stepGap.push_back( (DomainInterval[i][1]-DomainInterval[i][0])/(double)samples[i] );
                fSamples[i] = samples[i] + 1;
            }
            
            t[0] = DomainInterval[0][0];
            for (int d0 = 0; d0 < fSamples[0]; d0++)
            {
                std::ostringstream tempStream;
                tempStream << fileName << "_" << "t" << d0 << "_coord.txt";

                std::ofstream myfile;
                myfile.open(tempStream.str());

                t[1] = DomainInterval[1][0];
                for (int d1 = 0; d1 < fSamples[1]; d1++)
                {
                    t[2] = DomainInterval[2][0];
                    for (int d2 = 0; d2 < fSamples[2]; d2++)
                    {
                        t[3] = DomainInterval[3][0];
                        for (int d3 = 0; d3 < fSamples[3]; d3++)
                        {
                            pt = BezierMap(t);
                            for (int j = 0; j < TargetDimeion; j++)
                            {
                                myfile << std::setprecision(5) << pt.coord[j] << " ";
                            }
                            
                            if (DomainDimeion >= 4 && d3 != fSamples[3]-1)
                            {
                                myfile << std::endl;
                                t[3] = t[3] + stepGap[3];
                            }
                            
                        }

                        if (DomainDimeion >= 3 && d2 != fSamples[2]-1)
                        {
                            myfile << std::endl;
                            t[2] = t[2] + stepGap[2];
                        }

                    }

                    if (DomainDimeion >= 2 && d1 != fSamples[1]-1)
                    {
                        myfile << std::endl;
                        t[1] = t[1] + stepGap[1];
                    }

                }

                if (DomainDimeion >= 1)
                {
                    myfile << std::endl;
                    t[0] = t[0] + stepGap[0];
                }

                myfile.close();
            
                std::this_thread::sleep_for(std::chrono::seconds(1));
            }
        }


        /*  ==============================================
            Write the configuration data in the target space of the Bezier image into file with domain variables.
            ============================================== */
        void WriteConfigWithVariableFile(std::string fileName, std::vector<int> samples)
        {
            std::ofstream myfile;
            myfile.open (fileName);
            std::array<double, 4> t = {{0.0, 0.0, 0.0, 0.0}};
            std::vector<double> stepGap;
            std::array<int, 4> fSamples = {{1,1,1,1}};
            coordinate pt;

            for (int i = 0; i < samples.size(); i++)
            {
                stepGap.push_back( (DomainInterval[i][1]-DomainInterval[i][0])/(double)samples[i] );
                fSamples[i] = samples[i] + 1;
            }
            
            t[0] = DomainInterval[0][0];
            for (int d0 = 0; d0 < fSamples[0]; d0++)
            {
                t[1] = DomainInterval[1][0];
                for (int d1 = 0; d1 < fSamples[1]; d1++)
                {
                    t[2] = DomainInterval[2][0];
                    for (int d2 = 0; d2 < fSamples[2]; d2++)
                    {
                        t[3] = DomainInterval[3][0];
                        for (int d3 = 0; d3 < fSamples[3]; d3++)
                        {
                            pt = BezierMap(t);
                            
                            for (int j = 0; j < DomainDimeion; j++)
                            {
                                myfile << std::setprecision(5) << t[j] << " ";
                            }

                            for (int j = 0; j < TargetDimeion; j++)
                            {
                                myfile << std::setprecision(5) << pt.coord[j] << " ";
                            }
                            
                            if (DomainDimeion >= 4 && d3 != fSamples[3]-1)
                            {
                                myfile << std::endl;
                                t[3] = t[3] + stepGap[3];
                            }
                            
                        }

                        if (DomainDimeion >= 3 && d2 != fSamples[2]-1)
                        {
                            myfile << std::endl;
                            t[2] = t[2] + stepGap[2];
                        }

                    }

                    if (DomainDimeion >= 2 && d1 != fSamples[1]-1)
                    {
                        myfile << std::endl;
                        t[1] = t[1] + stepGap[1];
                    }

                }

                if (DomainDimeion >= 1)
                {
                    myfile << std::endl;
                    t[0] = t[0] + stepGap[0];
                }

            }

            myfile.close();
        }    
    
    };

