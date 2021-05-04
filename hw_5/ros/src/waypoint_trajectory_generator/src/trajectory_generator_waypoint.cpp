#include "trajectory_generator_waypoint.h"
#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h> 

using namespace std;    
using namespace Eigen;

TrajectoryGeneratorWaypoint::TrajectoryGeneratorWaypoint(){}
TrajectoryGeneratorWaypoint::~TrajectoryGeneratorWaypoint(){}

//define factorial function, input i, output i!
int TrajectoryGeneratorWaypoint::Factorial(int x)
{
    int fac = 1;
    for(int i = x; i > 0; i--)
        fac = fac * i;
    return fac;
}
/*

    STEP 2: Learn the "Closed-form solution to minimum snap" in L5, then finish this PolyQPGeneration function

    variable declaration: input       const int d_order,                    // the order of derivative
                                      const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
                                      const Eigen::MatrixXd &Vel,           // boundary velocity
                                      const Eigen::MatrixXd &Acc,           // boundary acceleration
                                      const Eigen::VectorXd &Time)          // time allocation in each segment
                          output      MatrixXd PolyCoeff(m, 3 * p_num1d);   // position(x,y,z), so we need (3 * p_num1d) coefficients

*/

Eigen::MatrixXd TrajectoryGeneratorWaypoint::PolyQPGeneration(
            const int d_order,                    // the order of derivative
            const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
            const Eigen::MatrixXd &Vel,           // boundary velocity
            const Eigen::MatrixXd &Acc,           // boundary acceleration
            const Eigen::VectorXd &Time)          // time allocation in each segment
{
    // enforce initial and final velocity and accleration, for higher order derivatives, just assume them be 0;
    int p_order   = 2 * d_order - 1;              // the order of polynomial
    int p_num1d   = p_order + 1;                  // the number of variables in each segment

    int m = Time.size();                          // the number of segments
    MatrixXd PolyCoeff = MatrixXd::Zero(m, 3 * p_num1d);           // position(x,y,z), so we need (3 * p_num1d) coefficients
    VectorXd Px(p_num1d * m), Py(p_num1d * m), Pz(p_num1d * m);

    MatrixXd Q = MatrixXd::Zero(p_num1d * m, p_num1d * m);
    get_Q(p_order, Time, Q);

    MatrixXd M = MatrixXd::Zero(p_num1d * m, p_num1d * m);
    get_M(p_order, Time, M);
    MatrixXd M_inv = M.inverse();

    int num_dF = 2*d_order+(m-1)*1;
    int num_dP = (m-1)*(d_order-1);
    MatrixXd Ct = MatrixXd::Zero(p_num1d * m, num_dF + num_dP);
    get_Ct(d_order, p_order, Time, Ct);
    MatrixXd C = Ct.transpose();
    
    MatrixXd R =  C * M_inv.transpose() * Q * M_inv * Ct;
    MatrixXd R_fp = R.topRightCorner(num_dF, num_dP);
    MatrixXd R_pp = R.bottomRightCorner(num_dP, num_dP);


    for (int dim_idx = 0; dim_idx < 3; dim_idx++)
    {
        VectorXd path = Path.col(dim_idx);

        VectorXd dF = VectorXd::Zero(num_dF);
        VectorXd dP = VectorXd::Zero(num_dP);

        dF(0) = path(0, dim_idx);
        dF(1) = Vel(0, dim_idx);
        dF(2) = Acc(0, dim_idx);
        dF(m+d_order-1) = path(m, dim_idx);
        dF(m+d_order) = Vel(1, dim_idx);
        dF(m+d_order+1) = Acc(1, dim_idx);

        for (int mid_idx = 1; mid_idx < m; mid_idx++)
        {
            dF(d_order-1+mid_idx) = path(mid_idx, dim_idx);
        }
        
        dP = -1.0 * R_pp.inverse() * R_fp.transpose() * dF;

        VectorXd d = VectorXd::Zero(num_dF+num_dP);
        d << dF, dP;

        VectorXd poly_coef = M_inv * Ct * d;
        MatrixXd poly_coef_row = poly_coef.transpose();

        for (int k = 0; k < m; k++)
        {
            PolyCoeff.block(k, i*p_num1d, 1, p_num1d) = poly_coef_row.block(0, k*p_num1d, 1, p_num1d);
        }
    }

    return PolyCoeff;
}

void get_Q(const int p_order,const Eigen::VectorXd &ts, Eigen::MatrixXd &Q) {
    int n_seg = ts.size();
    
    for (int k=0; k<n_seg; k++) {
        Eigen::MatrixXd Q_k = MatrixXd::Zero(p_order+1, p_order+1);
        for (int i=4; i <= p_order; i++) {
            for (int l=4; l <= p_order; l++) {
                Q_k(i+1,l+1) = i*(i-1)*(i-2)*(i-3)*l*(l-1)*(l-2)*(l-3) * pow(ts(k),(i + l - p_order)) / (i + l - p_order);
            }
        }
        Q.block(k * (p_order+1), k * (p_order+1), (p_order+1), (p_order+1)) = Q_k;    
    }
}


void get_M(int p_order, const Eigen::VectorXd &ts, Eigen::MatrixXd &M) {
    int n_seg = ts.size();
    
    for (int k = 0; k < n_seg; k++)
    {
        MatrixXd Mk = MatrixXd::Zero(p_order+1, p_order+1);
        for (int i = 0; i < d_order; i++) Mk(i, i) = Factorial(i);
        for (int i = d_order; i <= p_order; i++)
        {
            for (int j = i - d_order; j <= p_order; j++)
            {
                Mk(i, j) = Factorial(j) / Factorial(j-(i-d_order)) * pow(ts(k), j-(i-d_order));
            }
        }

        M.block(k * (p_order+1), k * (p_order+1), (p_order+1), (p_order+1)) = Mk;
    }
}

void get_Ct(int d_order, int p_order, const Eigen::VectorXd &ts, Eigen::MatrixXd &Ct){
    int m = ts.size();
    int p_num1d = p_order+1;
    
    MatrixXd C0 = MatrixXd::Identity(d_order, d_order);
    MatrixXd CT = MatrixXd::Identity(d_order, d_order);

    Ct.block(0, 0, d_order, d_order) = C0;
    Ct.block(d_order*(2*m-1), d_order+m-1, d_order, d_order) = CT;

    VectorXd V_f = VectorXd::Zero(d_order+1);
    V_f(0) = 1;
    V_f(d_order) = 1;
    for (int k = 0; k < m-1; k++)
        Ct.block(d_order+2*d_order*k, d_order+k, d_order+1, 1) = V_f;
    
    MatrixXd M_p = MatrixXd::Zero(p_num1d, d_order-1);
    M_p.block(1, 0, d_order-1, d_order-1) = MatrixXd::Identity(d_order-1, d_order-1);
    M_p.block(d_order+1, 0, d_order-1, d_order-1) = MatrixXd::Identity(d_order-1, d_order-1);

    for (int k = 0; k < m-1; k++)
        Ct.block(d_order+(p_num1d)*k, (d_order-1)*k+2*d_order+(m-1)*1, p_num1d, d_order-1) = M_p;
}
