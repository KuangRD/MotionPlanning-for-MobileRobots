#ifndef _TRAJECTORY_GENERATOR_WAYPOINT_H_
#define _TRAJECTORY_GENERATOR_WAYPOINT_H_

#include <Eigen/Eigen>
#include <vector>

class TrajectoryGeneratorWaypoint {
    private:
		double _qp_cost;
		Eigen::MatrixXd _Q;
		Eigen::VectorXd _Px, _Py, _Pz;
    public:
        TrajectoryGeneratorWaypoint();

        ~TrajectoryGeneratorWaypoint();

        Eigen::MatrixXd PolyQPGeneration(
            const int order,
            const Eigen::MatrixXd &Path,
            const Eigen::MatrixXd &Vel,
            const Eigen::MatrixXd &Acc,
            const Eigen::VectorXd &Time);
        
        int Factorial(int x);
        
        void get_Q(const int p_order,const Eigen::VectorXd &ts, Eigen::MatrixXd &Q) ;

        void get_M(int d_order, int p_order, const Eigen::VectorXd &ts, Eigen::MatrixXd &M);

        void get_Ct(int d_order, int p_order, const Eigen::VectorXd &ts, Eigen::MatrixXd &Ct);
};
        

#endif
