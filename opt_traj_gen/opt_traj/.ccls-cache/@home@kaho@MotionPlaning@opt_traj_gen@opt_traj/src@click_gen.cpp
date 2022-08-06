#include "lec5_hw/visualizer.hpp"
#include "lec5_hw/trajectory.hpp"

#include <ros/ros.h>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/PoseStamped.h>

#include <cmath>
#include <iostream>
#include <vector>
#include <Eigen/Core>

struct Config
{
    std::string targetTopic;
    double clickHeight;
    std::vector<double> initialVel;
    std::vector<double> initialAcc;
    std::vector<double> terminalVel;
    std::vector<double> terminalAcc;
    double allocationSpeed;
    double allocationAcc;
    int maxPieceNum;

    Config(const ros::NodeHandle &nh_priv)
    {
        nh_priv.getParam("TargetTopic", targetTopic);
        nh_priv.getParam("ClickHeight", clickHeight);
        nh_priv.getParam("InitialVel", initialVel);
        nh_priv.getParam("InitialAcc", initialAcc);
        nh_priv.getParam("TerminalVel", terminalVel);
        nh_priv.getParam("TerminalAcc", terminalAcc);
        nh_priv.getParam("AllocationSpeed", allocationSpeed);
        nh_priv.getParam("AllocationAcc", allocationAcc);
        nh_priv.getParam("MaxPieceNum", maxPieceNum);
    }
};

double timeTrapzVel(const double dist,
                    const double vel,
                    const double acc)
{
    const double t = vel / acc;
    const double d = 0.5 * acc * t * t;

    if (dist < d + d)
    {
        return 2.0 * sqrt(dist / acc);
    }
    else
    {
        return 2.0 * t + (dist - 2.0 * d) / vel;
    }
}

int fatorial(const int n)
{
    int result = 1;
    for (int i = 1; i <= n; i++)
    {
        result *= i;
    }
    return result;
}

Eigen::MatrixXd computeEi(const double time){
    Eigen::MatrixXd Ei = Eigen::MatrixXd::Zero(6,6);
    for(int i = 0; i < 6; ++i){
        for(int j = 0; j < 6; ++j){
            if(i == 0){
                Ei(i ,j) = pow(time, j);
                continue;
            }else if(i == 1) {
                Ei(i, j) = pow(time, j);
                continue;
            }else if(i - 1 <= j){
                Ei(i ,j) = fatorial(j) / fatorial(j - i + 1) * pow(time, j - i + 1);
            }
        }
    }
    return Ei;
}

void minimumJerkTrajGen(
    // Inputs:
    const int pieceNum,
    const Eigen::Vector3d &initialPos,
    const Eigen::Vector3d &initialVel,
    const Eigen::Vector3d &initialAcc,
    const Eigen::Vector3d &terminalPos,
    const Eigen::Vector3d &terminalVel,
    const Eigen::Vector3d &terminalAcc,
    const Eigen::Matrix3Xd &intermediatePositions,
    const Eigen::VectorXd &timeAllocationVector,
    // Outputs:
    Eigen::MatrixX3d &coefficientMatrix)
{
    // coefficientMatrix is a matrix with 6*piece num rows and 3 columes
    // As for a polynomial c0+c1*t+c2*t^2+c3*t^3+c4*t^4+c5*t^5,
    // each 6*3 sub-block of coefficientMatrix is
    // --              --
    // | c0_x c0_y c0_z |
    // | c1_x c1_y c1_z |
    // | c2_x c2_y c2_z |
    // | c3_x c3_y c3_z |
    // | c4_x c4_y c4_z |
    // | c5_x c5_y c5_z |
    // --              --
    double Tm = timeAllocationVector[pieceNum - 1];
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(6 * pieceNum , 6 * pieceNum);
    Eigen::MatrixXd b = Eigen::MatrixXd::Zero(6 * pieceNum , 3);
    Eigen::MatrixXd F0 = Eigen::MatrixXd::Zero(3, 6);
    Eigen::MatrixXd Fj = Eigen::MatrixXd::Zero(6, 6);
    Eigen::MatrixXd Em = Eigen::MatrixXd::Zero(3, 6);

    for(int i = 0; i < 6; ++i){
        if(i == 0) continue;
        for(int j = 0; j < 6; ++j){
            if( j == i - 1){
                Fj(i, j) = -1.0 * fatorial(j );
            }
        }
    }

    F0 << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 2.0, 0.0, 0.0, 0.0;

    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 6; ++j){
            if(i <= j){
                Em(i,j) = fatorial(j) / fatorial(j - i) * pow(Tm, j - i);
            }
        }
    }
//    std::cout << "Tm: " << Tm << std::endl;
//    std::cout << "Em: " << std::endl << Em << std::endl;

    for(int i = 0; i <= pieceNum; ++i){
        if(i == 0){
            Eigen::MatrixXd D0 = Eigen::MatrixXd::Zero(3, 3);
            D0.row(0) = initialPos;
            D0.row(1) = initialVel;
            D0.row(2) = initialAcc;
            b.block(0, 0, 3, 3) = D0;
            M.block(0, 0, 3, 6) = F0;
            continue;
        }else if (i == pieceNum){
            Eigen::MatrixXd Dn = Eigen::MatrixXd::Zero(3, 3);
            Dn.row(0) = terminalPos;
            Dn.row(1) = terminalVel;
            Dn.row(2) = terminalAcc;
            b.block(6 * i - 3, 0, 3, 3) = Dn;
            M.block(6 * i - 3, 6 * i - 6, 3, 6) = Em;
            continue;
        }else{
            Eigen::MatrixXd Di = Eigen::MatrixXd::Zero(6, 3);
//            std::cout << "waypoints: " << std::endl << intermediatePositions.transpose().row(i - 1) << std::endl;
            Di.row(0) = intermediatePositions.transpose().row(i - 1);
            b.block(6 * i - 3, 0, 6, 3) = Di;
        }
        for(int j = 0; j < pieceNum; ++j){
            if(i == j){
                M.block(6 * i - 3, 6 * i, 6, 6) = Fj;
            }else if(j == i - 1){
                Eigen::MatrixXd Ei = computeEi(timeAllocationVector(i - 1));
                M.block(6 * i - 3, 6 * j, 6, 6) = Ei;
            }
        }
    }
    ROS_WARN("BIVP Solver activated");
        std::cout << "----------------M Matrix---------------" << std::endl;
        std::cout << M << std::endl;
        std::cout << "----------------M Matrix---------------" << std::endl;
        std::cout << "----------------b Matrix---------------" << std::endl;
        std::cout << b << std::endl;
        std::cout << "----------------b Matrix---------------" << std::endl;
        Eigen::MatrixXd M_inv = M.inverse();
        coefficientMatrix = M_inv * b;

}

class ClickGen
{
private:
    Config config;

    ros::NodeHandle nh;
    ros::Subscriber targetSub;

    Visualizer visualizer;

    Eigen::Matrix3Xd positions;
    Eigen::VectorXd times;
    int positionNum;
    Trajectory<5> traj;

public:
    ClickGen(const Config &conf,
             ros::NodeHandle &nh_)
        : config(conf),
          nh(nh_),
          visualizer(nh),
          positions(3, config.maxPieceNum + 1),
          times(config.maxPieceNum),
          positionNum(0)
    {
        targetSub = nh.subscribe(config.targetTopic, 1,
                                 &ClickGen::targetCallBack, this,
                                 ros::TransportHints().tcpNoDelay());
    }

    void targetCallBack(const geometry_msgs::PoseStamped::ConstPtr &msg)
    {
        if (positionNum > config.maxPieceNum)
        {
            positionNum = 0;
            traj.clear();
        }

        positions(0, positionNum) = msg->pose.position.x;
        positions(1, positionNum) = msg->pose.position.y;
        positions(2, positionNum) = std::fabs(msg->pose.orientation.z) * config.clickHeight;

        if (positionNum > 0)
        {
            const double dist = (positions.col(positionNum) - positions.col(positionNum - 1)).norm();
            times(positionNum - 1) = timeTrapzVel(dist, config.allocationSpeed, config.allocationAcc);
        }

        ++positionNum;

        if (positionNum > 1)
        {
            const int pieceNum = positionNum - 1;
            const Eigen::Vector3d initialPos = positions.col(0);
            const Eigen::Vector3d initialVel(config.initialVel[0], config.initialVel[1], config.initialVel[2]);
            const Eigen::Vector3d initialAcc(config.initialAcc[0], config.initialAcc[1], config.initialAcc[2]);
            const Eigen::Vector3d terminalPos = positions.col(pieceNum);
            const Eigen::Vector3d terminalVel(config.terminalVel[0], config.terminalVel[1], config.terminalVel[2]);
            const Eigen::Vector3d terminalAcc(config.terminalAcc[0], config.terminalAcc[1], config.terminalAcc[2]);
            const Eigen::Matrix3Xd intermediatePositions = positions.middleCols(1, pieceNum - 1);
            const Eigen::VectorXd timeAllocationVector = times.head(pieceNum);

            Eigen::MatrixX3d coefficientMatrix = Eigen::MatrixXd::Zero(6 * pieceNum, 3);

            minimumJerkTrajGen(pieceNum,
                               initialPos, initialVel, initialAcc,
                               terminalPos, terminalVel, terminalAcc,
                               intermediatePositions,
                               timeAllocationVector,
                               coefficientMatrix);

            traj.clear();
            traj.reserve(pieceNum);
            for (int i = 0; i < pieceNum; i++)
            {
                traj.emplace_back(timeAllocationVector(i),
                                  coefficientMatrix.block<6, 3>(6 * i, 0).transpose().rowwise().reverse());
            }
        }

        visualizer.visualize(traj, positions.leftCols(positionNum));

        return;
    }
};

int main(int argc, char **argv)
{
    ros::init(argc, argv, "click_gen_node");
    ros::NodeHandle nh_;
    ClickGen clickGen(Config(ros::NodeHandle("~")), nh_);
    ros::spin();
    return 0;
}
