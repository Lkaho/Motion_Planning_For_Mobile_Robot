## Trajectory Optimization

### 1 算法流程

#### 1.1 入口函数：void targetCallBack(const geometry_msgs::PoseStamped::ConstPtr &msg)

订阅rivz中用户输入的目标点信息：

```c++
{
    /.../
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
    /.../
}
```

#### 1.2 时间分配函数 timeTrapzVel

根据梯形速度曲线对轨迹的时间进行分配

```c++
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
```

#### 1.3 根据BIVP的问题形式对轨迹的系数进行求解

与VJ Kumar论文minimum snap求解方法的不同在于，该种求解方法无需使用目标函数求解。因为由泛函极值条件可知若最小化Jerk，该曲线一定是一条5次多项式，那么根据初末p，v，a的约束，以及中间点直到snap的连续性约束(未知数个数与等式约束条件相等)，即可求解该条最优的曲线。

```c++
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
    // Please computed coefficientMatrix of the minimum-jerk trajectory
    // in this function

    // ------------------------ Put your solution below -----------------------
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

    // ------------------------ Put your solution above ------------------------
}
```

### 2 运行结果

![min_jerkpng](/home/kaho/mpCourseHw/opt_traj/src/lec5_hw/min_jerkpng.png)