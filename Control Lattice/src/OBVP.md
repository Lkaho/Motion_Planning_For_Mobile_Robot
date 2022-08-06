# $OBVP\ report$

main function see grid_path_searcher/src/demo_node.cpp and hw_tool.cpp

## 1 Forward integration

```c++
 // 前向积分，每个控制量在[umin, umax]区间分成三份，共27个motion promitive
  for (int i = 0; i <= _discretize_step; i++) { // acc_input_ax
    TraLibrary[i] = new TrajectoryStatePtr *[_discretize_step + 1];
    for (int j = 0; j <= _discretize_step; j++) { // acc_input_ay
      TraLibrary[i][j] = new TrajectoryStatePtr[_discretize_step + 1];
      for (int k = 0; k <= _discretize_step; k++) { // acc_input_az

        vector<Vector3d> Position;
        vector<Vector3d> Velocity;
        acc_input(0) = double(-_max_input_acc + i * (2 * _max_input_acc /
                                                     double(_discretize_step)));
        acc_input(1) = double(-_max_input_acc + j * (2 * _max_input_acc /
                                                     double(_discretize_step)));
        acc_input(2) =
            double(k * (2 * _max_input_acc / double(_discretize_step)) +
                   0.1); // acc_input_az >0.1

		delta_time = _time_interval / double(_time_step);
//每组控制量积分时长为1.25s，每次步长为delta_time
        for (int step = 0; step <= _time_step; step++) {
          pos += vel * delta_time;
          vel += acc_input * delta_time;
          Position.push_back(pos);
          Velocity.push_back(vel);
          double coord_x = pos(0);
          double coord_y = pos(1);
          double coord_z = pos(2);
          // check if if the trajectory face the obstacle
          if (_homework_tool->isObsFree(coord_x, coord_y, coord_z) != 1) {
            collision = true;
          }
        }
      }
```

## 2 OBVP Problem Solve

### 2.1 calculate the close form of derivative of J

用sympy计算J对T的导数的表达式，令其为0,乘以T^4后为4次多项式，得到T每一阶乘方系数的闭式表达，这里参考了助教的github对sympy库进行学习及运用

````python


from sympy import *
init_printing()
px0, py0, pz0, vx0, vy0, vz0, pxf, pyf, pzf, vxf, vyf, vzf, T= symbols('px0 py0 pz0 vx0 vy0 vz0 pxf pyf pzf vxf vyf vzf T')

delta_px = pxf - vx0 * T - px0
delta_py = pyf - vy0 * T - py0
delta_pz = pzf - vz0 * T - pz0
delta_vx = vxf - vx0
delta_vy = vyf - vy0
delta_vz = vzf - vz0

b = Matrix([delta_px, delta_py, delta_pz, delta_vx, delta_vy, delta_vz])

A = Matrix([
    [1/6 * T ** 3, 0, 0, 1/2 * T** 2, 0, 0],
    [0,  1/6 * T ** 3, 0, 0, 1/2 * T ** 2,  0],
    [0, 0, 1/6 * T ** 3, 0, 0, 1/2 * T ** 2],
    [1/2  * T ** 2, 0, 0, T, 0, 0],
    [0, 1/2 * T ** 2, 0, 0, T, 0],
    [0, 0, 1/2 * T ** 2, 0, 0, T]
])

x = simplify(A ** (-1) * b)

M = Matrix(
    [
        [(T**3)/3,        0,        0, (T**2)/2,        0,        0],
        [       0, (T**3)/3,        0,        0, (T**2)/2,        0],
        [       0,        0, (T**3)/3,        0,        0, (T**2)/2],
        [(T**2)/2,        0,        0,        T,        0,        0],
        [       0, (T**2)/2,        0,        0,        T,        0],
        [       0,        0, (T**2)/2,        0,        0,        T],
    ]
)

J = collect(
    expand(
        simplify(
            Transpose(x) * M * x
        )[0] + T
    ),
    syms=T
)

print(J)

dotJ = collect(
    expand(
        simplify(diff(J, T))
    ),
    syms = T
)
print(dotJ)

dotJ_new = collect(
    expand(
        simplify(dotJ * T**4)
    ),
    syms = T
)
````

### 2.2 hw_tool.cpp

利用多项式伴随矩阵求根求得$optimal \ T*$

```c++
const double coef_order3 = 0.0;
  const double coef_order2 =
      1.0 * (-4.0 * vx0 * vx0 - 4.0 * vx0 * vxf - 4.0 * vxf * vxf -
             4.0 * vy0 * vy0 - 4.0 * vy0 * vyf - 4.0 * vyf * vyf -
             4.0 * vz0 * vz0 - 4.0 * vz0 * vzf - 4.0 * vzf * vzf);
  const double coef_order1 =
      1.0 * (-24.0 * px0 * vx0 - 24.0 * px0 * vxf + 24.0 * pxf * vx0 +
             24.0 * pxf * vxf - 24.0 * py0 * vy0 - 24.0 * py0 * vyf +
             24.0 * pyf * vy0 + 24.0 * pyf * vyf - 24.0 * pz0 * vz0 -
             24.0 * pz0 * vzf + 24.0 * pzf * vz0 + 24.0 * pzf * vzf);
  const double coef_order0 =
      1.0 * (-36.0 * px0 * px0 + 72.0 * px0 * pxf - 36.0 * pxf * pxf -
             36.0 * py0 * py0 + 72.0 * py0 * pyf - 36.0 * pyf * pyf -
             36.0 * pz0 * pz0 + 72.0 * pz0 * pzf - 36.0 * pzf * pzf);
  Eigen::MatrixXd mat_adj = Eigen::MatrixXd::Zero(4, 4);

  Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic>
      matrix_eigenvalues;
  mat_adj(0, 3) = coef_order3;
  mat_adj(1, 3) = -1.0 * coef_order2;
  mat_adj(2, 3) = -1.0 * coef_order1;
  mat_adj(3, 3) = -1.0 * coef_order0;
  mat_adj(1, 0) = 1.0;
  mat_adj(2, 1) = 1.0;
  mat_adj(3, 2) = 1.0;

  matrix_eigenvalues = mat_adj.eigenvalues();
```

### 2.3  RVIZ

The green line is the optimal trajectory.

Red lines are invalid trajectory.

![img]([https://github.com/Lkaho/Motion_Planning_For_Mobile_Robot/blob/main/Control%20Lattice/src/OBVP.png)
