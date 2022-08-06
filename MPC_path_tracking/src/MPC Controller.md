# MPC based Trajectory Following  Controller

## 1、MPC Controller Design

### 1.1 System Model 

Due to low vehicle speed, we ultilze the kinematic model as the system model:
$$
\left\{\begin{aligned}
\dot{p}_{x} &=v \cos (\phi) \\
\dot{p}_{y} &=v \sin (\phi) \\
\dot{\phi} &=\frac{v}{L} \tan (\delta) \\
\dot{v} &=a
\end{aligned}\right.
$$
State: $X = [x\ y \ \phi\ v]^T$, Input: $U = [a\ \delta]$

- $Linearization$

  For improving the solving efficiency, we linearize the nonlinear system model. There are two different method.

  - Linearization along the state trajectory:

    We can linearize the system model based on the result of last run.

  - ###### Linearization along the reference trajectory

    we ca linearize the system mode based on the reference state.

Define the linearization point as follows:

$\bar{X} = [\bar{p_x}\ \bar{p_y}\ \bar{phi}\ \bar{v}]^T$ ,	$\bar{U} = [\bar{a}\ \bar{\delta}]^T$

The linear model:
$$
\left[\begin{array}{c}
\dot{p_{x}} \\
\dot{p_{y}} \\
\dot{\phi} \\
\dot{v}
\end{array}\right]=\left(\begin{array}{cccc}
0 & 0 & -\bar{v} \sin \bar{\phi} & \cos \bar{\phi} \\
0 & 0 & \bar{v} \cos \bar{\phi} & \sin \bar{\phi} \\
0 & 0 & 0 & \frac{\tan \bar{\delta}}{L} \\
0 & 0 & 0 & 0
\end{array}\right)\left[\begin{array}{c}
p_{x} \\
p_{y} \\
v \\
\phi
\end{array}\right]+\left(\begin{array}{cc}
0 & 0 \\
0 & 0 \\
0 & \frac{\bar{v}}{L} \frac{1}{\cos ^{2} \bar{\delta}} \\
1 & 0
\end{array}\right)\left[\begin{array}{c}
a \\
\delta
\end{array}\right]+\left(\begin{array}{c}
\bar{v} \sin \bar{\phi} \\
-\bar{v} \bar{\phi} \cos \bar{\phi} \\
-\bar{v} \frac{\bar{\delta}}{L \cos ^{2} \bar{\delta}} \\
0
\end{array}\right)
$$

- Discretizaton

  For computer calculation, we discretize the continuous system with forward Euler difference:
  $$
  \begin{aligned}
  A_{\text {discrete }} &=I+\Delta T \times A_{\text {linearized }} \\
  B_{\text {discrete }} &=\Delta T \times B_{\text {linearized }} \\
  g_{\text {discrete }} &=\Delta T \times g_{\text {linearized }}
  \end{aligned}
  $$

### 1.2 Define Cost Function

The objective is to follow the trajectory accuracy, the cost function can define as follow:
$$
\begin{aligned}
&\min \sum_{k=1}^{N}(x_k - \bar x_k)^TQ(x_k-\bar x_k)+ u_{k-1}^TRu_{k-1} \\
&\text { s.t. } x_{0}=x \text {, }\\
&x_{k+1}=f\left(x_{k}, u_{k}\right), \quad k=0, \ldots, N-1\\
&\underline{x} \leq x_{k} \leq \bar{x}, \quad k=1, \ldots, N\\
&\underline{u} \leq u_{k} \leq \bar{u}, \quad k=0, \ldots, N-1
\end{aligned}
$$

### 1.3 Prediction

$$
\eta = AA x_0 + BBU +G
$$

$\eta = [x_1\ x_2 \ ... x_N]^T$, $U = [u_0\ u_1 \ u_2 ... u_{N-1}]^T$
$$
AA=\left[\begin{array}{c}
A_{0} \\
A_{1} A_{0} \\
\vdots \\
\prod_{k=0}^{N-1} A_{k}
\end{array}\right]
$$

$$
BB=\left[\begin{array}{ccccc}
B_{0} & 0 & \ldots & 0 & 0 \\
A_{1} B_{0} & B_{1} & \ldots & 0 & 0 \\
\vdots & \vdots & \ddots & \vdots & \vdots \\
\prod_{k=1}^{N-1} A_{k} B_{0} & \prod_{k=2}^{N-1} A_{k} B_{1} & \ldots & A_{N-1} B_{N-2} & B_{N-1}
\end{array}\right]
$$

$$
G=\left[\begin{array}{c}
g_{0} \\
A_{1} g_{0}+g_{1} \\
\vdots \\
\sum_{n=0}^{N-2}\left(\prod_{k=n+1}^{N-1} A_{k}\right) g_{n}+g_{N-1}
\end{array}\right]
$$

Thus we can reomove the equility constrain of the optimize problem.

### 1.4 OSQP Solver

we adopt the QP solver OSQP to solve the QP problem.

OSQP solves convex quadratic programs (QPs) of the form
$$
\begin{array}{ll}
\operatorname{minimize} & \frac{1}{2} x^{T} P x+q^{T} x \\
\text { subject to } & l \leq A x \leq u
\end{array}
$$
where $x \in \mathbf{R}^{n}$ is the optimization variable. The objective function is defined by a positive semidefinite matrix $P \in \mathbf{S}_{+}^{n}$ and vector $q \in \mathbf{R}^{n}$. The linear constraints are defined by matrix $A \in \mathbf{R}^{m \times n}$ and vectors $l$ and $u$ so that $l_{i} \in \mathbf{R} \cup\{-\infty\}$ and $u_{i} \in \mathbf{R} \cup\{+\infty\}$ for all $i \in\{1, \ldots, m\}$.

Subtitutiing the predciton model into the cost function, we have
$$
U^T(BB^TQBB)U + 2(x_0AA^TQBB+GG^TQBB-q_xQBB)U
$$
$q_x = [\bar x_1,\ ...,\bar x_N]$

## 2 、Result

![out](/home/kaho/mpCourseHw/mpc_ws/src/out.gif)
