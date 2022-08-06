## 一、代码流程梳理，仅包含关键代码

### 1. $test/ planner.cpp$规划函数调用入口

```c++
TesterPathFinder构造函数中初始化了 rrt_star_ptr

shared_ptr<path_plan::RRTStar> rrt_star_ptr_;

TesterPathFinder(const ros::NodeHandle &nh) : nh_(nh)
    {
        env_ptr_.reset(new env::OccMap);
        env_ptr_->init(nh_);
		
    	//rrt* planner initilization
        rrt_star_ptr_.reset(new path_plan::RRTStar(nh_, env_ptr_));
        rrt_star_ptr_->setVisualizer(vis_ptr_);

        goal_sub_ = nh_.subscribe("/goal", 1, &TesterPathFinder::goalCallback, this);
    }

//goaL回调函数触发plan
void goalCallback(const geometry_msgs::PoseStamped::ConstPtr &goal_msg)
    {	
    //planer begin
        bool rrt_star_res = rrt_star_ptr_->plan(start_, goal_);
        if (rrt_star_res)
        {
		  /*..*/
        }
    }
```

### 2. $rrt/ star.h$

#### 2.1 $RRTStart$构造函数

````c++
    RRTStar(const ros::NodeHandle &nh, const env::OccMap::Ptr &mapPtr) : nh_(nh), map_ptr_(mapPtr)
    {
      //为最大采样点数开辟一块内存，之后所有的new node都从 nodes_pool_里边取
      valid_tree_node_nums_ = 0;
      nodes_pool_.resize(max_tree_node_nums_);
      for (int i = 0; i < max_tree_node_nums_; ++i)
      {
        nodes_pool_[i] = new TreeNode;
      }
    }
````

#### 2.2 $plan$函数

````c++
 bool rrt_star(const Eigen::Vector3d &s, const Eigen::Vector3d &g){
      /* kd tree init */
      kdtree *kd_tree = kd_create(3);
      //Add start and goal nodes to kd tree
      kd_insert3(kd_tree, start_node_->x[0], start_node_->x[1], start_node_->x[2], start_node_);
     // main loop
     int idx = 0;
     // 采样次数/采样时间到了就终止
      for (idx = 0; (ros::Time::now() - rrt_start_time).toSec() < search_time_ && valid_tree_node_nums_ < max_tree_node_nums_; ++idx)
      {
          //1、采样
        Eigen::Vector3d x_rand;
        sampler_.samplingOnce(x_rand);
          //2、找采样点最近的点
        struct kdres *p_nearest = kd_nearest3(kd_tree, x_rand[0], x_rand[1], x_rand[2]);
          //数据类型转换，该nearest_node就是初始的父节点
        RRTNode3DPtr nearest_node = (RRTNode3DPtr)kd_res_item_data(p_nearest);
        //生成x_new
        Eigen::Vector3d x_new = steer(nearest_node->x, x_rand, steer_length_);
          
        //3、找x_new某个range范围内的节点
        vector<RRTNode3DPtr> neighbour_nodes;
        struct kdres *nbr_set;
        nbr_set = kd_nearest_range3(kd_tree, x_new[0], x_new[1], x_new[2], search_radius_);
          
         //数据类型转换，压入neighborlist
       while (!kd_res_end(nbr_set))
        {
          RRTNode3DPtr curr_node = (RRTNode3DPtr)kd_res_item_data(nbr_set);
          neighbour_nodes.emplace_back(curr_node);
        }
          
        //4、对x_new的一些数据初始化：cost_from_p, cost_from_start
         double dist2nearest = calDist(nearest_node->x, x_new);
        double min_dist_from_start(nearest_node->cost_from_start + dist2nearest);
        double cost_from_p(dist2nearest);
        //min_node is default parent of x_new
        RRTNode3DPtr min_node(nearest_node); 
          
        //5、对临近节点集遍历，找到代价更小的cost_to_come
          //store_result保存了x_nee与临近节点的距离与碰撞查询结果，rewire的时候就不用再做一遍
        std::vector<std::pair<double, bool>> store_result;
        for (auto &curr_node : neighbour_nodes)
        {
          double cost_frome_curr_node = 0;
          bool isCollision = true;
          if(map_ptr_->isSegmentValid(curr_node->x, x_new))
          {
            isCollision = false;
            cost_frome_curr_node = calDist(curr_node->x, x_new);
            double cost_from_start_xnew  = cost_frome_curr_node + curr_node->cost_from_start;
            if(cost_from_start_xnew < min_dist_from_start){
              min_dist_from_start = cost_from_start_xnew;
              min_node = curr_node;
              cost_from_p  = cost_frome_curr_node;
            }
          }
          store_result.push_back(make_pair(cost_frome_curr_node, isCollision));
        }
          
        //6、得到x_new最优的解后，便可以将该点加入rrt树与kd tree中
        RRTNode3DPtr new_node(nullptr);
        new_node = addTreeNode(min_node, x_new, min_dist_from_start, cost_from_p);
        kd_insert3(kd_tree, x_new[0], x_new[1], x_new[2], new_node);
      }
     //7、尝试连接终点，判断
     /*
     	...
     */
     //8、剪枝
     //遍历邻近节点集，若其从x_new到其的cost_to_come值减小，则将其父节点设为x_new
     //这里加了一个判断：若f_new + h > curent_best_cost，则不对这个点重连
     int i = 0;
      for (auto &curr_node : neighbour_nodes)
        {
          double best_cost_before_rewire = goal_node_->cost_from_start;
          // ! -------------------------------------
          if(store_result[i].second) continue;
          double cost_to_come_cost_new = store_result[i].first + new_node->cost_from_start;
          double heur_val = heuristicFun(new_node->x, curr_node->x);
          if(cost_to_come_cost_new < curr_node->cost_from_start){
            if(cost_to_come_cost_new + heur_val > best_cost_before_rewire) continue;
            changeNodeParent(curr_node, new_node, store_result[i].first);
          }
          // ! -------------------------------------
          if (best_cost_before_rewire > goal_node_->cost_from_start)
          {
            vector<Eigen::Vector3d> curr_best_path;
            fillPath(goal_node_, curr_best_path);
            path_list_.emplace_back(curr_best_path);
            solution_cost_time_pair_list_.emplace_back(goal_node_->cost_from_start, (ros::Time::now() - rrt_start_time).toSec());
          }
          ++i;
        }
 }
````

## 二、结果显示

![rrt*](/home/kaho/sample_method/src/hw_picture_and_report/rrt*.png)

## 三、提升效果

### 3.1 根据Informed RRT* 论文改动采样函数

```c++
void informedSampleOnce(Eigen::Vector3d &sample, const Eigen::Vector3d &target, 
                          const Eigen::Vector3d &start, const double c_min, const double c_best)
  {
    sample[0] = normal_rand_(gen_);
    sample[1] = normal_rand_(gen_);
    sample[2] = normal_rand_(gen_);
    double r = pow(uniform_rand_(gen_), 1.0 / 3.0);
    sample = r * sample.normalized();
    Eigen::Vector3d x_center = (target + start) / 2;

  // transform Matrix: L
    Eigen::Matrix3d L;
    double L_diag1 = c_best / 2;
    double L_diag2 = sqrt(pow(c_best,2)-pow(c_min,2)) / 2;
    double L_diag3(L_diag2);
    
    L << L_diag1, 0, 0,
         0, L_diag2, 0,
         0, 0, L_diag3;
  
  // Rotation Matrix: C
    Eigen::MatrixXd M;
    Eigen::Vector3d I(1, 0, 0);
    Eigen::Vector3d a1 = (target - start).normalized();
    M = a1 * I.transpose();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV ); 
    Eigen::Matrix3d V = svd.matrixV(), U = svd.matrixU(); 
    Eigen::Matrix3d C, diag;
    double C_diag3 = U.determinant() * V.determinant();
    diag << 1, 0, 0,
            0, 1, 0,
            0, 0, C_diag3;
    C = U * diag * V.transpose();
    //get sample point in the hyperellipsoid:
    sample = C * L * sample + x_center;
  };
```

### 3.2 改进效果演示

![informedRRT*](/home/kaho/sample_method/src/hw_picture_and_report/informedRRT*.png)