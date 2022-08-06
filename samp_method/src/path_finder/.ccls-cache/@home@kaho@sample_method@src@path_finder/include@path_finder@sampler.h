/*
Copyright (C) 2022 Hongkai Ye (kyle_yeh@163.com)
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.
*/
#ifndef _BIAS_SAMPLER_
#define _BIAS_SAMPLER_

#include <ros/ros.h>
#include <Eigen/Eigen>
#include <Eigen/SVD>
#include <random>

 


class BiasSampler
{
public:
  BiasSampler()
  {
    std::random_device rd; //a random number engine used to seed the Mersenne twister engine
    gen_ = std::mt19937_64(rd());
    uniform_rand_ = std::uniform_real_distribution<double>(0.0, 1.0);
    normal_rand_ = std::normal_distribution<double>(0.0, 1.0);
    range_.setZero();
    origin_.setZero();
  };

  void setSamplingRange(const Eigen::Vector3d origin, const Eigen::Vector3d range)
  {
    origin_ = origin;
    range_ = range;
  }

  void samplingOnce(Eigen::Vector3d &sample)
  {
    sample[0] = uniform_rand_(gen_);
    sample[1] = uniform_rand_(gen_);
    sample[2] = uniform_rand_(gen_);
    //scale the sample point to the sampling range:
    //sample[0] = range_[0] * sample[0];
    //sample[1] = range_[1] * sample[1];
    //sample[2] = range_[2] * sample[2];
    sample.array() *= range_.array();  //scale the sample point to the sampling range: sma
    sample += origin_;
  };
  
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

  // (0.0 - 1.0)
  double getUniRandNum()
  {
    return uniform_rand_(gen_);
  }

private:
  Eigen::Vector3d range_, origin_;
  std::mt19937_64 gen_; //6
  std::uniform_real_distribution<double> uniform_rand_;
  std::normal_distribution<double> normal_rand_;
};

#endif
