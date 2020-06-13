#ifndef GMMREG_BASE_H_
#define GMMREG_BASE_H_

#ifdef WIN32
#include <windows.h>
#else
#include "port_ini.h"
#endif

#include <vector>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>

namespace gmmreg {

class Base {
 public:
  Base() {
    strcpy(common_section_, "FILES");
  }

  void Run(const char* f_config);
  virtual void PerformTransform(const vnl_vector<double>&) = 0;
  virtual double BendingEnergy() = 0;  // serving as a regularization term
  virtual void ComputeGradient(const double lambda,
      const vnl_matrix<double>& gradient, vnl_matrix<double>& grad_all) = 0;

  virtual ~Base() {}

 protected:
  /* m: # of points in model */
  /* s: # of points in scene */
  /* n: # of points in ctrl_pts */
  /* d: dimensionality, e.g. 2 for 2D points, 3 for 3D points */
  int m_, n_, s_, d_;

  // each row is a sample point
  vnl_matrix<double> model_, scene_, ctrl_pts_, transformed_model_;

  double sigma_, lambda_;
  vnl_matrix<double> kernel_;
  int b_normalize_;
  vnl_vector<double> model_centroid_, scene_centroid_;
  char section_[80], common_section_[80];

  unsigned int level_;
  std::vector<float> v_scale_;
  std::vector<int> v_func_evals_;

  // load input data from files
  virtual int PrepareInput(const char* input_config);
  int SetCtrlPts(const char* filename);
  void SaveTransformed(const char* filename,
      const vnl_vector<double>&, const char* f_config);
  void MultiScaleOptions(const char* f_config);
  void DenormalizeAll();
  int Initialize(const char* f_config);

  friend class RigidFunc;
  friend class ThinPlateSplineFunc;
  friend class GaussianRadialBasisFunc;

 private:
  double model_scale_, scene_scale_;
  void PrepareCommonOptions(const char* f_config);
  virtual void PrepareOwnOptions(const char* f_config) = 0;
  virtual void PrepareBasisKernel() = 0;
  virtual int SetInitParams(const char* filename) = 0;
  virtual void SaveResults(const char* filename,
      const vnl_vector<double>&) = 0;
  virtual void StartRegistration(vnl_vector<double>& params) = 0;
};

}  // namespace gmmreg

#endif  // GMMREG_BASE_H_
