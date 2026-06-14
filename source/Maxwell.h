//
// 文件名:  Maxwell.h
// 软件包:  JAUMIN applications
// 版权  :  (c) 2004-2015 北京应用物理与计算数学研究所
//          (c) 2013-2015 中物院高性能数值模拟软件中心
// 版本号:  $Revision$
// 修改  :  $Data$
// 描述  :  麦克斯维方程的网格片策略类.
//

#ifndef included_Maxwell
#define included_Maxwell

#include "fem/GridInfo.h"
#include "utils/JAUMIN_App.h"
#include "utils/JAUMIN_Macros.h"

using namespace std;
using namespace JAUMIN;

template <class TYPE> class UserParam {
public:
  enum { DIM = 3 };
  enum FunctionType { INVALID, SOURCE, BOUNDARY_1, BOUNDARY_3 };
  FunctionType ftype;

  int inner_face_set_id;
  int outer_face_set_id;
  double frequency; // so defined k_0
  double epsilon_r;
  double mu_r;
  double eta;
  double E_0_direction[DIM];
  double k_0_direction[DIM];
  double chi[DIM];

public:
  UserParam() {
    inner_face_set_id = 9;
    outer_face_set_id = 10;
    frequency = 3.0e8;
    epsilon_r = 1.0;
    mu_r = 1.0;
    eta = 1.0;
    E_0_direction[0] = 1;
    E_0_direction[1] = 0;
    E_0_direction[2] = 0;
    k_0_direction[0] = 0;
    k_0_direction[1] = 0;
    k_0_direction[2] = 1;
    chi[0] = 0;
    chi[1] = 0;
    chi[2] = 0;
    ftype = INVALID;
  }
  void getFromInput(tbox::Pointer<tbox::Database> db) {
    inner_face_set_id = db->getIntegerWithDefault("inner_face_set_id", inner_face_set_id);
    outer_face_set_id = db->getIntegerWithDefault("outer_face_set_id", outer_face_set_id);
    frequency = db->getDoubleWithDefault("frequency", frequency);
    epsilon_r = db->getDoubleWithDefault("epsilon_r", epsilon_r);
    mu_r = db->getDoubleWithDefault("mu_r", mu_r);
    eta = db->getDoubleWithDefault("eta", eta);
    if (db->keyExists("E_0_direction")) {
      db->getDoubleArray("E_0_direction", E_0_direction, DIM);
    }
    if (db->keyExists("k_0_direction")) {
      db->getDoubleArray("k_0_direction", k_0_direction, DIM);
    }
    if (db->keyExists("chi")) {
      db->getDoubleArray("chi", chi, DIM);
    }
  }
  double get_k_0() const {
    double pi = 3.14159265;
    double c = 3.0e8;
    return 2 * pi * frequency / c;
  }
  void source(double x, double y, double z, TYPE *value) const { value[0] = 0.0; }
  void boundary_1(double x, double y, double z, TYPE *value) const { value[0] = 0.0; }
  void boundary_3(double x, double y, double z, TYPE *value) const {
    double R[DIM], normal[DIM];
    R[0] = x;
    R[1] = y;
    R[2] = z;
    double sum = dotProduct(DIM, R, R);
    for (int i = 0; i < DIM; i++) {
      normal[i] = R[i] / sqrt(sum);
    }
    dcomplex coef = dcomplex(0, 1) * get_k_0() *
                    exp(dotProduct(DIM, k_0_direction, R) * get_k_0() * dcomplex(0, 1));

    double temp1_1[DIM], temp1_2[DIM];
    crossProduct(DIM, k_0_direction, E_0_direction, temp1_1);
    crossProduct(DIM, normal, temp1_1, temp1_2);

    double temp2_1[DIM], temp2_2[DIM];
    crossProduct(DIM, normal, E_0_direction, temp2_1);
    crossProduct(DIM, normal, temp2_1, temp2_2);

    for (int i = 0; i < DIM; i++) {
      value[i] = chi[i] + (coef / mu_r) * temp1_2[i] + (coef / eta) * temp2_2[i];
    }
  }
  void operator()(double x, double y, double z, TYPE *value) const {
    switch (ftype) {
    case SOURCE:
      source(x, y, z, value);
      break;
    case BOUNDARY_1:
      boundary_1(x, y, z, value);
      break;
    case BOUNDARY_3:
      boundary_3(x, y, z, value);
      break;
    default:
      assert(false);
    }
  }
};

/**
 * @brief 该类从标准构件网格片策略类 algs::StandardComponentPatchStrategy 派生,
 * 实现求解麦克斯维方程的数值计算子程序.
 *
 * The time-harmonic Maxwell equations:
 *	curl curl E - k^2 E = J in Omega,
 *	E x n = g x n at boundary.
 */
class Maxwell : public algs::StandardComponentPatchStrategy<NDIM> {
public:
  /*!
   * @brief 构造函数.
   * @param object_name 输入参数, 字符串, 表示对象名称.
   * @param input_db    输入参数, 指针,   指向输入数据库.
   *
   * @note
   * 该函数主要完成以下操作:
   *  -# 初始化内部数据成员;
   *  -# 定义变量和数据片.
   */
  Maxwell(const string &object_name, tbox::Pointer<tbox::Database> input_db);

  /*!
   * @brief 析构函数.
   */
  virtual ~Maxwell();

  /// @name 重载基类 algs::StandardComponentPatchStrategy<NDIM> 的函数:
  // @{

  /*!
   * @brief 初始化指定的积分构件.
   *
   * 注册待填充的数据片或待调度内存空间的数据片到积分构件.
   *
   * @param component 输入参数, 指针, 指向待初始化的积分构件对象.
   */
  void initializeComponent(algs::IntegratorComponent<NDIM> *component) const;

  /**
   * @brief 初始化数据片（支持初值构件）.
   *
   * @param patch          输入参数, 网格片类,     表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示初始化的时刻.
   * @param initial_time   输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
   * @param component_name 输入参数, 字符串, 当前调用该函数的初值构件的名称.
   */
  void initializePatchData(hier::Patch<NDIM> &patch, const double time, const bool initial_time,
                           const string &component_name);

  /**
   * @brief  在单个网格片上执行归约计算(支撑归约构件).
   *
   * @param vector         输入输出参数, 指针, 指向归约向量.
   * @param len            输入参数, 整型, 归约向量的长度.
   * @param patch          输入参数, 网格片类, 待归约的网格片.
   * @param time           输入参数, 双精度浮点型, 表示当前时刻.
   * @param dt             输入参数, 双精度浮点型, 表示当前时间步长.
   * @param component_name 输入参数, 字符串, 归约构件的名称.
   *
   * @note
   *  vector是输入输出参数: 输入是已遍历网格片的归约结果,
   *  输出是输入值和当前网格片计算结果的归约值.
   */
  void reduceOnPatch(double *vector, int len, hier::Patch<NDIM> &patch, const double time,
                     const double dt, const string &component_name);

  /*!
   * @brief 完成单个网格片上的数值计算（支持数值构件）.
   *
   * 该函数基于显式迎风格式，实现通量和守恒量的计算。
   *
   * @param patch          输入参数, 网格片类,     表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示当前时刻.
   * @param dt             输入参数, 双精度浮点型, 表示时间步长.
   * @param initial_time   输入参数, 逻辑型,       表示当前是否为初始时刻.
   * @param component_name 输入参数, 字符串, 当前调用该函数的数值构件的名称.
   */
  void computeOnPatch(hier::Patch<NDIM> &patch, const double time, const double dt,
                      const bool initial_time, const string &component_name);

  //@}

  ///@name 自定义函数
  //@{

  /*!
   * @brief 注册绘图量.
   * @param javis_writer 输入参数, 指针, 表示 JaVis 数据输出器.
   */
  void registerPlotData(tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_writer);

  //@}

private:
  /*! declare MaxwellLevelIntegrator as friend. */
  friend class MaxwellLevelIntegrator;

  /*!@brief 从输入数据库读入数据.  */
  void getFromInput(tbox::Pointer<tbox::Database> db);
  /*!@brief 注册变量和数据片.  */
  void registerModelVariables();

  /*! initialize grid geometry.  */
  void initGeometryOnPatch(hier::Patch<NDIM> &patch);
  /*! compute source term.  */
  void initJOnPatch(hier::Patch<NDIM> &patch);
  /*! initialize analytical solution.  */
  void initAnalyticOnPatch(hier::Patch<NDIM> &patch);
  /*! setup dof map on patch */
  void setupDofInfoOnPatch(hier::Patch<NDIM> &patch);
  /*! generate the auxillary variables for AMS */
  void buildAMSVarOnPatch(hier::Patch<NDIM> &patch);
  /*! assemble matrix and rhs */
  void buildLHSAndRHSOnPatch(hier::Patch<NDIM> &patch, const double time, const double dt);
  void buildEMMatrixOnPatch(hier::Patch<NDIM> &patch, const double time, const double dt);
  /*! accept solution from solution vector */
  void acceptSolutionOnPatch(hier::Patch<NDIM> &patch);
  /*! postprocess for plotting */
  void postprocessOnPatch(hier::Patch<NDIM> &patch);
  /*! compute edge order */
  void computeEdgeOrderOnPatch(hier::Patch<NDIM> &patch);
  /*! compute edge flag */
  void computeEdgeFlagOnPatch(hier::Patch<NDIM> &patch);
  /*! compute Hcurl error */
  void computeHcurlErrOnPatch(double *vector, int len, hier::Patch<NDIM> &patch);

  /// Added by Yin-Da Wang
  void setupWPModalDofInfoOnPatch(hier::Patch<NDIM> &patch);
  void setupModalMatrixOnPatch(hier::Patch<NDIM> &patch, double time, const double dt);
  void updateMaterialProperties(hier::Patch<NDIM> &patch);
  void calculateModalPower(double *vector, int len, hier::Patch<NDIM> &patch, const double time,
                           const double dt);

  string d_object_name; /*!@brief 对象名.  */

  UserParam<dcomplex> d_user_param; // 方程的所有可变参数。
  int entity_face_num = 8;
  int simulation_setup;
  tbox::Array<double> frequency_setup;
  /// @brief  材料数据库变量
  int num_total_material;
  tbox::Array<tbox::Array<int> > material_entity;
  class Material {
  private:
  public:
    int mat_id; // 材料编号
    double d_mu_r;
    double d_epsilon_r;
    double d_sigma;
    double d_losstan;
  };
  tbox::Array<tbox::Pointer<Material> > Material_list;
  void InitMaterialComponent(tbox::Pointer<tbox::Database> db);

  /// @brief 边界条件数据库变量
  tbox::Array<int> ABC_face;
  tbox::Array<int> PEC_face;
  void InitBoundaryInfo(tbox::Pointer<tbox::Database> db);

  int analytical_waveport_num;
  int numerical_waveport_num;
  int lumpport_num;
  class ANAL_WP {
  private:
  public:
    double a;
    double b;
    double epsilon_r;
    tbox::Array<double> coordinate;
    tbox::Array<int> face_list;
    tbox::Array<double> direction;
    double power;
    bool activated;
  };
  class NUM_WP {
  private:
  public:
    tbox::Array<int> face_list;
    tbox::Array<int> edge_list;
    tbox::Array<double> direction;
    double power;
    bool activated;
  };
  class LP {
  private:
  public:
    double width;
    double height;
    tbox::Array<int> face_list;
    double re_impedance;
    double im_impedance;
    double voltage;
  };
  tbox::Array<tbox::Pointer<ANAL_WP> > ANAL_WP_list;
  tbox::Array<tbox::Pointer<NUM_WP> > NUM_WP_list;
  tbox::Array<tbox::Pointer<LP> > LP_list;

  int eigen_num = 2;
  tbox::Array<tbox::Array<double> > eigen_value_pairs;
  int active_mode = 0;

  tbox::Array<dcomplex> port_power;

  tbox::Array<int> d_matA_modal_array; // 模式矩阵的大小与端口数一致
  tbox::Array<int> d_matB_modal_array; // 模式矩阵的大小与端口数一致
  tbox::Array<tbox::Array<int> >
      d_sol_modal_array; // 外层大小和端口数一致，内层大小与特征值数量一致

  tbox::Array<tbox::Pointer<solv::DOFInfo<NDIM> > > d_dof_info_modal_array;

  void InitPortInfo(tbox::Pointer<tbox::Database> db);

  int d_E_id;       // electric field on edge.
  int d_Ec_id;      // electric field on node(for plotting).
  int d_Ec_norm_id; // electric field on node(for plotting).
  int d_Ec_arg_id;  // electric field on node(for plotting).

  int d_Ea_id;       // electric field, analytic solution.
  int d_Eac_id;      // electric field, analytic solution, on node(for plotting).
  int d_J_id;        // displacement current, the source term.
  int d_jacobian_id; // cell(tetrahedron) jacobian.
  int d_volume_id;   // cell(tetrahedron) volume.

  int d_edge_order_id; // edge order, consistent across every patch.
  int d_edge_flag_id;  // edge flag, whether the edge is on constraint boundary.
  int d_edge_owner_id; // edge owner, whether the edge is owned by the patch.

  int d_port_edge_flag_id;
  /*
  int d_matrix_id;   // patch matrix.
  int d_rhs_id;      // patch vector, the right hand side.
  int d_solution_id; // patch vector, the solution.
  */

  // 材料参数分配编号
  int d_epsilonr_id;
  int d_mur_id;
  int d_sigma_id;
  int d_losstan_id;
  int d_mat_id;

  int d_mat_re_id;
  int d_mat_im_id;
  int d_rhs_re_id;
  int d_rhs_im_id;
  int d_sol_re_id;
  int d_sol_im_id;
  int d_premat_id;
  tbox::Pointer<solv::DOFInfo<NDIM> > d_dof_info; // DOF distribution.

  // AMS 方法中的辅助边量和顶点自由度分布信息
  int d_gmat_id;   // 边基函数和节点基函数的梯度的映射关系
  int d_xcoord_id; // x方向坐标向量
  int d_ycoord_id; // y方向坐标向量
  int d_zcoord_id; // z方向坐标向量
  tbox::Pointer<solv::DOFInfo<NDIM> > d_dof_info_node; // DOF distribution.
};

#endif
