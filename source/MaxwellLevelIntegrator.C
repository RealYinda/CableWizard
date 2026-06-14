//
// 文件名:  MaxwellLevelIntegrator.C
// 软件包:  JAUMIN applications
// 版权  :  (c) 2004-2015 北京应用物理与计算数学研究所
//          (c) 2013-2015 中物院高性能数值模拟软件中心
// 版本号:  $Revision$
// 修改  :  $Data$
// 描述  :  麦克斯维问题的网格层时间积分算法的实现.
//

#include "utils/JAUMIN_App.h"
#include <iostream>
#include <map>
#include <stdlib.h>
// #include "MatrixVectorOperator.h"
using namespace JAUMIN;

#include "BaseEigenSolver.h"
#include "MaxwellLevelIntegrator.h"

#include "JAUMINMVOps.h"

/*************************************************************************
 * 构造函数.
 *************************************************************************/
MaxwellLevelIntegrator::MaxwellLevelIntegrator(const string &object_name, Maxwell *patch_strategy,
                                               tbox::Pointer<tbox::Database> solver_db) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!object_name.empty());
  TBOX_ASSERT(patch_strategy != NULL);
#endif

  d_object_name = object_name;
  d_patch_strategy = patch_strategy;

  // 从数据库中获取解法器类型.
  d_solver_db = solver_db;
  //  d_solver = new JPSOL::JPLinearSolver<NDIM,dcomplex>(d_solver_db);

  d_solver = new JPSOL::ComplexLinearSolver<NDIM, double>(d_solver_db);
}

/*************************************************************************
 * 析构函数.
 ************************************************************************/
MaxwellLevelIntegrator::~MaxwellLevelIntegrator() {}

/*************************************************************************
 * 初始化该积分算法: 创建所有计算需要的积分构件.
 * 这些构件所操作的数据片, 由函数 d_patch_strategy->initializeComponent() 指定.
 *************************************************************************/
void MaxwellLevelIntegrator::initializeLevelIntegrator(
    tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager) {
  // 初值构件
  d_ini_intc["INIT_FIELD"] =
      new algs::InitializeIntegratorComponent<NDIM>("INIT_FIELD", d_patch_strategy, manager);
  d_ini_intc["INIT_DOFINFO"] =
      new algs::InitializeIntegratorComponent<NDIM>("INIT_DOFINFO", d_patch_strategy, manager);

  // 内存构件
  d_mem_intc["ALLOC_MATVEC"] =
      new algs::MemoryIntegratorComponent<NDIM>("ALLOC_MATVEC", d_patch_strategy, manager);
  d_mem_intc["ALLOC_MODAL_MATVEC"] =
      new algs::MemoryIntegratorComponent<NDIM>("ALLOC_MODAL_MATVEC", d_patch_strategy, manager);

  // 数值构件
  d_num_intc["DOF_MAP"] =
      new algs::NumericalIntegratorComponent<NDIM>("DOF_MAP", d_patch_strategy, manager);
  d_num_intc["INIT_AMSVAR"] =
      new algs::NumericalIntegratorComponent<NDIM>("INIT_AMSVAR", d_patch_strategy, manager);
  d_num_intc["COMPUTE_LHSRHS"] =
      new algs::NumericalIntegratorComponent<NDIM>("COMPUTE_LHSRHS", d_patch_strategy, manager);
  d_num_intc["INIT_GEOM"] =
      new algs::NumericalIntegratorComponent<NDIM>("INIT_GEOM", d_patch_strategy, manager);
  d_num_intc["INIT_ANAL"] =
      new algs::NumericalIntegratorComponent<NDIM>("INIT_ANAL", d_patch_strategy, manager);
  d_num_intc["INIT_J"] =
      new algs::NumericalIntegratorComponent<NDIM>("INIT_J", d_patch_strategy, manager);
  d_num_intc["ACCEPT_SOLU"] =
      new algs::NumericalIntegratorComponent<NDIM>("ACCEPT_SOLU", d_patch_strategy, manager);
  d_num_intc["POSTPROCESS"] =
      new algs::NumericalIntegratorComponent<NDIM>("POSTPROCESS", d_patch_strategy, manager);
  d_num_intc["EDGE_ORDER"] =
      new algs::NumericalIntegratorComponent<NDIM>("EDGE_ORDER", d_patch_strategy, manager);
  d_num_intc["EDGE_FLAG"] =
      new algs::NumericalIntegratorComponent<NDIM>("EDGE_FLAG", d_patch_strategy, manager);

  d_num_intc["MODAL_DOF_MAP"] =
      new algs::NumericalIntegratorComponent<NDIM>("MODAL_DOF_MAP", d_patch_strategy, manager);
  d_num_intc["MODAL_MAT"] =
      new algs::NumericalIntegratorComponent<NDIM>("MODAL_MAT", d_patch_strategy, manager);

  // 数值(通信)构件
  d_com_intc["COMM_EDGE_FLAG"] =
      new algs::NumericalIntegratorComponent<NDIM>("COMM_EDGE_FLAG", d_patch_strategy, manager);
  d_com_intc["COMM_WP_POWER"] =
      new algs::NumericalIntegratorComponent<NDIM>("COMM_WP_POWER", d_patch_strategy, manager);

  // 规约构件：计算Hcurl误差.
  d_red_intc["SUM_ERR"] =
      new algs::ReductionIntegratorComponent<NDIM>("SUM_ERR", MPI_SUM, d_patch_strategy, manager);

  // 规约构件：计算端口的功率.
  d_red_intc["WP_POWER"] =
      new algs::ReductionIntegratorComponent<NDIM>("WP_POWER", MPI_SUM, d_patch_strategy, manager);
}

/*************************************************************************
 * 初始化指定网格层的数据片.
 ************************************************************************/
void MaxwellLevelIntegrator::initializeLevelData(
    const tbox::Pointer<hier::BasePatchLevel<NDIM> > level, const double init_data_time,
    const bool initial_time) {
  tbox::Pointer<hier::PatchLevel<NDIM> > patch_level = level;

  if (initial_time) {
    // Allocate space for field data.
    tbox::pout << "Initialize on level......" << endl;
    d_ini_intc["INIT_FIELD"]->initializeLevelData(level, init_data_time, initial_time);

    // Setup DofInfo for each patch.
    d_ini_intc["INIT_DOFINFO"]->initializeLevelData(level, init_data_time, initial_time);

    // Compute globaly consistent edge order, required by Nedelec element.
    d_num_intc["EDGE_ORDER"]->computing(patch_level, init_data_time, 0.);
    // Compute cell volume, jacobian, etc.
    d_num_intc["INIT_GEOM"]->computing(patch_level, init_data_time, 0.);
    // Initialize the analytic solution patch data.
    d_num_intc["INIT_ANAL"]->computing(patch_level, init_data_time, 0.);
    // Initialize the source term(displacement current) patch data.
    d_num_intc["INIT_J"]->computing(patch_level, init_data_time, 0.);

    // Mark boundary edges.
    d_num_intc["EDGE_FLAG"]->computing(patch_level, init_data_time, 0.);
    // Fill ghost region of edge flag.
    d_com_intc["COMM_EDGE_FLAG"]->computing(patch_level, init_data_time, 0.);
  }
}

/*************************************************************************
 * 计算时间步长.
 * 稳态问题，只计算一个时间步，时间步长为1.
 ************************************************************************/
double MaxwellLevelIntegrator::getLevelDt(const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
                                          const double dt_time, const bool initial_time,
                                          const int flag_last_dt, const double last_dt) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!level.isNull());
#endif
  return 1.0;
}

/*************************************************************************
 * 向前积分一个时间步.
 ************************************************************************/
int MaxwellLevelIntegrator::advanceLevel(const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
                                         const double current_time, const double predict_dt,
                                         const double max_dt, const double min_dt,
                                         const bool first_step, const int step_number,
                                         double &actual_dt) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!level.isNull());
#endif

  tbox::Pointer<hier::PatchLevel<NDIM> > patch_level = level;
  tbox::pout << "assemble matrix:  " << std::endl;
  // Set dof map
  d_num_intc["DOF_MAP"]->computing(patch_level, current_time, actual_dt);
  tbox::pout << "setup modal:  " << std::endl;
  d_num_intc["MODAL_DOF_MAP"]->computing(patch_level, current_time, actual_dt);
  d_mem_intc["ALLOC_MODAL_MATVEC"]->allocatePatchData(patch_level, current_time + predict_dt);
  d_num_intc["MODAL_MAT"]->computing(patch_level, current_time, actual_dt);
  const int port_num = d_patch_strategy->numerical_waveport_num;
  for (int pp = 0; pp < port_num; pp++) {
    int modal_A_mat = d_patch_strategy->d_matA_modal_array[pp];
    int modal_B_mat = d_patch_strategy->d_matB_modal_array[pp];
    tbox::Array<int> modal_sol_id = d_patch_strategy->d_sol_modal_array[pp];
    tbox::Pointer<JPSOL::JMatrix<NDIM, double> > d_A_matrix =
        new JPSOL::JMatrix<NDIM, double>(patch_level, modal_A_mat);
    tbox::Pointer<JPSOL::JMatrix<NDIM, double> > d_B_matrix =
        new JPSOL::JMatrix<NDIM, double>(patch_level, modal_B_mat);
    tbox::Array<tbox::Pointer<JPSOL::JVector<NDIM, double> > > d_port_sol(modal_sol_id.getSize());
    for (int i = 0; i < modal_sol_id.getSize(); i++) {
      d_port_sol[i] = new JPSOL::JVector<NDIM, double>(patch_level, modal_sol_id[i]);
    }
    tbox::Array<double> eigen_value(d_patch_strategy->eigen_num);
    tbox::Pointer<JPSOL::JPEigenSolver<NDIM, double> > PortModalSolver =
        new JPSOL::JPEigenSolver<NDIM, double>(d_solver_db->getDatabase("SolverModal"));
    PortModalSolver->set_A_Matrix(d_A_matrix);
    PortModalSolver->set_B_Matrix(d_B_matrix);
    tbox::pout << "Solving numerical waveport " << pp + 1 << endl;
    PortModalSolver->solve(d_port_sol, eigen_value);
    for (int e_n = 0; e_n < d_patch_strategy->eigen_num; e_n++) {
      d_patch_strategy->eigen_value_pairs[pp][e_n] = eigen_value[e_n];
    }
  }
  d_com_intc["COMM_WP_POWER"]->computing(patch_level, current_time, actual_dt);
  // 端口数量*模式数量*2
  const int length_sol = d_patch_strategy->eigen_num * port_num * 2;
  double modalpower[length_sol];
  tbox::pout << "Compute numerical waveport power" << endl;
  d_red_intc["WP_POWER"]->reduction(&modalpower[0], length_sol, patch_level, current_time,
                                    actual_dt);
  for (int p_n = 0; p_n < d_patch_strategy->numerical_waveport_num; p_n++) {
    tbox::pout<< "Port "<<p_n+1<<": "<<endl;
    for (int e_n = 0; e_n < d_patch_strategy->eigen_num; e_n++) {
      tbox::pout<<"Eigen number "<<e_n+1<<": ";
      tbox::pout<<modalpower[2 * (p_n * d_patch_strategy->eigen_num + e_n)]<<"+j"<<
      modalpower[2 * (p_n * d_patch_strategy->eigen_num + e_n)+1]<<endl;
    }
  }

  d_mem_intc["ALLOC_MODAL_MATVEC"]->deallocatePatchData(level);
  // 开辟矩阵、向量、有限元解数据片内存
  d_mem_intc["ALLOC_MATVEC"]->allocatePatchData(patch_level, current_time + predict_dt);

  // 逐个eigen值去求解，如果是线缆的话，两个特征根作为共模和差模
  for (int m = 0; m < d_patch_strategy->eigen_num; m++) {
    d_patch_strategy->active_mode = m;
    // 组装矩阵 组装右端项 加载约束
    d_num_intc["COMPUTE_LHSRHS"]->computing(patch_level, current_time, actual_dt);
    tbox::pout << "assemble matrix ok " << std::endl;
    /*
    // Solve
    int mat_id  = d_patch_strategy->d_matrix_id;
    int rhs_id  = d_patch_strategy->d_rhs_id;
    int sol_id  = d_patch_strategy->d_solution_id;

    if (first_step) {
      d_rhs    = new solv::JAUMINVector<NDIM, dcomplex>(patch_level,rhs_id);
      d_sol    = new solv::JAUMINVector<NDIM, dcomplex>(patch_level,sol_id);
      d_matrix = new solv::JAUMINMatrix<NDIM, dcomplex>(patch_level,mat_id);
    }
    d_solver->setMatrix(d_matrix);
    d_solver->setRHS(d_rhs);
    d_solver->solve(d_sol);
    */
#if 1
    d_num_intc["INIT_AMSVAR"]->computing(patch_level, current_time, actual_dt);

    int pre_mat_id = d_patch_strategy->d_premat_id;
    int gmat_id = d_patch_strategy->d_gmat_id;
    int xcoord_id = d_patch_strategy->d_xcoord_id;
    int ycoord_id = d_patch_strategy->d_ycoord_id;
    int zcoord_id = d_patch_strategy->d_zcoord_id;
    d_gmat = new solv::JAUMINMatrix<NDIM, double>(patch_level, gmat_id);
    d_xcoord = new solv::JAUMINVector<NDIM, double>(patch_level, xcoord_id);
    d_ycoord = new solv::JAUMINVector<NDIM, double>(patch_level, ycoord_id);
    d_zcoord = new solv::JAUMINVector<NDIM, double>(patch_level, zcoord_id);
    d_pre_matrix = new solv::JAUMINMatrix<NDIM, double>(patch_level, pre_mat_id);

    // 加载AMS方法的辅助变量
    d_solver->setGradMatrix(d_gmat);
    d_solver->setXCoord(d_xcoord);
    d_solver->setYCoord(d_ycoord);
    d_solver->setZCoord(d_zcoord);

    // 设置预条件矩阵
    d_solver->setPreMatrix(d_pre_matrix);

#endif

    // new  Solve
    int mat_re_id = d_patch_strategy->d_mat_re_id;
    int mat_im_id = d_patch_strategy->d_mat_im_id;
    int rhs_re_id = d_patch_strategy->d_rhs_re_id;
    int rhs_im_id = d_patch_strategy->d_rhs_im_id;
    int sol_re_id = d_patch_strategy->d_sol_re_id;
    int sol_im_id = d_patch_strategy->d_sol_im_id;

    if (first_step) {
      d_rhs_re = new solv::JAUMINVector<NDIM, double>(patch_level, rhs_re_id);
      d_rhs_im = new solv::JAUMINVector<NDIM, double>(patch_level, rhs_im_id);
      d_sol_re = new solv::JAUMINVector<NDIM, double>(patch_level, sol_re_id);
      d_sol_im = new solv::JAUMINVector<NDIM, double>(patch_level, sol_im_id);
      d_mat_re = new solv::JAUMINMatrix<NDIM, double>(patch_level, mat_re_id);
      d_mat_im = new solv::JAUMINMatrix<NDIM, double>(patch_level, mat_im_id);
    }

    d_solver->setPreMatrix(d_pre_matrix);
    d_solver->setMatrix(d_mat_re, d_mat_im);
    d_solver->setRHS(d_rhs_re, d_rhs_im);
    d_solver->solve(d_sol_re, d_sol_im);

    /*
    d_solver->setMatrix(d_pre_matrix);
    d_solver->setRHS(d_rhs_re);
    d_solver->setPreAlgorithm(d_pre_algorithm);
    d_solver->solve(d_sol_re);
    tbox::pout<< "the second amg solve" << endl;
    d_solver->setRHS(d_rhs_im);
    d_solver->solve(d_sol_im);
    */

    //    Move solution data to field data.
    d_num_intc["ACCEPT_SOLU"]->computing(patch_level, current_time, actual_dt);

    // Convert EdgeData to CellData for plotting,
    // since TeraVap can't plot edge data.
    d_num_intc["POSTPROCESS"]->computing(patch_level, current_time, actual_dt);
  }

#if 1
  // 计算数值解与解析解误差.
  double error = 0.;
  d_red_intc["SUM_ERR"]->reduction(&error, 1, patch_level, current_time, actual_dt);
  if (tbox::MPI::getRank() == 0) {
    printf("\n** Hcurl error squre: %0.6le\n", error);
    printf("\n** Hcurl error: %0.6le\n", sqrt(error));
  }
#endif
#if 0
  if(tbox::MPI::getRank() == 0) {
    std::cout << "\nreference value:\n";
    std::cout << "----------------------------------------------------------------\n";
    std::cout << "file_name = cube1o_r0.k          PHG: Hcurl error = 2.466566e+00\n"
              << "file_name = cube1o_r1.k          PHG: Hcurl error = 2.367742e+00\n"
              << "file_name = cube1o_r2.k          PHG: Hcurl error = 1.478600e+00\n"
              << "file_name = cube1o_r3.k          PHG: Hcurl error = 1.257240e+00\n"
              << "file_name = cube1o_r4.k          PHG: Hcurl error = 1.000217e+00\n"
              << "file_name = cube1o_r5.k          PHG: Hcurl error = 6.097758e-01\n"
              << "file_name = cube1o_r6.k          PHG: Hcurl error = 7.138828e-01\n"
              << std::endl;
  }
#endif
  // 释放临时数据片内存.
  d_mem_intc["ALLOC_MATVEC"]->deallocatePatchData(level);

  actual_dt = 1.0;
  return (0);
}

/*************************************************************************
 * 接收数值解.
 ************************************************************************/
void MaxwellLevelIntegrator::acceptTimeDependentSolution(
    const tbox::Pointer<hier::BasePatchLevel<NDIM> > level, const double new_time,
    const bool last_step) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!level.isNull());
#endif
}
