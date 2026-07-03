//
// 文件名:  Maxwell.C
// 软件包:  JAUMIN applications
// 版权  :  (c) 2004-2015 北京应用物理与计算数学研究所
//          (c) 2013-2015 中物院高性能数值模拟软件中心
// 版本号:  $Revision$
// 修改  :  $Data$
// 描述  :  Maxwell方程的网格片策略类的实现.
//

#include <algorithm>
#include <cstdlib>
#include <stdlib.h>

#include "JAUMIN_Macros.h"
#include "utils/JAUMIN_App.h"

#include "DOFInfo.h"
#include "Maxwell.h"
#include "fem/Nedelec.h"
#include "fem/TetGeom.h"
#include "fem/TetQuad.h"
#include "fem/TriGeom.h"
#include "fem/triNedelec.h"
typedef dcomplex NTYPE;

/*************************************************************************
 * 构造函数.
 *************************************************************************/
Maxwell::Maxwell(const string &object_name, tbox::Pointer<tbox::Database> input_db) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!object_name.empty());
  TBOX_ASSERT(!input_db.isNull());
#endif

  d_object_name = object_name;
  d_user_param.getFromInput(input_db);
  tbox::pout << "Reading input database......" << endl;
  getFromInput(input_db);

  // 注册变量和数据片.
  registerModelVariables();
}

/*************************************************************************
 * 构造函数.
 ************************************************************************/
Maxwell::~Maxwell() {}

/*************************************************************************
 * 注册变量和数据片.
 ************************************************************************/
void Maxwell::registerModelVariables() {
  // 棱边的电场量
  tbox::Pointer<pdat::EdgeVariable<NDIM, NTYPE> > E =
      new pdat::EdgeVariable<NDIM, NTYPE>("E", 1, 1);
  // 单元中心的电场量，用于输出和转换
  tbox::Pointer<pdat::CellVariable<NDIM, NTYPE> > Ec =
      new pdat::CellVariable<NDIM, NTYPE>("Ec", NDIM, 1);
  tbox::Pointer<pdat::CellVariable<NDIM, double> > Ec_norm =
      new pdat::CellVariable<NDIM, double>("Ec_norm", NDIM, 1);
  tbox::Pointer<pdat::CellVariable<NDIM, double> > Ec_arg =
      new pdat::CellVariable<NDIM, double>("Ec_arg", NDIM, 1);
  // 解析解在棱边上的值，用于比较
  tbox::Pointer<pdat::EdgeVariable<NDIM, NTYPE> > Ea =
      new pdat::EdgeVariable<NDIM, NTYPE>("Ea", 1, 1);
  // 解析解在单元中心的值
  tbox::Pointer<pdat::CellVariable<NDIM, double> > Eac =
      new pdat::CellVariable<NDIM, double>("Eac", NDIM, 1);
  // 右端项
  tbox::Pointer<pdat::EdgeVariable<NDIM, double> > J =
      new pdat::EdgeVariable<NDIM, double>("J", 1, 1);
  // 每一个单元的jacobi矩阵
  tbox::Pointer<pdat::CellVariable<NDIM, double> > jacobian =
      new pdat::CellVariable<NDIM, double>("jacobian", (NDIM + 1) * (NDIM + 1), 1);
  // 每个单元的体积
  tbox::Pointer<pdat::CellVariable<NDIM, double> > volume =
      new pdat::CellVariable<NDIM, double>("volume", 1, 1);
  // 边的方向。取值为1或者-1。因为边没有全局编号，所以不能像串行那样比较边的两个顶点编号来决定边的正向。
  // 第i条边取值为1表示edge_node_idx[edge_node_ext[i]]-->edge_node_idx[edge_node_ext[i]+1]为正向。
  // 其计算方式参见computeEdgeOrderOnPatch
  tbox::Pointer<pdat::EdgeVariable<NDIM, int> > edge_order =
      new pdat::EdgeVariable<NDIM, int>("edge_order", 1, 1);
  // 取值为1或0，1表示该边位于边界上。其计算方式参见computeEdgeFlagOnPatch
  tbox::Pointer<pdat::EdgeVariable<NDIM, int> > edge_flag =
      new pdat::EdgeVariable<NDIM, int>("edge_flag", 1, 1);
  // 数值波端口的edge flag，和端口号一致
  tbox::Pointer<pdat::EdgeVariable<NDIM, int> > port_edge_flag =
      new pdat::EdgeVariable<NDIM, int>("port_edge_flag", 1, 1);

  hier::VariableDatabase<NDIM> *variable_db = hier::VariableDatabase<NDIM>::getDatabase();
  tbox::Pointer<hier::VariableContext> current = variable_db->getContext("current");

  d_E_id = variable_db->registerVariableAndContext(E, current, 0);
  d_Ec_id = variable_db->registerVariableAndContext(Ec, current, 0);
  d_Ec_norm_id = variable_db->registerVariableAndContext(Ec_norm, current, 0);
  d_Ec_arg_id = variable_db->registerVariableAndContext(Ec_arg, current, 0);
  d_Ea_id = variable_db->registerVariableAndContext(Ea, current, 1);
  d_Eac_id = variable_db->registerVariableAndContext(Eac, current, 0);
  d_J_id = variable_db->registerVariableAndContext(J, current, 1);
  d_edge_order_id = variable_db->registerVariableAndContext(edge_order, current, 1);
  d_edge_flag_id = variable_db->registerVariableAndContext(edge_flag, current, 1);
  d_port_edge_flag_id = variable_db->registerVariableAndContext(port_edge_flag, current, 1);
  d_jacobian_id = variable_db->registerVariableAndContext(jacobian, current, 1);
  d_volume_id = variable_db->registerVariableAndContext(volume, current, 1);

  // 材料参数数据库
  tbox::Pointer<pdat::CellVariable<NDIM, double> > epsilonr =
      new pdat::CellVariable<NDIM, double>("epsilonr", 1, 1);
  tbox::Pointer<pdat::CellVariable<NDIM, double> > mur =
      new pdat::CellVariable<NDIM, double>("mur", 1, 1);
  tbox::Pointer<pdat::CellVariable<NDIM, double> > losstan =
      new pdat::CellVariable<NDIM, double>("losstan", 1, 1);
  tbox::Pointer<pdat::CellVariable<NDIM, double> > sigma =
      new pdat::CellVariable<NDIM, double>("sigma", 1, 1);
  tbox::Pointer<pdat::CellVariable<NDIM, int> > mat_id =
      new pdat::CellVariable<NDIM, int>("mat_id", 1, 1);
  d_epsilonr_id = variable_db->registerVariableAndContext(epsilonr, current, 1);
  d_mur_id = variable_db->registerVariableAndContext(mur, current, 1);
  d_losstan_id = variable_db->registerVariableAndContext(losstan, current, 1);
  d_sigma_id = variable_db->registerVariableAndContext(sigma, current, 1);
  d_mat_id = variable_db->registerVariableAndContext(mat_id, current, 1);

  // 自由度信息(四个逻辑型表示点，边，面，体上自由度是否非零)
  d_dof_info = new solv::DOFInfo<NDIM>(false, true, false, false);
  d_dof_info_node = new solv::DOFInfo<NDIM>(true, false, false, false);

  /*
  tbox::Pointer<pdat::VectorVariable<NDIM,NTYPE> > solution = new
  pdat::VectorVariable<NDIM,NTYPE>("solution", d_dof_info); // 向量变量， 解向量
  tbox::Pointer<pdat::VectorVariable<NDIM,NTYPE> > rhs = new
  pdat::VectorVariable<NDIM,NTYPE>("rhs", d_dof_info); // 向量变量， 右端项
  tbox::Pointer<pdat::CSRMatrixVariable<NDIM,NTYPE> > matrix = new
  pdat::CSRMatrixVariable<NDIM,NTYPE>("matrix", d_dof_info); // 矩阵变量， 矩阵

  d_solution_id = variable_db->registerVariableAndContext(solution, current, 1);
  d_rhs_id = variable_db->registerVariableAndContext(rhs, current, 1);
  d_matrix_id = variable_db->registerVariableAndContext(matrix, current, 1);
  */
  // 生成解，右端项，和系数矩阵的实部和虚部变量
  tbox::Pointer<pdat::VectorVariable<NDIM, double> > sol_re =
      new pdat::VectorVariable<NDIM, double>("solution_real",
                                             d_dof_info); // 向量变量， 解向量
  tbox::Pointer<pdat::VectorVariable<NDIM, double> > sol_im =
      new pdat::VectorVariable<NDIM, double>("solution_image",
                                             d_dof_info); // 向量变量， 解向量
  tbox::Pointer<pdat::VectorVariable<NDIM, double> > rhs_re =
      new pdat::VectorVariable<NDIM, double>("rhs_real",
                                             d_dof_info); // 向量变量， 右端项
  tbox::Pointer<pdat::VectorVariable<NDIM, double> > rhs_im =
      new pdat::VectorVariable<NDIM, double>("rhs_image",
                                             d_dof_info); // 向量变量， 右端项
  tbox::Pointer<pdat::CSRMatrixVariable<NDIM, double> > mat_re =
      new pdat::CSRMatrixVariable<NDIM, double>("mat_real",
                                                d_dof_info); // 矩阵变量， 矩阵
  tbox::Pointer<pdat::CSRMatrixVariable<NDIM, double> > mat_im =
      new pdat::CSRMatrixVariable<NDIM, double>("mat_image",
                                                d_dof_info); // 矩阵变量， 矩阵

  tbox::Pointer<pdat::CSRMatrixVariable<NDIM, double> > premat =
      new pdat::CSRMatrixVariable<NDIM, double>("pre_mat",
                                                d_dof_info); // 矩阵变量， 矩阵

  /// AMS 方法中的辅助边量
  tbox::Pointer<pdat::CSRMatrixVariable<NDIM, double> > gmat =
      new pdat::CSRMatrixVariable<NDIM, double>("gmat",
                                                d_dof_info); // 矩阵变量， 矩阵
  tbox::Pointer<pdat::VectorVariable<NDIM, double> > xcoord =
      new pdat::VectorVariable<NDIM, double>("xcoord", d_dof_info_node); // 向量变量， 解向量
  tbox::Pointer<pdat::VectorVariable<NDIM, double> > ycoord =
      new pdat::VectorVariable<NDIM, double>("ycoord", d_dof_info_node); // 向量变量， 解向量
  tbox::Pointer<pdat::VectorVariable<NDIM, double> > zcoord =
      new pdat::VectorVariable<NDIM, double>("zcoord", d_dof_info_node); // 向量变量， 解向量

  d_sol_re_id = variable_db->registerVariableAndContext(sol_re, current, 1);
  d_sol_im_id = variable_db->registerVariableAndContext(sol_im, current, 1);
  d_rhs_re_id = variable_db->registerVariableAndContext(rhs_re, current, 1);
  d_rhs_im_id = variable_db->registerVariableAndContext(rhs_im, current, 1);
  d_mat_re_id = variable_db->registerVariableAndContext(mat_re, current, 1);
  d_mat_im_id = variable_db->registerVariableAndContext(mat_im, current, 1);

  d_premat_id = variable_db->registerVariableAndContext(premat, current, 1);

  d_gmat_id = variable_db->registerVariableAndContext(gmat, current, 1);
  d_xcoord_id = variable_db->registerVariableAndContext(xcoord, current, 1);
  d_ycoord_id = variable_db->registerVariableAndContext(ycoord, current, 1);
  d_zcoord_id = variable_db->registerVariableAndContext(zcoord, current, 1);

  // 模式场求解的变量注册
  d_matA_modal_array.resizeArray(numerical_waveport_num);
  d_matB_modal_array.resizeArray(numerical_waveport_num);
  d_sol_modal_array.resizeArray(numerical_waveport_num);
  for (int pp = 0; pp < numerical_waveport_num; pp++) {
    stringstream matAname;
    matAname << "modal_matA_" << pp + 1;
    stringstream matBname;
    matAname << "modal_matB_" << pp + 1;
    tbox::Pointer<pdat::CSRMatrixVariable<NDIM, double> > Modal_A_matrix =
        new pdat::CSRMatrixVariable<NDIM, double>(matAname.str(),
                                                  d_dof_info_modal_array[pp]); //
    tbox::Pointer<pdat::CSRMatrixVariable<NDIM, double> > Modal_B_matrix =
        new pdat::CSRMatrixVariable<NDIM, double>(matBname.str(),
                                                  d_dof_info_modal_array[pp]); //
    REGISTER_VARIABLE(d_matA_modal_array[pp], Modal_A_matrix, current, 1);
    REGISTER_VARIABLE(d_matB_modal_array[pp], Modal_B_matrix, current, 1);
    d_sol_modal_array[pp].resizeArray(eigen_num);
    for (int nn = 0; nn < eigen_num; nn++) {
      stringstream solname;
      solname << "modal_sol_" << pp + 1 << "_" << nn;
      tbox::Pointer<pdat::VectorVariable<NDIM, double> > Modal_sol =
          new pdat::VectorVariable<NDIM, double>(solname.str(), d_dof_info_modal_array[pp]);
      d_sol_modal_array[pp][nn] = variable_db->registerVariableAndContext(Modal_sol, current, 1);
    }
  }
}

/*************************************************************************
 *  从输入数据库读入数据.
 ************************************************************************/
void Maxwell::getFromInput(tbox::Pointer<tbox::Database> db) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!db.isNull());
#endif
  entity_face_num = db->getInteger("entity_face_num");
  simulation_setup = db->getInteger("simulation_setup");
  frequency_setup = db->getDoubleArray("frequency_setup");
  tbox::pout << "Reading material properties......" << endl;
  InitMaterialComponent(db);
  tbox::pout << "Reading boundary properties......" << endl;
  InitBoundaryInfo(db);
  tbox::pout << "Reading port properties......" << endl;
  InitPortInfo(db);
}

/*************************************************************************
 * 注册绘图量.
 *************************************************************************/
void Maxwell::registerPlotData(tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_writer) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!(javis_writer.isNull()));
#endif
  javis_writer->registerPlotQuantity("E_h_norm", "VECTOR", d_Ec_norm_id);
  javis_writer->registerPlotQuantity("E_h_arg", "VECTOR", d_Ec_arg_id);
  // javis_writer->registerPlotQuantity("E_h", "VECTOR", d_Ec_id);
  // javis_writer->registerPlotQuantity("E_a", "VECTOR", d_Eac_id);
}

/*************************************************************************
 *  初始化指定的积分构件.
 ************************************************************************/
void Maxwell::initializeComponent(algs::IntegratorComponent<NDIM> *component) const {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(component);
#endif
  const string &component_name = component->getName();
  if (component_name == "ALLOC_MATVEC") {
    //    component->registerPatchData(d_matrix_id);
    //    component->registerPatchData(d_solution_id);
    //    component->registerPatchData(d_rhs_id);

    component->registerPatchData(d_mat_re_id);
    component->registerPatchData(d_mat_im_id);
    component->registerPatchData(d_rhs_re_id);
    component->registerPatchData(d_rhs_im_id);
    component->registerPatchData(d_sol_re_id);
    component->registerPatchData(d_sol_im_id);
    component->registerPatchData(d_premat_id);

    // AMS方法中变量
    component->registerPatchData(d_gmat_id);
    component->registerPatchData(d_xcoord_id);
    component->registerPatchData(d_ycoord_id);
    component->registerPatchData(d_zcoord_id);
  }
  if (component_name == "ALLOC_MODAL_MATVEC") {
    for (int num_port = 0; num_port < numerical_waveport_num; num_port++) {
      component->registerPatchData(d_matA_modal_array[num_port]);
      component->registerPatchData(d_matB_modal_array[num_port]);
      for (int num_eigen = 0; num_eigen < eigen_num; num_eigen++) {
        component->registerPatchData(d_sol_modal_array[num_port][num_eigen]);
      }
    }

  } else if (component_name == "INIT_FIELD") {
    component->registerInitPatchData(d_E_id);
    component->registerInitPatchData(d_Ea_id);
    component->registerInitPatchData(d_Ec_id);
    component->registerInitPatchData(d_Ec_norm_id);
    component->registerInitPatchData(d_Ec_arg_id);
    component->registerInitPatchData(d_Eac_id);
    component->registerInitPatchData(d_J_id);
    component->registerInitPatchData(d_edge_order_id);
    component->registerInitPatchData(d_edge_flag_id);
    component->registerInitPatchData(d_port_edge_flag_id);
    component->registerInitPatchData(d_jacobian_id);
    component->registerInitPatchData(d_volume_id);

    component->registerInitPatchData(d_epsilonr_id);
    component->registerInitPatchData(d_mur_id);
    component->registerInitPatchData(d_losstan_id);
    component->registerInitPatchData(d_sigma_id);
    component->registerInitPatchData(d_mat_id);

  } else if (component_name == "COMM_EDGE_FLAG") {
    // 注册待填充影像区的数据片
    component->registerCommunicationPatchData(d_edge_flag_id, d_edge_flag_id);
    component->registerCommunicationPatchData(d_port_edge_flag_id, d_port_edge_flag_id);
  } else if (component_name == "COMM_WP_POWER") {
    // 注册待填充影像区的数据片
    for (int num_port = 0; num_port < numerical_waveport_num; num_port++)
      for (int num_mod = 0; num_mod < eigen_num; num_mod++)
        component->registerCommunicationPatchData(d_sol_modal_array[num_port][num_mod],
                                                  d_sol_modal_array[num_port][num_mod]);

  } else if (component_name == "INIT_DOFINFO") {
    // 将自由度信息中的若干数据片注册到初始化构件。
    d_dof_info->registerToInitComponent(component);
    d_dof_info_node->registerToInitComponent(component);
    for (int i = 0; i < d_dof_info_modal_array.getSize(); ++i)
      d_dof_info_modal_array[i]->registerToInitComponent(component);
  }
}

/*************************************************************************
 *  初始化数据片（支持初值构件）.
 ************************************************************************/
void Maxwell::initializePatchData(hier::Patch<NDIM> &patch, const double time,
                                  const bool initial_time, const string &component_name) {
  // 本示例程序中，
  // 初值构件只用来开辟内存，不用来初始化数据。
  // 用数值构件来初始化数据。
}

/*************************************************************************
 * 完成单个网格片上的数值计算（支持数值构件）.
 ************************************************************************/
void Maxwell::computeOnPatch(hier::Patch<NDIM> &patch, const double time, const double dt,
                             const bool initial_time, const string &component_name) {

  if (component_name == "INIT_ANAL")
    return initAnalyticOnPatch(patch);
  if (component_name == "INIT_J")
    return initJOnPatch(patch);
  if (component_name == "INIT_GEOM")
    return initGeometryOnPatch(patch);
  if (component_name == "DOF_MAP")
    return setupDofInfoOnPatch(patch);
  if (component_name == "MODAL_DOF_MAP")
    return setupWPModalDofInfoOnPatch(patch);
  if (component_name == "INIT_AMSVAR")
    return buildAMSVarOnPatch(patch);
  if (component_name == "COMPUTE_LHSRHS")
    return buildLHSAndRHSOnPatch(patch, time, dt);
  if (component_name == "ACCEPT_SOLU")
    return acceptSolutionOnPatch(patch);
  if (component_name == "POSTPROCESS")
    return postprocessOnPatch(patch);
  if (component_name == "EDGE_ORDER")
    return computeEdgeOrderOnPatch(patch);
  if (component_name == "EDGE_FLAG")
    return computeEdgeFlagOnPatch(patch);
  if (component_name == "MODAL_MAT")
    return setupModalMatrixOnPatch(patch, time, dt);
  if (component_name == "COMM_EDGE_FLAG")
    return;
  if (component_name == "COMM_WP_POWER")
    return;
  if (component_name == "COMM_CELL_FLAG")
    return;
  TBOX_ERROR("Component \"" << component_name << "\" not matched.\n");
}

/*************************************************************************
 * 完成单个网格片上的规约计算（支持规约构件）.
 ************************************************************************/
void Maxwell::reduceOnPatch(double *vector, int len, hier::Patch<NDIM> &patch, const double time,
                            const double dt, const string &component_name) {
  if (component_name == "SUM_ERR")
    return computeHcurlErrOnPatch(vector, len, patch);
  if (component_name == "WP_POWER")
    return calculateModalPower(vector, len, patch, time, dt);
  TBOX_ERROR("Component \"" << component_name << "\" not matched.\n");
}

/*************************************************************************
 *  计算单元体积、Jacobian矩阵
 ************************************************************************/
void Maxwell::initGeometryOnPatch(hier::Patch<NDIM> &patch) {
  int num_cell_ghost = patch.getNumberOfCells(1);
  tbox::Pointer<pdat::CellData<NDIM, double> > jacobian = patch.getPatchData(d_jacobian_id);
  tbox::Pointer<pdat::CellData<NDIM, double> > volume = patch.getPatchData(d_volume_id);

  DECLARE_ADJACENCY(patch, cell, node, Cell, Node);
  appu::TetGeom tetrahedron(patch);

  for (int cell = 0; cell < num_cell_ghost; cell++) {
    if (cell_node_ext[cell + 1] - cell_node_ext[cell] > 4) {
    }
    (*volume)(0, cell) = tetrahedron.volume(cell);
    tetrahedron.jacobian(cell, &((*jacobian)(0, cell)));
  }
}

/*************************************************************************
 *  计算源项（直接由解析函数获得）
 ************************************************************************/
void Maxwell::initJOnPatch(hier::Patch<NDIM> &patch) {
  int num_edge_ghost = patch.getNumberOfEdges(1);
  tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
      patch.getPatchGeometry()->getNodeCoordinates();
  tbox::Pointer<pdat::EdgeData<NDIM, int> > edge_order = patch.getPatchData(d_edge_order_id);
  tbox::Pointer<pdat::EdgeData<NDIM, double> > J = patch.getPatchData(d_J_id);
  tbox::Array<int> edge_node_ext, edge_node_idx;
  patch.getPatchTopology()->getEdgeAdjacencyNodes(edge_node_ext, edge_node_idx);

  for (int edge = 0; edge < num_edge_ghost; edge++) {
    int v0 = edge_node_idx[edge_node_ext[edge] + 0];
    int v1 = edge_node_idx[edge_node_ext[edge] + 1];
    int inorder = (*edge_order)(0, edge);
    if (!inorder)
      std::swap(v0, v1);
    double x0, y0, z0, x1, y1, z1;
    x0 = (*node_coord)(0, v0);
    y0 = (*node_coord)(1, v0);
    z0 = (*node_coord)(2, v0);
    x1 = (*node_coord)(0, v1);
    y1 = (*node_coord)(1, v1);
    z1 = (*node_coord)(2, v1);
    double funcvalue[NDIM] = {0.0};
    // Analytic::func_f((x0 + x1) * 0.5, (y0 + y1) * 0.5, (z0 + z1) * 0.5,
    // funcvalue);
    x1 -= x0;
    y1 -= y0;
    z1 -= z0;
    double dof = x1 * funcvalue[0] + y1 * funcvalue[1] + z1 * funcvalue[2];
    (*J)(0, edge) = dof;
  }
}

/*************************************************************************
 *  初始化解析解
 ************************************************************************/
void Maxwell::initAnalyticOnPatch(hier::Patch<NDIM> &patch) {
  int num_edge_ghost = patch.getNumberOfEdges(1);
  int num_cell = patch.getNumberOfCells(0);
  tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
      patch.getPatchGeometry()->getNodeCoordinates();
  tbox::Pointer<pdat::CellData<NDIM, double> > cell_coord =
      patch.getPatchGeometry()->getCellCoordinates();
  tbox::Pointer<pdat::EdgeData<NDIM, int> > edge_order = patch.getPatchData(d_edge_order_id);
  tbox::Pointer<pdat::EdgeData<NDIM, NTYPE> > Ea = patch.getPatchData(d_Ea_id);
  tbox::Pointer<pdat::CellData<NDIM, double> > Eac = patch.getPatchData(d_Eac_id);
  tbox::Array<int> edge_node_ext, edge_node_idx;
  patch.getPatchTopology()->getEdgeAdjacencyNodes(edge_node_ext, edge_node_idx);
  tbox::Array<int> cell_edge_ext, cell_edge_idx;
  patch.getPatchTopology()->getCellAdjacencyEdges(cell_edge_ext, cell_edge_idx);

  for (int edge = 0; edge < num_edge_ghost; edge++) {
    int v0 = edge_node_idx[edge_node_ext[edge] + 0];
    int v1 = edge_node_idx[edge_node_ext[edge] + 1];
    int inorder = (*edge_order)(0, edge);
    if (!inorder)
      std::swap(v0, v1);
    double x0, y0, z0, x1, y1, z1;
    x0 = (*node_coord)(0, v0);
    y0 = (*node_coord)(1, v0);
    z0 = (*node_coord)(2, v0);
    x1 = (*node_coord)(0, v1);
    y1 = (*node_coord)(1, v1);
    z1 = (*node_coord)(2, v1);
    double funcvalue[NDIM] = {0.0};
    // Analytic::func_u((x0 + x1) * 0.5, (y0 + y1) * 0.5, (z0 + z1) * 0.5,
    // funcvalue);
    x1 -= x0;
    y1 -= y0;
    z1 -= z0;
    double dof = x1 * funcvalue[0] + y1 * funcvalue[1] + z1 * funcvalue[2];
    (*Ea)(0, edge) = dof;
  }

  for (int cell = 0; cell < num_cell; cell++) {
    //    double x, y, z;
    double funcvalue[NDIM] = {0.0};
    //    x = (*cell_coord)(0, cell);
    //    y = (*cell_coord)(1, cell);
    //    z = (*cell_coord)(2, cell);
    // Analytic::func_u(x, y, z, funcvalue);
    (*Eac)(0, cell) = funcvalue[0];
    (*Eac)(1, cell) = funcvalue[1];
    (*Eac)(2, cell) = funcvalue[2];
  }
}

/*************************************************************************
 *  建立网格片上的左端项.
 ************************************************************************/
void Maxwell::setupDofInfoOnPatch(hier::Patch<NDIM> &patch) {
  int num_edge_ghost = patch.getNumberOfEdges(1);
  int *dis_ptr = d_dof_info->getDOFDistribution(patch, hier::EntityUtilities::EDGE);
  for (int i = 0; i < num_edge_ghost; ++i)
    dis_ptr[i] = 1;
  d_dof_info->buildPatchDOFMapping(patch);

  int num_node_ghost = patch.getNumberOfNodes(1);
  int *dis_ptr_node = d_dof_info_node->getDOFDistribution(patch, hier::EntityUtilities::NODE);
  for (int i = 0; i < num_node_ghost; ++i)
    dis_ptr_node[i] = 1;
  d_dof_info_node->buildPatchDOFMapping(patch);
}

/*************************************************************************
 * 计算AMS的辅助边量：
 * G： gmat 节点级函数的梯度到棱基函数的插值映射关系
 * xcoord x方向的坐标向量
 * ycoord y方向的坐标向量
 * zcoord z方向的坐标向量
 ************************************************************************/
void Maxwell::buildAMSVarOnPatch(hier::Patch<NDIM> &patch) {

  tbox::Pointer<pdat::VectorData<NDIM, double> > xcoord = patch.getPatchData(d_xcoord_id);
  tbox::Pointer<pdat::VectorData<NDIM, double> > ycoord = patch.getPatchData(d_ycoord_id);
  tbox::Pointer<pdat::VectorData<NDIM, double> > zcoord = patch.getPatchData(d_zcoord_id);

  tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
      patch.getPatchGeometry()->getNodeCoordinates();

  int num_node_ghost = patch.getNumberOfNodes(1);
  /// 节点自由度信息中的映射信息。
  int *node_dof_map = d_dof_info_node->getDOFMapping(patch, hier::EntityUtilities::NODE);

  for (int i = 0; i < num_node_ghost; ++i) {

    int index = node_dof_map[i];
    (*xcoord)(index) = (*node_coord)(0, i);
    (*ycoord)(index) = (*node_coord)(1, i);
    (*zcoord)(index) = (*node_coord)(2, i);
  }

  tbox::Pointer<pdat::EdgeData<NDIM, int> > edge_order = patch.getPatchData(d_edge_order_id);
  tbox::Pointer<pdat::CSRMatrixData<NDIM, double> > gmat = patch.getPatchData(d_gmat_id);

  int num_edge_ghost = patch.getNumberOfEdges(1);
  /// 边自由度信息中的映射信息。
  int *edge_dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::EDGE);

  tbox::Array<int> edge_node_ext, edge_node_idx;
  patch.getPatchTopology()->getEdgeAdjacencyNodes(edge_node_ext, edge_node_idx);

  for (int edge = 0; edge < num_edge_ghost; edge++) {
    int v0 = edge_node_idx[edge_node_ext[edge] + 0];
    int v1 = edge_node_idx[edge_node_ext[edge] + 1];
    int inorder = (*edge_order)(0, edge);
    if (!inorder)
      std::swap(v0, v1);
    int index = edge_dof_map[edge];
    int nd1 = node_dof_map[v0];
    int nd2 = node_dof_map[v1];

    gmat->addMatrixValue(index, nd1, -1.0);
    gmat->addMatrixValue(index, nd2, 1.0);
  }
  gmat->assemble();
}

/*************************************************************************
 *  建立网格片上的左端项、右端项，并施加边界条件.
 ************************************************************************/
void Maxwell::buildLHSAndRHSOnPatch(hier::Patch<NDIM> &patch, const double time, const double dt) {
  int num_cell_ghost = patch.getNumberOfCells(1);
  int num_edge = patch.getNumberOfEdges(0);

  // 自由度映射信息
  int *dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::EDGE);
  //  // 向量数据片
  //  tbox::Pointer< pdat::VectorData<NDIM,NTYPE> > vec =
  //  patch.getPatchData(d_rhs_id);
  //  // 矩阵数据片
  //  tbox::Pointer< pdat::CSRMatrixData<NDIM,NTYPE> > mat =
  //  patch.getPatchData(d_matrix_id);

  // 向量数据片 实部和虚部
  tbox::Pointer<pdat::VectorData<NDIM, double> > vec_re = patch.getPatchData(d_rhs_re_id);
  tbox::Pointer<pdat::VectorData<NDIM, double> > vec_im = patch.getPatchData(d_rhs_im_id);
  // 矩阵数据片 实部和虚部
  tbox::Pointer<pdat::CSRMatrixData<NDIM, double> > mat_re = patch.getPatchData(d_mat_re_id);
  tbox::Pointer<pdat::CSRMatrixData<NDIM, double> > mat_im = patch.getPatchData(d_mat_im_id);

  tbox::Pointer<pdat::CSRMatrixData<NDIM, double> > pre_mat = patch.getPatchData(d_premat_id);

  // Shape Function
  appu::Nedelec shapefunc(patch, patch.getPatchData(d_edge_order_id),
                          patch.getPatchData(d_jacobian_id));
  appu::triNedelec trishapefunc(patch, patch.getPatchData(d_edge_order_id),
                                patch.getPatchData(d_jacobian_id));

  // Quadrature Rule
  appu::TetQuad quad(patch, patch.getPatchData(d_volume_id), patch.getPatchData(d_jacobian_id));

  tbox::Pointer<pdat::EdgeData<NDIM, double> > J = patch.getPatchData(d_J_id);
  tbox::Pointer<pdat::EdgeData<NDIM, NTYPE> > Ea = patch.getPatchData(d_Ea_id);
  tbox::Pointer<pdat::EdgeData<NDIM, int> > edge_flag = patch.getPatchData(d_edge_flag_id);

  DECLARE_ADJACENCY(patch, cell, face, Cell, Face);
  DECLARE_ADJACENCY(patch, cell, edge, Cell, Edge);
  DECLARE_ADJACENCY(patch, face, edge, Face, Edge);
  DECLARE_ADJACENCY(patch, face, node, Face, Node);
  // 参数定义
  int step = time / dt + 0.5;
  double freq = frequency_setup[step];
  double omega = 2 * M_PI * freq;
  double mu_0 = 4 * M_PI * 1e-7;
  double epsr_max = 10.;
  double epsilon_0 = 8.8542 * 1e-12;
  double k_0_square = (omega) * (omega)*mu_0 * epsilon_0;
  double zita = k_0_square * epsr_max;
  // 这里先把系数提出来
  double input_power = NUM_WP_list[0]->power;
  NTYPE coefofunit = 1.0 * sqrt(2. * input_power) / (sqrt(port_power[active_mode]));

  tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
      patch.getPatchGeometry()->getNodeCoordinates();

#if 0
  {
  // 计算任意一点point在35435号单元上的重心坐标。
  tbox::Pointer< pdat::CellData<NDIM,double> > jacobian_all = patch.getPatchData(d_jacobian_id);

  double point [] = {0.22961, 0.456723, 1.17694, 1.0};
  for (int cell = 0; cell < num_cell_ghost; cell++) {
      const int NBASIS = NDIM+1;
      if (cell == 35435) {
        double lambda[NBASIS] = {0,0,0,0};
        double(*jacobian)[NBASIS] = (double(*)[NBASIS]) (&((*jacobian_all)(0, cell)));
        cout << "lambda is " << endl;
        for (int i=0; i<NBASIS; i++) {
            for (int j=0; j<NBASIS; j++) {
                lambda[i] += jacobian[i][j] * point[j];
            }
            cout << lambda[i] << endl;
        }
      }
  }
  }
#endif

  // CELLPART: cell_part_id[cell] 用于保存cell单元的part id。
  tbox::Array<int> cell_part_id(num_cell_ghost);
  // 初始化为-1, 也就是不属于任何part
  for (int i = 0; i < num_cell_ghost; i++) {
    cell_part_id[i] = -1;
  }

  // 经过以下两个循环，每一个单元都应该属于一个part。否则程序有问题。
  int setid;
  int gw_width = 1;
  setid = 1;
  if (patch.hasEntitySet(setid, hier::EntityUtilities::CELL, gw_width)) {
    tbox::Array<int> indices =
        patch.getEntityIndicesInSet(setid, hier::EntityUtilities::CELL, gw_width);
    for (int i = 0; i < indices.size(); i++) {
      cell_part_id[indices[i]] = setid;
    }
  }
  setid = 2;
  if (patch.hasEntitySet(setid, hier::EntityUtilities::CELL, gw_width)) {
    tbox::Array<int> indices =
        patch.getEntityIndicesInSet(setid, hier::EntityUtilities::CELL, gw_width);
    for (int i = 0; i < indices.size(); i++) {
      cell_part_id[indices[i]] = setid;
    }
  }

  // 遍历单元(包括影像单元)
  for (int cell = 0; cell < num_cell_ghost; cell++) {
    const int nedge = (NDIM == 2 ? 3 : 6);
    double curl_x_curl[nedge][nedge];
    double phi_x_phi[nedge][nedge];
    double face_phi_x_phi[nedge][nedge];
    NTYPE fun_x_phi[nedge];
    NTYPE face_fun_x_phi[nedge];
    NTYPE A[nedge][nedge], b[nedge];
    double A_aux[nedge][nedge];
    int mapping[nedge];

    // compute \int \curl\phi_i \cdot \curl\phi_j
    quad.quadCurlBasDotCurlBas<appu::Nedelec>(cell, &shapefunc, 0, &(curl_x_curl[0][0]));

    // compute \int \phi_i \cdot \phi_j
    quad.quadBasDotBas<appu::Nedelec>(cell, &shapefunc, 2, &(phi_x_phi[0][0]));

    // right hand side: \int J * phi_i
    d_user_param.ftype = d_user_param.SOURCE;
    quad.quadFunctionDotBas<appu::Nedelec, NTYPE>(cell, d_user_param, &shapefunc, 2,
                                                  &(fun_x_phi[0]));

    // face quadrature
    for (int i = 0; i < nedge; i++) {
      face_fun_x_phi[i] = 0.0;
      for (int j = 0; j < nedge; j++) {
        face_phi_x_phi[i][j] = 0.0;
      }
    }
    /// 数值波端口的边界条件
    for (int num_port = 0; num_port < numerical_waveport_num; num_port++) {
      // 取出目前的wave port
      tbox::Pointer<NUM_WP> now_num_wp = NUM_WP_list[num_port];
      // 判断该端口是否为激励端口 (activated==true: 入射激励; false: 仅匹配BC)
      const bool is_input_port = now_num_wp->activated;

      // 取出模式dof映射
      int *modal_dof_map_edge =
          d_dof_info_modal_array[num_port]->getDOFMapping(patch, hier::EntityUtilities::EDGE);
      int *modal_dof_map_node =
          d_dof_info_modal_array[num_port]->getDOFMapping(patch, hier::EntityUtilities::NODE);
      double port_beta = sqrt(zita - zita / eigen_value_pairs[num_port][active_mode]);

      // 吸收边界条件：对所有模式求和 gama，确保端口面正确吸收所有传播模式
      dcomplex gama_sum(0, 0);
      for (int m = 0; m < eigen_num; m++) {
        double beta_m = sqrt(zita - zita / eigen_value_pairs[num_port][m]);
        gama_sum += dcomplex(0, 1) * beta_m;
      }
      for (int face_list_num = 0; face_list_num < now_num_wp->face_list.getSize();
           face_list_num++) {
        if (HAS_ENTITY_SET(patch, now_num_wp->face_list[face_list_num], FACE, 1)) {
          DECLARE_ENTITY_SET(patch, port_i_list, now_num_wp->face_list[face_list_num], FACE, 1);
          std::sort(port_i_list.getPointer(), port_i_list.getPointer() + port_i_list.size());
          int nface = 4;
          // 对本单元的所有的面进行循环
          for (int loc_ff = 0; loc_ff < nface; loc_ff++) {
            int glo_ff = cell_face_idx[cell_face_ext[cell] + loc_ff];
            bool is_on_port_face = std::binary_search(
                port_i_list.getPointer(), port_i_list.getPointer() + port_i_list.size(), glo_ff);
            if (!is_on_port_face)
              continue;
            // BC 项：所有端口（激励+匹配）都需要
            quad.faceQuadBasDotBas(cell, loc_ff, &shapefunc, 2, &(face_phi_x_phi[0][0]));

            // 激励项：仅输入端口需要
            if (is_input_port) {
              GET_PATCH_DATA(patch, modal_sol, d_sol_modal_array[num_port][active_mode], Vector, double);
              double *modal_vec_pointer = modal_sol->getPointer();
              NTYPE edge_coef[3] = {0, 0, 0}; // 面上边的插值系数
              NTYPE node_coef[3] = {0, 0, 0}; // 面上结点的插值系数
              for (int i = 0; i < 3; i++) {
                int edge = face_edge_idx[face_edge_ext[glo_ff] + i];
                int edge_map = modal_dof_map_edge[edge];
                if (edge_map == -1)
                  edge_coef[i] = 0; // 赋值PEC的边
                else
                  edge_coef[i] = coefofunit * (modal_vec_pointer)[edge_map] / port_beta;
              }

              for (int i = 0; i < 3; i++) {
                int node = face_node_idx[face_node_ext[glo_ff] + i];
                int node_map = modal_dof_map_node[node];
                if (node_map == -1)
                  node_coef[i] = 0;
                else
                  node_coef[i] = coefofunit * (modal_vec_pointer)[node_map] * dcomplex(0, 1);
              }
              const int nnode = 3;
              // 记得在老版本的DoubleVector中把.h文件拷贝到新版本的库里
              tbox::Array<hier::DoubleVector<NDIM> > vertex_3d(nnode);
              for (int i = 0; i < nnode; ++i) {
                int node_id = face_node_idx[face_node_ext[glo_ff] + i];
                for (int k = 0; k < NDIM; ++k) {
                  vertex_3d[i][k] = (*node_coord)(k, node_id);
                }
              }
              tbox::Array<hier::DoubleVector<NDIM - 1> > vertex_2d(nnode);
              // 转到端口面的二维局部坐标
              Transfer3DCoordTo2D(vertex_3d, now_num_wp->direction, vertex_2d);
              double area = fabs(0.5 * (AREA2DQUAD(vertex_2d[0], vertex_2d[1], vertex_2d[2])));
              bool raw_order = AREA2DQUAD(vertex_2d[0], vertex_2d[1], vertex_2d[2]) > 0;
              if (!raw_order) {
                std::swap(vertex_2d[1], vertex_2d[2]);
              }
              tbox::Array<tbox::Array<double> > BasGrad(3);
              GradientOn2DCoord(vertex_2d, BasGrad);
              int faceorder = 0;
              if (raw_order)
                faceorder = 1;
              quad.faceQuadfunctionDotBas_modal<appu::Nedelec, appu::triNedelec, dcomplex, dcomplex>(
                  cell, loc_ff, &edge_coef[0], &node_coef[0], &shapefunc, &trishapefunc, BasGrad,
                  faceorder, 0, &(face_fun_x_phi[0]), now_num_wp->direction);
            } // end of is_input_port
          }
        }
      }
    }

    if (patch.getPatchGeometry()->hasEntitySet(d_user_param.outer_face_set_id,
                                               hier::EntityUtilities::FACE,
                                               patch.getNumberOfFaces(1))) {
      // 获取本网格片的所有边界面的编号（包含影像区）
      tbox::Array<int> outer_face = patch.getPatchGeometry()->getEntityIndicesInSet(
          d_user_param.outer_face_set_id, hier::EntityUtilities::FACE, patch.getNumberOfFaces(1));
      // 为方便查找，首先对边界面编号进行排序。
      std::sort(outer_face.getPointer(), outer_face.getPointer() + outer_face.size());

      int nface = 4;
      // 对本单元的所有的面进行循环
      for (int face = 0; face < nface; face++) {
        double tmp_face_phi_x_phi[nedge][nedge];
        NTYPE tmp_face_fun_x_phi[nedge];

        // 获取本单元的第face个面在网格片的编号
        int patch_face_id = cell_face_idx[cell_face_ext[cell] + face];
        // 判断该面是否属于边界面。
        bool is_in_outer_face = std::binary_search(
            outer_face.getPointer(), outer_face.getPointer() + outer_face.size(), patch_face_id);
        // 如果面不属于边界面，直接跳过该面的积分。
        if (!is_in_outer_face) {
          continue;
        }
        // 计算单刚中的面积分项
        quad.faceQuadBasDotBas(cell, face, &shapefunc, 2, &(tmp_face_phi_x_phi[0][0]));
        // 计算单荷中的面积分项
        d_user_param.ftype = d_user_param.BOUNDARY_3;
        quad.faceQuadFunctionDotBas<appu::Nedelec, NTYPE>(cell, face, d_user_param, &shapefunc, 2,
                                                          &(tmp_face_fun_x_phi[0]));

        // 累加到边界积分的单刚和单荷上。
        for (int i = 0; i < nedge; i++) {
          face_fun_x_phi[i] += tmp_face_fun_x_phi[i];
          for (int j = 0; j < nedge; j++) {
            face_phi_x_phi[i][j] += tmp_face_phi_x_phi[i][j];
          }
        }

#define ANSWER_CHECK 0
#if ANSWER_CHECK
        int bdry_edge[] = {1232, 1246, 1911};
        tbox::Array<int> face_node_ext, face_node_idx;
        patch.getPatchTopology()->getFaceAdjacencyNodes(face_node_ext, face_node_idx);
        bool equal = true;
        for (int i = face_node_ext[patch_face_id]; i < face_node_ext[patch_face_id]; i++) {
          if (!std::binary_search(bdry_edge, bdry_edge + 3, face_node_idx[i])) {
            equal = false;
          }
        }
        if (equal) {
          // print tmp_face_phi_x_phi, and tmp_face_fun_x_phi
          NTYPE A3[nedge][nedge], B1[nedge];
          dcomplex coef2 = dcomplex(0, 1) * d_user_param.get_k_0() / d_user_param.eta;
          for (int i = 0; i < nedge; i++) {
            for (int j = 0; j < nedge; j++)
              A3[i][j] = coef2 * face_phi_x_phi[i][j];
          }
          for (int i = 0; i < nedge; i++) {
            B1[i] = face_fun_x_phi[i];
          }
          double t;
          t = sin(4);
        }
#endif
      }
    }

    // CELLPART 计算单刚和单荷
    if (cell_part_id[cell] == 1) {
      dcomplex coef1 = -pow(d_user_param.get_k_0(), 2.0) * d_user_param.epsilon_r;
      dcomplex coef2 = dcomplex(0, 1) * d_user_param.get_k_0() / d_user_param.eta;
      for (int i = 0; i < nedge; i++) {
        b[i] = fun_x_phi[i] - face_fun_x_phi[i];
        for (int j = 0; j < nedge; j++) {
          A[i][j] = curl_x_curl[i][j] + coef1 * phi_x_phi[i][j] + coef2 * face_phi_x_phi[i][j];
          A_aux[i][j] =
              curl_x_curl[i][j] + abs(coef1) * phi_x_phi[i][j] + abs(coef2) * face_phi_x_phi[i][j];
        }
      }
    } else if (cell_part_id[cell] == 2) {
      dcomplex coef1 = -pow(d_user_param.get_k_0(), 2.0) * d_user_param.epsilon_r;
      dcomplex coef2 = dcomplex(0, 1) * d_user_param.get_k_0() / d_user_param.eta;
      for (int i = 0; i < nedge; i++) {
        b[i] = fun_x_phi[i] - face_fun_x_phi[i];
        for (int j = 0; j < nedge; j++) {
          A[i][j] = curl_x_curl[i][j] + coef1 * phi_x_phi[i][j] + coef2 * face_phi_x_phi[i][j];
          A_aux[i][j] =
              curl_x_curl[i][j] + abs(coef1) * phi_x_phi[i][j] + abs(coef2) * face_phi_x_phi[i][j];
        }
      }
    } else {
      TBOX_WARNING("The cell " << cell << " does not belong to any part, please check!");
    }

#if ANSWER_CHECK
    NTYPE A1[nedge][nedge], A2[nedge][nedge], A3[nedge][nedge];
    for (int i = 0; i < nedge; i++) {
      for (int j = 0; j < nedge; j++)
        A1[i][j] = curl_x_curl[i][j];
    }

    for (int i = 0; i < nedge; i++) {
      for (int j = 0; j < nedge; j++)
        A2[i][j] = coef1 * phi_x_phi[i][j];
    }

    for (int i = 0; i < nedge; i++) {
      for (int j = 0; j < nedge; j++)
        A3[i][j] = coef2 * face_phi_x_phi[i][j];
    }
#endif
    // Apply constraint boundary condition.
    // We apply Dirichlet BC for Maxwell equations: n x u = n x g
    // here by letting: u_i = t . g, where t is the i-th edge vector.
    // FIXME: Is there better formula?
    for (int i = 0; i < nedge; i++) {
      int eidx = cell_edge_idx[cell_edge_ext[cell] + i];
      if ((*edge_flag)(0, eidx)) {
        // cout << "dirichlet edge is " <<  eidx << endl;
        //    a11 a12 a13  b1     a11  0  a12  b1-v2*a12
        //    a21 a22 a33  b2 ==>  0   1   0       v2
        //    a31 a32 a33  b3     a31  0  a33  b3-v2*a32
        // NTYPE v = (*Ea)(0, eidx);
        NTYPE v = 0; // 第一类边界条件为零。
        for (int j = 0; j < nedge; j++)
          b[j] = (j == i ? v : (b[j] - A[j][i] * v));
        for (int j = 0; j < nedge; j++) {
          A[i][j] = (i == j ? 1. : 0.);
          A_aux[i][j] = (i == j ? 1. : 0.);
        }
        for (int j = 0; j < nedge; j++) {
          A[j][i] = (i == j ? 1. : 0.);
          A_aux[j][i] = (i == j ? 1. : 0.);
        }
      }
    }

    // element index to patch matrix index mapping
    for (int i = 0; i < nedge; i++) {
      int edge = cell_edge_idx[cell_edge_ext[cell] + i];
      mapping[i] = dof_map[edge];
    }

    /*
    // add to patch matrix
    for(int i = 0; i < nedge; i++) {
      int edge = cell_edge_idx[cell_edge_ext[cell] + i];
      if(edge >= num_edge) // only non-ghost row will assemble
        continue;
      for(int j = 0; j < nedge; j++)
        mat->addMatrixValue(mapping[i], mapping[j],A[i][j]);
    }

    // add to patch vector
    for(int i = 0; i < nedge; i++) {
      int edge = cell_edge_idx[cell_edge_ext[cell] + i];
      if(edge >= num_edge) // only non-ghost row will assemble
        continue;
      vec->addVectorValue(mapping[i], b[i]);
    }
    */

    // complex 类操作  ->real();   ->image(); conj(complex  )
    // add to patch matrix(real and image)
    for (int i = 0; i < nedge; i++) {
      int edge = cell_edge_idx[cell_edge_ext[cell] + i];
      if (edge >= num_edge) // only non-ghost row will assemble
        continue;
      for (int j = 0; j < nedge; j++) {
        double val_re = (A[i][j]).real();
        double val_im = (A[i][j]).imag();
        mat_re->addMatrixValue(mapping[i], mapping[j], val_re);
        mat_im->addMatrixValue(mapping[i], mapping[j], val_im);
        pre_mat->addMatrixValue(mapping[i], mapping[j], A_aux[i][j]);
      }
    }

    // add to patch vector(real and image)
    for (int i = 0; i < nedge; i++) {
      int edge = cell_edge_idx[cell_edge_ext[cell] + i];
      if (edge >= num_edge) // only non-ghost row will assemble
        continue;
      double val_re = (b[i]).real();
      double val_im = (b[i]).imag();
      vec_re->addVectorValue(mapping[i], val_re);
      vec_im->addVectorValue(mapping[i], val_im);
    }
  }

  /*  // 矩阵组装。
  mat->assemble();
  */
  mat_re->assemble();
  mat_im->assemble();
  pre_mat->assemble();
}

/*************************************************************************
 *  建立网格片上的电磁系统离散化正定方程.
 ************************************************************************/
void Maxwell::buildEMMatrixOnPatch(hier::Patch<NDIM> &patch, const double time, const double dt) {
  int num_cell_ghost = patch.getNumberOfCells(1);
  int num_edge = patch.getNumberOfEdges(0);
  // 自由度映射信息
  int *dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::EDGE);
  // 向量数据片 实部和虚部
  tbox::Pointer<pdat::VectorData<NDIM, double> > vec_re = patch.getPatchData(d_rhs_re_id);
  tbox::Pointer<pdat::VectorData<NDIM, double> > vec_im = patch.getPatchData(d_rhs_im_id);

  // 矩阵数据片 实部和虚部
  tbox::Pointer<pdat::CSRMatrixData<NDIM, double> > mat_re = patch.getPatchData(d_mat_re_id);
  tbox::Pointer<pdat::CSRMatrixData<NDIM, double> > mat_im = patch.getPatchData(d_mat_im_id);
  // 预处理矩阵
  tbox::Pointer<pdat::CSRMatrixData<NDIM, double> > pre_mat = patch.getPatchData(d_premat_id);

  // 四面体Nedelec基函数
  appu::Nedelec shapefunc(patch, patch.getPatchData(d_edge_order_id),
                          patch.getPatchData(d_jacobian_id));
  // 三角形Nedelec基函数
  appu::triNedelec trishapefunc(patch, patch.getPatchData(d_edge_order_id),
                                patch.getPatchData(d_jacobian_id));
  // 积分信息
  appu::TetQuad quad(patch, patch.getPatchData(d_volume_id), patch.getPatchData(d_jacobian_id));
  GET_PATCH_DATA(patch, edge_flag, d_edge_flag_id, Edge, int);
  DECLARE_ADJACENCY(patch, cell, face, Cell, Face);
  DECLARE_ADJACENCY(patch, cell, edge, Cell, Edge);
  DECLARE_ADJACENCY(patch, face, edge, Face, Edge);
  DECLARE_ADJACENCY(patch, face, node, Face, Node);
  int step = time / dt + 0.5;
  double freq = frequency_setup[step];
  double omega = 2 * M_PI * freq;
  double mu_0 = 4 * M_PI * 1e-7;
  double epsr_max = 10.;
  double epsilon_0 = 8.8542 * 1e-12;
  double k_0_square = (omega) * (omega)*mu_0 * epsilon_0;
  double zita = k_0_square * epsr_max;
  // 先提取功率参数
  double input_power = NUM_WP_list[0]->power;
  NTYPE coefofunit = 1.0 * sqrt(2. * input_power) / (sqrt(port_power[active_mode]));

  GET_PATCH_DATA(patch, epsilonr_data, d_epsilonr_id, Cell, double);
  GET_PATCH_DATA(patch, mur_data, d_mur_id, Cell, double);
  GET_PATCH_DATA(patch, sigma_data, d_sigma_id, Cell, double);
  GET_PATCH_DATA(patch, losstan_data, d_losstan_id, Cell, double);

  tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
      patch.getPatchGeometry()->getNodeCoordinates();
  for (int cell = 0; cell < num_cell_ghost; cell++) {
    const int nedge = (NDIM == 2 ? 3 : 6);
    double curl_x_curl[nedge][nedge], phi_x_phi[nedge][nedge], face_phi_x_phi[nedge][nedge];
    NTYPE fun_x_phi[nedge], face_fun_x_phi[nedge];
    NTYPE A[nedge][nedge], b[nedge];
    double A_aux[nedge][nedge];
    int mapping[nedge];
    // 旋度项
    quad.quadCurlBasDotCurlBas<appu::Nedelec>(cell, &shapefunc, 0, &(curl_x_curl[0][0]));
    // 位移项
    quad.quadBasDotBas<appu::Nedelec>(cell, &shapefunc, 2, &(phi_x_phi[0][0]));
    double alpha = (*sigma_data)(0, cell) * omega * mu_0;
    NTYPE pure_imag(0, 1);
    NTYPE alpha_imag = pure_imag * alpha;
    NTYPE beta_imag = k_0_square * ((*epsilonr_data)(0, cell) * (1 + (*losstan_data)(0, cell)));
    NTYPE coef = k_0_square * beta_imag - alpha_imag; // 矩阵方程系数
    for (int i = 0; i < nedge; i++) {
      b[i] = 0;
      for (int j = 0; j < nedge; j++) {
        A[i][j] = curl_x_curl[i][j] - coef * phi_x_phi[i][j];
        A_aux[i][j] = curl_x_curl[i][j] + abs(coef) * phi_x_phi[i][j];
      }
    }
    /*************************************************************************
     *  数值波端口的边界条件
     ************************************************************************/

    // 数值波端口部分部分
    for (int num_port = 0; num_port < numerical_waveport_num; num_port++) {
      // 取出目前的wave port
      tbox::Pointer<NUM_WP> now_num_wp = NUM_WP_list[num_port];
      // 判断该端口是否为激励端口 (activated==true: 入射激励; false: 仅匹配BC)
      const bool is_input_port = now_num_wp->activated;

      // 取出模式dof映射
      int *modal_dof_map_edge =
          d_dof_info_modal_array[num_port]->getDOFMapping(patch, hier::EntityUtilities::EDGE);
      int *modal_dof_map_node =
          d_dof_info_modal_array[num_port]->getDOFMapping(patch, hier::EntityUtilities::NODE);
      double port_beta = sqrt(zita - zita / eigen_value_pairs[num_port][active_mode]);

      // 吸收边界条件：对所有模式求和 gama，确保端口面正确吸收所有传播模式
      dcomplex gama_sum(0, 0);
      for (int m = 0; m < eigen_num; m++) {
        double beta_m = sqrt(zita - zita / eigen_value_pairs[num_port][m]);
        gama_sum += dcomplex(0, 1) * beta_m;
      }
      for (int face_list_num = 0; face_list_num < now_num_wp->face_list.getSize();
           face_list_num++) {
        if (HAS_ENTITY_SET(patch, now_num_wp->face_list[face_list_num], FACE, 1)) {
          DECLARE_ENTITY_SET(patch, port_i_list, now_num_wp->face_list[face_list_num], FACE, 1);
          std::sort(port_i_list.getPointer(), port_i_list.getPointer() + port_i_list.size());
          int nface = 4;
          // 对本单元的所有的面进行循环
          for (int loc_ff = 0; loc_ff < nface; loc_ff++) {
            int glo_ff = cell_face_idx[cell_face_ext[cell] + loc_ff];
            bool is_on_port_face = std::binary_search(
                port_i_list.getPointer(), port_i_list.getPointer() + port_i_list.size(), glo_ff);
            if (!is_on_port_face)
              continue;
            // BC 项：所有端口（激励+匹配）都需要
            quad.faceQuadBasDotBas(cell, loc_ff, &shapefunc, 2, &(face_phi_x_phi[0][0]));

            // 激励项：仅输入端口需要
            if (is_input_port) {
              GET_PATCH_DATA(patch, modal_sol, d_sol_modal_array[num_port][active_mode], Vector, double);
              double *modal_vec_pointer = modal_sol->getPointer();
              NTYPE edge_coef[3] = {0, 0, 0}; // 面上边的插值系数
              NTYPE node_coef[3] = {0, 0, 0}; // 面上结点的插值系数
              for (int i = 0; i < 3; i++) {
                int edge = face_edge_idx[face_edge_ext[glo_ff] + i];
                int edge_map = modal_dof_map_edge[edge];
                if (edge_map == -1)
                  edge_coef[i] = 0; // 赋值PEC的边
                else
                  edge_coef[i] = coefofunit * (modal_vec_pointer)[edge_map] / port_beta;
              }

              for (int i = 0; i < 3; i++) {
                int node = face_node_idx[face_node_ext[glo_ff] + i];
                int node_map = modal_dof_map_node[node];
                if (node_map == -1)
                  node_coef[i] = 0;
                else
                  node_coef[i] = coefofunit * (modal_vec_pointer)[node_map] * dcomplex(0, 1);
              }
              const int nnode = 3;
              // 记得在老版本的DoubleVector中把.h文件拷贝到新版本的库里
              tbox::Array<hier::DoubleVector<NDIM> > vertex_3d(nnode);
              for (int i = 0; i < nnode; ++i) {
                int node_id = face_node_idx[face_node_ext[glo_ff] + i];
                for (int k = 0; k < NDIM; ++k) {
                  vertex_3d[i][k] = (*node_coord)(k, node_id);
                }
              }
              tbox::Array<hier::DoubleVector<NDIM - 1> > vertex_2d(nnode);
              // 转到端口面的二维局部坐标
              Transfer3DCoordTo2D(vertex_3d, now_num_wp->direction, vertex_2d);
              bool raw_order = AREA2DQUAD(vertex_2d[0], vertex_2d[1], vertex_2d[2]) > 0;
              if (!raw_order) {
                std::swap(vertex_2d[1], vertex_2d[2]);
              }
              tbox::Array<tbox::Array<double> > BasGrad(3);
              GradientOn2DCoord(vertex_2d, BasGrad);
              int faceorder = 0;
              if (raw_order)
                faceorder = 1;
              quad.faceQuadfunctionDotBas_modal<appu::Nedelec, appu::triNedelec, dcomplex,
                                                dcomplex>(
                  cell, loc_ff, &edge_coef[0], &node_coef[0], &shapefunc, &trishapefunc, BasGrad,
                  faceorder, 0, &(face_fun_x_phi[0]), now_num_wp->direction);
              dcomplex gama_in = dcomplex(0, 1) * port_beta;
              for (int i = 0; i < nedge; i++) {
                b[i] += ((-2.0 * gama_in * face_fun_x_phi[i]));
              }
            } // end of is_input_port

            double abs_gama_sum = sqrt(gama_sum.real() * gama_sum.real() +
                                       gama_sum.imag() * gama_sum.imag());
            for (int i = 0; i < nedge; i++) {
              for (int j = 0; j < nedge; j++) {
                A[i][j] += (gama_sum * face_phi_x_phi[i][j]);
                A_aux[i][j] += (abs_gama_sum * face_phi_x_phi[i][j]);
              }
            }
          } // end of face
        } // end of have entity set
      } // end of face list
    } // end of port

    /*************************************************************************
     *  最外层吸收边界条件
     ************************************************************************/

    // 吸收边界部分
    for (int abc_ff = 0; abc_ff < ABC_face.getSize(); abc_ff++) {
      if (HAS_ENTITY_SET(patch, ABC_face[abc_ff], FACE, 1)) {
        DECLARE_ENTITY_SET(patch, ABC_i_face, ABC_face[abc_ff], FACE, 1);
        std::sort(ABC_i_face.getPointer(), ABC_i_face.getPointer() + ABC_i_face.size());
        int nface = 4;
        for (int loc_ff = 0; loc_ff < nface; loc_ff++) {
        }
      }
    }
  }
}

/*************************************************************************
 * 计算边的方向：
 * 通过坐标确定全局一致的边的方向。
 * 注意：这种做法可能会因为坐标在网格片间存在偏差
 *       导致方向无法全局一致（对网格坐标不须更新
 *       的应用程序，可以认为不会发生这种情况）。
 *       将来JAUMIN框架须通过更为稳健的算法向用户
 *       提供边的全局一致方向。
 ************************************************************************/
void Maxwell::computeEdgeOrderOnPatch(hier::Patch<NDIM> &patch) {
  int num_edge_ghost = patch.getNumberOfEdges(1);
  tbox::Array<int> edge_node_ext, edge_node_idx;
  patch.getPatchTopology()->getEdgeAdjacencyNodes(edge_node_ext, edge_node_idx);
  tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
      patch.getPatchGeometry()->getNodeCoordinates();
  tbox::Pointer<pdat::EdgeData<NDIM, int> > edge_order = patch.getPatchData(d_edge_order_id);

  for (int edge = 0; edge < num_edge_ghost; edge++) {
    int v0 = edge_node_idx[edge_node_ext[edge] + 0];
    int v1 = edge_node_idx[edge_node_ext[edge] + 1];
    (*edge_order)(0, edge) = -1;
    for (int d = 0; d < NDIM; d++) {
      if ((*node_coord)(d, v0) == (*node_coord)(d, v1))
        continue;
      (*edge_order)(0, edge) = ((*node_coord)(d, v0) < (*node_coord)(d, v1) ? 1 : 0);
      break;
    }
    TBOX_ASSERT((*edge_order)(0, edge) != -1);
  }
}

/*************************************************************************
 * 边的边界标识
 ************************************************************************/
void Maxwell::computeEdgeFlagOnPatch(hier::Patch<NDIM> &patch) {
  tbox::Pointer<pdat::EdgeData<NDIM, int> > edge_flag = patch.getPatchData(d_edge_flag_id);
  GET_PATCH_DATA(patch, port_edge_flag, d_port_edge_flag_id, Edge, int);
  edge_flag->fillAll(0);
  port_edge_flag->fillAll(0);

  if (!patch.getPatchGeometry()->hasEntitySet(
          d_user_param.inner_face_set_id, hier::EntityUtilities::FACE, patch.getNumberOfFaces(1)))
    return;

  int num_edge = patch.getNumberOfEdges(0);
  tbox::Array<int> face_edge_ext, face_edge_idx;
  patch.getPatchTopology()->getFaceAdjacencyEdges(face_edge_ext, face_edge_idx);
  const tbox::Array<int> &bdry_face =
      patch.getEntityIndicesInSet(d_user_param.inner_face_set_id, hier::EntityUtilities::FACE, 1);

  for (int f = 0; f < bdry_face.size(); f++) {
    int nedge = 3;
    int face = bdry_face[f];
    for (int e = 0; e < nedge; e++) {
      int edge = face_edge_idx[face_edge_ext[face] + e];
      if (edge >= num_edge)
        continue;
      (*edge_flag)(0, edge) = 1;
    }
  }
  for (int num_wp_num = 0; num_wp_num < NUM_WP_list.getSize(); num_wp_num++) {
    tbox::Pointer<NUM_WP> now_num_wp = NUM_WP_list[num_wp_num];
    for (int edge_list_num = 0; edge_list_num < now_num_wp->edge_list.getSize(); edge_list_num++) {
      if (HAS_ENTITY_SET(patch, now_num_wp->edge_list[edge_list_num] + entity_face_num, EDGE, 1)) {
        DECLARE_ENTITY_SET(patch, port_i_list,
                           now_num_wp->edge_list[edge_list_num] + entity_face_num, EDGE, 1);
        for (int ee = 0; ee < port_i_list.size(); ee++) {
          (*port_edge_flag)(0, ee) = num_wp_num * 2;
        }
      }
    }
    for (int face_list_num = 0; face_list_num < now_num_wp->face_list.getSize(); face_list_num++) {
      if (HAS_ENTITY_SET(patch, now_num_wp->face_list[face_list_num], FACE, 1)) {
        DECLARE_ENTITY_SET(patch, port_i_list, now_num_wp->face_list[face_list_num], FACE, 1);
        for (int ff = 0; ff < port_i_list.size(); ff++) {
          for (int loc_ee = face_edge_ext[ff]; loc_ee < face_edge_ext[ff + 1]; loc_ee++) {
            int glo_ee = face_edge_idx[loc_ee];
            (*port_edge_flag)(0, glo_ee) = num_wp_num * 2 + 1;
          }
        }
      }
    }
  }
}

/*************************************************************************
 * 将解向量数据片复制到棱心量数据片.
 ************************************************************************/
void Maxwell::acceptSolutionOnPatch(hier::Patch<NDIM> &patch) {
  int num_edge = patch.getNumberOfEdges(0);
  tbox::Pointer<pdat::EdgeData<NDIM, NTYPE> > E = patch.getPatchData(d_E_id);
  //  tbox::Pointer< pdat::VectorData<NDIM,NTYPE> > solution =
  //  patch.getPatchData(d_solution_id); tbox::Pointer<
  //  pdat::EdgeData<NDIM,double> > edge_coord =
  //  patch.getPatchGeometry()->getEdgeCoordinates(); tbox::Pointer<
  //  pdat::EdgeData<NDIM,NTYPE> > Ea = patch.getPatchData(d_Ea_id);
  //  tbox::Pointer< pdat::EdgeData<NDIM,int> > edge_flag =
  //  patch.getPatchData(d_edge_flag_id);

  // 获取数值解的实部和虚部
  tbox::Pointer<pdat::VectorData<NDIM, double> > sol_re = patch.getPatchData(d_sol_re_id);
  tbox::Pointer<pdat::VectorData<NDIM, double> > sol_im = patch.getPatchData(d_sol_im_id);

  // 自由度信息中的映射信息。
  int *dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::EDGE);

  for (int edge = 0; edge < num_edge; edge++) {
    int dof = dof_map[edge];
    //    (*E)(0, edge) = (solution->getPointer())[dof];

    double val_re = (sol_re->getPointer())[dof];
    double val_im = (sol_im->getPointer())[dof];
    (*E)(0, edge) = dcomplex(val_re, val_im);
  }
}

/*************************************************************************
 * 将棱心量数据片转为单元中心量数据片用于可视化输出.
 ************************************************************************/
void Maxwell::postprocessOnPatch(hier::Patch<NDIM> &patch) {
  int num_cell = patch.getNumberOfCells(0);
  tbox::Pointer<pdat::EdgeData<NDIM, NTYPE> > E = patch.getPatchData(d_E_id);
  tbox::Pointer<pdat::CellData<NDIM, NTYPE> > Ec = patch.getPatchData(d_Ec_id);
  tbox::Pointer<pdat::CellData<NDIM, double> > Ec_norm = patch.getPatchData(d_Ec_norm_id);
  tbox::Pointer<pdat::CellData<NDIM, double> > Ec_arg = patch.getPatchData(d_Ec_arg_id);
  tbox::Array<int> cell_edge_ext, cell_edge_idx;
  patch.getPatchTopology()->getCellAdjacencyEdges(cell_edge_ext, cell_edge_idx);

  // 处理基函数的值和curl值
  appu::Nedelec shapefunc(patch, patch.getPatchData(d_edge_order_id),
                          patch.getPatchData(d_jacobian_id));

  for (int cell = 0; cell < num_cell; cell++) {
    int nedge = (NDIM == 2 ? 3 : 6);
    NTYPE dof[nedge];

    for (int edge = 0; edge < nedge; edge++)
      dof[edge] = (*E)(0, cell_edge_idx[cell_edge_ext[cell] + edge]);

    double lambda[] = {0.25, 0.25, 0.25, 0.25};
    double phi[nedge][NDIM];
    shapefunc.basis(cell, lambda, &(phi[0][0]));

    NTYPE v[] = {0., 0., 0.};
    for (int edge = 0; edge < nedge; edge++) {
      v[0] += dof[edge] * phi[edge][0];
      v[1] += dof[edge] * phi[edge][1];
      v[2] += dof[edge] * phi[edge][2];
    }

    for (int i = 0; i < NDIM; i++) {
      (*Ec)(i, cell) = NTYPE(v[i]);
      (*Ec_norm)(i, cell) = norm(NTYPE(v[i]));
      (*Ec_arg)(i, cell) = arg(NTYPE(v[i]));
    }
  }
}

/*************************************************************************
 * 计算Hcurl误差: || E - E_h ||_hcurl.
 ************************************************************************/
void Maxwell::computeHcurlErrOnPatch(double *vector, int len, hier::Patch<NDIM> &patch) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(len == 1);
#endif

  int num_cell = patch.getNumberOfCells(0);
  tbox::Array<int> cell_edge_ext, cell_edge_idx;
  patch.getPatchTopology()->getCellAdjacencyEdges(cell_edge_ext, cell_edge_idx);
  tbox::Pointer<pdat::EdgeData<NDIM, NTYPE> > Ea = patch.getPatchData(d_Ea_id);
  tbox::Pointer<pdat::EdgeData<NDIM, NTYPE> > E = patch.getPatchData(d_E_id);
  tbox::Pointer<pdat::CellData<NDIM, double> > cell_coord =
      patch.getPatchGeometry()->getCellCoordinates();

  // 处理基函数的值和curl值
  // Shape Function
  appu::Nedelec shapefunc(patch, patch.getPatchData(d_edge_order_id),
                          patch.getPatchData(d_jacobian_id));

  // 处理高斯积分
  //  Quadrature Rule
  appu::TetQuad quad(patch, patch.getPatchData(d_volume_id), patch.getPatchData(d_jacobian_id));

  double err = 0.;
  for (int cell = 0; cell < num_cell; cell++) {
    int nedge = (NDIM == 2 ? 3 : 6);
    NTYPE dof_err[nedge];
    double err_hcurl = 0.;
    for (int edge = 0; edge < nedge; edge++) {
      int eidx = cell_edge_idx[cell_edge_ext[cell] + edge];
      dof_err[edge] = (*E)(0, eidx) - (*Ea)(0, eidx);
    }
    // compute \int \curl(E - Ea) \cdot \curl(E - Ea)
    double curlbas_x_curlbas[nedge][nedge];
    quad.quadCurlBasDotCurlBas(cell, &shapefunc, 0, &(curlbas_x_curlbas[0][0]));
    for (int i = 0; i < nedge; i++) {
      for (int j = 0; j < nedge; j++) {
        err_hcurl += (conj(dcomplex(dof_err[i])) * curlbas_x_curlbas[i][j] * dof_err[j]).real();
      }
    }
    // quad.quadCurlDofDotCurlDof<appu::Nedelec>(cell, &shapefunc, 0, dof_err,
    // dof_err, &err_hcurl);
    err += err_hcurl;
  }

  *vector = err;
}

/**
 * @brief PatchStrategy::InitMaterialComponent
 * @param patch
 * 材料初始化模块的功能
 */
void Maxwell::InitMaterialComponent(tbox::Pointer<tbox::Database> db) {
  if (db->isDatabase("Material")) {
    tbox::Pointer<tbox::Database> material_db = db->getDatabase("Material");
    num_total_material = material_db->getInteger("material_num");
    material_entity.resizeArray(num_total_material);
    Material_list.resizeArray(num_total_material);
    for (int mat_num = 0; mat_num < num_total_material; mat_num++) {
      stringstream materialname;
      materialname << "material" << mat_num + 1;
      stringstream materialentity;
      materialentity << "material_" << mat_num + 1;
      material_entity[mat_num] = material_db->getIntegerArray(materialentity.str());
      tbox::Pointer<tbox::Database> now_material_db = material_db->getDatabase(materialname.str());
      tbox::Pointer<Material> material_i = new Material;
      material_i->mat_id = mat_num + 1;
      material_i->d_mu_r = now_material_db->getDouble("mu_r");
      material_i->d_epsilon_r = now_material_db->getDouble("epsilon_r");
      material_i->d_sigma = now_material_db->getDouble("sigma");
      material_i->d_losstan = now_material_db->getDouble("losstan");
      Material_list[mat_num] = material_i;
    }
  } else {
    TBOX_ERROR(d_object_name << ": " << " No key `Material' found in input." << endl);
  }
}

/**
 * @brief PatchStrategy::InitMaterialComponent
 * @param patch
 * 边界条件初始化模块的功能
 */
void Maxwell::InitBoundaryInfo(tbox::Pointer<tbox::Database> db) {
  if (db->isDatabase("Boundary")) {
    tbox::Pointer<tbox::Database> boundary_db = db->getDatabase("Boundary");
    ABC_face = boundary_db->getIntegerArray("ABC_face");
    PEC_face = boundary_db->getIntegerArray("PEC_face");
  } else {
    TBOX_ERROR(d_object_name << ": " << " No key `Boundary' found in input." << endl);
  }
}

void Maxwell::InitPortInfo(tbox::Pointer<tbox::Database> db) {
  if (db->isDatabase("Port")) {
    tbox::Pointer<tbox::Database> port_db = db->getDatabase("Port");
    /// 解析波端口数量
    analytical_waveport_num = port_db->getInteger("analytical_waveport_num");
    /// 数值波端口数量
    numerical_waveport_num = port_db->getInteger("numerical_waveport_num");
    /// 集总端口数量
    lumpport_num = port_db->getInteger("lumpport_num");

    ANAL_WP_list.resizeArray(analytical_waveport_num);
    NUM_WP_list.resizeArray(numerical_waveport_num);
    LP_list.resizeArray(lumpport_num);
    d_dof_info_modal_array.resizeArray(numerical_waveport_num);
    eigen_value_pairs.resizeArray(numerical_waveport_num);
    port_power.resizeArray(numerical_waveport_num * eigen_num);
    for (int num_port = 0; num_port < numerical_waveport_num; num_port++) {
      d_dof_info_modal_array[num_port] = new solv::DOFInfo<NDIM>(true, true, false, false);
      eigen_value_pairs[num_port].resizeArray(eigen_num);
    }
    for (int anal_wp_num = 0; anal_wp_num < analytical_waveport_num; anal_wp_num++) {
    }

    for (int num_wp_num = 0; num_wp_num < numerical_waveport_num; num_wp_num++) {
      stringstream portname;
      portname << "NUM_WP" << num_wp_num + 1;
      tbox::Pointer<tbox::Database> now_port_db = port_db->getDatabase(portname.str());
      tbox::Pointer<NUM_WP> num_wp_i = new NUM_WP;
      num_wp_i->face_list = now_port_db->getIntegerArray("face_list");
      num_wp_i->edge_list = now_port_db->getIntegerArray("edge_list");
      num_wp_i->direction = now_port_db->getDoubleArray("direction");
      num_wp_i->power = now_port_db->getDouble("power");
      num_wp_i->activated = now_port_db->getBool("activated");
      NUM_WP_list[num_wp_num] = num_wp_i;
    }
    for (int lp_num = 0; lp_num < lumpport_num; lp_num++) {
    }
  } else {
    TBOX_ERROR(d_object_name << ": " << " No key `Port' found in input." << endl);
  }
}

void Maxwell::setupWPModalDofInfoOnPatch(hier::Patch<NDIM> &patch) {
  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
  int num_nodes = patch.getNumberOfNodes(1);
  int num_edges = patch.getNumberOfEdges(1);
  DECLARE_ADJACENCY(patch, edge, node, Edge, Node);
  DECLARE_ADJACENCY(patch, face, node, Face, Node);
  GET_PATCH_DATA(patch, port_edge_flag, d_port_edge_flag_id, Edge, int);
  for (int num_port = 0; num_port < numerical_waveport_num; num_port++) {
    tbox::Pointer<solv::DOFInfo<NDIM> > this_port_dof = d_dof_info_modal_array[num_port];
    int *dis_ptr_modal_edge = this_port_dof->getDOFDistribution(patch, hier::EntityUtilities::EDGE);
    int *dis_ptr_modal_node = this_port_dof->getDOFDistribution(patch, hier::EntityUtilities::NODE);
    // 首先所有节点赋为0自由度
    for (int i = 0; i < num_nodes; ++i) {
      dis_ptr_modal_node[i] = 0;
    }
    for (int face_list_num = 0; face_list_num < NUM_WP_list[num_port]->face_list.getSize();
         face_list_num++) {
      if (HAS_ENTITY_SET(patch, face_list_num < NUM_WP_list[num_port]->face_list[face_list_num],
                         FACE, 1)) {
        DECLARE_ENTITY_SET(patch, wp_face_list,
                           face_list_num < NUM_WP_list[num_port]->face_list[face_list_num], FACE,
                           1);
        // 首先所有节点赋为0自由度
        for (int i = 0; i < num_nodes; ++i) {
          dis_ptr_modal_node[i] = 0;
        }
        // 将所有端口面赋予自由度
        for (int ff = 0; ff < wp_face_list.getSize(); ff++) {
          /// 遍历所有面上的结点
          for (int loc_nn = face_node_ext[ff]; loc_nn < face_node_ext[ff + 1]; loc_nn++) {
            int glo_nn = face_node_idx[loc_nn];
            dis_ptr_modal_node[glo_nn] = 1;
          }
        }
      }
    }

    for (int glo_ee = 0; glo_ee < num_edges; ++glo_ee) {
      if ((*port_edge_flag)(0, glo_ee) == 2 * num_port) {
        dis_ptr_modal_edge[glo_ee] = 1;
      } else {
        dis_ptr_modal_edge[glo_ee] = 0;
      }
      if ((*port_edge_flag)(0, glo_ee) == 2 * num_port + 1) {
        int node1 = edge_node_idx[edge_node_ext[glo_ee]];
        int node2 = edge_node_idx[edge_node_ext[glo_ee] + 1];
        dis_ptr_modal_node[node1] = 0;
        dis_ptr_modal_node[node2] = 0;
      }
    }
    this_port_dof->buildPatchDOFMapping(patch);
    d_dof_info_modal_array[num_port] = this_port_dof;
  }
}

void Maxwell::setupModalMatrixOnPatch(hier::Patch<NDIM> &patch, double time, const double dt) {
  int num_edge = patch.getNumberOfEdges(0);
  /// 取出本地PatchGeometry.
  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
  /// 取出本地PatchTopology.
  tbox::Pointer<hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();
  /// 取出本地Patch的结点坐标数组.
  tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord = patch_geo->getNodeCoordinates();
  /// 取出本地Patch的结点坐标数组.
  tbox::Pointer<pdat::EdgeData<NDIM, double> > edge_coord = patch_geo->getEdgeCoordinates();
  DECLARE_ADJACENCY(patch, cell, node, Cell, Node);
  DECLARE_ADJACENCY(patch, cell, edge, Cell, Edge);
  DECLARE_ADJACENCY(patch, cell, face, Cell, Face);

  DECLARE_ADJACENCY(patch, face, node, Face, Node);
  DECLARE_ADJACENCY(patch, face, edge, Face, Edge);
  DECLARE_ADJACENCY(patch, face, cell, Face, Cell);

  DECLARE_ADJACENCY(patch, edge, node, Edge, Node);
  DECLARE_ADJACENCY(patch, edge, cell, Edge, Cell);
  // 参数定义
  int step = time / dt + 0.5;
  double freq = frequency_setup[step];
  double omega = 2 * M_PI * freq;
  double mu_0 = 4 * M_PI * 1e-7;
  double epsr_max = 10.;
  double epsilon_0 = 8.8542 * 1e-12;
  double k_0_square = (omega) * (omega)*mu_0 * epsilon_0;
  double zita = k_0_square * epsr_max;
  GET_PATCH_DATA(patch, epsilonr_data, d_epsilonr_id, Cell, double);
  GET_PATCH_DATA(patch, mur_data, d_mur_id, Cell, double);
  GET_PATCH_DATA(patch, sigma_data, d_sigma_id, Cell, double);
  GET_PATCH_DATA(patch, losstan_data, d_losstan_id, Cell, double);
  GET_PATCH_DATA(patch, matid_data, d_mat_id, Cell, int);

  appu::triNedelec shapefunc(patch, patch.getPatchData(d_edge_order_id),
                             patch.getPatchData(d_jacobian_id));
  appu::TriQuad quad(patch, patch.getPatchData(d_volume_id), patch.getPatchData(d_jacobian_id));

  for (int num_port = 0; num_port < numerical_waveport_num; num_port++) {
    // 遍历所有的端口
    int *dof_map_modal_edge =
        d_dof_info_modal_array[num_port]->getDOFMapping(patch, hier::EntityUtilities::EDGE);
    int *dof_map_modal_node =
        d_dof_info_modal_array[num_port]->getDOFMapping(patch, hier::EntityUtilities::NODE);
    GET_PATCH_DATA(patch, mat_A, d_matA_modal_array[num_port], CSRMatrix, double);
    GET_PATCH_DATA(patch, mat_B, d_matB_modal_array[num_port], CSRMatrix, double);
    mat_A->fill(0.0);
    mat_B->fill(0.0);
    tbox::Pointer<NUM_WP> now_num_wp = NUM_WP_list[num_port];
    for (int face_list_num = 0; face_list_num < now_num_wp->face_list.getSize(); face_list_num++) {
      if (HAS_ENTITY_SET(patch, now_num_wp->face_list[face_list_num], FACE, 1)) {
        DECLARE_ENTITY_SET(patch, port_i_list, now_num_wp->face_list[face_list_num], FACE, 1);
        for (int ff = 0; ff < port_i_list.getSize(); ff++) {
          const int nnode = 3;
          const int nedge = 3;

          double NODEBasdotNODEBas[nnode][nnode] = {0};
          double NODEGraddotNODEGrad[nnode][nnode] = {0};
          double EDGEBasdotEDGEBas[nedge][nedge] = {0};
          double EDGECurldotEDGECurl[nedge][nedge] = {0};
          double EDGECurldotNODEGrad[nedge][nnode] = {0};
          double NODEGraddotEDGECurl[nnode][nedge] = {0};
          int mapping_edge[nedge];
          int mapping_node[nnode];

          // 获取面所属的三棱柱体单元
          int glo_ff = port_i_list[ff];
          int glo_cc = face_cell_idx[face_cell_ext[glo_ff]];
          double epsr = (*epsilonr_data)(0, glo_cc);
          double mur = 1.0;
          double Att[nedge][nedge], Btt[nedge][nedge];
          double Btz[nedge][nnode], Bzt[nnode][nedge], Bzz[nnode][nnode];
          // 记得在老版本的DoubleVector中把.h文件拷贝到新版本的库里
          tbox::Array<hier::DoubleVector<NDIM> > vertex_3d(nnode);
          for (int i = 0; i < nnode; ++i) {
            int node_id = face_node_idx[face_node_ext[glo_ff] + i];
            for (int k = 0; k < NDIM; ++k) {
              vertex_3d[i][k] = (*node_coord)(k, node_id);
            }
          }
          tbox::Array<hier::DoubleVector<NDIM - 1> > vertex_2d(nnode);
          // 转到端口面的二维局部坐标
          Transfer3DCoordTo2D(vertex_3d, now_num_wp->direction, vertex_2d);
          double area = fabs(0.5 * (AREA2DQUAD(vertex_2d[0], vertex_2d[1], vertex_2d[2])));
          bool raw_order = AREA2DQUAD(vertex_2d[0], vertex_2d[1], vertex_2d[2]) > 0;
          if (!raw_order) {
            std::swap(vertex_2d[1], vertex_2d[2]);
          }
          tbox::Array<tbox::Array<double> > BasGrad(3);
          GradientOn2DCoord(vertex_2d, BasGrad);
          quad.nodeBasdotBas<double>(3, raw_order, area, &NODEBasdotNODEBas[0][0]);
          quad.nodeGraddotGrad<double>(3, raw_order, area, BasGrad, &NODEGraddotNODEGrad[0][0]);
          quad.edgeBasDotedgeBas<appu::triNedelec>(glo_ff, 3, raw_order, area, &shapefunc,
                                                   &EDGEBasdotEDGEBas[0][0], BasGrad);
          quad.edgeCurlDotedgeCurl<appu::triNedelec>(glo_ff, 3, raw_order, area, &shapefunc,
                                                     &EDGECurldotEDGECurl[0][0], BasGrad);
          quad.edgeBasDotNodeGrad<appu::triNedelec>(glo_ff, 3, raw_order, area, &shapefunc,
                                                    &EDGECurldotNODEGrad[0][0], BasGrad);
          quad.NodeGradDotedgeBas<appu::triNedelec>(glo_ff, 3, raw_order, area, &shapefunc,
                                                    &NODEGraddotEDGECurl[0][0], BasGrad);
          for (int i = 0; i < nedge; i++) {
            for (int j = 0; j < nedge; j++) {
              Att[i][j] = (1.0 / mur) * EDGECurldotEDGECurl[i][j] -
                          k_0_square * epsr * EDGEBasdotEDGEBas[i][j];
              Btt[i][j] = (1.0 / mur) * EDGEBasdotEDGEBas[i][j];
            }
          }
          for (int i = 0; i < nnode; i++) {
            for (int j = 0; j < nnode; j++) {
              Bzz[i][j] = (1.0 / mur) * NODEGraddotNODEGrad[i][j] -
                          k_0_square * epsr * NODEBasdotNODEBas[i][j];
            }
          }
          for (int i = 0; i < nedge; i++) {
            for (int j = 0; j < nnode; j++) {
              Btz[i][j] = (1.0 / mur) * EDGECurldotNODEGrad[i][j];
            }
          }
          for (int i = 0; i < nnode; i++) {
            for (int j = 0; j < nedge; j++) {
              Bzt[i][j] = (1.0 / mur) * NODEGraddotEDGECurl[i][j];
            }
          }
          for (int i = 0; i < nedge; i++) {
            int edge = face_edge_idx[face_edge_ext[glo_ff] + i];
            mapping_edge[i] = dof_map_modal_edge[edge];
          }
          for (int i = 0; i < nnode; i++) {
            int node = face_node_idx[face_node_ext[glo_ff] + i];
            mapping_node[i] = dof_map_modal_node[node];
          }

          // 边-边映射关系
          for (int i = 0; i < nedge; i++) {
            for (int j = 0; j < nedge; j++) {
              if (mapping_edge[i] == -1 || mapping_edge[j] == -1)
                continue;
              mat_B->addMatrixValue(mapping_edge[i], mapping_edge[j], Btt[i][j] + Att[i][j] / zita);
              mat_A->addMatrixValue(mapping_edge[i], mapping_edge[j], Btt[i][j]);
            }
          }

          // 点-边映射关系
          for (int i = 0; i < nnode; i++) {
            for (int j = 0; j < nedge; j++) {
              if (mapping_node[i] == -1 || mapping_edge[j] == -1)
                continue;
              mat_B->addMatrixValue(mapping_node[i], mapping_edge[j], Bzt[i][j]);
              mat_A->addMatrixValue(mapping_node[i], mapping_edge[j], Bzt[i][j]);
            }
          }

          // 边-点映射关系
          for (int i = 0; i < nedge; i++) {
            for (int j = 0; j < nnode; j++) {
              if (mapping_edge[i] == -1 || mapping_node[j] == -1)
                continue;
              mat_B->addMatrixValue(mapping_edge[i], mapping_node[j], Btz[i][j]);
              mat_A->addMatrixValue(mapping_edge[i], mapping_node[j], Btz[i][j]);
            }
          }

          // 点-点映射关系
          for (int i = 0; i < nnode; i++) {
            for (int j = 0; j < nnode; j++) {
              if (mapping_node[i] == -1 || mapping_node[j] == -1)
                continue;
              mat_B->addMatrixValue(mapping_node[i], mapping_node[j], Bzz[i][j]);
              mat_A->addMatrixValue(mapping_node[i], mapping_node[j], Bzz[i][j]);
            }
          }
        }
      }
    }
    mat_A->assemble();
    mat_B->assemble();
  }
}

void Maxwell::updateMaterialProperties(hier::Patch<NDIM> &patch) {
  GET_PATCH_DATA(patch, epsilonr_data, d_epsilonr_id, Cell, double);
  GET_PATCH_DATA(patch, mur_data, d_mur_id, Cell, double);
  GET_PATCH_DATA(patch, sigma_data, d_sigma_id, Cell, double);
  GET_PATCH_DATA(patch, losstan_data, d_losstan_id, Cell, double);
  GET_PATCH_DATA(patch, matid_data, d_mat_id, Cell, int);
  for (int imaterial = 0; imaterial < Material_list.getSize(); imaterial++) {
    tbox::Array<int> material_i_list = material_entity[imaterial];
    if (HAS_ENTITY_SET(patch, material_i_list[imaterial], CELL, 1)) {
      DECLARE_ENTITY_SET(patch, material_i_cell_list, material_i_list[imaterial], CELL, 1);
      for (int jcount = 0; jcount < material_i_cell_list.getSize(); jcount++) {
        int glo_cc = material_i_cell_list[jcount];
        (*epsilonr_data)(0, glo_cc) = Material_list[imaterial]->d_epsilon_r;
        (*mur_data)(0, glo_cc) = Material_list[imaterial]->d_mu_r;
        (*losstan_data)(0, glo_cc) = Material_list[imaterial]->d_losstan;
        (*sigma_data)(0, glo_cc) = Material_list[imaterial]->d_sigma;
        (*matid_data)(0, glo_cc) = Material_list[imaterial]->mat_id;
      }
    }
  }
}

void Maxwell::calculateModalPower(double *vector, int len, hier::Patch<NDIM> &patch,
                                  const double time, const double dt) {
  /// 取出本地PatchGeometry.
  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
  /// 取出本地PatchTopology.
  tbox::Pointer<hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();
  /// 取出本地Patch的结点坐标数组.
  tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord = patch_geo->getNodeCoordinates();
  DECLARE_ADJACENCY(patch, cell, node, Cell, Node);
  DECLARE_ADJACENCY(patch, cell, edge, Cell, Edge);
  DECLARE_ADJACENCY(patch, cell, face, Cell, Face);

  DECLARE_ADJACENCY(patch, face, node, Face, Node);
  DECLARE_ADJACENCY(patch, face, edge, Face, Edge);
  DECLARE_ADJACENCY(patch, face, cell, Face, Cell);

  DECLARE_ADJACENCY(patch, edge, node, Edge, Node);
  DECLARE_ADJACENCY(patch, edge, cell, Edge, Cell);
  // 参数定义
  int step = time / dt + 0.5;
  double freq = frequency_setup[step];
  double omega = 2 * M_PI * freq;
  double mu_0 = 4 * M_PI * 1e-7;
  double epsr_max = 10.;
  double epsilon_0 = 8.8542 * 1e-12;
  double k_0_square = (omega) * (omega)*mu_0 * epsilon_0;
  double zita = k_0_square * epsr_max;
  GET_PATCH_DATA(patch, epsilonr_data, d_epsilonr_id, Cell, double);
  GET_PATCH_DATA(patch, mur_data, d_mur_id, Cell, double);
  GET_PATCH_DATA(patch, sigma_data, d_sigma_id, Cell, double);
  GET_PATCH_DATA(patch, losstan_data, d_losstan_id, Cell, double);
  for (int num_port = 0; num_port < numerical_waveport_num; num_port++) {
    const int nnode = 3;
    const int nedge = 3;
    tbox::Pointer<NUM_WP> now_num_wp = NUM_WP_list[num_port];
    // 遍历所有的端口
    int *dof_map_modal_edge =
        d_dof_info_modal_array[num_port]->getDOFMapping(patch, hier::EntityUtilities::EDGE);
    int *dof_map_modal_node =
        d_dof_info_modal_array[num_port]->getDOFMapping(patch, hier::EntityUtilities::NODE);
    for (int num_eigen = 0; num_eigen < eigen_num; num_eigen++) {
      GET_PATCH_DATA(patch, modal_sol, d_sol_modal_array[num_port][num_eigen], Vector, double);
      appu::TriQuad quad(patch, patch.getPatchData(d_volume_id), patch.getPatchData(d_jacobian_id));
      appu::triNedelec shapefunc(patch, patch.getPatchData(d_edge_order_id),
                                 patch.getPatchData(d_jacobian_id));

      double *modal_vec_pointer = modal_sol->getPointer();
      double beta = sqrt(zita - zita / eigen_value_pairs[num_port][num_eigen]);
      for (int face_list_num = 0; face_list_num < now_num_wp->face_list.getSize();
           face_list_num++) {
        if (HAS_ENTITY_SET(patch, now_num_wp->face_list[face_list_num], FACE, 1)) {
          DECLARE_ENTITY_SET(patch, port_i_list, now_num_wp->face_list[face_list_num], FACE, 1);
          for (int ff = 0; ff < port_i_list.getSize(); ff++) {
            int glo_ff = port_i_list[ff];
            double edge_e_coef[3] = {0, 0, 0};
            double node_e_coef[3] = {0, 0, 0};
            // 轮流对场赋值
            for (int i = 0; i < 3; i++) {
              int glo_ee = face_edge_idx[face_edge_ext[glo_ff] + i];
              int edge_mapping = dof_map_modal_edge[glo_ee];
              // 如果这条边没有自由度的话就是PEC边
              if (edge_mapping < 0)
                edge_e_coef[i] = 0;
              else
                edge_e_coef[i] = modal_vec_pointer[edge_mapping];
            }

            for (int i = 0; i < 3; i++) {
              int glo_nn = face_node_idx[face_node_ext[glo_ff] + i];
              int node_mapping = dof_map_modal_node[glo_nn];
              // 如果这条边没有自由度的话就是PEC边
              if (node_mapping < 0)
                node_e_coef[i] = 0;
              else
                node_e_coef[i] = modal_vec_pointer[node_mapping];
            }
            // 记得在老版本的DoubleVector中把.h文件拷贝到新版本的库里
            tbox::Array<hier::DoubleVector<NDIM> > vertex_3d(nnode);
            for (int i = 0; i < nnode; ++i) {
              int node_id = face_node_idx[face_node_ext[glo_ff] + i];
              for (int k = 0; k < NDIM; ++k) {
                vertex_3d[i][k] = (*node_coord)(k, node_id);
              }
            }
            tbox::Array<hier::DoubleVector<NDIM - 1> > vertex_2d(nnode);
            // 转到端口面的二维局部坐标
            Transfer3DCoordTo2D(vertex_3d, now_num_wp->direction, vertex_2d);
            double area = fabs(0.5 * (AREA2DQUAD(vertex_2d[0], vertex_2d[1], vertex_2d[2])));
            bool raw_order = AREA2DQUAD(vertex_2d[0], vertex_2d[1], vertex_2d[2]) > 0;
            tbox::Array<tbox::Array<double> > BasGrad(3);
            if (!raw_order) {
              std::swap(vertex_2d[1], vertex_2d[2]);
            }
            GradientOn2DCoord(vertex_2d, BasGrad);
            NTYPE facepower = 0.;
            quad.faceEPower(glo_ff, 3, raw_order, area, &shapefunc, &edge_e_coef[0],
                            &node_e_coef[0], &facepower, BasGrad, beta, omega, mu_0);
            vector[2 * (num_port * eigen_num + num_eigen)] += facepower.real();
            vector[2 * (num_port * eigen_num + num_eigen) + 1] += facepower.imag();
          }
        }
      }
    }
  }
}