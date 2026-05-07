//
// 文件名:	MaxwellLevelIntegrator.h
// 软件包:	JAUMIN applications
// 版权  :      (c) 2004-2015 北京应用物理与计算数学研究所
//              (c) 2013-2015 中物院高性能数值模拟软件中心
// 版本号:	$Revision:$
// 修改  :	$Date:$
// 描述  :	麦克斯维网格层时间积分算法
//

#ifndef included_MaxwellLevelIntegrator
#define included_MaxwellLevelIntegrator
 
#include <map>
#include <string>
#include "utils/JAUMIN_App.h"
#include "JAUMINVector.h"
#include "JAUMINMatrix.h"
#include "ComplexLinearSolver.h"

using namespace JAUMIN;

#include "Maxwell.h"

/**
 * @brief 该类从网格层时间积分算法策略类 algs::TimeIntegratorLevelStrategy 派生,
 * 实现麦克斯维方程在网格层上的求解流程.
 */
class MaxwellLevelIntegrator 
: public algs::TimeIntegratorLevelStrategy<NDIM>
{
public:
   /**
    * @brief 构造函数.
    * @param object_name     输入参数, 字符串, 表示对象名称.
    * @param patch_strategy  输入参数, 指针, 麦克斯维方程网格片积分算法.
    * @param solver_db       输入参数, 数据库指针, 指向解法器的数据库. 
    */     
   MaxwellLevelIntegrator(const string& object_name, 
                          Maxwell* patch_strategy,
                          tbox::Pointer<tbox::Database> solver_db);

   /**
    * @brief 析构函数.
    */
   virtual ~MaxwellLevelIntegrator();

   ///@name 重载基类algs::TimeIntegratorLevelStrategy<NDIM>的函数
   //@{

   /**
    * @brief 初始化该积分算法: 创建所有计算需要的积分构件.
    * @param manager 输入参数, 指针, 指向积分构件管理器.
    */
   void initializeLevelIntegrator(
           tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager);

   /**
    * @brief 初始化指定网格层的数据片.
    *
    * @param level          输入参数, 指针,         指向待初始化网格层.
    * @param init_data_time 输入参数, 双精度浮点型, 表示初始化的时刻.
    * @param initial_time   输入参数, 逻辑型,       真值表示当前时刻为计算的初始时刻.
    */
   void initializeLevelData(const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                            const double init_data_time,
                            const bool   initial_time = true);

   /**
    * @brief 返回指定网格层的时间步长. 
    *
    * @param level        输入参数, 指针,         指向网格层.
    * @param dt_time      输入参数, 双精度浮点型, 表示计算时间步长的当前时刻.
    * @param initial_time 输入参数, 逻辑型,       真值表示当前时刻为计算的初始时刻.
    * @param flag_last_dt 输入参数, 整型,         表示上个时间步积分返回的状态.
    * @param last_dt      输入参数, 双精度浮点型, 表示上个时间步长.
    *
    * @return 双精度浮点型, 表示网格层的时间步长.
    */
   double getLevelDt(const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                     const double dt_time, 
                     const bool initial_time,
                     const int  flag_last_dt,
                     const double last_dt);

   /**
    * @brief 网格层向前积分一个时间步. 
    *
    * @param level        输入参数, 指针,         指向待积分的网格层.
    * @param current_time 输入参数, 双精度浮点型, 表示时间步的起始时刻.
    * @param predict_dt   输入参数, 双精度浮点型, 表示为该时间步预测的时间步长.
    * @param max_dt       输入参数, 双精度浮点型, 表示时间步允许的最大时间步长.
    * @param min_dt       输入参数, 双精度浮点型, 表示时间步允许的最小时间步长.
    * @param first_step   输入参数, 逻辑型,       真值当前为重构后或时间步序列的第1步.
    * @param step_number  输入参数, 整型,         表示积分步数.
    * @param actual_dt    输出参数, 双精度浮点型, 表示时间步实际采用的时间步长.
    *
    * @return 整型, 表示该时间步积分的状态. 
    *
    * @note
    * 相关的网格片计算函数见 Maxwell::computeOnPatch().
    */ 
   int advanceLevel(const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                    const double current_time, 
                    const double predict_dt,
                    const double max_dt,
                    const double min_dt,
                    const bool   first_step, 
                    const int    step_number,
                    double &     actual_dt);

   /**
    * @brief 更新网格层的状态到新的时刻.
    * 
    * @param level           输入参数, 指针, 指向网格层.
    * @param new_time        输入参数, 双精度浮点型, 表示新的时刻.
    * @param deallocate_data 输入参数, 逻辑型, 真值表示接收数值解后, 释放新值数据片的内存空间.
    */
   void acceptTimeDependentSolution(const tbox::Pointer< hier::BasePatchLevel<NDIM> > level, 
                                    const double new_time, 
                                    const bool deallocate_data);

   //@}
 private:
   /*!@brief 对象名称. */
   string d_object_name; 

   /*!@brief 麦克斯维的网格片积分算法. */
   Maxwell* d_patch_strategy;

   /*!@brief 解法器: 求解矩阵系统 */
   /* old solve
   //   tbox::Pointer<tbox::Database>                    d_solver_db; 
   //   tbox::Pointer<JPSOL::JPLinearSolver<NDIM, dcomplex> >     d_solver;
   //   tbox::Pointer<solv::JAUMINVector<NDIM, dcomplex> > d_sol;
   //   tbox::Pointer<solv::JAUMINMatrix<NDIM, dcomplex> > d_matrix;

   //   tbox::Pointer<solv::JAUMINVector<NDIM, dcomplex> > d_rhs;
   */
   
   /*!@brief 解法器: 求解矩阵系统 */
   tbox::Pointer<tbox::Database>                    d_solver_db; 

   tbox::Pointer<solv::JAUMINMatrix<NDIM, double> > d_mat_re;
   tbox::Pointer<solv::JAUMINMatrix<NDIM, double> > d_mat_im;      
   tbox::Pointer<solv::JAUMINVector<NDIM, double> > d_sol_re;
   tbox::Pointer<solv::JAUMINVector<NDIM, double> > d_sol_im;
   tbox::Pointer<solv::JAUMINVector<NDIM, double> > d_rhs_re;
   tbox::Pointer<solv::JAUMINVector<NDIM, double> > d_rhs_im;

   //预条件矩阵
   tbox::Pointer<solv::JAUMINMatrix<NDIM, double> > d_pre_matrix;

   //   AMS 方法辅助变量
   tbox::Pointer<solv::JAUMINVector<NDIM, double> > d_xcoord;
   tbox::Pointer<solv::JAUMINVector<NDIM, double> > d_ycoord;
   tbox::Pointer<solv::JAUMINVector<NDIM, double> > d_zcoord;
   tbox::Pointer<solv::JAUMINMatrix<NDIM, double> > d_gmat;   
   tbox::Pointer<JPSOL::ComplexLinearSolver<NDIM, double> >     d_solver;

   
   /*!@brief 内存构件 */
   std::map<std::string, tbox::Pointer<algs::MemoryIntegratorComponent<NDIM> > >     d_mem_intc;
   /*!@brief 内存构件 */
   std::map<std::string, tbox::Pointer<algs::InitializeIntegratorComponent<NDIM> > > d_ini_intc;
   /*!@brief 数值构件(专门用来遍历网格片进行计算) */
   std::map<std::string, tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > >  d_num_intc;
   /*!@brief 数值构件(专门用来填充影像区) */
   std::map<std::string, tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > >  d_com_intc;
   /*!@brief 规约构件 */
   std::map<std::string, tbox::Pointer<algs::ReductionIntegratorComponent<NDIM> > >  d_red_intc;
};

#endif
