//
// 文件名: main.C
// 软件包: JAUMIN application
// 版权  : (c) 2004-2015 北京应用物理与计算数学研究所
//         (c) 2013-2015 中物院高性能数值模拟软件中心
// 版本号: $Revision$
// 修改  : $Data$
// 描述  : 主控程序(四面体非结构网格上, 求解时谐麦克斯维方程问题).
//

#include "utils/JAUMIN_App.h"
using namespace JAUMIN;

#include "MaxwellLevelIntegrator.h"
#include "Maxwell.h"

/*!
*************************************************************************
* @brief   在网格文件名前面追加input文件的路径
************************************************************************
*/
static void 
prefixInputDirName(const string& input_filename, 
                   tbox::Pointer<tbox::Database> input_db)
{
  string path_name = "";
  string::size_type slash_pos = input_filename.find_last_of( '/' );
  if(slash_pos != string::npos)
    path_name = input_filename.substr(0, slash_pos + 1);      

  string mesh_file = input_db->getDatabase("GridGeometry")
    ->getDatabase("MeshImportationParameter")
    ->getString("file_name");

  slash_pos = mesh_file.find_first_of( '/' );
  if(slash_pos != 0){
    input_db->getDatabase("GridGeometry")
      ->getDatabase("MeshImportationParameter")
      ->putString("file_name", path_name + mesh_file);
  }
}


/*!
*************************************************************************
*                                                                     
* @brief 基于JAUMIN框架的非结构网格, 求解时谐麦克斯维方程问题.         
*
* 该程序分以下几个步骤:
* -# 预处理: 初始化MPI和JAUMIN环境, 解析输入文件, 读取主程序控制参数;
* -# 创建网格片层次结构和时间积分算法类对象, 主要包括:
*       - 网格几何       hier::GridGeometry<NDIM>
*       - 网格拓扑       hier::GridTopology<NDIM>
*       - 网格片层次结构 hier::PatchHierarchy<NDIM>
*    -# 网格片积分算法 Maxwell 
*       - 应用级: 提供求解时谐麦克斯维方程方程的数值计算子程序
*    -# 网格层积分算法 MaxwellLevelIntegrator
*       - 应用级: 提供基于网格层的时谐麦克斯维方程方程求解流程
*    -# 网格层时间积分算法 algs::HierarchyTimeIntegrator<NDIM>
* -# 初始化网格片层次结构和物理量数据片;
* -# 循环: 时间步进;
* -# 后处理: 释放应用类对象, 释放JAUMIN和MPI内部资源.
*                                                                      
************************************************************************
*/
int main( int argc, char *argv[]) {
  // 初始化MPI和JAUMIN环境.
  tbox::MPI::init(&argc, &argv);
  tbox::JAUMINManager::startup();

  {
    /*******************************************************************************
     *                               预  处  理                                    *
     *******************************************************************************/
    // 解析命令行参数:
    bool is_from_restart = false;         // 此次运行是否为重启动.
    string input_filename;                // 输入文件名.
    string restart_read_dirname;          // 重启动文件所在路径
    int restore_num = -1;                 // 重启动时间序号
    if ( (argc != 2) && (argc != 4) ) {
      tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
      << "<restart dir> <restore number> [options]\n"
      << "  options:\n"
      << "  none at this time"
      << endl;
      tbox::MPI::abort();
      return (-1);
    } else {
      input_filename = argv[1];

      if (argc == 4) {  // 重启动情形
        is_from_restart = true;
        restart_read_dirname = argv[2];    // 重启动文件所在路径
        restore_num = atoi(argv[3]);       // 重启动时间序号
      }
    }

    /// 把信息输出到log文件
    tbox::plog << "input_filename = " << input_filename << endl;
    if ( is_from_restart ) {
      tbox::plog << "restart_read_dirname = " << restart_read_dirname << endl;
      tbox::plog << "restore_num = " << restore_num << endl;
    }

    // 解析输入文件的计算参数到输入数据库, 称之为根数据库.
    tbox::Pointer<tbox::Database> input_db = new tbox::InputDatabase("input_db");
    tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

    // 在网格文件名前面追加input文件的路径
    prefixInputDirName(input_filename, input_db);

    // 从根数据库中获得名称为"Main"的子数据库.
    tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

    // 从"Main"子数据库中获取日志文件控制参数.
    string log_file_name = main_db->getString("log_file_name");
    bool log_all_nodes   = main_db->getBoolWithDefault("log_all_nodes",false);
    if (log_all_nodes) {
      tbox::PIO::logAllNodes(log_file_name);
    } else {
      tbox::PIO::logOnlyNodeZero(log_file_name);
    }

    // 从"Main"子数据库中获取可视化输出控制参数.
    int javis_dump_interval = main_db->getIntegerWithDefault("javis_dump_interval",0);
    string javis_dump_dirname;
    int javis_number_procs_per_file = -1;
    if ( javis_dump_interval > 0 ) {
      javis_dump_dirname = main_db->getString("javis_dump_dirname");
      javis_number_procs_per_file =
        main_db->getIntegerWithDefault("javis_number_procs_per_file",1);
    }

    // 从"Main"子数据库中获取重启动输出控制参数.
    int restart_dump_interval = main_db->getIntegerWithDefault("restart_dump_interval",0);
    string restart_dump_dirname;
    if ( restart_dump_interval > 0 ) {
      restart_dump_dirname = main_db->getString("restart_dump_dirname");
    }

    // 打开重启动输入文件
    tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();
    if ( is_from_restart ) {
      restart_manager->openRestartFile(restart_read_dirname,
                                       restore_num,
                                       tbox::MPI::getNodes() );
    }

    // 创建计时器.
    tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
    tbox::Pointer<tbox::Timer> t_write_javis_data =
      tbox::TimerManager::getManager()->getTimer("apps::main::write_javis_data");
    tbox::Pointer<tbox::Timer> t_write_restart =
      tbox::TimerManager::getManager()->getTimer("apps::main::write_restart_data");
    /*******************************************************************************
     *                     创建网格片层次结构和积分算法类对象                      *
     *******************************************************************************/
    //(1) 创建非结构网格几何.
    tbox::Pointer<hier::GridGeometry<NDIM> > grid_geometry =
      new hier::GridGeometry<NDIM>("GridGeometry",
                                   input_db->getDatabase("GridGeometry"));

    //(2) 创建非结构网格拓扑.
    tbox::Pointer<hier::GridTopology<NDIM> > grid_topology =
      new hier::GridTopology<NDIM>("GridTopology",
                                   input_db->getDatabase("GridTopology"));

    //(3) 创建网格片层次结构.
    tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy =
      new hier::PatchHierarchy<NDIM>("PatchHierarchy",
                                     grid_geometry,
                                     grid_topology,
                                     true);

    //(4) 创建网格片时间积分算法类（应用级: 提供求解时谐麦克斯维方程方程的数值计算子程序）.
    Maxwell* maxwell_model = new Maxwell("Maxwell",
                                         input_db->getDatabase("Maxwell"));

    //(5) 创建网格层时间积分算法类（应用级: 提供基于网格层的时谐麦克斯维方程方程求解流程）.
    MaxwellLevelIntegrator* maxwell_level_integrator =
      new MaxwellLevelIntegrator("MaxwellLevelIntegrator",
                                 maxwell_model,
				 input_db->getDatabase("Solver"));

#ifdef DEBUG_RESTART
        cout << "++++++In Maxwellection main(), before new algs::HierarchyTimeIntegrator<NDIM>()." << endl;
#endif

    //(6) 创建网格层时间积分算法类.
    tbox::Pointer<algs::HierarchyTimeIntegrator<NDIM> > time_integrator =
      new algs::HierarchyTimeIntegrator<NDIM>("HierarchyTimeIntegrator",
                                              input_db->getDatabase("HierarchyTimeIntegrator"),
                                              patch_hierarchy,
                                              maxwell_level_integrator,
                                              true);

    /*******************************************************************************
     *  可视化输出器类.
     *******************************************************************************/
    tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_data_writer;
    if (javis_dump_interval > 0) {
      javis_data_writer = new appu::JaVisDataWriter<NDIM>("Maxwell_JaVis_Writer",
                                                          javis_dump_dirname,
                                                          javis_number_procs_per_file);

      maxwell_model->registerPlotData(javis_data_writer); //注册待输出的绘图量.
    }

    /*******************************************************************************
     *                初 始 化 网 格 层 和 物 理 量                    *
     *******************************************************************************/
    // 初始化网格片层次结构和物理量.
    time_integrator->initializeHierarchy();

#ifdef DEBUG_RESTART
        cout << "++++++In Maxwell main(), after time_integrator->initializeHierarchy()." << endl;
#endif

    // 输出初始化信息到日志文件.
    tbox::plog << "\nCheck input data and variables before simulation:" << endl;
    tbox::plog << "Input database..." << endl;
    input_db->printClassData(tbox::plog);
    tbox::plog << "\nVariable database..." << endl;
    hier::VariableDatabase<NDIM>::getDatabase()->printClassData(tbox::plog);

    // 输出可视化数据场到可视化输出器.
    if (javis_dump_interval > 0) {
      t_write_javis_data->start();
      javis_data_writer->writePlotData(patch_hierarchy,
                                       time_integrator->getIntegratorStep(),
                                       time_integrator->getIntegratorTime());
      t_write_javis_data->stop();
    }

    // 关闭重启动输入文件
    if (is_from_restart) {
      restart_manager->closeRestartFile();
    }

    /************************************************************************************
     *                              时  间  步  进                                      *
     ************************************************************************************/

    double loop_time     = time_integrator->getIntegratorTime(); // 当前模拟时间
    double loop_time_end = time_integrator->getEndTime();        // 模拟结束时间

    int iteration_num    = time_integrator->getIntegratorStep(); // 当前时间步数

#ifdef DEBUG_RESTART
        cout << "++++++In Maxwell main(), before while()" << endl
                    << "loop_time = " << loop_time << endl
                    << "loop_time_end = " << loop_time_end << endl
                    << "iteration_num = " << iteration_num << endl
                    << "time_integrator->stepsRemaining() = " << time_integrator->stepsRemaining() << endl;
#endif

    while ((loop_time < loop_time_end) && time_integrator->stepsRemaining()) {
      iteration_num = time_integrator->getIntegratorStep() + 1;

      tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
      tbox::pout << "At begining of timestep # " <<  iteration_num - 1 << endl;
      tbox::pout << "Simulation time is " << loop_time << endl;
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;

#ifdef DEBUG_RESTART
        cout << "++++++In Maxwell main(), before time_integrator->advanceHierarchy()." << endl;
#endif

      // 将数值解推进一个时间步, 返回其时间步长.
      double dt_actual = time_integrator->advanceHierarchy();

#ifdef DEBUG_RESTART
        cout << "++++++In Maxwell main(), after time_integrator->advanceHierarchy()." << endl;
#endif

      loop_time += dt_actual;

      tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
      tbox::pout << "At end of timestep # " <<  iteration_num - 1 << endl;
      tbox::pout << "Dt = " << dt_actual << ", Simulation time is " << loop_time << endl;
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;

      if (iteration_num % 10 == 0)
        tbox::TimerManager::getManager()->print(tbox::plog);
      // 输出重启动数据.
      if ((restart_dump_interval>0) && ((iteration_num % restart_dump_interval) == 0)) {
        t_write_restart->start();
        restart_manager->writeRestartFile(restart_dump_dirname,iteration_num);
        t_write_restart->stop();
      }

      // 输出可视化数据.
      if ((javis_dump_interval > 0) && ((iteration_num % javis_dump_interval) == 0 ||
          !((loop_time < loop_time_end) && time_integrator->stepsRemaining()))) {
        t_write_javis_data->start();
        javis_data_writer->writePlotData(patch_hierarchy,
                                         iteration_num,
                                         loop_time);
        t_write_javis_data->stop();
      }
    }

    /************************************************************************************
     *                               模  拟  结  束                                     *
     ************************************************************************************/
    // 释放对象.
    if (maxwell_level_integrator) {
      delete maxwell_level_integrator;
    }

    if (maxwell_model) {
      delete maxwell_model;
    }

    // 输出计时器统计的时间数据.
    tbox::TimerManager::getManager()->print(tbox::plog);
    tbox::MemoryUtilities::printMemoryInfo(tbox::plog);
	tbox::MemoryUtilities::printMaxMemory(tbox::plog);
  }

  // 释放JAUMIN和MPI内部资源.
  tbox::JAUMINManager::shutdown();
  tbox::MPI::finalize();

  return(0);
}
