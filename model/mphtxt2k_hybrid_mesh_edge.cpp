//
// Created by liekkas on 2026/1/26.
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip> // I/O流控制头文件（类似于格式化输出）

using namespace std;

typedef unsigned int uint;

// const uint N = 13421772;

const uint MAX_NUM_INLINE = 8; // 描述单元-结点关系中单元内结点的最大数量（六面体为8）

const uint NUM_NODE_EDG = 2;
const uint NUM_NODE_TRI = 3;
const uint NUM_NODE_QUAD = 4;
const uint NUM_NODE_TET = 4;
const uint NUM_NODE_PYR = 5;
const uint NUM_NODE_PRISM = 6;
const uint NUM_NODE_HEX = 8;

string STR_NUM_VTX("# number of mesh vertices");        // 顶点数量
string STR_COORDS("# Mesh vertex coordinates");         // 顶点坐标
string STR_VTX("1 # number of vertices per element");   // 顶点
string STR_EDG("2 # number of vertices per element");   // 线段
string STR_TRI("3 # number of vertices per element");   // 三角形
string STR_QUAD("4 # number of vertices per element");  // 四边形
string STR_TET("4 # number of vertices per element");   // 四面体
string STR_PYR("5 # number of vertices per element");   // 金字塔
string STR_PRISM("6 # number of vertices per element"); // 三棱柱
string STR_HEX("8 # number of vertices per element");   // 六面体

enum class ReadState : uint
{
    Standby,
    ReadCoordinates,
    CheckType,
    ReadVtx,
    ReadEdg,
    ReadQuad,
    ReadHex,
    ReadTri,
    ReadTet,
    ReadPyr,
    ReadPrism,
    End
};

struct NumOfModel
{
    uint n_nodes = 0; // number of nodes
    uint n_element_types = 0;
    uint n_boby_elements = 0;
    // number of all kinds of elements
    uint n_vtx = 0;   // 顶点 1 # number of nodes per element
    uint n_edg = 0;   // 棱边 2 # number of nodes per element
    uint n_tri = 0;   // 三角形 3 # number of nodes per element
    uint n_quad = 0;  // 四边形 4 # number of nodes per element
    uint n_tet = 0;   // 四面体 4 # number of nodes per element
    uint n_pyr = 0;   // 金字塔 5 # number of nodes per element
    uint n_prism = 0; // 三棱柱 6 # number of nodes per element
    uint n_hex = 0;   // 六面体 8 # number of nodes per element
};

struct InputInfo
{
    int BoundaryNum = 0;    // 边界条件数量
    double len_scale = 1.0; // 模型尺寸系数
    string in_filename;     // name of input file
    string out_filename;    // name of output file
    vector<uint> BoundaryNCounts;
    vector<vector<uint>> BoundaryNArr;
};

InputInfo read_input_file(const string &filename = "input")
{
    InputInfo info;
    ifstream in(filename);
    if (!in)
    {
        throw runtime_error("cannot open input file: " + filename);
    }

    string line;

    // 跳过空行和以 '%' 开头的注释行
    auto skip_comment_and_empty = [&](string &l) -> bool
    {
        while (std::getline(in, l))
        {
            // 去掉前导空白
            size_t pos = l.find_first_not_of(" \t\r\n");
            if (pos == string::npos)
                continue; // 空行
            if (l[pos] == '%')
                continue; // 注释行
            l = l.substr(pos);
            return true;
        }
        return false; // EOF
    };

    // 1) #scale
    if (!skip_comment_and_empty(line))
        throw runtime_error("unexpected EOF when reading #scale");
    if (line != "#scale")
        throw runtime_error("expect #scale, got: " + line);

    if (!skip_comment_and_empty(line))
        throw runtime_error("unexpected EOF when reading scale value");
    {
        istringstream iss(line);
        iss >> info.len_scale;
    }

    // 2) #in_filename
    if (!skip_comment_and_empty(line))
        throw runtime_error("unexpected EOF when reading #in_filename");
    if (line != "#in_filename")
        throw runtime_error("expect #in_filename, got: " + line);

    if (!skip_comment_and_empty(line))
        throw runtime_error("unexpected EOF when reading in_filename");
    info.in_filename = line; // 整行作为文件名

    // 3) #out_filename
    if (!skip_comment_and_empty(line))
        throw runtime_error("unexpected EOF when reading #out_filename");
    if (line != "#out_filename")
        throw runtime_error("expect #out_filename, got: " + line);

    if (!skip_comment_and_empty(line))
        throw runtime_error("unexpected EOF when reading out_filename");
    info.out_filename = line;

    // 4) #BoundaryNum
    if (!skip_comment_and_empty(line))
        throw runtime_error("unexpected EOF when reading #BoundaryNum");
    if (line != "#BoundaryNum")
        throw runtime_error("expect #BoundaryNum, got: " + line);

    if (!skip_comment_and_empty(line))
        throw runtime_error("unexpected EOF when reading BoundaryNum");
    {
        istringstream iss(line);
        iss >> info.BoundaryNum;
    }
    if (info.BoundaryNum < 0)
    {
        throw runtime_error("BoundaryNum must be non-negative: " +
                            std::to_string(info.BoundaryNum));
    }

    // 分配数组
    info.BoundaryNCounts.assign(info.BoundaryNum, 0);
    info.BoundaryNArr.assign(info.BoundaryNum, vector<uint>());

    // 5) #Boudary_Imformation（注意拼写就是文件里的）
    if (!skip_comment_and_empty(line))
        throw runtime_error("unexpected EOF when reading #Boudary_Imformation");
    if (line != "#Boudary_Imformation")
        throw runtime_error("expect #Boudary_Imformation, got: " + line);

    // 6) 逐个读取 Boundary1, Boundary2, ... BoundaryN
    for (int b = 0; b < info.BoundaryNum; ++b)
    {
        // 6.1 读取 "#BoundaryX"
        if (!skip_comment_and_empty(line))
            throw runtime_error("unexpected EOF when reading #Boundary" + to_string(b + 1));
        // 简单检查一下前缀
        if (line.rfind("#Boundary", 0) != 0)
        {
            throw runtime_error("expect #BoundaryX, got: " + line);
        }

        // 6.2 读取该边界下面的面数量
        if (!skip_comment_and_empty(line))
            throw runtime_error("unexpected EOF when reading face count for Boundary" + to_string(b + 1));
        {
            istringstream iss(line);
            if (!(iss >> info.BoundaryNCounts[b]))
            {
                throw runtime_error("failed to parse Boundary face count from line: " + line);
            }
        }
        // 如果该边界的面数量为0，则其arr数据也为空，直接跳到下一个边界
        if (info.BoundaryNCounts[b] == 0)
        {
            // 需要消耗掉 #BoundaryX_arr 和它下面的一行（通常是0）
            if (!skip_comment_and_empty(line))
                throw runtime_error("unexpected EOF when reading #Boundary" + to_string(b + 1) + "_arr for empty boundary");
            if (!skip_comment_and_empty(line)) // 消耗掉那个 0
                throw runtime_error("unexpected EOF when reading '0' for empty Boundary" + to_string(b + 1));
            continue; // 继续下一个 for 循环
        }

        // 6.3 读取 "#BoundaryX_arr"
        if (!skip_comment_and_empty(line))
            throw runtime_error("unexpected EOF when reading #Boundary" + to_string(b + 1) + "_arr");
        if (line.rfind("#Boundary", 0) != 0 ||
            line.find("_arr") == string::npos)
        {
            throw runtime_error("expect #BoundaryX_arr, got: " + line);
        }

        // 6.4 根据已知的数量 info.BoundaryNCounts[b]，读取所有面的ID
        // 每行最多 4 个无符号整数
        // 不足 4 个用 0 补足（文件中可以直接写 0）
        uint faces_written = 0;
        uint faces_to_read = info.BoundaryNCounts[b];
        info.BoundaryNArr[b].resize(faces_to_read);

        while (faces_written < faces_to_read)
        {
            if (!skip_comment_and_empty(line))
            {
                throw runtime_error("unexpected EOF for Boundary" + to_string(b + 1) + "_arr");
            }
            istringstream iss(line); // 读一行
            for (int k = 0; k < 4; ++k)
            {
                uint val;
                if (!(iss >> val))
                {
                    throw runtime_error("format error for Boundary" + to_string(b + 1) + "_arr: line is too short.");
                }
                if (faces_written < faces_to_read)
                {
                    info.BoundaryNArr[b][faces_written] = val;
                    faces_written++;
                }
                else
                {
                    break;
                }
            }
        }
    }
    return info;
}

// the end of the line of a file:
// in windows '/n'-->'/r''/n'
// in linux '/n'-->'/n'
// when use getline() read a line data from the file
// if this file is from windows, the end of the data we read will have a char '\r'
// but there is no '\r' when file from linux
// this function is deal with the word '\r'
void deal_s(string &s1)
{
    string::size_type s_length = 0;
    s_length = s1.size();
    if (s_length > 0 && s1[s_length - 1] == '\r')
        s1 = s1.substr(0, s_length - 1);
}

// 读取坐标点函数
vector<array<double, 3>> read_nodes(ifstream &in, uint n_nodes, double scale)
{
    vector<array<double, 3>> nodes;
    nodes.reserve(n_nodes);
    string line;
    for (uint i = 0; i < n_nodes; ++i)
    {
        array<double, 3> p{};
        // 读一行里的三个 double
        if (!getline(in, line))
        {
            throw runtime_error("read_nodes: failed to read 3 doubles for node " + to_string(i));
        }
        istringstream iss(line);
        if (!(iss >> p[0] >> p[1] >> p[2]))
        {
            throw runtime_error("read_nodes: failed to read 3 doubles for node " + to_string(i) + ", line: " + line);
        }
        p[0] *= scale;
        p[1] *= scale;
        p[2] *= scale;
        nodes.push_back(p);
    }
    // 把最后一个空行丢掉
    string dummy;
    getline(in, dummy);
    return nodes;
}

// 判断分类读取类型
void change_state(string &type_name, ReadState &type_state)
{
    // std :: cout << type_name << ": " << endl;
    if (type_name.find("vtx") != string::npos)
        type_state = ReadState::ReadVtx;
    else if (type_name.find("quad") != string::npos)
        type_state = ReadState::ReadQuad;
    else if (type_name.find("edg") != string::npos)
        type_state = ReadState::ReadEdg;
    else if (type_name.find("hex") != string::npos)
        type_state = ReadState::ReadHex;
    else if (type_name.find("tri") != string::npos)
        type_state = ReadState::ReadTri;
    else if (type_name.find("tet") != string::npos)
        type_state = ReadState::ReadTet;
    else if (type_name.find("pyr") != string::npos)
        type_state = ReadState::ReadPyr;
    else if (type_name.find("prism") != string::npos)
        type_state = ReadState::ReadPrism;
}

// 读取vtx
vector<uint> read_vtx(ifstream &in, uint n_vtx)
{
    string line;
    // 1) 读 "# Elements" 这一行
    getline(in, line);
    // 2) 读下面 n_vtx 行的元素值
    vector<uint> elements;
    elements.reserve(n_vtx);
    for (uint i = 0; i < n_vtx; ++i)
    {
        uint val;
        getline(in, line);
        istringstream iss(line);
        if (!(iss >> val))
        {
            throw runtime_error("read_nodes: failed to read  " + to_string(i) + ", line: " + line);
        }
        elements.push_back(val);
    }
    string dummy; // 把一个空行丢掉
    getline(in, dummy);
    // 3) 读 索引
    getline(in, dummy); // 跳掉<n_vtx> # number of geometric entity indices
    getline(in, dummy); // 跳掉# Geometric entity indices
    vector<uint> indices;
    indices.reserve(n_vtx);
    for (uint i = 0; i < n_vtx; ++i)
    {
        uint idx;
        getline(in, line);
        istringstream iss(line);
        if (!(iss >> idx))
        {
            throw runtime_error("read_nodes: failed to read  " + to_string(i) + ", line: " + line);
        }
        indices.push_back(idx);
    }
    // 4) 用 indices 对 elements 做重排
    // 约定：indices[k] 是原始 elements 中的下标，重排后：
    //   reordered[k] = elements[ indices[k] ]
    vector<uint> reordered;
    reordered.resize(elements.size());
    for (uint k = 0; k < indices.size(); ++k)
    {
        uint pos = indices[k];
        if (pos >= elements.size())
        {
            throw runtime_error("index out of range in geometric entity indices");
        }
        reordered[k] = elements[pos];
    }
    return reordered;
}

// 读 "# Elements" 下面的元素-
template <size_t K>
vector<array<uint, K>> read_elements_block(ifstream &in, uint n_elem)
{
    vector<array<uint, K>> elems;
    elems.reserve(n_elem);
    string line;
    // 1) 读 "# Elements" 这一行，丢掉内容
    getline(in, line);
    // 2) 读下面 n_elem 行，每行 K 个整数
    for (uint i = 0; i < n_elem; ++i)
    {
        array<uint, K> a{};
        if (!getline(in, line))
        {
            throw runtime_error("read_nodes: failed to read 3 doubles for node " + to_string(i));
        }
        istringstream iss(line);
        for (size_t j = 0; j < K; ++j)
        {
            if (!(iss >> a[j]))
            {
                throw runtime_error("read_nodes: failed to read 3 doubles for node " + to_string(i) + ", line: " + line);
            }
            a[j] += 1;
        }
        elems.push_back(a);
    }
    // 3) 吃掉最后一个数字行后面的换行
    getline(in, line);
    return elems;
}

// 读 "xxx # number of geometric entity indices" + "# Geometric entity indices" + n 行索引
vector<uint> read_geometric_indices_block(ifstream &in, uint n_elem, bool is_solid = true)
{
    vector<uint> indices;
    indices.reserve(n_elem);
    string line;
    // 1) 读 "XXX # number of geometric entity indices" 这一行
    if (!getline(in, line))
    {
        return indices;
    }
    // 2) 读 "# Geometric entity indices" 这一行
    getline(in, line);
    // 3) 读下面 n_elem 行，每行一个整数
    for (uint i = 0; i < n_elem; ++i)
    {
        uint idx = 0;
        if (!(in >> idx))
        {
            return indices;
        }
        // is_soild 为mphtxt当中域量从1开始,而边和面从0开始,需要加1
        if (is_solid)
            indices.push_back(idx);
        else
            indices.push_back(idx + 1);
    }
    // 4) 吃掉最后一个数字行后面的换行
    getline(in, line);
    return indices;
}

int main()
{
    cout << "start read..." << endl;
    // 面\体\边的总数
    uint NumOfFace = 0, NumOfSolid = 0, NumOfEdge = 0;
    // 存储容器
    vector<array<double, 3>> coordinates;
    vector<uint> inf_of_vtx;
    vector<array<uint, NUM_NODE_EDG>> inf_of_edg;
    vector<array<uint, NUM_NODE_TRI>> inf_of_tri;
    vector<array<uint, NUM_NODE_QUAD>> inf_of_quad;
    vector<array<uint, NUM_NODE_TET>> inf_of_tet;
    vector<array<uint, NUM_NODE_PYR>> inf_of_pyr;
    vector<array<uint, NUM_NODE_PRISM>> inf_of_prism;
    vector<array<uint, NUM_NODE_HEX>> inf_of_hex;

    vector<uint> entity_idcs_edg;
    vector<uint> entity_idcs_tri;
    vector<uint> entity_idcs_quad;
    vector<uint> entity_idcs_tet;
    vector<uint> entity_idcs_pyr;
    vector<uint> entity_idcs_prism;
    vector<uint> entity_idcs_hex;

    NumOfModel Model_in;
    string str_line;
    InputInfo InputInfo_data = read_input_file("input");

    ifstream in;
    in.open(InputInfo_data.in_filename);

    ReadState read_state = ReadState::Standby;
    if (in.fail() || (!in.is_open()))
        cout << "error(main()): Fail to open input file." << endl;
    else
    {
        while (getline(in, str_line))
        {
            deal_s(str_line);
            switch (read_state)
            {
            case ReadState::Standby:
            {
                if (str_line.find(STR_NUM_VTX) != string::npos)
                {
                    read_state = ReadState::ReadCoordinates;
                    istringstream iss(str_line);
                    iss >> Model_in.n_nodes; // 读取坐标点数量
                    // cout<<"start read Standby. Nodes : "<<Model_in.n_nodes<<endl;
                }
                break;
            }
            case ReadState::ReadCoordinates:
            {
                if (str_line == STR_COORDS)
                {

                    coordinates = read_nodes(in, Model_in.n_nodes, InputInfo_data.len_scale);

                    if (!getline(in, str_line))
                    {
                        throw runtime_error("unexpected EOF when reading number of element types");
                    }
                    istringstream iss(str_line);
                    if (!(iss >> Model_in.n_element_types))
                    {
                        throw runtime_error("failed to parse \"number of element types\" from line: " + str_line);
                    }
                    read_state = ReadState::CheckType;
                }
                break;
            }
            case ReadState::CheckType:
            {
                if (str_line.find("# type name") != string::npos)
                    change_state(str_line, read_state);
                break;
            }
            case ReadState::ReadVtx:
            {
                if (str_line == STR_VTX)
                {
                    cout << "*************" << endl;
                    getline(in, str_line);
                    istringstream iss(str_line);
                    if (!(iss >> Model_in.n_vtx))
                    {
                        throw runtime_error("failed to parse \"number of elements\" from line: " + str_line);
                    }
                    inf_of_vtx = read_vtx(in, Model_in.n_vtx);
                    read_state = ReadState::CheckType;
                }
                break;
            }
            case ReadState::ReadEdg:
            {
                if (str_line == STR_EDG)
                {
                    getline(in, str_line);
                    istringstream iss(str_line);
                    if (!(iss >> Model_in.n_edg))
                    {
                        throw runtime_error("failed to parse \"number of elements\" from line: " + str_line);
                    }
                    inf_of_edg = read_elements_block<NUM_NODE_EDG>(in, Model_in.n_edg);
                    entity_idcs_edg = read_geometric_indices_block(in, Model_in.n_edg, false);
                    read_state = ReadState::CheckType;
                }
                break;
            }
            case ReadState::ReadTri:
            {
                if (str_line == STR_TRI)
                {
                    getline(in, str_line);
                    istringstream iss(str_line);
                    if (!(iss >> Model_in.n_tri))
                    {
                        throw runtime_error("failed to parse \"number of elements\" from line: " + str_line);
                    }
                    inf_of_tri = read_elements_block<NUM_NODE_TRI>(in, Model_in.n_tri);
                    entity_idcs_tri = read_geometric_indices_block(in, Model_in.n_tri, false);
                    read_state = ReadState::CheckType;
                }
                break;
            }
            case ReadState::ReadQuad:
            {
                if (str_line == STR_QUAD)
                {
                    getline(in, str_line);
                    istringstream iss(str_line);
                    if (!(iss >> Model_in.n_quad))
                    {
                        throw runtime_error(
                            "failed to parse \"number of elements\" from line: " + str_line);
                    }
                    inf_of_quad = read_elements_block<NUM_NODE_QUAD>(in, Model_in.n_quad);
                    entity_idcs_quad = read_geometric_indices_block(in, Model_in.n_quad, false);
                    read_state = ReadState::CheckType;
                }
                break;
            }
            case ReadState::ReadTet:
            {
                if (str_line == STR_TET)
                {
                    getline(in, str_line);
                    istringstream iss(str_line);
                    if (!(iss >> Model_in.n_tet))
                    {
                        throw runtime_error(
                            "failed to parse \"number of elements\" from line: " + str_line);
                    }
                    inf_of_tet = read_elements_block<NUM_NODE_TET>(in, Model_in.n_tet);
                    entity_idcs_tet = read_geometric_indices_block(in, Model_in.n_tet);
                    read_state = ReadState::CheckType;
                }
                break;
            }
            case ReadState::ReadPyr:
            {
                if (str_line == STR_PYR)
                {
                    getline(in, str_line);
                    istringstream iss(str_line);
                    if (!(iss >> Model_in.n_pyr))
                    {
                        throw runtime_error("failed to parse \"number of elements\" from line: " + str_line);
                    }
                    inf_of_pyr = read_elements_block<NUM_NODE_PYR>(in, Model_in.n_pyr);
                    entity_idcs_pyr = read_geometric_indices_block(in, Model_in.n_pyr);
                    read_state = ReadState::CheckType;
                }
                break;
            }
            case ReadState::ReadPrism:
            {
                if (str_line == STR_PRISM)
                {
                    getline(in, str_line);
                    istringstream iss(str_line);
                    if (!(iss >> Model_in.n_prism))
                    {
                        throw runtime_error("failed to parse \"number of elements\" from line: " + str_line);
                    }
                    inf_of_prism = read_elements_block<NUM_NODE_PRISM>(in, Model_in.n_prism);
                    entity_idcs_prism = read_geometric_indices_block(in, Model_in.n_prism);
                    read_state = ReadState::CheckType;
                }
                break;
            }
            case ReadState::ReadHex:
            {
                if (str_line == STR_HEX)
                {
                    getline(in, str_line);
                    istringstream iss(str_line);
                    if (!(iss >> Model_in.n_hex))
                    {
                        throw runtime_error("failed to parse \"number of elements\" from line: " + str_line);
                    }
                    inf_of_hex = read_elements_block<NUM_NODE_HEX>(in, Model_in.n_hex);
                    entity_idcs_hex = read_geometric_indices_block(in, Model_in.n_hex);
                    read_state = ReadState::CheckType;
                }
                break;
            }
            default:
                break;
            }
        }
        cout << "The total number of nodes is " << Model_in.n_nodes << endl;
        cout << "The number of vtx is " << Model_in.n_vtx << endl;
        cout << "The number of two nodes element at edges is " << Model_in.n_edg << endl;

        if (Model_in.n_tri > 0)
            cout << "The number of triangle elements at faces is " << Model_in.n_tri << endl;
        if (Model_in.n_quad > 0)
            cout << "The number of quadrangle elements at faces is " << Model_in.n_quad << endl;
        Model_in.n_boby_elements = Model_in.n_tet + Model_in.n_pyr +
                                   Model_in.n_prism + Model_in.n_hex;
        cout << "The total number of 3D elements is " << Model_in.n_boby_elements << endl;
        if (Model_in.n_tet > 0)
            cout << "  Number of tetrahedrons is " << Model_in.n_tet << endl;
        if (Model_in.n_pyr > 0)
            cout << "  Number of pyramids is " << Model_in.n_pyr << endl;
        if (Model_in.n_prism > 0)
            cout << "  Number of prisms is " << Model_in.n_prism << endl;
        if (Model_in.n_hex > 0)
            cout << "  Number of hexahedrons is " << Model_in.n_hex << endl;

        in.close();
        in.clear();
    }
    uint max_tri = entity_idcs_tri.empty() ? 0 : *max_element(entity_idcs_tri.begin(), entity_idcs_tri.end());
    uint max_quad = entity_idcs_quad.empty() ? 0 : *max_element(entity_idcs_quad.begin(), entity_idcs_quad.end());
    uint max_pyr = entity_idcs_pyr.empty() ? 0 : *max_element(entity_idcs_pyr.begin(), entity_idcs_pyr.end());
    uint max_tet = entity_idcs_tet.empty() ? 0 : *max_element(entity_idcs_tet.begin(), entity_idcs_tet.end());
    uint max_hex = entity_idcs_hex.empty() ? 0 : *max_element(entity_idcs_hex.begin(), entity_idcs_hex.end());
    uint max_prism = entity_idcs_prism.empty() ? 0 : *max_element(entity_idcs_prism.begin(), entity_idcs_prism.end());

    NumOfFace = max(max_tri, max_quad);
    NumOfSolid = max({max_pyr, max_tet, max_hex, max_prism});
    NumOfEdge = *max_element(entity_idcs_edg.begin(), entity_idcs_edg.end());

    cout << "The number of edges is " << NumOfEdge << endl;
    cout << "The number of faces is " << NumOfFace << endl;
    cout << "The number of entities is " << NumOfSolid << endl;

    /////////////////////////////////////write file/////////////////////////////////////

    ofstream wr_file;
    // cout << "Start write" << endl;
    wr_file.open(InputInfo_data.out_filename);

    //////////////////////////写入文件头///////////////////////////
    wr_file << "*KEYWORD\n"
            << "*TITLE\n";
    for (uint k = 0; k < NumOfSolid; k++)
    {
        wr_file << "*PART\n"
                << "LSHELL1\n";
        wr_file << "$#" << setw(8) << "pid" << setw(10) << "secid" << setw(10) << "mid"
                << setw(10) << "eosid" << setw(10) << "hgid" << setw(10) << "grav" << setw(10) << "adpopt" << setw(10) << "tmid" << '\n';
        wr_file << setw(10) << k + 1 << setw(10) << "0" << setw(10) << "0"
                << setw(10) << "0" << setw(10) << "0" << setw(10) << "0" << setw(10) << "0" << setw(10) << "0" << '\n';
        // setw(10)是针对紧接的输出，其格式作设置，如果字符宽度小于10则用空格补齐，大于等于10则不需要空格补齐，对于紧接的输出后面的数据格式没有影响
    }

    ////////////////////////写入边界条件信息////////////////////////
    for (uint k = 0; k < InputInfo_data.BoundaryNum; k++)
    {
        wr_file << "*SET_NODE_LIST\n";
        wr_file << "$#" << setw(8) << "sid" << setw(10) << "da1" << setw(10) << "da2"
                << setw(10) << "da3" << setw(10) << "da4" << setw(10) << "solver" << '\n';

        wr_file << setw(10) << k + 1 << setw(10) << "0.0" << setw(10) << "0.0"
                << setw(10) << "0.0" << setw(10) << "0.0" << "MECH" << '\n';

        wr_file << "$#" << setw(8) << "nid1" << setw(10) << "nid2" << setw(10) << "nid3"
                << setw(10) << "nid4" << setw(10) << "nid5" << setw(10) << "nid6" << setw(10) << "nid7" << setw(10) << "nid8" << '\n';
        for (int ii = 0; ii < InputInfo_data.BoundaryNCounts[k]; ii++)
        {
            // 先写三角形
            for (int kk = 0; kk < Model_in.n_tri; kk++)
            {
                if (InputInfo_data.BoundaryNArr[k][ii] == entity_idcs_tri[kk])
                {
                    wr_file << setw(10) << inf_of_tri[kk][0];
                    wr_file << setw(10) << inf_of_tri[kk][1];
                    wr_file << setw(10) << inf_of_tri[kk][2];
                    wr_file << setw(10) << inf_of_tri[kk][2];
                    for (int kk1 = 0; kk1 < 4; kk1++)
                        wr_file << setw(10) << 0;
                    wr_file << '\n';
                }
            }

            // 再写四边形
            for (int kk = 0; kk < Model_in.n_quad; kk++)
            {
                if (InputInfo_data.BoundaryNArr[k][ii] == entity_idcs_quad[kk])
                {
                    wr_file << setw(10) << inf_of_quad[kk][0];
                    wr_file << setw(10) << inf_of_quad[kk][2];
                    wr_file << setw(10) << inf_of_quad[kk][3];
                    wr_file << setw(10) << inf_of_quad[kk][1];
                    for (int kk1 = 0; kk1 < 4; kk1++)
                        wr_file << setw(10) << 0;
                    wr_file << '\n';
                }
            }
        }
    }

    // 这里写入模型几何面上的面单元，k代表第几个面，每个*SET_SEGMENT下面的数据为共面的四边行单元节点标号
    ////////////////////////////write information of face/////////////////////////////////////////////
    for (uint k = 1; k < NumOfFace + 1; k++)
    {
        wr_file << "*SET_SEGMENT\n";
        wr_file << setw(10) << k << setw(10) << "0.0" << setw(10) << "0.0"
                << setw(10) << "0.0" << setw(10) << "0.0" << "MECH" << '\n';
        // 先写三角形
        for (uint kk = 0; kk < Model_in.n_tri; kk++)
        {
            // 遍历所有几何面上的三角形单元
            array<uint, MAX_NUM_INLINE> map = {0, 1, 2, 2, 2, 2, 2, 2};
            if (k == entity_idcs_tri[kk])
            {
                for (int kk1 = 0; kk1 < MAX_NUM_INLINE; kk1++)
                {
                    uint kkk = map[kk1];
                    wr_file << setw(10) << inf_of_tri[kk][kkk];
                }
                wr_file << '\n';
            }
        }
        // 再写四边形
        for (uint kk = 0; kk < Model_in.n_quad; kk++)
        {
            // 遍历所有几何面上的四边形单元
            array<uint, MAX_NUM_INLINE> map = {0, 2, 3, 1, 3, 3, 3, 3};
            if (k == entity_idcs_quad[kk])
            {
                for (int kk1 = 0; kk1 < MAX_NUM_INLINE; kk1++)
                {
                    uint kkk = map[kk1];
                    wr_file << setw(10) << inf_of_quad[kk][kkk];
                }
                wr_file << '\n';
            }
        }
    }

    // 这里写入线单元
    for (int k = 1; k < NumOfEdge + 1; k++)
    {
        wr_file << "*SET_SEGMENT\n";
        wr_file << std::setw(10) << k + NumOfFace << std::setw(10) << "0.0" << std::setw(10) << "0.0"
                << std::setw(10) << "0.0" << std::setw(10) << "0.0" << "MECH" << '\n';
        // std::cout<<NumOfFace<<" "<<k<<std::endl;;
        for (unsigned long ii = 0; ii < Model_in.n_edg; ii++)
        {
            if (entity_idcs_edg[ii] == k)
            {
                // std::cout<<ii<<"  k"<<k<<std::endl;
                wr_file << std::setw(10) << inf_of_edg[ii][0];
                wr_file << std::setw(10) << inf_of_edg[ii][1];
                wr_file << std::setw(10) << inf_of_edg[ii][1];
                wr_file << std::setw(10) << inf_of_edg[ii][1];
                for (int kk = 0; kk < 4; kk++)
                    wr_file << std::setw(10) << "0.000";
                wr_file << '\n';
            }
        }
    }

    ///////////////////////////////write information of node//////////////////////////////////////////
    wr_file << "*NODE\n";
    wr_file << "$#" << setw(8) << "nid" << setw(8) << "x" << setw(8) << "y" << setw(10) << "z" << setw(10) << "tc"
            << setw(10) << "rc" << '\n';

    for (uint k = 0; k < Model_in.n_nodes; k++)
    {
        wr_file << setw(8) << k + 1;
        for (int kk = 0; kk < 3; kk++)
            wr_file << setw(16) << setprecision(10) << coordinates[k][kk];
        wr_file << setw(8) << "0" << setw(8) << "0" << '\n';
    }

    /////////////////////////write information of element_solid///////////////////////////////////////
    wr_file << "*ELEMENT_SOLID\n";
    wr_file << "$#" << setw(8) << "eid" << setw(8) << "pid" << setw(8) << "n1" << setw(8) << "n2"
            << setw(8) << "n3" << setw(8) << "n4" << setw(8) << "n5" << setw(8) << "n6" << setw(8) << "n7" << setw(8) << "n8" << '\n';
    for (uint k = 0; k < Model_in.n_tet; k++)
    {
        // 遍历所有的四面体单元
        array<uint, MAX_NUM_INLINE> map = {0, 1, 2, 3, 3, 3, 3, 3};
        wr_file << setw(8) << k + 1;
        wr_file << setw(8) << entity_idcs_tet[k];
        for (uint kk = 0; kk < MAX_NUM_INLINE; kk++)
        {
            uint kkk = map[kk];
            wr_file << setw(8) << inf_of_tet[k][kkk];
        }
        wr_file << '\n';
    }
    for (uint k = 0; k < Model_in.n_pyr; k++)
    {
        // 遍历所有的金字塔单元
        array<uint, MAX_NUM_INLINE> map = {0, 1, 3, 2, 4, 4, 4, 4};
        wr_file << setw(8) << k + 1 + Model_in.n_tet;
        wr_file << setw(8) << entity_idcs_pyr[k];
        for (uint kk = 0; kk < MAX_NUM_INLINE; kk++)
        {
            uint kkk = map[kk];
            wr_file << setw(8) << inf_of_pyr[k][kkk];
        }
        wr_file << '\n';
    }
    for (uint k = 0; k < Model_in.n_prism; k++)
    {
        // 遍历所有的三棱柱单元
        array<uint, MAX_NUM_INLINE> map = {0, 1, 2, 2, 3, 4, 5, 5};
        wr_file << setw(8) << k + 1 + Model_in.n_tet + Model_in.n_pyr;
        wr_file << setw(8) << entity_idcs_prism[k];
        for (uint kk = 0; kk < MAX_NUM_INLINE; kk++)
        {
            uint kkk = map[kk];
            wr_file << setw(8) << inf_of_prism[k][kkk];
        }
        wr_file << '\n';
    }
    for (uint k = 0; k < Model_in.n_hex; k++)
    {
        // 遍历所有的六面体单元
        array<uint, MAX_NUM_INLINE> map = {0, 2, 3, 1, 4, 6, 7, 5};
        wr_file << setw(8) << k + 1 + Model_in.n_tet + Model_in.n_pyr + Model_in.n_prism;
        wr_file << setw(8) << entity_idcs_hex[k];
        for (uint kk = 0; kk < MAX_NUM_INLINE; kk++)
        {
            uint kkk = map[kk];
            wr_file << setw(8) << inf_of_hex[k][kkk];
        }
        wr_file << '\n';
    }

    wr_file << "*END";
    wr_file.close();
    cout << "The output file is already written." << endl;

    return 0;
}
