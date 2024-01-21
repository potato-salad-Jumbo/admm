#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
// #include<winnt.h>
// #include <windows.h>
#define INF 1.0e13
#define ZERO 1.0e-13
#define INVALID -1      /* Represents an invalid value. */
#define WAS_IN_QUEUE -7 /* Shows that the node was in the queue before. (7 is for luck.) */

using namespace std;
using namespace std::chrono;


class CNode //节点
{
public:
  int ID;
  string Name;
  int Origin_ID = -1;
  vector<int> IncomingLink;
  vector<int> OutgoingLink;
};

class CLink // LinkID --> Fft，Capacity, BlockID, Alpha, Power
{
public:
  int ID;
  CNode *pInNode;
  CNode *pOutNode;
  double FreeFlowTravelTime;
  double Capacity;
  double Alpha = 0.15;
  double Power = 4.0;
  int BlockID;
  double OptimalSolution;
};

class CODPair // ODPairID-->
{
public:
  int ori, des, desID, oriID;
  double ODdemand;
};

class COrigin
{
public:
  int ID;
  CNode *pOriginNode;
  vector<int> DestinationNode;
  vector<double> ODDemand;
  vector<int> ODPairIndex;
};


void ReadData(string fileName);
void ADMM();
void NetworkInitialization();
void ADMM_Initialization();
void ADMM_DualUpdate();
void ADMM_PrimalUpdate();
void ADMM_PrimalUpdateStage2();
double nodelinksum();
void linkSubproblem_v3(int lId);
void linkSubproblem_v4(int lId);
void linkSubproblem_dummy(int lId);
double GetUEGap();
void GetShortestPath(int orinode, double *ShortestPathCost, int *ShortestPathParent);


vector<CNode> m_Node;
vector<CLink> m_Link;
vector<int> *BlockSet;
int max_BlockId;
// vector<int> OIdSet;
vector<COrigin> m_Origin;
double **m_ONDem;
double **m_ONDual;
double **m_OLPrimal;
double **pre_Dual;
///////////////////////////
double *LinkFlow;
double *LinkCost;
double *LinkCostDiff;
double *LinkFreeTravelTime;
double *LinkCostCoef;
double *LinkCostDiffCoef;
double *LinkOptimalSolution;
double *prelinkflow;

///////////////////////////
double penalty = 0.08; // Small:1       SF:0.008;0.005会更快,但是存在强烈震荡;     An:0.026/全部OD时0.08,耗时7s左右，0.4耗时1秒多；
                       // CS：设置1，耗时3s多，不需要设置内循环，而AN网络需要设置内循环
double rg = 10;
int counter;
//***SP****************************************
int FIRST_THRU_NODE;
double step;
double curobjective;
double preobjective;
double curgradientdirection;
double pregradientdirection;
///////////////////////////
double p_dual = 0, p_primal = 0, p_rg = 0, s_dual, e_dual, s_primal, e_primal, s_rg, e_rg, s_primal_eff, e_primal_eff, max_period, p_primal_eff = 0, pprg = 0;
/////////////////////////////


void ReadData(string fileName)
{
  // cout  << fileName << endl;
  int count;
  string lineStr;   //读取每行的字符串
  string inlineStr; //每行断开
  stringstream stream;
  int intData;
  double doubleData;
  // read node_set//////////////////////////////////
  ifstream nodeFile("ADMM/data/" + fileName + "_node.csv");
  if (nodeFile.is_open())
  {
    count = 0;
    while (getline(nodeFile, lineStr)) //当没有读取到文件末尾时循环继续
    {
      CNode Node;
      stringstream ss(lineStr);
      getline(ss, inlineStr, '\t');
      Node.Name = inlineStr; // Node的name是连续的，从1开始，ID从0开始，间隔-1；
      Node.ID = count;
      count++;
      m_Node.push_back(Node);
    }
    nodeFile.close();
    // printf("Read node data finish!\n");
    // test
    // for (int i = 0; i < m_Node.size(); i++)
    // {
    //     cout  << m_Node[i].ID << ','<< m_Node[i].Name<< endl;
    // }
  }
  else
  {
    printf("Read node data wrong!\n");
  }

  // read link_set//////////////////////////////////
  ifstream linkFile("data/" + fileName + "_link.csv");
  count = 0;
  if (linkFile.is_open())
  {
    while (getline(linkFile, lineStr))
    {
      CLink Link;
      stringstream ss(lineStr);
      Link.ID = count;
      count++;

      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> intData;
      m_Node[intData - 1].OutgoingLink.push_back(Link.ID);
      Link.pInNode = &m_Node[intData - 1];
      stream.clear();

      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> intData;
      m_Node[intData - 1].IncomingLink.push_back(Link.ID);
      Link.pOutNode = &m_Node[intData - 1];
      stream.clear();

      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> doubleData;
      Link.FreeFlowTravelTime = doubleData;
      stream.clear();

      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> doubleData;
      Link.Capacity = doubleData;
      stream.clear();

      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> intData;
      Link.BlockID = intData - 1; //要保证Block Name从1开始依次累加，与Block ID有固定的-1关系
      stream.clear();

      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> doubleData;
      Link.OptimalSolution = doubleData;
      stream.clear();

      m_Link.push_back(Link);
    }
    linkFile.close();
    // printf("Read link data finish!\n");
    // for (int i=0; i< m_Link.size();i++)
    //     cout<<m_Link[i].ID<<", "<<m_Link[i].pInNode->Name<<", "<<m_Link[i].pOutNode->Name<<", "<<m_Link[i].BlockID<<endl;
  }
  else
  {
    printf("Read link data wrong!\n");
  }
  ////////////////////////////////////////////////////////////////////
  // for (int i = 0; i < m_Node.size(); i++)
  // {
  //     cout <<i<<", " << m_Node[i].ID << ',' << m_Node[i].Name<< endl;
  //     for (int j = 0; j < m_Node[i].IncomingLink.size(); j++)
  //         printf("  %d\n",m_Node[i].IncomingLink[j]);
  //     // for (int j = 0; j < m_Node[i].OutgoingLink.size(); j++)
  //     //     printf("  %d\n",m_Node[i].OutgoingLink[j]);
  // }

  // establish Block_set//////////////////////////////////
  count = 0;
  max_BlockId = 0;
  for (int i = 0; i < m_Link.size(); i++)
    if (m_Link[i].BlockID > max_BlockId)
      max_BlockId = m_Link[i].BlockID;
  // cout<<max_BlockId<<endl;
  BlockSet = new vector<int>[max_BlockId + 1]; //要保证Block Name与Block ID有固定的-1关系
  for (int i = 0; i < max_BlockId + 1; i++)
  {
    vector<int> LinkIdSet;
    for (int j = 0; j < m_Link.size(); j++)
      if (i == m_Link[j].BlockID)
        LinkIdSet.push_back(j);
    BlockSet[i] = LinkIdSet;
  }
  // for (int i=0; i< max_BlockId+1;i++)
  // int i=9;
  //     for (int j=0; j< BlockSet[i].size();j++)
  //         cout<<i<<","<<BlockSet[i][j]<<endl;

  // read OD demand and origin set//////////////////////////////////
  // read fileName_od.csv
  ifstream odFile("data/" + fileName + "_od.csv");
  CNode *pNode;
  COrigin *pOrigin;
  int line = -1;
  if (odFile.is_open())
  {
    while (getline(odFile, lineStr)) //当没有读取到文件末尾时循环继续
    {
      line++;
      stringstream ss(lineStr);
      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> intData;
      pNode = &m_Node[intData - 1]; //要保证nodeName与nodeID有固定的-1关系
      stream.clear();

      if (pNode->Origin_ID == -1)
      {
        // pOrigin = new COrigin();
        COrigin Origin;
        Origin.ID = m_Origin.size(); //也是要保持originID，这个很有问题。
        Origin.pOriginNode = pNode;
        pNode->Origin_ID = Origin.ID; // Origin.ID表示m_Origin中的Index，这个Index记录在了pNode中。
        m_Origin.push_back(Origin);   // node->m_ori  m_ori->node都有可以获取，双向键值对。
        pOrigin = &m_Origin[pNode->Origin_ID];
      }
      else
      {
        pOrigin = &m_Origin[pNode->Origin_ID];
      }
      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> intData;
      pNode = &m_Node[intData - 1]; //这里获取的是des的信息，也是Node的ID与Name有固定关系。
      stream.clear();

      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> doubleData;
      // cout<<doubleData<<endl;
      if (doubleData > 0)
      {
        pOrigin->DestinationNode.push_back(pNode->ID); //这里只是获得了Node的ID信息。
        pOrigin->ODDemand.push_back(doubleData);
        pOrigin->ODPairIndex.push_back(line); //这里也是为了做双向索引。
      }
      stream.clear();
    }
    odFile.close();
    // printf("Read origin data finish!\n");
  }
  else
  {
    printf("Read origin data wrong!\n");
  }

  // Read O_N based Demand table
  m_ONDem = new double *[m_Origin.size()];
  for (int oId = 0; oId < m_Origin.size(); oId++)
  {
    m_ONDem[oId] = new double[m_Node.size()];
  }
  ifstream onFile("data/" + fileName + "_on.csv");
  int oriId, desId;
  if (onFile.is_open())
  {
    while (getline(onFile, lineStr))
    {
      stringstream ss(lineStr);

      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> intData;
      oriId = intData - 1; //获取的也是NodeID，前提是默认Node ID与Name相差-1
      stream.clear();

      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> intData;
      desId = intData - 1;
      stream.clear();

      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> doubleData;
      m_ONDem[oriId][desId] = doubleData;
      stream.clear();
      // cout<<oriId<<","<<desId<<endl;
    }
    onFile.close();
    // printf("Read on data finish!\n");
    // for(int oi=0;oi<m_Origin.size();oi++)
    //     for (int ni=0;ni<m_Node.size();ni++)
    //         cout<<oi+1<<","<<ni+1<<", "<<m_ONDem[oi][ni]<<endl;
  }
  else
  {
    printf("Read on data wrong!\n");
  }
}


void NetworkInitialization()
{
  CLink *pLink;
  double temp;
  //********************************************
  LinkCostDiff = new double[m_Link.size()];
  LinkFlow = new double[m_Link.size()];
  LinkCost = new double[m_Link.size()];
  LinkFreeTravelTime = new double[m_Link.size()];
  LinkCostCoef = new double[m_Link.size()];
  LinkCostDiffCoef = new double[m_Link.size()];
  LinkOptimalSolution = new double[m_Link.size()];
  prelinkflow = new double[m_Link.size()];
  //********************************************
  for (int link = 0; link < m_Link.size(); link++)
  {
    pLink = &m_Link.at(link);
    // pLink->pInNode->OutgoingLink.push_back(link);
    // pLink->pOutNode->IncomingLink.push_back(link);
    LinkFreeTravelTime[link] = pLink->FreeFlowTravelTime;
    temp = pow(pLink->Capacity, pLink->Power); //传统的BPR函数，即4次方
    // temp = pow (pLink->Capacity, 2);                   //修改后的BPR函数，即2或1次方

    LinkCostCoef[link] = LinkFreeTravelTime[link] * pLink->Alpha / temp;
    LinkCostDiffCoef[link] = LinkFreeTravelTime[link] * pLink->Alpha * pLink->Power / temp; //传统的BPR函数，即4次方
    // LinkCostDiffCoef[link] = LinkFreeTravelTime[link] * pLink->Alpha * 2 / temp;      //修改后的BPR函数，即2或1次方
    LinkOptimalSolution[link] = pLink->OptimalSolution;
  }
  //***************************************************
}


void ADMM_Initialization()
{
  // dual_variable
  m_ONDual = new double *[m_Origin.size()];
  pre_Dual = new double *[m_Origin.size()];
  for (int oId = 0; oId < m_Origin.size(); oId++)
  {
    m_ONDual[oId] = new double[m_Node.size()];
    pre_Dual[oId] = new double[m_Node.size()];
    for (int nId = 0; nId < m_Node.size(); nId++)
    {
      m_ONDual[oId][nId] = 0.0;
      pre_Dual[oId][nId] = 0.0;
      // cout<<oId<<","<<nId<<","<<m_ONDual[oId][nId]<<endl;
    }
  }
  // primal_variable
  m_OLPrimal = new double *[m_Origin.size()];
  for (int oId = 0; oId < m_Origin.size(); oId++)
  {
    m_OLPrimal[oId] = new double[m_Link.size()];
    for (int lId = 0; lId < m_Link.size(); lId++)
    {
      m_OLPrimal[oId][lId] = 0.0;
    }
  }
  for (int link = 0; link < m_Link.size(); link++)
  {
    LinkFlow[link] = 0.0;
    prelinkflow[link] = 0.0;
  }
}

#include <chrono> // 确保包含了必要的头文件

void ADMM_DualUpdate()
{
  double oId_nId_outgoingLfCons = 0;
  double oId_nId_IncomingLfCons = 0;
  double on_conscons = 0;
  CNode *pNode;
  double max_dual = 0;

  for (int oId = 0; oId < m_Origin.size(); oId++)
  {
    for (int nId = 0; nId < m_Node.size(); nId++)
    {
      pNode = &m_Node.at(nId);
      oId_nId_outgoingLfCons = 0;
      for (int lId = 0; lId < pNode->OutgoingLink.size(); lId++)
        oId_nId_outgoingLfCons += m_OLPrimal[oId][pNode->OutgoingLink[lId]];
      oId_nId_IncomingLfCons = 0;
      for (int lId = 0; lId < pNode->IncomingLink.size(); lId++)
        oId_nId_IncomingLfCons += m_OLPrimal[oId][pNode->IncomingLink[lId]];
      on_conscons = oId_nId_outgoingLfCons - oId_nId_IncomingLfCons - m_ONDem[oId][nId];
      m_ONDual[oId][nId] += penalty * on_conscons;
    }
  }
}

void ADMM_PrimalUpdate()
{
  double max_primal = 0;
  for (int bId = 0; bId <= max_BlockId; bId++)
  {
    max_primal = 0;
    for (int lId = 0; lId < BlockSet[bId].size(); lId++)
    {
      auto s_primal_eff = std::chrono::high_resolution_clock::now(); // 开始计时
      linkSubproblem_v3(BlockSet[bId][lId]); // 执行原本的操作
      auto e_primal_eff = std::chrono::high_resolution_clock::now(); // 结束计时
      auto duration = std::chrono::duration_cast<std::chrono::microseconds>(e_primal_eff - s_primal_eff).count();
      if (max_primal < duration)
        max_primal = duration;
    }
    p_primal_eff += max_primal;
  }
}

void ADMM_PrimalUpdate_NOBlock()
{
  // s_primal = omp_get_wtime();
  for (int link = 0; link < m_Link.size(); link++)
  {
    linkSubproblem_v3(link);
  }
}


double ObjofLinkSubproblem(int lId)
{
  double objective = 0.0;
  double tail_o_succecons = 0, head_o_succecons = 0;
  double tail_o_predecons = 0, head_o_predecons = 0;
  double tail_o_cons = 0, head_o_cons = 0, conscons = 0;
  double FirstTerm = 0, SecondTerm = 0, ThirdTerm = 0;
  double A1 = 0, A2 = 0;
  CLink *pLink;
  pLink = &m_Link[lId];
  //**********计算子问题目标函数的第一项******
  FirstTerm = LinkFreeTravelTime[lId] * LinkFlow[lId] + LinkCostCoef[lId] * pow(LinkFlow[lId], 5) / 5;

  for (int oId = 0; oId < m_Origin.size(); oId++)
  {
    tail_o_succecons = 0;
    for (int lId2 = 0; lId2 < pLink->pInNode->OutgoingLink.size(); lId2++)
      tail_o_succecons += m_OLPrimal[oId][pLink->pInNode->OutgoingLink[lId2]]; //累加除待解决变量之外的origin based link flow（outgoing）
    tail_o_predecons = 0;
    for (int lId2 = 0; lId2 < pLink->pInNode->IncomingLink.size(); lId2++)
      tail_o_predecons += m_OLPrimal[oId][pLink->pInNode->IncomingLink[lId2]]; //累加除待解决变量之外的origin based link flow（incoming）
    tail_o_cons = tail_o_succecons - tail_o_predecons;
    tail_o_cons -= m_ONDem[oId][pLink->pInNode->ID];
    A1 = (m_ONDual[oId][pLink->pInNode->ID] * tail_o_cons) + (penalty / 2) * (tail_o_cons) * (tail_o_cons);
    SecondTerm += A1;
    // cout<<" pLink_InNode_o_Outgoing - pLink_InNode_o_Incoming-ONDem: "<<tail_o_cons<<endl;
    // head node=out node
    head_o_succecons = 0;
    for (int lId2 = 0; lId2 < pLink->pOutNode->OutgoingLink.size(); lId2++)
    {
      // cout<<" head_outgonging_links   "<<lId2<<", "<<pLink->pOutNode->OutgoingLink[lId2]<<": ("<< m_Link[pLink->pOutNode->OutgoingLink[lId2]].pInNode->Name <<", "<<m_Link[pLink->pOutNode->OutgoingLink[lId2]].pOutNode->Name<<")"<<endl;
      head_o_succecons += m_OLPrimal[oId][pLink->pOutNode->OutgoingLink[lId2]];
      // cout<<" -> "<<m_OLPrimal[oId][pLink->pOutNode->OutgoingLink[lId2]]<<endl;
    }
    // cout<<" pLink_OutNode_o_Outgoing Link Flow Sum WO itself: "<<head_o_succecons<<endl;
    head_o_predecons = 0;
    for (int lId2 = 0; lId2 < pLink->pOutNode->IncomingLink.size(); lId2++)
      head_o_predecons += m_OLPrimal[oId][pLink->pOutNode->IncomingLink[lId2]];
    // cout<<" pLink_OutNode_o_Incoming Link Flow Sum WO itself: "<<head_o_predecons<<endl;
    head_o_cons = head_o_succecons - head_o_predecons;
    head_o_cons -= m_ONDem[oId][pLink->pOutNode->ID];
    // cout<<" pLink_OutNode_o_Outgoing - pLink_OutNode_o_Incoming-ONDem: "<<head_o_cons<<endl;
    A2 = (m_ONDual[oId][pLink->pOutNode->ID] * head_o_cons) + (penalty / 2) * (head_o_cons) * (head_o_cons);
    ThirdTerm += A2;
  }

  objective = FirstTerm + SecondTerm + ThirdTerm;

  return objective;
}


void linkSubproblem_v3(int lId)
{ 
  CLink *pLink;
  pLink = &m_Link[lId];
  double beta1, alpha1, ared1, pred1, m1; //应用Armijo Rule确定步长的参数
  double temp_olf;
  double *update_OLPrimal_1st;
  double *update_OLPrimal_1st_cons;
  update_OLPrimal_1st = new double[m_Origin.size()];
  update_OLPrimal_1st_cons = new double[m_Origin.size()];
  double templinkFlow = LinkFlow[lId];
  double temppow = templinkFlow * templinkFlow * templinkFlow;
  double snd = LinkCostDiffCoef[lId] * temppow + 2 * penalty; //传统的BPR函数，即4次方
  double l_BPR = LinkFreeTravelTime[lId] + LinkCostCoef[lId] * temppow * templinkFlow; //传统的BPR函数，即4次方
  double o_Lp1std_linear, gap, dv;
  double tail_o_succecons = 0, head_o_succecons = 0;
  double tail_o_predecons = 0, head_o_predecons = 0;
  double tail_o_cons = 0, head_o_cons = 0, conscons = 0;
  for (int oId = 0; oId < m_Origin.size(); oId++)
  {
    o_Lp1std_linear = m_OLPrimal[oId][lId] * penalty * 2;
    tail_o_succecons = 0;
    for (int lId2 = 0; lId2 < pLink->pInNode->OutgoingLink.size(); lId2++)
      if (pLink->pInNode->OutgoingLink[lId2] != lId)                             //
        tail_o_succecons += m_OLPrimal[oId][pLink->pInNode->OutgoingLink[lId2]]; //累加除待解决变量之外的origin based link flow（outgoing）
    tail_o_predecons = 0;
    for (int lId2 = 0; lId2 < pLink->pInNode->IncomingLink.size(); lId2++)
      if (pLink->pInNode->IncomingLink[lId2] != lId)
        tail_o_predecons += m_OLPrimal[oId][pLink->pInNode->IncomingLink[lId2]]; //累加除待解决变量之外的origin based link flow（incoming）
    tail_o_cons = tail_o_succecons - tail_o_predecons;
    tail_o_cons -= m_ONDem[oId][pLink->pInNode->ID];
    head_o_succecons = 0;
    for (int lId2 = 0; lId2 < pLink->pOutNode->OutgoingLink.size(); lId2++)
    {
     if (pLink->pOutNode->OutgoingLink[lId2] != lId)
      {
        head_o_succecons += m_OLPrimal[oId][pLink->pOutNode->OutgoingLink[lId2]];;
      }
    }
    head_o_predecons = 0;
    for (int lId2 = 0; lId2 < pLink->pOutNode->IncomingLink.size(); lId2++)
      if (pLink->pOutNode->IncomingLink[lId2] != lId)
        head_o_predecons += m_OLPrimal[oId][pLink->pOutNode->IncomingLink[lId2]];
    head_o_cons = head_o_succecons - head_o_predecons;
    head_o_cons -= m_ONDem[oId][pLink->pOutNode->ID];
    conscons = tail_o_cons - head_o_cons;
    update_OLPrimal_1st_cons[oId] = penalty * conscons + m_ONDual[oId][pLink->pInNode->ID] - m_ONDual[oId][pLink->pOutNode->ID];
    update_OLPrimal_1st[oId] = l_BPR + o_Lp1std_linear + update_OLPrimal_1st_cons[oId];
    m_OLPrimal[oId][lId] = m_OLPrimal[oId][lId] - (update_OLPrimal_1st[oId]) / snd;
    if (m_OLPrimal[oId][lId] < 0)
      m_OLPrimal[oId][lId] = 0;
  }
  double terminate, snd_15, d_15, terminategap = 10, terminate_0 = 10; // termination criteria
  int subiter = 0;
  double gama = 0;
  terminate = 0;
  for (int oId = 0; oId < m_Origin.size(); oId++)
    terminate += fabs(m_OLPrimal[oId][lId] * update_OLPrimal_1st[oId]); // flow(k+1)*gradient(k)
  templinkFlow = 0;
  for (int oId = 0; oId < m_Origin.size(); oId++)
    templinkFlow += m_OLPrimal[oId][lId];
  LinkFlow[lId] = templinkFlow;
  {
    subiter = subiter + 1;
    //**************Fixed Step Size*****************
    step = 1.0;
    for (int oId = 0; oId < m_Origin.size(); oId++)
    {
      temppow = templinkFlow * templinkFlow * templinkFlow; //这里可以改成link flow[lid]
      //传统的BPR函数，即4次方
      snd = LinkCostDiffCoef[lId] * temppow + 2 * penalty;                          //由于link flow是变化的，snd也会变化，
      l_BPR = LinkFreeTravelTime[lId] + LinkCostCoef[lId] * temppow * templinkFlow; // 1st中的BPR time也会随着link flow的变化而变化。
      update_OLPrimal_1st[oId] = l_BPR + m_OLPrimal[oId][lId] * penalty * 2 + update_OLPrimal_1st_cons[oId]; // origin-based link flow会变化，常数项不变了；
      temp_olf = m_OLPrimal[oId][lId] - step * update_OLPrimal_1st[oId] / snd; ///马上update，接着做投影；
      if (temp_olf < 0)
        temp_olf = 0;
      templinkFlow -= m_OLPrimal[oId][lId];
      templinkFlow += temp_olf; //这里只更新一个，很聪明的做法；
      m_OLPrimal[oId][lId] = temp_olf;
    }
    terminate = 0;
    for (int oId = 0; oId < m_Origin.size(); oId++)
      terminate += fabs(update_OLPrimal_1st[oId] * m_OLPrimal[oId][lId]);
    gap = 0;
    for (int oId = 0; oId < m_Origin.size(); oId++)
      gap += update_OLPrimal_1st[oId];

    terminategap = fabs(terminate - terminate_0);
    terminate_0 = terminate;
    LinkFlow[lId] = templinkFlow;
  }
  LinkFlow[lId] = templinkFlow;
  delete[] update_OLPrimal_1st;
  delete[] update_OLPrimal_1st_cons;
}


void GetShortestPath(int orinode, double *ShortestPathCost, int *ShortestPathParent)
{
  CNode *pNode;
  CLink *pLink;
  int m_nNode = m_Node.size();
  int now, NewNode;
  double NewCost;
  int QueueFirst, QueueLast;
  int *QueueNext = new int[m_nNode](); //存需要挨个检查的节点列表，每个元素是下一个需要检查的节点编号，或者注明本节点加入过队伍
  for (int node = 0; node < m_nNode; node++)
  {
    QueueNext[node] = INVALID;
    ShortestPathCost[node] = INF;
    ShortestPathParent[node] = INVALID;
  }
  now = orinode;
  QueueNext[now] = WAS_IN_QUEUE;
  ShortestPathParent[now] = INVALID;
  ShortestPathCost[now] = 0.0;
  QueueFirst = QueueLast = INVALID;
  while ((now != INVALID) && (now != WAS_IN_QUEUE))
  {
    pNode = &m_Node[now]; //对应于网络中的点
    if (now >= FIRST_THRU_NODE - 1 || now == orinode)
    {
      for (int index = 0; index < pNode->OutgoingLink.size(); index++)
      {
        /* For every link that terminate at "now": */
        pLink = &m_Link[pNode->OutgoingLink[index]]; //找对应node点与index点之间的路段编号
        NewNode = pLink->pOutNode->ID;
        NewCost = ShortestPathCost[now] + LinkCost[pLink->ID];
        if (ShortestPathCost[NewNode] > NewCost)
        {
          /* If the new lable is better than the old one, correct it, and make sure that the new node to the queue. */
          ShortestPathCost[NewNode] = NewCost;
          ShortestPathParent[NewNode] = pLink->ID;
          /* If the new node was in the queue before, add it as the first in the queue. */
          if (QueueNext[NewNode] == WAS_IN_QUEUE)
          {
            QueueNext[NewNode] = QueueFirst;
            QueueFirst = NewNode;
            if (QueueLast == INVALID)
              QueueLast = NewNode;
          }
          /* If the new node is not in the queue, and wasn't there before, add it at the end of the queue. */
          else if (QueueNext[NewNode] == INVALID && NewNode != QueueLast)
          {
            if (QueueLast != INVALID)
            {
              /*Usually*/
              QueueNext[QueueLast] = NewNode;
              QueueLast = NewNode;
            }
            else
            {
              /* If the queue is empty, initialize it. */
              QueueFirst = QueueLast = NewNode;
              QueueNext[QueueLast] = INVALID;
            }
          }
          /* If the new node is in the queue, just leave it there. (Do nothing) */
        }
      }
    }
    /* Get the first node out of the queue, and use it as the current node. */
    now = QueueFirst;
    if ((now == INVALID) || (now == WAS_IN_QUEUE))
      break;
    QueueFirst = QueueNext[now];
    QueueNext[now] = WAS_IN_QUEUE;
    if (QueueLast == now)
      QueueLast = INVALID;
  }
  delete[] QueueNext;
}


double GetUEGap()
{
  // double start1 = omp_get_wtime();
  double link_systemcost = 0.0, od_systemcost = 0.0, temppow;
  int ori, link, des;
  for (link = 0; link < m_Link.size(); link++)
  {
    double templf = LinkFlow[link];
    temppow = templf * templf * templf;
    LinkCost[link] = LinkFreeTravelTime[link] + LinkCostCoef[link] * temppow * templf; //传统的BPR函数，即4次方
    // LinkCost[link]=LinkFreeTravelTime[link] + LinkCostCoef[link]*templf*templf;     //传统的BPR函数，即2次方
    // LinkCost[link]=LinkFreeTravelTime[link] + LinkCostCoef[link]*templf;     //传统的BPR函数，即1次方
    link_systemcost += LinkFlow[link] * LinkCost[link];
  }
  // double end1 = omp_get_wtime();
  // pprg = (end1 - start1);
  double max = 0.0;
  for (int ori = 0; ori < m_Origin.size(); ori++)
  {
    s_rg = omp_get_wtime();
    double *ShortestPathCost = new double[m_Node.size()];
    int *ShortestPathParent = new int[m_Node.size()];
    COrigin *pOrigin = &m_Origin.at(ori);
    GetShortestPath(pOrigin->pOriginNode->ID, ShortestPathCost, ShortestPathParent);
    for (int des = 0; des < pOrigin->DestinationNode.size(); des++)
      od_systemcost += pOrigin->ODDemand.at(des) * ShortestPathCost[pOrigin->DestinationNode.at(des)];
    delete[] ShortestPathParent;
    delete[] ShortestPathCost;
    e_rg = omp_get_wtime();
    if (ori == 0)
      max = (e_rg - s_rg);
    // if (counter == 40)
    // cout << "Ori, " << ori << " , CPU Time, " << setprecision(20)   << e_rg-s_rg << endl;
  }
  // cout << endl;
  p_rg += (max + pprg);
  return 1 - od_systemcost / link_systemcost;
}


double ObjectiveValue()
{
  double obj = 0.0;
  double link_integral = 0.0;
  CLink *pLink;
  for (int link = 0; link < m_Link.size(); link++)
  {
    pLink = &m_Link.at(link);
    link_integral = LinkFreeTravelTime[link] * LinkFlow[link] + LinkCostCoef[link] * pow(LinkFlow[link], 5) / 5;
    // link_integral = LinkFreeTravelTime[link] * LinkFlow[link] + LinkCostCoef[link] * pow(LinkFlow[link], 3) / 3;       //修改后的BPR函数，即2次方
    // link_integral = LinkFreeTravelTime[link] * LinkFlow[link] + LinkCostCoef[link] * pow(LinkFlow[link], 2) / 2;       //修改后的BPR函数，即1次方
    obj += link_integral;
  }
  return obj;
}


void ADMM()
{
  counter = 0;
  int InnerLoop = 1;
  CLink *pLink;
  NetworkInitialization();
  ADMM_Initialization();
  double inter;
  double ECT = 0;
  double prerg = 1000;
  // double start= clock();
  double start = omp_get_wtime();
  // ofstream outfile;
  // outfile.open("D:/!!!ParallelComputingGroup/Code/ADMM/Result_An.csv");
  // outfile << "Counter\tInter-Start\tSum \tRg" << endl;
  while (rg > 1.0e-5) // rg > 1.0e-10  // rg>1e-6  //counter<2066;
  // while (counter < 46512)
  {
    // s_dual = omp_get_wtime();
    // s_dual = clock();
    ADMM_DualUpdate();
    // e_dual = omp_get_wtime();
    // e_dual = clock();
    // p_dual += (e_dual-s_dual);
    // double s_primal = clock();
    ADMM_PrimalUpdate();
    // double e_primal = clock();
    // p_primal += (e_primal - s_primal)/CLOCKS_PER_SEC;
    // ADMM_PrimalUpdate_NOBlock();
    // for (int ori =0; ori < m_Origin.size(); ori++)
    // {
    //     for (int link =0; link < m_Link.size(); link++)
    //     {
    //         for (int block =0; block < BlockSet[1].size(); block++)
    //         {
    //             pLink = &m_Link[BlockSet[1][block]];
    //             if (link == pLink->ID)
    //             {
    //                 printf("1. Primal[%d][%d]:%f \n", ori, link, m_OLPrimal[ori][link]);
    //             }
    //         }
    //     }
    //     for (int node =0; node < m_Node.size(); node++)
    //     {
    //         printf("Dual[%d][%d]:%f \n", ori, node, m_ONDual[ori][node]);
    //     }
    // }
    // s_rg = omp_get_wtime();
    // s_rg = clock();
    rg = GetUEGap();
    // e_rg = omp_get_wtime();
    // e_rg = clock();
    // p_rg += (e_rg - s_rg);
    counter += 1;
    // cout<<counter<<", RG is: "<<rg<<endl;
    // inter=omp_get_wtime();
    // double obj = ObjectiveValue();
    rg = fabs(rg);

    while (rg > prerg)
    {
      // counter += 1;
      // s_dual = omp_get_wtime();
      ADMM_DualUpdate();
      // e_dual = omp_get_wtime();
      // e_dual = clock();
      // p_dual += (e_dual-s_dual);
      // double s_primal = clock();
      ADMM_PrimalUpdate();
      // double e_primal = clock();
      // p_primal += (e_primal - s_primal)/CLOCKS_PER_SEC;
      // s_rg = omp_get_wtime();
      // s_rg = clock();
      rg = fabs(GetUEGap());
      // e_rg = omp_get_wtime();
      // e_rg = clock();
      // p_rg += (e_rg - s_rg);
      // if (counter == 46512 || counter == 46511)
      // {
      //     cout << counter << " , " << rg << " , "<< inter - start << " , " << p_primal_eff <<endl;  //p_dual+p_rg+p_primal_eff
      //     break;
      // }
      // else
      cout << "xx" << endl;
    }
    prerg = rg;
    inter = omp_get_wtime();
    // ECT = p_dual+p_primal_eff +p_rg;
    cout << counter << " , " << rg << " , " << inter - start << " , " << p_primal_eff << endl; // p_dual+p_rg+p_primal_eff
  }

  // double end = clock();

  // cout<<counter<<", RG is: "<<rg<<endl;
  // outfile.close();
  //
  // cout  << end-start << endl;
}


int main()
{
  FIRST_THRU_NODE = 1;

  ReadData("AnAll"); // Small  //SF  //C:/Users/xych/Qsync/working/ADMM2006/Code/ADMM
  // printf("stage1\n");
  ADMM();

  // system("pause");
  return 0;
}
