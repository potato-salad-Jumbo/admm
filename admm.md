# RAW

```c++
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
    // first oder derivative, link and origin based, second part, with oblf variable;
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
        head_o_succecons += m_OLPrimal[oId][pLink->pOutNode->OutgoingLink[lId2]];
      }
    }
    head_o_predecons = 0;
    for (int lId2 = 0; lId2 < pLink->pOutNode->IncomingLink.size(); lId2++)
      if (pLink->pOutNode->IncomingLink[lId2] != lId)
        head_o_predecons += m_OLPrimal[oId][pLink->pOutNode->IncomingLink[lId2]];
    // cout<<" pLink_OutNode_o_Incoming Link Flow Sum WO itself: "<<head_o_predecons<<endl;
    head_o_cons = head_o_succecons - head_o_predecons;
    head_o_cons -= m_ONDem[oId][pLink->pOutNode->ID];
    // cout<<" pLink_OutNode_o_Outgoing - pLink_OutNode_o_Incoming-ONDem: "<<head_o_cons<<endl;
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
  while (terminate > 1e-5 && terminategap > 1e-3) // //terminate > 1e-9   //counter<1  //
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

```

# kernal

```
__global__ void updateOLPrimal1st_kernel(double* update_OLPrimal_1st, double* m_OLPrimal, int m_OriginSize, double snd, double l_BPR, double penalty, double* tail_cons, double* head_cons, double* m_ONDual) {

  int oId = blockIdx.x * blockDim.x + threadIdx.x;

  if (oId < m_OriginSize) {
    // 并行计算update_OLPrimal_1st
    double o_Lp1std_linear = m_OLPrimal[oId] * penalty * 2;
    double conscons = tail_cons[oId] - head_cons[oId];
    update_OLPrimal_1st[oId] = l_BPR + o_Lp1std_linear + conscons * penalty + m_ONDual[oId][tail] - m_ONDual[oId][head]; 
  }
}

__global__ void updateTerminate_kernel(double* terminate, double* update_OLPrimal_1st, double* m_OLPrimal, int m_OriginSize) {
  
  int oId = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (oId < m_OriginSize) {
    // 并行计算terminate
    atomicAdd(terminate, fabs(m_OLPrimal[oId] * update_OLPrimal_1st[oId])); 
  }
}

__global__ void updateOLPrimal_kernel(double* m_OLPrimal, double* update_OLPrimal_1st, int m_OriginSize, double step, double snd) {

  int oId = blockIdx.x * blockDim.x + threadIdx.x;

  if (oId < m_OriginSize) {
    // 并行更新m_OLPrimal
    double temp = m_OLPrimal[oId] - step * update_OLPrimal_1st[oId] / snd;
    if (temp < 0) temp = 0;
    m_OLPrimal[oId] = temp;
  }

}
```

# host

```c++

void linkSubproblem_v3(int lId) {
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
    
    
    double *d_update_OLPrimal_1st, *d_m_OLPrimal, *d_tail_cons, *d_head_cons, *d_ONDual;
    double *d_terminate;
    double terminate_host; // This will be used to copy the terminate variable back to the host

    int m_OriginSize = m_Origin.size();
    size_t size = m_OriginSize * sizeof(double);

    // Allocate memory on the device
    cudaMalloc((void **)&d_update_OLPrimal_1st, size);
    cudaMalloc((void **)&d_m_OLPrimal, size);
    cudaMalloc((void **)&d_tail_cons, size);
    cudaMalloc((void **)&d_head_cons, size);
    cudaMalloc((void **)&d_ONDual, size * m_NodeSize); // Assuming m_NodeSize is the size of the second dimension for m_ONDual
    cudaMalloc((void **)&d_terminate, sizeof(double));

    // Initialize terminate on device to 0
    cudaMemset(d_terminate, 0, sizeof(double));

    // Copy data from host to device
    cudaMemcpy(d_m_OLPrimal, m_OLPrimal, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_tail_cons, tail_o_cons, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_head_cons, head_o_cons, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_ONDual, m_ONDual, size * m_NodeSize, cudaMemcpyHostToDevice);

    // Determine the block size and grid size
    int blockSize = 256; // This is an example block size, you may need to tune it based on your GPU
    int gridSize = (m_OriginSize + blockSize - 1) / blockSize;

    // Launch kernels
    updateOLPrimal1st_kernel<<<gridSize, blockSize>>>(d_update_OLPrimal_1st, d_m_OLPrimal, m_OriginSize, snd, l_BPR, penalty, d_tail_cons, d_head_cons, d_ONDual);
    updateTerminate_kernel<<<gridSize, blockSize>>>(d_terminate, d_update_OLPrimal_1st, d_m_OLPrimal, m_OriginSize); 
    updateOLPrimal_kernel<<<gridSize, blockSize>>>(d_m_OLPrimal, d_update_OLPrimal_1st, m_OriginSize, step, snd);

    // Copy results back to host
    cudaMemcpy(m_OLPrimal, d_m_OLPrimal, size, cudaMemcpyDeviceToHost);
    cudaMemcpy(&terminate_host, d_terminate, sizeof(double), cudaMemcpyDeviceToHost);

    // Check for errors
    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        // Handle error...
    }

    // Use the terminate_host variable as needed...

    // Free device memory
    cudaFree(d_update_OLPrimal_1st);
    cudaFree(d_m_OLPrimal);
    cudaFree(d_tail_cons);
    cudaFree(d_head_cons);
    cudaFree(d_ONDual);
    cudaFree(d_terminate);

    // ... [The rest of your existing code after CUDA parts] ...
}

```

