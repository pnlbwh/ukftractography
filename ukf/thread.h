/**
 * \file thread.h
 * \brief Contains the functions that deal with the distribution of the threads for mutli-threaded execution.
 * \author Yinpeng Li (mousquetaires@unc.edu)
*/

#ifndef THREAD_H_
#define THREAD_H_

#include <vector>
#include "tractography.h"

typedef std::vector<int>      WorkList;
typedef std::vector<WorkList> WorkDistribution;

// class ProgressThread: public Thread
// {
// protected:
// const std::vector<TractographyThread *>& work_threads_ ;
// public:
// ProgressThread(const std::vector<TractographyThread *>& work_threads) ;
// void operator()() ;
// ~ProgressThread() ;
// } ;

// const int PROGRESS_REPORT_INTERVAL = 2 ;	//Report progress once per 2 seconds

WorkDistribution GenerateWorkDistribution(const int num_threads, const int total_num_works);

void * ThreadCallback(void *arg);

struct thread_struct
  {
  Tractography *tractography_;
  WorkDistribution* work_distribution;
  std::vector<SeedPointInfo>* seed_infos_;
  bool branching_;
  int num_tensors_;
  std::vector<Fiber>* output_fiber_group_;
  std::vector<std::vector<SeedPointInfo> >* branching_seed_info_vec;
  std::vector<std::vector<BranchingSeedAffiliation> >* branching_seed_affiliation_vec;
  };

#endif
