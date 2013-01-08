/**
 * \file thread.h
 * \brief Contains the functions that deal with the distribution of the threads for mutli-threaded execution.
 * \author Yinpeng Li (mousquetaires@unc.edu)
*/

#ifndef THREAD_H_
#define THREAD_H_

#include <vector>
#include <boost/thread.hpp>
#include "tractography.h"

typedef std::vector<int> WorkList ;
typedef std::vector<WorkList> WorkDistribution ;

/**
 * \class Thread
 * \brief Generic class for threads
*/
class Thread
{
public:
  Thread() {}
  virtual void operator()() = 0 ;
  virtual ~Thread() = 0 ;
} ;


/**
 * \class TractographyThread 
 * \brief One thread of the tractography 
*/
class TractographyThread: public Thread
{
protected:
  const WorkList& work_list_ ;
  int index_ ;
  mutable boost::mutex index_mutex_ ;
  int id_ ;
  Tractography *tractography_ ;	// Pointer to the algorithm implementation
  std::vector<SeedPointInfo>& seed_infos_ ;
  bool branching_ ;
  int num_tensors_ ;
  std::vector<Fiber>& output_fiber_group_ ;
public:
  TractographyThread(	const WorkList& work_list,
                      const int id,
                      Tractography *tractography,
                      std::vector<SeedPointInfo>& seed_infos,
                      bool branching,
                      int num_tensors,
                      std::vector<Fiber>& output_fiber_group
                    ) ;
  int get_num_works() const ;
  int get_index() const ;
  int get_id() const ;
  void lock_index() const ;
  void operator()() ;
  ~TractographyThread() ;
public:
  std::vector<SeedPointInfo> branching_seed_info_ ;
  std::vector<BranchingSeedAffiliation> branching_seed_affiliation_ ;
} ;

/**
 * \class ProgressThread 
 * \brief Thread for monitoring the progress of tractogrpahy 
*/
class ProgressThread: public Thread
{
protected:
  const std::vector<TractographyThread *>& work_threads_ ;
public:
  ProgressThread(const std::vector<TractographyThread *>& work_threads) ;
  void operator()() ;
  ~ProgressThread() ;
} ;

const int PROGRESS_REPORT_INTERVAL = 2 ;	//Report progress once per 2 seconds

WorkDistribution GenerateWorkDistribution(const int num_threads, const int total_num_works) ;

#endif
