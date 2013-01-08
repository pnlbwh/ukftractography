/**
 * \file thread.cc
 * \brief implementation of thread.h
 * \author Yinpeng Li (mousquetaires@unc.edu)
*/

#include "thread.h"
#include <boost/progress.hpp>
#include <cassert>

WorkDistribution GenerateWorkDistribution(const int num_threads, const int total_num_works)
{
  WorkDistribution distribution ;
  distribution.resize(num_threads) ;

  for (int i = 0; i < num_threads; i++) {
    int iter = i ;
    while (iter < total_num_works) {
      distribution[i].push_back(iter) ;
      iter += num_threads ;
    }
  }

  return distribution ;
}

Thread::~Thread() {}



ProgressThread::ProgressThread(const std::vector<TractographyThread *>& work_threads): work_threads_(work_threads)
{
}

void ProgressThread::operator()()
{
  const int num_work_threads = static_cast<int>(work_threads_.size()) ;

  int total_num_works = 0 ;
  for (size_t i = 0; i < work_threads_.size(); i++) {
    total_num_works += work_threads_[i]->get_num_works() ;
  }

  std::cout << "Number of threads: " << num_work_threads << std::endl ;
  std::cout << "Work distributed: " << total_num_works << std::endl ;

  boost::progress_display disp(static_cast<unsigned long>(total_num_works)) ;

  int num_works_completed = 0 ;

  while (num_works_completed < total_num_works) {
    boost::this_thread::sleep(boost::posix_time::seconds(PROGRESS_REPORT_INTERVAL)) ;

    int temp = 0 ;
    for (size_t i = 0; i < work_threads_.size(); i++) {
      temp += work_threads_[i]->get_index() ;
    }

    disp += temp - num_works_completed ;

    num_works_completed = temp ;
  }
}

ProgressThread::~ProgressThread() {}






TractographyThread::TractographyThread(	const WorkList& work_list,
                                        const int id,
                                        Tractography *tractography,
                                        std::vector<SeedPointInfo>& seed_infos,
                                        bool branching,
                                        int num_tensors,
                                        std::vector<Fiber>& output_fiber_group
                                      ):
  work_list_(work_list),
  index_(0),
  id_(id),
  tractography_(tractography),
  seed_infos_(seed_infos),
  branching_(branching),
  num_tensors_(num_tensors),
  output_fiber_group_(output_fiber_group)
{
}

int TractographyThread::get_num_works() const
{
  return static_cast<int>(work_list_.size()) ;
}

int TractographyThread::get_index() const
{
  boost::lock_guard<boost::mutex> lock(index_mutex_) ;
  return index_ ;
}

int TractographyThread::get_id() const
{
  return id_ ;
}

void TractographyThread::lock_index() const
{
  index_mutex_.lock() ;
}

void TractographyThread::operator()()
{
  //assert(id_ != 0) ;
  for (WorkList::const_iterator it = work_list_.begin(); it != work_list_.end(); it++) {

    if (num_tensors_ == 3) {
      tractography_->Follow3T(id_, *it, seed_infos_[*it], output_fiber_group_[*it], branching_, branching_seed_info_, branching_seed_affiliation_) ;
    } else if (num_tensors_ == 2) {
      tractography_->Follow2T(id_, *it, seed_infos_[*it], output_fiber_group_[*it], branching_, branching_seed_info_, branching_seed_affiliation_) ;
    } else {
      tractography_->Follow1T(id_, seed_infos_[*it], output_fiber_group_[*it]);
    }

    assert(branching_seed_info_.size() == branching_seed_affiliation_.size()) ;

    {
      // Increase the index under the protection of mutex
      boost::lock_guard<boost::mutex> lock(index_mutex_) ;
      index_++ ;
    }

  }
}

TractographyThread::~TractographyThread() {}
