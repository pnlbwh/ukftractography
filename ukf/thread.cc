/**
 * \file thread.cc
 * \brief implementation of thread.h
 * \author Yinpeng Li (mousquetaires@unc.edu)
*/

#include "thread.h"
#include "itkMultiThreader.h"
#include <cassert>

WorkDistribution GenerateWorkDistribution(const int num_threads, const int total_num_works)
{
  WorkDistribution distribution;

  distribution.resize(num_threads);
  for( int i = 0; i < num_threads; i++ )
    {
    int iter = i;
    while( iter < total_num_works )
      {
      distribution[i].push_back(iter);
      iter += num_threads;
      }
    }

  return distribution;
}

// void ProgressThread::operator()()
// {
// const int num_work_threads = static_cast<int>(work_threads_.size()) ;

// int total_num_works = 0 ;
// for (size_t i = 0; i < work_threads_.size(); i++) {
// total_num_works += work_threads_[i]->get_num_works() ;
// }

// std::cout << "Number of threads: " << num_work_threads << std::endl ;
// std::cout << "Work distributed: " << total_num_works << std::endl ;

// boost::progress_display disp(static_cast<unsigned long>(total_num_works)) ;

// int num_works_completed = 0 ;

// while (num_works_completed < total_num_works) {
// boost::this_thread::sleep(boost::posix_time::seconds(PROGRESS_REPORT_INTERVAL)) ;

// int temp = 0 ;
// for (size_t i = 0; i < work_threads_.size(); i++) {
// temp += work_threads_[i]->get_index() ;
// }

// disp += temp - num_works_completed ;

// num_works_completed = temp ;
// }
// }

ITK_THREAD_RETURN_TYPE ThreadCallback(void *arg)
{
  int id_ =
    ( (itk::MultiThreader::ThreadInfoStruct *)( arg ) )->ThreadID;
  thread_struct *str =
    (thread_struct *)( ( (itk::MultiThreader::ThreadInfoStruct *)( arg ) )->UserData );
  WorkDistribution                                     work_distribution = *str->work_distribution;
  WorkList &                                           work_list_ = work_distribution[id_];
  std::vector<Fiber>&                                  output_fiber_group_ = *str->output_fiber_group_;
  std::vector<SeedPointInfo>&                          seed_infos_ = *str->seed_infos_;
  std::vector<std::vector<SeedPointInfo> >&            branching_seed_info_vec = *str->branching_seed_info_vec;
  std::vector<std::vector<BranchingSeedAffiliation> >& branching_seed_affiliation_vec =
    *str->branching_seed_affiliation_vec;

  for( WorkList::const_iterator it = work_list_.begin(); it != work_list_.end(); it++ )
    {

    if( str->num_tensors_ == 3 )
      {
      str->tractography_->Follow3T(id_, *it, seed_infos_[*it], output_fiber_group_[*it],
                                   str->branching_, branching_seed_info_vec[id_], branching_seed_affiliation_vec[id_]);
      }
    else if( str->num_tensors_ == 2 )
      {
      str->tractography_->Follow2T(id_, *it, seed_infos_[*it], output_fiber_group_[*it],
                                   str->branching_, branching_seed_info_vec[id_], branching_seed_affiliation_vec[id_]);
      }
    else
      {
      str->tractography_->Follow1T(id_, seed_infos_[*it], output_fiber_group_[*it]);
      }

    assert(branching_seed_info_vec[id_].size() == branching_seed_affiliation_vec[id_].size() );
    }

  return ITK_THREAD_RETURN_VALUE;
}
