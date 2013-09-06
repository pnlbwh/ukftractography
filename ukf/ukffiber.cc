/*
 * \file fiber.cc
 * \brief Implements functions defined in fiber.h
 * \author Yinpeng Li (mousquetaires@unc.edu)
*/

#include "ukffiber.h"
#include <iostream>

void PostProcessFibers( const std::vector<UKFFiber>& raw_primary,
                        const std::vector<UKFFiber>& raw_branch,
                        const std::vector<BranchingSeedAffiliation>& branching_seed_affiliation,
                        const bool branches_only,
                        std::vector<UKFFiber>& fibers)
{
  assert(fibers.empty() );
  const int num_half_fibers = static_cast<int>(raw_primary.size() );
  assert( (num_half_fibers > 0) && (num_half_fibers % 2 == 0) );

  const bool record_fa = !raw_primary[0].fa.empty();
  const bool record_fa2 = !raw_primary[0].fa2.empty();
  const bool record_trace = !raw_primary[0].trace.empty();
  const bool record_trace2 = !raw_primary[0].trace2.empty();
  const bool record_free_water = !raw_primary[0].free_water.empty();
  const bool record_normMSE = !raw_primary[0].normMSE.empty();
  const bool record_cov = !raw_primary[0].covariance.empty();

  const int num_primary_fibers = branches_only ? 0 : num_half_fibers / 2;
  const int num_branches = static_cast<int>(raw_branch.size() );

  assert(num_branches == static_cast<int>(branching_seed_affiliation.size() ) );

  std::vector<int> num_points_on_primary_fiber(num_primary_fibers); // Number of points on each full primary fiber
  std::vector<int> num_points_on_branch(num_branches);              // Number of points on each full branch
  // Compute the numbers of points on full primary fibers
  for( int i = 0; i < num_primary_fibers; i++ )
    {
    num_points_on_primary_fiber[i] =
      static_cast<int>(raw_primary[2 * i].position.size() + raw_primary[2 * i + 1].position.size() ) - 1;
    // The two fibers share the same seed point
    }
  // Compute the numbers of points on full branches
  for( int i = 0; i < num_branches; i++ )
    {
    const size_t fiber_index = branching_seed_affiliation[i].fiber_index_;
    const int    position_on_fiber = branching_seed_affiliation[i].position_on_fiber_;

    if( fiber_index % 2 == 0 )
      {
      num_points_on_branch[i] = static_cast<int>(raw_primary[fiber_index + 1].position.size()
                                                 + raw_branch[i].position.size() )
        + position_on_fiber - 1;
      }
    else
      {
      num_points_on_branch[i] = static_cast<int>(raw_primary[fiber_index - 1].position.size()
                                                 + raw_branch[i].position.size() )
        + position_on_fiber - 1;
      }
    }

  // Compute the number of valid full primary fibers and the number of valid full branches
  int num_valid_primary_fibers = 0;
  int num_valid_branch = 0;
  for( int i = 0; i < num_primary_fibers; i++ )
    {
    if( num_points_on_primary_fiber[i] >= MINIMUM_NUM_POINTS_ON_FIBER )
      {
      num_valid_primary_fibers++;
      }
    }
  for( int i = 0; i < num_branches; i++ )
    {
    const size_t fiber_index = branching_seed_affiliation[i].fiber_index_;
    const int    position_on_fiber = branching_seed_affiliation[i].position_on_fiber_;

    if( (num_points_on_branch[i] >= MINIMUM_NUM_POINTS_ON_FIBER) &&
        (position_on_fiber + FIBER_TAIL_THRESHOLD < static_cast<int>(raw_primary[fiber_index].position.size() ) )
        // NOTE that when a branch originates near the end of a primary fiber, it's very probable that this branch
        // will contain lots of error, so this kind of branch is deemed as invalid
        )
      {
      num_valid_branch++;
      }
    }

  fibers.resize(num_valid_primary_fibers + num_valid_branch);

  // Join half fibers to full fibers
  int counter = 0;
  for( int i = 0; i < num_primary_fibers; i++ )
    {

    if( num_points_on_primary_fiber[i] < MINIMUM_NUM_POINTS_ON_FIBER )
      {
      continue;
      }

    const UKFFiber& first_half = raw_primary[2 * i];
    const UKFFiber& second_half = raw_primary[2 * i + 1];

    fibers[counter].position.resize(num_points_on_primary_fiber[i]);
    if( record_fa )
      {
      fibers[counter].fa.resize(num_points_on_primary_fiber[i]);
      }
    if( record_fa2 )
      {
      fibers[counter].fa2.resize(num_points_on_primary_fiber[i]);
      }
    if( record_trace )
      {
      fibers[counter].trace.resize(num_points_on_primary_fiber[i]);
      }
    if( record_trace2 )
      {
      fibers[counter].trace2.resize(num_points_on_primary_fiber[i]);
      }
    if( record_free_water )
      {
      fibers[counter].free_water.resize(num_points_on_primary_fiber[i]);
      }
    if( record_normMSE )
      {
      fibers[counter].normMSE.resize(num_points_on_primary_fiber[i]);
      }
    fibers[counter].norm.resize(num_points_on_primary_fiber[i]);
    fibers[counter].state.resize(num_points_on_primary_fiber[i]);
    if( record_cov )
      {
      fibers[counter].covariance.resize(num_points_on_primary_fiber[i]);
      }

    int k = 0;
    // The first point in the first_half, namely the seed point in the first half, is excluded
    for( int j = static_cast<int>(first_half.position.size() ) - 1; j > 0; j-- )
      {
      fibers[counter].position[k] = first_half.position[j];
      if( record_fa )
        {
        fibers[counter].fa[k] = first_half.fa[j];
        }
      if( record_fa2 )
        {
        fibers[counter].fa2[k] = first_half.fa2[j];
        }
      if( record_trace )
        {
        fibers[counter].trace[k] = first_half.trace[j];
        }
      if( record_trace2 )
        {
        fibers[counter].trace2[k] = first_half.trace2[j];
        }
      if( record_free_water )
        {
        fibers[counter].free_water[k] = first_half.free_water[j];
        }
      if( record_normMSE )
        {
        fibers[counter].normMSE[k] = first_half.normMSE[j];
        }
      fibers[counter].norm[k] = first_half.norm[j];
      fibers[counter].state[k] = first_half.state[j];
      if( record_cov )
        {
        fibers[counter].covariance[k] = first_half.covariance[j];
        }
      k++;
      }
    for( int j = 0; j < static_cast<int>(second_half.position.size() ); j++ )
      {
      fibers[counter].position[k] = second_half.position[j];
      if( record_fa )
        {
        fibers[counter].fa[k] = second_half.fa[j];
        }
      if( record_fa2 )
        {
        fibers[counter].fa2[k] = second_half.fa2[j];
        }
      if( record_trace )
        {
        fibers[counter].trace[k] = second_half.trace[j];
        }
      if( record_trace2 )
        {
        fibers[counter].trace2[k] = second_half.trace2[j];
        }
      if( record_free_water )
        {
        fibers[counter].free_water[k] = second_half.free_water[j];
        }
      if( record_normMSE )
        {
        fibers[counter].normMSE[k] = second_half.normMSE[j];
        }
      fibers[counter].norm[k] = second_half.norm[j];
      fibers[counter].state[k] = second_half.state[j];
      if( record_cov )
        {
        fibers[counter].covariance[k] = second_half.covariance[j];
        }
      k++;
      }

    assert(k == num_points_on_primary_fiber[i]);

    counter++;
    }

  assert(counter == num_valid_primary_fibers);
  // Backtrace the branches
  for( int i = 0; i < num_branches; i++ )
    {

    if( num_points_on_branch[i] < MINIMUM_NUM_POINTS_ON_FIBER )
      {
      continue;
      }

    const size_t fiber_index = branching_seed_affiliation[i].fiber_index_;
    const int    position_on_fiber = branching_seed_affiliation[i].position_on_fiber_;
    if( position_on_fiber + FIBER_TAIL_THRESHOLD >= static_cast<int>(raw_primary[fiber_index].position.size() ) )
      {
      continue;
      }

    int first_half_index = -1;                             // The first half is the reversed one
    int second_half_index = static_cast<int>(fiber_index); // The second half is the primary fiber from which the branch
                                                           // originates

    if( fiber_index % 2 == 0 )
      {
      first_half_index = static_cast<int>(fiber_index) + 1;
      }
    else
      {
      first_half_index = static_cast<int>(fiber_index) - 1;
      }

    const UKFFiber& first_half = raw_primary[first_half_index];
    const UKFFiber& second_half = raw_primary[second_half_index];
    const UKFFiber& branch = raw_branch[i];  // This is the un-back-traced branch

    fibers[counter].position.resize(num_points_on_branch[i]);
    if( record_fa )
      {
      fibers[counter].fa.resize(num_points_on_branch[i]);
      }
    if( record_fa2 )
      {
      fibers[counter].fa2.resize(num_points_on_branch[i]);
      }
    if( record_trace )
      {
      fibers[counter].trace.resize(num_points_on_branch[i]);
      }
    if( record_trace2 )
      {
      fibers[counter].trace2.resize(num_points_on_branch[i]);
      }
    if( record_free_water )
      {
      fibers[counter].free_water.resize(num_points_on_branch[i]);
      }
    if( record_normMSE )
      {
      fibers[counter].normMSE.resize(num_points_on_branch[i]);
      }
    fibers[counter].norm.resize(num_points_on_branch[i]);
    fibers[counter].state.resize(num_points_on_branch[i]);
    if( record_cov )
      {
      fibers[counter].covariance.resize(num_points_on_branch[i]);
      }

    int k = 0;
    // The first point in the first_half, namely the seed point in the first half, is excluded
    for( int j = static_cast<int>(first_half.position.size() ) - 1; j > 0; j-- )
      {
      fibers[counter].position[k] = first_half.position[j];
      if( record_fa )
        {
        fibers[counter].fa[k] = first_half.fa[j];
        }
      if( record_fa2 )
        {
        fibers[counter].fa2[k] = first_half.fa2[j];
        }
      if( record_trace )
        {
        fibers[counter].trace[k] = first_half.trace[j];
        }
      if( record_trace2 )
        {
        fibers[counter].trace2[k] = first_half.trace2[j];
        }
      if( record_free_water )
        {
        fibers[counter].free_water[k] = first_half.free_water[j];
        }
      if( record_normMSE )
        {
        fibers[counter].normMSE[k] = first_half.normMSE[j];
        }
      fibers[counter].norm[k] = first_half.norm[j];
      fibers[counter].state[k] = first_half.state[j];
      if( record_cov )
        {
        fibers[counter].covariance[k] = first_half.covariance[j];
        }
      k++;
      }
    // The point in the second_half where the branch originates is also excluded
    for( int j = 0; j < position_on_fiber; j++ )
      {
      fibers[counter].position[k] = second_half.position[j];
      if( record_fa )
        {
        fibers[counter].fa[k] = second_half.fa[j];
        }
      if( record_fa2 )
        {
        fibers[counter].fa2[k] = second_half.fa2[j];
        }
      if( record_trace )
        {
        fibers[counter].trace[k] = second_half.trace[j];
        }
      if( record_trace2 )
        {
        fibers[counter].trace2[k] = second_half.trace2[j];
        }
      if( record_free_water )
        {
        fibers[counter].free_water[k] = second_half.free_water[j];
        }
      if( record_normMSE )
        {
        fibers[counter].normMSE[k] = second_half.normMSE[j];
        }
      fibers[counter].norm[k] = second_half.norm[j];
      fibers[counter].state[k] = second_half.state[j];
      if( record_cov )
        {
        fibers[counter].covariance[k] = second_half.covariance[j];
        }
      k++;
      }
    for( int j = 0; j < static_cast<int>(branch.position.size() ); j++ )
      {
      fibers[counter].position[k] = branch.position[j];
      if( record_fa )
        {
        fibers[counter].fa[k] = branch.fa[j];
        }
      if( record_fa2 )
        {
        fibers[counter].fa2[k] = branch.fa2[j];
        }
      if( record_trace )
        {
        fibers[counter].trace[k] = branch.trace[j];
        }
      if( record_trace2 )
        {
        fibers[counter].trace2[k] = branch.trace2[j];
        }
      if( record_free_water )
        {
        fibers[counter].free_water[k] = branch.free_water[j];
        }
      if( record_normMSE )
        {
        fibers[counter].normMSE[k] = branch.normMSE[j];
        }
      fibers[counter].norm[k] = branch.norm[j];
      fibers[counter].state[k] = branch.state[j];
      if( record_cov )
        {
        fibers[counter].covariance[k] = branch.covariance[j];
        }
      k++;
      }

    assert(k == num_points_on_branch[i]);

    counter++;
    }

  assert(counter == num_valid_primary_fibers + num_valid_branch);
}
