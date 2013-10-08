/**
 * \file seed.cc
 * \brief implementation of seed.h
*/

#include "seed.h"
#include <iostream>

void PrintSeedInfo(const std::vector<SeedPointInfo>& seed_infos)
{
  for( size_t i = 0; i < seed_infos.size(); ++i )
    {
    const SeedPointInfo& inf = seed_infos[i];
    std::cout << "Seed point " << i << ":" << std::endl;
    std::cout << "  Pos: (" << inf.point[0] << ", " << inf.point[1] << ", "
              << inf.point[2] << ")" << std::endl;
    std::cout << "  State:";
    for( size_t j = 0; j < inf.state.size(); ++j )
      {
      std::cout << " " << inf.state[j];
      }
    std::cout << std::endl;
    std::cout << "  Dir: (" << inf.start_dir[0] << ", " << inf.start_dir[1]
              << ", " << inf.start_dir[2] << ")" << std::endl;
    std::cout << "  FA:" << inf.fa << std::endl;
    }
}
