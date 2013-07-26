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
    std::cout << "  Pos: (" << inf.point._[0] << ", " << inf.point._[1] << ", "
              << inf.point._[2] << ")" << std::endl;
    std::cout << "  State:";
    for( size_t j = 0; j < inf.state.size(); ++j )
      {
      std::cout << " " << inf.state[j];
      }
    std::cout << std::endl;
    std::cout << "  Dir: (" << inf.start_dir._[0] << ", " << inf.start_dir._[1]
              << ", " << inf.start_dir._[2] << ")" << std::endl;
    std::cout << "  FA:" << inf.fa << std::endl;
    }
}
