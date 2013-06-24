#include "vtkReader.h"
#include "CompareFiber.h"
#include <iostream>

#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif

//DELETEME static int RegressionTestImage (const char *, const char *, int, bool);

int main(int argc, char *argv[])
{
   if(argc < 3){
      std::cerr << "Usage:" << std::endl;
      std::cerr << argv[0] << "testFiber compareFiber" << std::endl;
      return -1;
   }

   std::vector<Fiber> test_fibers;
   std::vector<Fiber> compare_fibers;

   vtkReader * const reader = new vtkReader();

   std::cout << "** Reading test VTK file... " << argv[1] << "\n";
   reader->SetInputPath(argv[1]);
   reader->SetOutputFibers(test_fibers);
   reader->SetReadFieldData(true);
   if (reader->Run()) {
     std::cerr << "Error reading the test file" << std::endl;
     return EXIT_FAILURE;
   }

   std::cout << "** Reading compare VTK file..." << argv[2] << "\n";
   reader->SetInputPath(argv[2]);
   reader->SetOutputFibers(compare_fibers);
   reader->SetReadFieldData(true);
   if (reader->Run()) {
     std::cerr << "Error reading the compare file" << std::endl;
     return EXIT_FAILURE;
   }

   if (test_fibers.size() != compare_fibers.size()) {
     std::cerr << "Test and compare fibers are not equal in amount." << std::endl;
     return EXIT_FAILURE;
   }

  for (unsigned int i=0; i<test_fibers.size(); ++i) {
    if (test_fibers[i].Points.size() != compare_fibers[i].Points.size()) {
     std::cerr << "One of the test fibers is different in length compared to the compare fiber" << std::endl;
     return EXIT_FAILURE;
    }
    // add a test to compare point to field size
    for (unsigned int j=0; j<test_fibers[i].Points.size(); ++j) {
      double difference = ((test_fibers[i].Points[j]._[0] + test_fibers[i].Points[j]._[1] + test_fibers[i].Points[j]._[2]) -
			  (compare_fibers[i].Points[j]._[0] + compare_fibers[i].Points[j]._[1] + compare_fibers[i].Points[j]._[2]));

      difference = (difference < 0) ? -difference : difference; // fix sign

     if( difference>1e-9 ){
      std::cerr << "Difference in Points is above tolerance (1e-9): " << difference << std::endl;
      return EXIT_FAILURE;
    }
    // do the same for NMSE field


    }

  }
  std::cout << "Test succeded!\n";
  return EXIT_SUCCESS;
}
