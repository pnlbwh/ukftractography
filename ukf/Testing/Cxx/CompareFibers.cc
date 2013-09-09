#if 0
#include "vtkReader.h"
#else
#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkCellArray.h"
#include "vtkSmartPointer.h"
#include "math.h"
#endif
#include "fiber.h"
#include <iostream>

#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif


int main(int argc, char *argv[])
{
  if( argc < 3 )
    {
    std::cerr << "Usage:" << std::endl;
    std::cerr << argv[0] << "testFiber compareFiber" << std::endl;
    return -1;
    }

#if 0
  std::vector<Fiber> test_fibers;
  std::vector<Fiber> compare_fibers;

  vtkReader<Fiber> * const reader = new vtkReader<Fiber>();

  std::cout << "** Reading test VTK file... " << argv[1] << "\n";
  reader->SetInputPath(argv[1]);
  reader->SetOutputFibers(test_fibers);
  reader->SetReadFieldData(true);
  if( reader->Run() )
    {
    std::cerr << "Error reading the test file" << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "** Reading compare VTK file..." << argv[2] << "\n";
  reader->SetInputPath(argv[2]);
  reader->SetOutputFibers(compare_fibers);
  reader->SetReadFieldData(true);
  if( reader->Run() )
    {
    std::cerr << "Error reading the compare file" << std::endl;
    return EXIT_FAILURE;
    }
  if( test_fibers.size() != compare_fibers.size() )
    {
    std::cerr << "Test and compare fibers are not equal in amount." << std::endl;
    return EXIT_FAILURE;
    }
  for( unsigned int i = 0; i < test_fibers.size(); ++i )
    {
    if( test_fibers[i].Points.size() != compare_fibers[i].Points.size() )
      {
      std::cerr << "One of the test fibers is different in length compared to the compare fiber" << std::endl;
      return EXIT_FAILURE;
      }
    // add a test to compare point to field size
    for( unsigned int j = 0; j < test_fibers[i].Points.size(); ++j )
      {
      double difference =
        ( (test_fibers[i].Points[j]._[0] + test_fibers[i].Points[j]._[1] + test_fibers[i].Points[j]._[2])
        - (compare_fibers[i].Points[j]._[0] + compare_fibers[i].Points[j]._[1]
           + compare_fibers[i].Points[j]._[2]) );

      difference = (difference < 0) ? -difference : difference; // fix sign

      if( difference > 1e-9 )
        {
        std::cerr << "Difference in Points is above tolerance (1e-9): " << difference << std::endl;
        return EXIT_FAILURE;
        }
      // do the same for NMSE field

      }

    }
#else
  vtkSmartPointer<vtkPolyDataReader> reader1 =
    vtkSmartPointer<vtkPolyDataReader>::New();
  vtkSmartPointer<vtkPolyDataReader> reader2 =
    vtkSmartPointer<vtkPolyDataReader>::New();
  reader1->SetFileName(argv[1]);
  reader2->SetFileName(argv[2]);
  reader1->Update();
  reader2->Update();
  vtkSmartPointer<vtkPolyData> input1 = reader1->GetOutput();
  vtkSmartPointer<vtkPolyData> input2 = reader2->GetOutput();
  std::cerr << "Input 1" << std::endl
            << "Verts " << input1->GetNumberOfVerts()
            << " Lines " << input1->GetNumberOfLines()
            << " Polys " << input1->GetNumberOfPolys()
            << " Strips " << input1->GetNumberOfStrips()
            << std::endl;
  std::cerr << "Input 2" << std::endl
            << "Verts " << input2->GetNumberOfVerts()
            << " Lines " << input2->GetNumberOfLines()
            << " Polys " << input2->GetNumberOfPolys()
            << " Strips " << input2->GetNumberOfStrips()
            << std::endl;

  vtkIdType size1 = input1->GetNumberOfPoints();
  vtkIdType size2 = input2->GetNumberOfPoints();
  if(size1 != size2)
    {
    std::cerr << "first file fiber has " << size1
              << " points, second file has " << size2
              << std::endl;
    return EXIT_FAILURE;
    }

  for(vtkIdType i = 0; i < size1; ++i)
    {
    const double *pt1 = input1->GetPoint(i);
    const double *pt2 = input2->GetPoint(i);
    double distance = sqrt( ((pt1[0] - pt2[0]) * (pt1[0] - pt2[0]))
                            + ((pt1[1] - pt2[1]) * (pt1[1] - pt2[1]))
                            + ((pt1[2] - pt2[2]) * (pt1[2] - pt2[2])));
    if(distance > 1.0E-3)
      {
      std::cerr << "Difference in Points is above tolerance (1e-9): " << distance << std::endl;
      return EXIT_FAILURE;
      }
    }

#endif
  std::cout << "Test succeded!\n";
  return EXIT_SUCCESS;
}
