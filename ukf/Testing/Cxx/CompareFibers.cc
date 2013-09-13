#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkCellArray.h"
#include "vtkSmartPointer.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "math.h"
#include "fiber.h"
#include <iostream>
#include "itksys/SystemTools.hxx"

#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif

vtkSmartPointer<vtkPolyData>
ReadPolyData(const char *filename)
{
  const std::string ext(itksys::SystemTools::GetFilenameExtension(filename));
  if(ext == ".vtp")
    {
    vtkSmartPointer<vtkXMLPolyDataReader> reader =
      vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName(filename);
    reader->Update();
    return vtkSmartPointer<vtkPolyData>(reader->GetOutput());
    }
  else
    {
    vtkSmartPointer<vtkPolyDataReader> reader =
      vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(filename);
    reader->Update();
    return vtkSmartPointer<vtkPolyData>(reader->GetOutput());
    }
}
int main(int argc, char *argv[])
{
  int rval = EXIT_SUCCESS;
  if( argc < 3 )
    {
    std::cerr << "Usage:" << std::endl;
    std::cerr << argv[0] << "testFiber compareFiber" << std::endl;
    return EXIT_FAILURE;
    }

  vtkSmartPointer<vtkPolyData> input1 = ReadPolyData(argv[1]);
  vtkSmartPointer<vtkPolyData> input2 = ReadPolyData(argv[2]);
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

  vtkPointData *pd1 = input1->GetPointData();
  vtkPointData *pd2 = input2->GetPointData();
  vtkPointData *masterPD;
  // pd1->Print(std::cerr);
  // pd2->Print(std::cerr);
  int numComponents1 = pd1->GetNumberOfComponents();
  int numComponents2 = pd2->GetNumberOfComponents();
  int maxNumComponents = numComponents1 > numComponents2 ? numComponents1 : numComponents2;
  masterPD = pd1;
  if(numComponents1 != numComponents2)
    {
    std::cerr << "Number of components mismatch:" << std::endl
              << "Components1 = " << numComponents1
              << " Components2 = " << numComponents2
              << std::endl;

    masterPD = numComponents1 > numComponents2 ? pd1 : pd2;
    }
  for(int i = 0; i < maxNumComponents; ++i)
    {
    const char *curname = masterPD->GetArrayName(i);
    std::cerr << curname << std::endl;
    vtkDataArray *array1 = pd1->GetArray(curname);
    vtkDataArray *array2 = pd2->GetArray(curname);
    if(array1 == 0)
      {
      std::cerr << "component " << curname << " missing in first file" << std::endl;
      rval = EXIT_FAILURE;
      continue;
      }
    if(array2 == 0)
      {
      std::cerr << "component " << curname << " missing in second file" << std::endl;
      rval = EXIT_FAILURE;
      continue;
      }
    vtkFloatArray *farray1 = vtkFloatArray::SafeDownCast(array1);
    vtkFloatArray *farray2 = vtkFloatArray::SafeDownCast(array1);
    vtkIdType farray1size = farray1->GetNumberOfTuples();
    vtkIdType farray2size = farray2->GetNumberOfTuples();
    if(farray1size != farray2size)
      {
      std::cerr << "Array size mismatch:" << std::endl
                << "array1 = " << farray1size
                << " array2 = " << farray2size
                << std::endl;
      return EXIT_FAILURE;
      }
    for(vtkIdType j = 0; j < farray1size; ++j)
      {
      float f1 = farray1->GetValue(j);
      float f2 = farray2->GetValue(j);
      if(fabs(f1 - f2) > 1.0e-3)
        {
        std::cerr << "Difference is aove tolerance "
                  << f1 << " " << f2 << std::endl;
        return EXIT_FAILURE;
        }
      }
    }

  if(rval == EXIT_SUCCESS)
    {
    std::cout << "Test succeded!\n";
    }
  return rval;
}
