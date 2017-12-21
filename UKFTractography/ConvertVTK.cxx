#include "ConvertVTKCLP.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPolyDataIO.h"

int main(int argc, char *argv[])
{
  PARSE_ARGS ;

  if(inputFile.empty())
    {
    std::cerr << "Missing input filename" << std::endl;
    return EXIT_FAILURE;
    }
  if(outputFile.empty())
    {
    std::cerr << "Missing output filename" << std::endl;
    return EXIT_FAILURE;
    }
  vtkSmartPointer<vtkPolyData> pd;
  try
    {
    pd = vtkPolyDataIO::Read(inputFile.c_str());
    }
  catch(...)
    {
    std::cerr << "Failed reading " << inputFile << std::endl;
    return EXIT_FAILURE;
    }
  try
    {
    vtkPolyDataIO::Write(pd, outputFile.c_str(), !writeAscii, !writeUnCompressed);
    }
  catch(...)
    {
    std::cerr << "Failed writing " << outputFile << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
