#include "ConvertVTKCLP.h"
#include "itksys/SystemTools.hxx"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkXMLPolyDataWriter.h"

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

void
WritePolyData(const vtkPolyData *pd, const char *filename, bool binary, bool compressed)
{
  const std::string ext(itksys::SystemTools::GetFilenameExtension(filename));
  // all the casting is because vtk SetInput isn't const-correct.
  vtkDataObject *dataObject =
    const_cast<vtkDataObject *>(static_cast<const vtkDataObject *>(pd));
  if(ext == ".vtp")
    {
    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
      vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    if(binary)
      {
      writer->SetDataModeToBinary();
      }
    else
      {
      writer->SetDataModeToAscii();
      }
    if(compressed)
      {
      writer->SetCompressorTypeToZLib();
      }
    else
      {
      writer->SetCompressorTypeToNone();
      }
    writer->SetInput(dataObject);
    writer->SetFileName(filename);
    writer->Write();
    }
  else
    {
    vtkSmartPointer<vtkPolyDataWriter> writer =
      vtkSmartPointer<vtkPolyDataWriter>::New();
    if(binary)
      {
      writer->SetFileTypeToBinary();
      }
    writer->SetInput(dataObject);
    writer->SetFileName(filename);
    writer->Write();
    }
}

int main(int argc, char *argv[])
{
  // PARSE_ARGS ;
  GENERATE_LOGO;
  GENERATE_XML;
  GENERATE_TCLAP;
  GENERATE_ECHOARGS;
  GENERATE_ProcessInformationAddressDecoding;

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
    pd = ReadPolyData(inputFile.c_str());
    }
  catch(...)
    {
    std::cerr << "Failed reading " << inputFile << std::endl;
    return EXIT_FAILURE;
    }
  try
    {
    WritePolyData(pd, outputFile.c_str(), writeBinary, writeCompressed);
    }
  catch(...)
    {
    std::cerr << "Failed writing " << outputFile << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
