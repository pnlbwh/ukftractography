#ifndef __vtkReadWrite_h
#define __vtkReadWrite_h
#include "itksys/SystemTools.hxx"
#include "vtkVersionMacros.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkXMLPolyDataWriter.h"

class vtkPolyDataIO
{
public:
  static vtkSmartPointer<vtkPolyData>
  Read(const char *filename)
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

  static void
  Write(const vtkPolyData *pd, const char *filename, bool binary, bool compressed)
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
#if (VTK_MAJOR_VERSION < 6)
        writer->SetInput(dataObject);
#else
        writer->SetInputData(dataObject);
#endif
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
#if (VTK_MAJOR_VERSION < 6)
        writer->SetInput(dataObject);
#else
        writer->SetInputData(dataObject);
#endif
        writer->SetFileName(filename);
        writer->Write();
        }
    }

};

#endif // __vtkReadWrite_h

