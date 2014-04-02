#ifndef __fiberbundle_h
#define __fiberbundle_h
#include "fiber.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkDataArray.h"
#include <map>
#include <vector>
class fiberbundle
{
public:
  typedef std::vector<Fiber> FiberVector;

  void ReadFibers( std::string inputFibersFileName );

  void WriteFibers( std::string outputFibersFileName, bool writeAscii, bool writeUnCompressed );

  FiberVector &GetFibers() { return m_FiberBundle; }

  void Print();

private:
  template <typename TArray>
  bool AddArray(Fiber &curFiber,
                vtkIdType nPoints,
                const vtkIdType *pts,
                const std::string &curname,
                vtkDataArray *curArray)
    {
      TArray *array = TArray::SafeDownCast(curArray);
      if(array == 0)
        {
        return false;
        }
      if(array->GetNumberOfComponents() == 1)
        {
        std::vector<float> curVec;
        for(vtkIdType curPoint = 0; curPoint < nPoints; ++curPoint)
          {
          curVec.push_back(array->GetValue(pts[curPoint]));
          }
        curFiber.Fields[curname] = curVec;
        }
      else if(array->GetNumberOfComponents() == 9)
        {
        stdMat_t curVec;
        for(vtkIdType curPoint = 0; curPoint < nPoints; ++curPoint)
          {
          const double *curval = array->GetTuple(pts[curPoint]);
          mat33_t mat;
          for(unsigned int i = 0, v = 0; i < 3; ++i)
            {
            for(unsigned int j = 0; j < 3; ++j, ++v)
              {
              mat(i,j) = curval[v];
              }
            }
          curVec.push_back(mat);
          }
        curFiber.Tensors[curname] = curVec;
        }
      return true;
    }
  FiberVector m_FiberBundle;
  vtkSmartPointer<vtkPolyData> m_PolyData;
  std::string m_InputFibersFileName;
};

#endif // __fiberbundle_h
