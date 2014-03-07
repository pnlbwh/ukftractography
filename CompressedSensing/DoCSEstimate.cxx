#include "DoCSEstimate.h"
#include "BuildRidges.h"
#include "BuildSensor.h"
#include "MultiSample.h"
#include "BPDN_homotopy.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include "itkSubtractImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkAbsImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkVectorImageToImageAdaptor.h"
#include "itkCastImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
// #include "itkMedianImageFilter.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
template <typename TImage>
void SubtractImage(const TImage *t1, const TImage *t2, typename TImage::Pointer &rval)
{
  typedef itk::SubtractImageFilter<TImage,TImage,TImage>  SubtractFilterType;
  typename SubtractFilterType::Pointer sub1 = SubtractFilterType::New();
  sub1->SetInput1(t1);
  sub1->SetInput2(t2);
  sub1->Update();
  rval = sub1->GetOutput();
}

template <typename TImage>
void AddImage(const TImage *t1, const TImage *t2, typename TImage::Pointer &rval)
{
  typedef itk::AddImageFilter<TImage,TImage,TImage>  AddFilterType;
  typename AddFilterType::Pointer sub1 = AddFilterType::New();
  sub1->SetInput1(t1);
  sub1->SetInput2(t2);
  sub1->Update();
  rval = sub1->GetOutput();
}

template <typename TImage>
void MultiplyImage(const TImage *t1, const TImage *t2,typename TImage::Pointer &rval)
{
  typedef itk::MultiplyImageFilter<TImage,TImage,TImage>  MultiplyFilterType;
  typename MultiplyFilterType::Pointer sub1 = MultiplyFilterType::New();
  sub1->SetInput1(t1);
  sub1->SetInput2(t2);
  sub1->Update();
  rval = sub1->GetOutput();
}

template <typename TImage>
void MultiplyImage(const TImage *t1, typename TImage::PixelType constant, typename TImage::Pointer &rval)
{
  typedef itk::MultiplyImageFilter<TImage,TImage,TImage>  MultiplyFilterType;
  typename MultiplyFilterType::Pointer sub1 = MultiplyFilterType::New();
  sub1->SetInput1(t1);
  sub1->SetConstant(constant);
  sub1->Update();
  rval = sub1->GetOutput();
}

template <typename TImage>
void AllocImage(const itk::ImageBase<TImage::ImageDimension> *templateImage, typename TImage::Pointer &outp)
{
  outp = TImage::New();
  outp->CopyInformation(templateImage);
  outp->SetRegions(templateImage->GetLargestPossibleRegion());
  outp->Allocate();
}

template <typename TImage>
void AllocImage(const itk::ImageBase<TImage::ImageDimension> *templateImage,
                typename TImage::PixelType initval,
                typename TImage::Pointer &outp)
{
  AllocImage<TImage>(templateImage,outp);
  outp->FillBuffer(initval);
}

template <typename TImage>
void AllocVecImage(const itk::ImageBase<TImage::ImageDimension> *templateImage,
                   unsigned long vecSize,
                   typename TImage::Pointer &outp)
{
  outp = TImage::New();
  outp->CopyInformation(templateImage);
  outp->SetRegions(templateImage->GetLargestPossibleRegion());
  outp->SetNumberOfComponentsPerPixel(vecSize);
  outp->Allocate();
}

template <typename TImage>
double PNormImageDiff(const TImage *im1, const TImage *im2)
{
  itk::ImageRegionConstIterator<TImage> im1It(im1,im1->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<TImage> im2It(im2,im2->GetLargestPossibleRegion());
  double maxPixel = 0;
  for(im1It.GoToBegin(), im2It.GoToBegin(); !im1It.IsAtEnd(); ++im1It,++im2It)
    {
    const double current = std::abs(im1It.Get() - im2It.Get());
    if(current > maxPixel)
      {
      maxPixel = current;
      }
    }
  return maxPixel;
}

typedef std::vector<unsigned long> ShiftVectorType;
ShiftVectorType MakeShiftVector(unsigned long dim, signed short direction)
{
  ShiftVectorType rval;
  if(direction < 0)
    {
    rval.push_back(0);
    for(unsigned int i = 0; i < dim - 1; ++i)
      {
      rval.push_back(i);
      }
    }
  else
    {
    for(unsigned int i = 1; i < dim; ++i)
      {
      rval.push_back(i);
      }
    rval.push_back(dim - 1);
    }
  return rval;
}

template<typename TImage>
void SubtractShift(const TImage *inImage,unsigned axis,ShiftVectorType &shift, typename TImage::Pointer &rval)
{
  AllocImage<TImage>(inImage,rval);
  itk::ImageRegionIterator<TImage> rvalIt(rval,rval->GetLargestPossibleRegion());
  itk::ImageRegionConstIteratorWithIndex<TImage> inIt(inImage,inImage->GetLargestPossibleRegion());

  for(rvalIt.GoToBegin(), inIt.GoToBegin(); !rvalIt.IsAtEnd(); ++rvalIt, ++inIt)
    {
    typename TImage::IndexType index = inIt.GetIndex();
    index[axis] = shift[index[axis]];
    rvalIt.Set(inImage->GetPixel(index) - inIt.Get());
    }
}

template <class ImageType>
void
WriteImage(const ImageType *image,
           const std::string & filename)
{
  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename  WriterType::Pointer writer = WriterType::New();
  writer->UseCompressionOn();
  writer->UseInputMetaDataDictionaryOn();
  writer->SetFileName( filename.c_str() );
  writer->SetInput(image);
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "Exception Object caught: " << std::endl;
    std::cout << err << std::endl;
    throw;
    }
}

} // anon namespace

DoCSEstimate
::DoCSEstimate(NrrdFile &nrrdFile,MaskImageType::Pointer &maskImage,MatrixType &newGradients) :
  m_NrrdFile(nrrdFile), m_MaskImage(maskImage), m_NewGradients(newGradients),
  m_Pi(std::atan(static_cast<double>(1.0))*4)
{
  this->m_SpaceOrigin = nrrdFile.GetImage()->GetOrigin();
  this->m_BValue = nrrdFile.GetBValue();
  DWIVectorImageType::DirectionType dir = nrrdFile.GetImage()->GetDirection();
  this->m_SpaceDirections.conservativeResize(3,3);
  for(unsigned int i = 0; i < 3; ++i)
    {
    for(unsigned int j = 0; j < 3; ++j)
      {
      this->m_SpaceDirections(i,j) = dir[i][j];
      }
    }
}

/** Go from 3D Vector image to 2D Eigen Matrix.
 *  Just to keep things interesting, Eigen Matrices are column major,
 *  and ITK Images are row major. And the ITK image is [X Y X] image
 *  of Vectors with length numGradients, and the 2D Eigen matrix
 * is organized as [ row grad ] with row = X * Y * Z and grad as the numGradient.
 */
void
DoCSEstimate
::Reshape(const DWIVectorImageType::Pointer &src,unsigned rows, unsigned columns, MatrixType &outMatrix)
{
  outMatrix.resize(rows,columns);
  // S=reshape(S,[n d]);
  DWIVectorImageType::IndexType index = {{0,0,0}};
  DWIVectorImageType::SizeType size = src->GetLargestPossibleRegion().GetSize();
  for(unsigned int z = 0; z < size[2]; ++z)
    {
    index[2] = z;
    for(unsigned int x = 0; x < size[0]; ++x)
        {
        index[0] = x;
        for(unsigned int y = 0; y < size[1]; ++y)
          {
          index[1] = y;
          unsigned long row = (z * size[1] * size[0]) + (y * size[0]) + x;
          const DWIVectorImageType::PixelType &current = this->m_IntensityData->GetPixel(index);
          for(unsigned int grad = 0; grad < columns; ++grad)
            {
            outMatrix(row,grad) = current[grad];
            }
          }
        }
    }
}

/** Go from 2D Eigen Matrix to 3D Vector Image
 *  Reverses the Reshape operation above.
 *
 * TODO: both Eigen::Matrix and itk::Image access performance is
 * enhanced by locality of reference, but given they difference in
 * organization, one or the other accesses is going to have to skip
 * around in memeory. No idea which is better.
 */
void
DoCSEstimate
::Reshape(const MatrixType &src, DWIVectorImageType::Pointer templateImage, DWIVectorImageType::Pointer &outImage)
{
  unsigned int numComponents = templateImage->GetNumberOfComponentsPerPixel();
  AllocVecImage<DWIVectorImageType>(templateImage,numComponents,outImage);

  DWIVectorImageType::SizeType size = templateImage->GetLargestPossibleRegion().GetSize();
  DWIVectorImageType::IndexType index;
  for(unsigned int z = 0; z < size[2]; ++z)
    {
    index[2] = z;
      for(unsigned int x = 0; x < size[0]; ++x)
        {
        index[0] = x;
        for(unsigned int y = 0; y < size[1]; ++y)
          {
          index[1] = y;
          DWIVectorImageType::PixelType curPixel(numComponents);
          unsigned int row = (z * index[1] * index[0]) + (y * index[0]) + x;
          for(unsigned int grad = 0; grad < numComponents; ++grad)
            {
            curPixel[grad] = src(row,grad);
            }
        outImage->SetPixel(index,curPixel);
        }
      }
    }
}

bool
DoCSEstimate
::AverageB0AndExtractIntensity()
{

  MatrixType gradients = this->m_NrrdFile.GetGradients();
  std::vector<unsigned> b0Indices;
  std::vector<unsigned> gradientIndices;

  for(unsigned i = 0; i < gradients.rows(); ++i)
    {
    double norm = gradients.row(i).norm();
    if(norm < Eps)
      {
      b0Indices.push_back(i);
      }
    else
      {
      gradientIndices.push_back(i);
      }
    }

  const unsigned int numB0(b0Indices.size());
  if(numB0 == 0)
    {
    std::cerr << "input data doesn't have a b0 image" << std::endl;
    return false;
    }

  // allocate b0 average
  DWIVectorImageType::Pointer originalNrrd = this->m_NrrdFile.GetImage();
  DWIVectorImageType::RegionType region = originalNrrd->GetLargestPossibleRegion();
  AllocImage<B0AvgImageType>(originalNrrd,0.0,this->m_AverageB0);

  itk::ImageRegionConstIterator<DWIVectorImageType> b0It(originalNrrd,region);
  itk::ImageRegionIterator<B0AvgImageType> avgIt(this->m_AverageB0,region);

  for(b0It.GoToBegin(), avgIt.GoToBegin(); !b0It.IsAtEnd() && !avgIt.IsAtEnd(); ++b0It, ++avgIt)
    {
    const DWIVectorImageType::PixelType val = b0It.Get();
    B0AvgImageType::PixelType avg(0.0);
    for(unsigned int i = 0; i < numB0; ++i)
      {
      avg += val[b0Indices[i]];
      }
    avg /= static_cast<double>(numB0);
    if(std::fabs(avg) <= Eps)
      {
      avg = 1.0;
      }
    avgIt.Set(avg);
    }
  //
  // separate out gradient volumes
  AllocVecImage<DWIVectorImageType>(originalNrrd,gradientIndices.size(),this->m_IntensityData);

  const unsigned int numGradientIndices = gradientIndices.size();
  //
  // make new vectorimage with just the gradient values
  itk::ImageRegionIterator<DWIVectorImageType> newIt(this->m_IntensityData,region);
  for(b0It.GoToBegin(), newIt.GoToBegin(); !b0It.IsAtEnd() && !newIt.IsAtEnd(); ++b0It,++newIt)
    {
    DWIVectorImageType::PixelType newVal(numGradientIndices);
    const DWIVectorImageType::PixelType oldVal = b0It.Get();
    for(unsigned int i = 0; i < numGradientIndices; ++i)
      {
      newVal[i] = oldVal[gradientIndices[i]];
      }
    newIt.Set(newVal);
    }

  // get just the non-B0 gradients
  for(unsigned int i = 0; i < numGradientIndices; ++i)
    {
    this->m_VoxelLatticeAlignedGradientDirections.conservativeResize(i+1,3);
    this->m_VoxelLatticeAlignedGradientDirections.row(i) = gradients.row(gradientIndices[i]);
    }
  return true;
}


bool
DoCSEstimate
::Compute()
{
  if(!this->AverageB0AndExtractIntensity())
    {
    return false;
    }
  // MeasurementFrame was a parameter in the original computation
  // but isn't needed here -- gradients are converted to world
  // frame when the file gets loaded.

  //
  // MATLAB CODE IS
  // DWIIntensityData = single(DWIIntensityData) ./ averagedB0(:,:,:,ones(1,numGradientDirections));
  // %remove negative values
  // DWIIntensityData(DWIIntensityData<0)=eps;
  itk::ImageRegionIterator<DWIVectorImageType>
    ivIt(this->m_IntensityData,this->m_IntensityData->GetLargestPossibleRegion());
  unsigned int numGrad = this->m_IntensityData->GetNumberOfComponentsPerPixel();

  {
  itk::ImageRegionConstIterator<B0AvgImageType>
    b0It(this->m_AverageB0,this->m_AverageB0->GetLargestPossibleRegion());
  for(ivIt.Begin(),b0It.Begin(); !ivIt.IsAtEnd() && !b0It.IsAtEnd(); ++ivIt,++b0It)
    {
    DWIVectorImageType::PixelType vec = ivIt.Get();
    double curB0 = b0It.Get();
    if(std::fabs(curB0) > Eps)
      {
      for(unsigned int i = 0; i < numGrad; ++i)
        {
        vec[i] /= curB0;
        if(vec[i] < 0)
          {
          vec[i] = Eps;
          }
        }
      }
    else
      {
      for(unsigned int i = 0; i < numGrad; ++i)
        {
        vec[i] = Eps;
        }
      }
    ivIt.Set(vec);
    }
  }

  MatrixType estimatedGradients = this->m_NewGradients;
  // this matlab code doesn't seem to do anything:
  // n0 = size(new_gradients,1);
  // new_gradients=new_gradients(1:n0,:);

  DWIVectorImageType::SizeType DWIVecImageSize =
    this->m_IntensityData->GetLargestPossibleRegion().GetSize();
  unsigned numGradientVoxels =
    DWIVecImageSize[0]
    * DWIVecImageSize[1]
    * DWIVecImageSize[2];

  MatrixType icos3;
  Icosahedron3(icos3);
  MatrixType D0(3,3);
  D0 = MatrixType::Zero(3,3);
  D0(0,0) = 1.e-6 * 300;
  D0(1,1) = 1.e-6 * 300;
  D0(2,2) = 1.e-6 * 1700;

  double rho;
  double p;

  const unsigned J(2);
  this->OptimBand(m_BValue,D0,icos3,rho,p);
  // rho=rho*2^(J*p);
  rho = rho * std::pow(2,(J * p));
  // psi=buildridges(J,rho,1,p,0);
  MatrixType psi = BuildRidges(J,rho,1,p,0);

  std::vector<MatrixType> v;
  unsigned long M;
  // [v,M]=multisample(J);
  MultiSample(J,v,M);
  // A=buildsensor(gradientDirections,v,psi);
  MatrixType A = BuildSensor(this->m_VoxelLatticeAlignedGradientDirections,v,psi);
  // A0=buildsensor(new_gradients,v,psi);
  MatrixType A0 = BuildSensor(this->m_NewGradients,v,psi);
  //
  // %parameter setup
  // lmd=0.06;                   % Lagrangian L1
  const double lmd(0.06);
  // myu=0.05;                   % Lagrangian TV
  const double myu(0.05);
  // gama=0.5;                   % Bregman parameter
  const double gama(0.5);
  // NIT=2000;                   % Max number of FISTA iterations
  const unsigned NIT(2000);
  // tol=1e-3;                   % Relative change in Chambolle
  const double tol(1.0e-3);
  //
  // id=find(mask~=0);
  std::vector<unsigned char> id;
  MaskImageType::SizeType maskSize = this->m_MaskImage->GetLargestPossibleRegion().GetSize();
  MaskImageType::IndexType maskIndex;
  for(unsigned int z=0; z < maskSize[2]; ++z)
    {
    maskIndex[2] = z;
    for(unsigned int x=0; x < maskSize[0]; ++x)
      {
      maskIndex[0] = x;
      for(unsigned int y=0; y < maskSize[1]; ++y)
        {
        maskIndex[1] = y;
        id.push_back(this->m_MaskImage->GetPixel(maskIndex));
        }
      }
    }

  //
  // matlabpool('open');
  //
  // u=step2(DWIIntensityData,myu,tol);
  DWIVectorImageType::Pointer u;
  this->step2(this->m_IntensityData,myu,tol,u);
  u->SetMetaDataDictionary(this->m_NrrdFile.GetImage()->GetMetaDataDictionary());
  WriteImage<DWIVectorImageType>(u,"/scratch/kent/ukf/build/UKFTractography-build/CompressedSensing/SmoothedImage.nrrd");

  // c=step1(DWIIntensityData,A,lmd,NIT,id);   % initialization of ridgelet coefficients
  MatrixType c = step1(this->m_IntensityData,A,lmd,NIT,id);
  //
  // Ac=reshape(reshape(c,[numGradientVoxels M])*A',[nx ny nz numGradientDirections]);
  // NOTE: The inner 'reshape' call makes no sense, since the returned
  // Matrix already has the same shape.
  DWIVectorImageType::Pointer Ac;
  {
  MatrixType ctmp = c * A.transpose();
  this->Reshape(ctmp,this->m_IntensityData,Ac);
  }
  // p=Ac-u;
  typedef itk::SubtractImageFilter<DWIVectorImageType,DWIVectorImageType,DWIVectorImageType> SubtractFilterType;
  SubtractFilterType::Pointer subFilter = SubtractFilterType::New();
  subFilter->SetInput1(Ac);
  subFilter->SetInput2(u);
  subFilter->Update();
  DWIVectorImageType::Pointer P = subFilter->GetOutput();
  //
  //
  // TNIT=3;                     % number of outer iterations
  // for itr=1:TNIT,
  const int TNIT(3);
  DWIVectorImageType::Pointer estimatedSignal;
  for(unsigned int itr = 0; itr < TNIT; ++itr)
    {
    typedef itk::AddImageFilter<DWIVectorImageType,DWIVectorImageType,DWIVectorImageType> AddImageFilterType;
    typedef itk::SubtractImageFilter<DWIVectorImageType,DWIVectorImageType,DWIVectorImageType> SubtractImageFilterType;
    // fprintf(1,'Iteration %d of %d\t',[itr TNIT]);
    std::cerr << "Iteration "
              << itr+1
              << " of "
              << TNIT << std::endl;
    //
    // t=u-p;
    subFilter->SetInput1(u);
    subFilter->SetInput2(P);
    subFilter->Update();
    DWIVectorImageType::Pointer t = subFilter->GetOutput();
    // c=step1(t,A,lmd/gama,NIT,id);
    c = this->step1(t,A,lmd/gama,NIT,id);
    // Ac=reshape(reshape(c,[numGradientVoxels M])*A',[nx ny nz numGradientDirections]);
    // NOTE - the inner reshape makes no sense, since c is already
    // shaped that way.
    {
    MatrixType ctmp = c * A.transpose();
    this->Reshape(ctmp,this->m_IntensityData,Ac);
    }
    //
    // t=(1/(1+gama))*(DWIIntensityData+gama*(Ac+p));
    {
    typedef itk::Image<DWIVectorImageType::InternalPixelType,
                       DWIVectorImageType::ImageDimension> ScalarImageType;
    typedef itk::MultiplyImageFilter<DWIVectorImageType,ScalarImageType,DWIVectorImageType>
      MultScalarFilterType;
    // Ac + p
    AddImageFilterType::Pointer add1 = AddImageFilterType::New();
    add1->SetInput1(Ac);
    add1->SetInput2(P);
    // gama * (Ac + p)
    MultScalarFilterType::Pointer mult1 = MultScalarFilterType::New();
    mult1->SetInput1(add1->GetOutput());
    mult1->SetConstant(gama);
    // DWIIntensityData + gama * (Ac + p)
    AddImageFilterType::Pointer add2 = AddImageFilterType::New();
    add2->SetInput1(this->m_IntensityData);
    add2->SetInput2(mult1->GetOutput());
    add2->Update();
    // t = (1 / (1 + gama)) * (DWIIntensityData + gama * (Ac + p))
    const double invGama = 1.0 / (1.0 * gama);
    MultScalarFilterType::Pointer mult2 = MultScalarFilterType::New();
    mult2->SetInput(add2->GetOutput());
    mult2->SetConstant(invGama);
    mult2->Update();
    t = mult2->GetOutput();
    }
    // u = step2(t,myu/(1+gama),tol);
    this->step2(this->m_IntensityData,myu/(1.0 + gama),tol,u);
    // p=p+(Ac-u);
    {
    SubtractFilterType::Pointer subFilter = SubtractFilterType::New();
    AddImageFilterType::Pointer addFilter = AddImageFilterType::New();
    subFilter->SetInput1(Ac);
    subFilter->SetInput2(u);
    addFilter->SetInput1(P);
    addFilter->SetInput2(subFilter->GetOutput());
    addFilter->Update();
    P = addFilter->GetOutput();
    }
    //
    // tv=sum(sum(sum(sqrt((Ac(:,[2:ny,ny],:,:)-Ac).^2+(Ac([2:nx,nx],:,:,:)-Ac).^2+(Ac(:,:,[2:nz,nz],:)-Ac).^2))));
    // f=0.5*sum((Ac(:)-DWIIntensityData(:)).^2)+lmd*sum(abs(c(:)))+myu*sum(tv);
    //
    // fprintf(1,'Cost = %f\n',f);
    // end
    //
    // matlabpool('close');
    //
    // estimatedSignal=reshape(reshape(c,[numGradientVoxels M])*A0',[nx ny nz n0]);   %estimated signal
    {
    MatrixType ctmp = c * A0.transpose();
    this->Reshape(ctmp,this->m_IntensityData,estimatedSignal);
    }
    //
    // % Scale up  gradient direction data by averagedB0 data
    // estimatedSignal = estimatedSignal.* averagedB0(:,:,:,ones(1,n0)); %multiply by the B0 image
    {
    // can't multiply vector image by scalar image, so do it the hard
    // way.
    typedef itk::MultiplyImageFilter<DWIVectorImageType,B0AvgImageType,DWIVectorImageType>
      MultScalarFilterType;
    MultScalarFilterType::Pointer mult = MultScalarFilterType::New();
    mult->SetInput1(estimatedSignal);
    mult->SetInput2(this->m_AverageB0);
    mult->Update();
    estimatedSignal = mult->GetOutput();
    }
    }
  // %% Insert the B0 back in
  // estimatedSignal = cat(4,averagedB0,estimatedSignal); %add B0 to the data
  // estimatedGradients=[0 0 0;new_gradients];
  const unsigned newGradCount(estimatedSignal->GetNumberOfComponentsPerPixel() + 1);
  DWIVectorImageType::Pointer newImage;
  AllocVecImage<DWIVectorImageType>(estimatedSignal,newGradCount,newImage);

  itk::ImageRegionConstIterator<DWIVectorImageType> estimIt(estimatedSignal,estimatedSignal->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<B0AvgImageType> b0It(this->m_AverageB0,this->m_AverageB0->GetLargestPossibleRegion());
  itk::ImageRegionIterator<DWIVectorImageType> newIt(newImage,newImage->GetLargestPossibleRegion());
  for(estimIt.GoToBegin(), b0It.GoToBegin(), newIt.GoToBegin();
      !estimIt.IsAtEnd() && !b0It.IsAtEnd() && newIt.IsAtEnd();
      ++estimIt, ++b0It, ++newIt)
    {
    DWIVectorImageType::PixelType curEstim = estimIt.Get();
    DWIVectorImageType::PixelType newVec;
    newVec.SetSize(curEstim.GetSize()+1);
    newVec[0] = b0It.Get();
    for(unsigned int i = 1; i < newGradCount; ++i)
      {
      newVec[i] = curEstim[i-1];
      }
    }
  this->m_NrrdFile.SetImage(newImage);
  MatrixType newGrad(estimatedGradients.rows()+1,3);
  newGrad(0,0) = newGrad(0,1) = newGrad(0,2) = 0.0;;
  newGrad.block(1,0,estimatedGradients.rows(), estimatedGradients.cols()) = estimatedGradients;
  this->m_NrrdFile.SetGradients(newGrad);
  return true;
}

void DoCSEstimate
::OptimBand(double B,const MatrixType &D,const MatrixType &ico3,double  &rho, double  &p)
{
  MatrixType u = ico3;

  MatrixType u2 = u*D;
  MatrixType u3(u.rows(),u.cols());

  MatrixType S(u.rows(),1);

  // S=exp(-b*sum((u*D).*u,2));
  for(unsigned i = 0; i < u.rows(); ++i)
    {
    for(unsigned j = 0; j < u.cols(); ++j)
      {
      u3(i,j) = u2(i,j) * u(i,j);
      }
    }
  for(unsigned i = 0; i < u.rows(); ++i)
    {
    S(i) = std::exp(-B * u3.row(i).sum());
    }

  const unsigned int m(22);
  // Y=spheriharms(u,m);
  MatrixType Y = SphericalHarmonics(u,22); // KW checked
  // DumpToFile(Y,"/scratch/kent/ukf/build/UKFTractography-build/CompressedSensing/YMatrix.csv");

  // C=(2*(0:m)'+1)/(4*pi);
  MatrixType C(m+1,1);
  for(unsigned int i = 0; i < m+1; ++i)
    {
    C(i) = (2.0 * static_cast<double>(i) + 1.0);
    C(i) /= 4.0 * this->m_Pi;
    }
  // KW Checked C

  //Lmd=ones(m+1,1);
  MatrixType Lmd(m+1,1);
  Lmd = MatrixType::Ones(m+1,1);
  //Lmd(2:2:m+1)=0;
  for(unsigned int i = 1; i < m+1; i += 2)
    {
    Lmd(i) = 0;
    }
  // NOTE: below computation uses loop index
  // K so when used as index needs to be adjusted down.
  //for k=2:2:m,  KW Checked
  for(unsigned int k = 2; k <= m; k += 2)
    {
    // Lmd(k+1)=-Lmd(k-1)*(k-1)/k;
    Lmd(k) = -Lmd(k - 2) * (k - 1)/k;
    }

  std::vector<unsigned> ind(m+1);

  //ind=[1; ((1:m).^2)'+(2:m+1)'];

  ind[0] = 0; // KW Ind checked
  for(unsigned int i = 1; i < m+1; ++i)
    {
    ind[i] = ((i * i) + ( i + 1)) - 1;
    }
  // s=(Y'*Y)\(Y'*S);
  MatrixType s = Y.transpose() * Y;
  MatrixType s2 = Y.transpose() * S;
  s = s.ldlt().solve(s2);
//  MatrixType s = (Y.transpose() * Y).ldlt().solve(Y.transpose() * S);
  // h=(s(ind).*sqrt(1./C));
  MatrixType h(m+1,1);
  for(unsigned int i = 0; i < m+1; ++i)
    {
    h(i,0) = s(ind[i]) * sqrt(1.0/C(i));
    }

  // phi=zeros(m+1,1);
  MatrixType phi;
  phi = MatrixType::Zero(m+1,1);
  // phi(1:2:end)=h(1:2:end)./Lmd(1:2:end);
  for(unsigned i = 0; i < m+1; i += 2)
    {
    phi(i) = h(i) / Lmd(i);
    }
  // phi=phi./phi(1);
  phi /= phi(0);

  MatrixType x(m/2,1);
  // x=(2:2:m)';
  for(unsigned i = 0; i < m/2; ++i)
    {
    x(i) = 2 * (i+1);
    }
  // yy=log(-log(phi(3:2:end)));
  MatrixType yy(m/2,1);
  for(unsigned int i = 0, j = 2; i < m / 2; ++i, j += 2)
    {
    yy(i,0) = std::log(-std::log(phi(j)));
    }
  // X=[ones(length(x),1) log(x)];
  MatrixType X;
  X = MatrixType::Ones(m/2,2);
  for(unsigned int i = 0; i < m/2; ++i)
    {
    X(i,0) = 1.0;
    X(i,1) = log(x(i));
    }
  // alpha=(X'*X)\(X'*yy);
  MatrixType alpha1 = X.transpose() * X;
  MatrixType alpha2 = X.transpose() * yy;
  Eigen::FullPivLU<MatrixType> lu2(alpha1);
  MatrixType alpha = lu2.solve(alpha2);
  // rho=exp(alpha(1));
  rho = std::exp(alpha(0));
  // p=alpha(2);
  p = alpha(1);
}

/** use an array of indices to shift voxels in a particular direction */
template <typename ImageType>
ImageType *shiftImage(const typename ImageType::Pointer &input,
                      unsigned int dimension,
                      const std::vector<unsigned int> &indices)
{
  typename ImageType::Pointer rval = AllocImage<ImageType>(input);
  itk::ImageRegionConstIteratorWithIndex<ImageType> fromIt(input,input->GetLargestPossibleRegion());
  for(fromIt.GoToBegin(); !fromIt.IsAtEnd(); ++fromIt)
    {
    typename ImageType::IndexType index = fromIt.GetIndex();
    index[dimension] = indices[index[dimension]];
    rval->SetPixel(index,fromIt.Get());
    }
  rval.GetPointer()->Register();
  return rval.GetPointer();
}

void
DoCSEstimate
::ToVecImage(const B0AvgImageType *fromImage, unsigned int gradientIndex, DWIVectorImageType::Pointer &target)
{
  itk::ImageRegionConstIterator<B0AvgImageType>
    uIt(fromImage, fromImage->GetLargestPossibleRegion());
  itk::ImageRegionIterator<DWIVectorImageType>
    dwiIt(target,target->GetLargestPossibleRegion());
  for(uIt.GoToBegin(), dwiIt.GoToBegin(); !uIt.IsAtEnd() && !dwiIt.IsAtEnd(); ++uIt, ++dwiIt)
    {
    DWIVectorImageType::PixelType pixel = dwiIt.Get();
    pixel[gradientIndex] = uIt.Get();
    dwiIt.Set(pixel);
    }
}

DoCSEstimate::B0AvgImageType *
DoCSEstimate
::FromVecImage(DWIVectorImageType::Pointer inputImage,unsigned gradientIndex)
{
  B0AvgImageType::Pointer f;
  AllocImage<B0AvgImageType>(inputImage,f);

  itk::ImageRegionConstIterator<DWIVectorImageType>
    dwiIt(inputImage,inputImage->GetLargestPossibleRegion());
  itk::ImageRegionIterator<B0AvgImageType> fIt(f,f->GetLargestPossibleRegion());
  for(dwiIt.GoToBegin(), fIt.GoToBegin(); !dwiIt.IsAtEnd() && !fIt.IsAtEnd(); ++dwiIt, ++fIt)
    {
    DWIVectorImageType::PixelType curVec = dwiIt.Get();
    fIt.Set(curVec[gradientIndex]);
    }
  f.GetPointer()->Register();
  return f.GetPointer();
}



/** Denoise image
 *  this method has to behave very differently than matlab.
 *  for one thing, Eigen only does 2D matrices, so we use ITK images.
 *  for another, matlab expects a sequence of volumes, one per
 *  gradient, and the itk gradient image is a vector image, and each
 *  voxel has one scalar per gradient.
 */
void
DoCSEstimate
::tvdenoise3(DWIVectorImageType::Pointer &inputImage,
             unsigned int gradientIndex,
             double lambda,double tol,DWIVectorImageType::Pointer &target)
{
  typedef B0AvgImageType TImage;
  TImage::Pointer f = this->FromVecImage(inputImage,gradientIndex);
#if 0
  typedef itk::MedianImageFilter<TImage,TImage> MedianFilterType;
  MedianFilterType::Pointer medianFilter = MedianFilterType::New();
  MedianFilterType::InputSizeType radius;
  radius.Fill(1);
  medianFilter->SetRadius(radius);
  medianFilter->SetInput(f);
  medianFilter->Update();
  this->ToVecImage(medianFilter->GetOutput(),gradientIndex,target);
#else
// if (nargin < 3),
//     Tol = 1e-2;
// end
//
// if (lambda < 0),
//     error('Parameter lambda must be nonnegative.');
// end
//
// dt = 0.25/2;
  double dt = 0.25/2;
//
// N = size(f);
  const TImage::SizeType N =
    inputImage->GetLargestPossibleRegion().GetSize();
// id = [2:N(1),N(1)];
  ShiftVectorType id = MakeShiftVector(N[0],1);
// iu = [1,1:N(1)-1];
  ShiftVectorType iu = MakeShiftVector(N[0],-1);
// ir = [2:N(2),N(2)];
  ShiftVectorType ir = MakeShiftVector(N[1],1);
// il = [1,1:N(2)-1];
  ShiftVectorType il = MakeShiftVector(N[1],-1);
// ib = [2:N(3),N(3)];
  ShiftVectorType ib = MakeShiftVector(N[2],1);
// ifr = [1,1:N(3)-1];
  ShiftVectorType ifr = MakeShiftVector(N[2],-1);
//
// p1 = zeros(size(f));
  TImage::Pointer p1;
  AllocImage<TImage>(f,0.0,p1);
// p2 = zeros(size(f));
  TImage::Pointer p2;
  AllocImage<TImage>(f,0.0,p2);
// p3 = zeros(size(f));
  TImage::Pointer p3;
  AllocImage<TImage>(f,0.0,p3);
//
// divp = zeros(size(f));
  TImage::Pointer divp;
  AllocImage<TImage>(f,0.0,divp);
// lastdivp = ones(size(f));
  TImage::Pointer lastdivp;
  AllocImage<TImage>(f,1.0,lastdivp);
//
  TImage::Pointer fLambda;
  MultiplyImage<TImage>(f,lambda,fLambda);
// if (length(N) == 3),
//     while (norm(divp(:) - lastdivp(:),inf) > Tol),
  for(double normdivp = 1.0; normdivp > tol; )
    {
//         lastdivp = divp;
//         DO below to avod repeated allocations
//         z = divp - f*lambda;
    TImage::Pointer z;
    SubtractImage<TImage>(divp,fLambda,z);
    //
    // avoid re-allocating images
    TImage::Pointer tmpdivp = lastdivp;
    lastdivp = divp;
    divp = tmpdivp;

    itk::ImageRegionIterator<TImage> p1It(p1,p1->GetLargestPossibleRegion());
    itk::ImageRegionIterator<TImage> p2It(p2,p2->GetLargestPossibleRegion());
    itk::ImageRegionIterator<TImage> p3It(p3,p3->GetLargestPossibleRegion());
    itk::ImageRegionConstIteratorWithIndex<TImage>
      zIt(z,z->GetLargestPossibleRegion());

    for(p1It.GoToBegin(), p2It.GoToBegin(), p3It.GoToBegin(),zIt.GoToBegin();
        !p1It.IsAtEnd();
        ++p1It, ++p2It, ++p3It, ++zIt)
      {
//         z1 = z(:,ir,:) - z;
//         z2 = z(id,:,:) - z;
//         z3 = z(:,:,ib) - z;
      const double zVal = zIt.Get();
      TImage::IndexType zIndex;
      zIndex = zIt.GetIndex();
      zIndex[1] = ir[zIndex[1]];
      double z1val = z->GetPixel(zIndex) - zVal;
      zIndex = zIt.GetIndex();
      zIndex[0] = id[zIndex[0]];
      double z2val = z->GetPixel(zIndex) - zVal;
      zIndex = zIt.GetIndex();
      zIndex[2] = ib[zIndex[2]];
      double z3val = z->GetPixel(zIndex) - zVal;
//    denom = 1 + dt*sqrt(z1.^2 + z2.^2 + z3.^2);
      const double denom = 1 + dt * std::sqrt( (z1val*z1val) + (z2val*z2val) + (z3val*z3val));
      const double p1val(p1It.Get()), p2val(p2It.Get()), p3val(p3It.Get());
//         p1 = (p1 + dt*z1)./denom;
      p1It.Set((p1val + dt * z1val) / denom);
//         p2 = (p2 + dt*z2)./denom;
      p2It.Set((p2val + dt * z2val) / denom);
//         p3 = (p3 + dt*z3)./denom;
      p3It.Set((p3val + dt * z3val) / denom);
      }
//         divp = p1 - p1(:,il,:) + p2 - p2(iu,:,:) + p3 - p3(:,:,ifr);
    //  image allocation happens once
    //  divp = AllocImage<TImage>(f);
    itk::ImageRegionIteratorWithIndex<TImage> divpIt(divp,divp->GetLargestPossibleRegion());
    for(divpIt.GoToBegin(), p1It.GoToBegin(), p2It.GoToBegin(), p3It.GoToBegin();
        !p1It.IsAtEnd();
        ++divpIt, ++p1It, ++p2It, ++p3It)
      {
      const TImage::IndexType curIndex = divpIt.GetIndex();
      TImage::IndexType shiftIndex;
      shiftIndex = curIndex;
      shiftIndex[2] = il[curIndex[2]];
      const double p1val = p1It.Get() - p1->GetPixel(shiftIndex);
      shiftIndex = curIndex;
      shiftIndex[0] = iu[curIndex[0]];
      const double p2val = p2It.Get() - p2->GetPixel(shiftIndex);
      shiftIndex = curIndex;
      shiftIndex[2] = ifr[curIndex[2]];
      const double p3val = p3It.Get() - p3->GetPixel(shiftIndex);
      divpIt.Set(p1val + p2val + p3val);
      }
//     end
    normdivp = PNormImageDiff<TImage>(divp,lastdivp);
    }
// end
//
// u = f - divp/lambda;
  typedef itk::MultiplyImageFilter<TImage,TImage,TImage>  MultiplyFilterType;
  typedef itk::SubtractImageFilter<TImage,TImage,TImage>  SubtractFilterType;
  MultiplyFilterType::Pointer mult = MultiplyFilterType::New();
  SubtractFilterType::Pointer sub = SubtractFilterType::New();
  mult->SetInput1(divp);
  mult->SetConstant(1.0/lambda);
  sub->SetInput1(f);
  sub->SetInput2(mult->GetOutput());
  sub->Update();
  TImage::Pointer u = sub->GetOutput();
  this->ToVecImage(u,gradientIndex,target);
#endif
}

// function u = tvdenoise3(f,lambda,Tol)
// %This function extends the chambolle's algo to 3D images
// %TVDENOISE  Total variation grayscale and color image denoising
// %   u = TVDENOISE(f,lambda) denoises the input image f.  The smaller
// %   the parameter lambda, the stronger the denoising.
// %
// %   The output u approximately minimizes the Rudin-Osher-Fatemi (ROF)
// %   denoising model
// %
// %       Min  TV(u) + lambda/2 || f - u ||^2_2,
// %        u
// %
// %   where TV(u) is the total variation of u.  If f is a color image (or any
// %   array where size(f,3) > 1), the vectorial TV model is used,
// %
// %       Min  VTV(u) + lambda/2 || f - u ||^2_2.
// %        u
// %
// %   TVDENOISE(...,Tol) specifies the stopping tolerance (default 1e-2).
// %
// %   The minimization is solved using Chambolle's method,
// %      A. Chambolle, "An Algorithm for Total Variation Minimization and
// %      Applications," J. Math. Imaging and Vision 20 (1-2): 89-97, 2004.
// %   When f is a color image, the minimization is solved by a generalization
// %   of Chambolle's method,
// %      X. Bresson and T.F. Chan,  "Fast Minimization of the Vectorial Total
// %      Variation Norm and Applications to Color Image Processing", UCLA CAM
// %      Report 07-25.
// %
// %   Example:
// %   f = double(imread('barbara-color.png'))/255;
// %   f = f + randn(size(f))*16/255;
// %   u = tvdenoise(f,12);
// %   subplot(1,2,1); imshow(f); title Input
// %   subplot(1,2,2); imshow(u); title Denoised
//
// % Pascal Getreuer 2007-2008
//
// if (nargin < 3),
//     Tol = 1e-2;
// end
//
// if (lambda < 0),
//     error('Parameter lambda must be nonnegative.');
// end
//
// dt = 0.25/2;
//
// N = size(f);
// id = [2:N(1),N(1)];
// iu = [1,1:N(1)-1];
// ir = [2:N(2),N(2)];
// il = [1,1:N(2)-1];
// ib = [2:N(3),N(3)];
// ifr = [1,1:N(3)-1];
//
// p1 = zeros(size(f));
// p2 = zeros(size(f));
// p3 = zeros(size(f));
//
// divp = zeros(size(f));
// lastdivp = ones(size(f));
//
// if (length(N) == 3),
//     while (norm(divp(:) - lastdivp(:),inf) > Tol),
//         lastdivp = divp;
//         z = divp - f*lambda;
//         z1 = z(:,ir,:) - z;
//         z2 = z(id,:,:) - z;
//         z3 = z(:,:,ib) - z;
//         denom = 1 + dt*sqrt(z1.^2 + z2.^2 + z3.^2);
//         p1 = (p1 + dt*z1)./denom;
//         p2 = (p2 + dt*z2)./denom;
//         p3 = (p3 + dt*z3)./denom;
//         divp = p1 - p1(:,il,:) + p2 - p2(iu,:,:) + p3 - p3(:,:,ifr);
//     end
// end
//
// u = f - divp/lambda;

void
DoCSEstimate
::step2(DWIVectorImageType::Pointer &inputImage,double myu, double tol,DWIVectorImageType::Pointer &rval)
{
  std::cerr << "Step 2" << std::endl;
  unsigned int numGradients = inputImage->GetNumberOfComponentsPerPixel();
 AllocVecImage<DWIVectorImageType>(inputImage,numGradients,rval);

  typedef DWIVectorImageType::SizeValueType sizetype;
  for(unsigned int k = 0; k < numGradients; ++k)
    {
    std::cerr << "Gradient " << k << std::endl;
    this->tvdenoise3(inputImage,k,1/myu,tol,rval);
    }
}

// function [u] = step2(u,myu,tol)
//
// r=1/myu;
//
// parfor k=1:size(u,4),
//    u(:,:,:,k)=tvdenoise3(u(:,:,:,k),r,tol);
// end

// function [X] = step1(S,A,lmd,NIT,id)
// NOTE the S parameter is the intensity data.
MatrixType
DoCSEstimate
::step1(DWIVectorImageType::Pointer &inputImage,MatrixType &A,double lmd,unsigned NIT,std::vector<unsigned char> &id)
{
  std::cerr << "Step 1" << std::endl;
  // [nx, ny, nz, d]=size(S);
  const ImageSizeType size = inputImage->GetLargestPossibleRegion().GetSize();
  // n=nx*ny*nz;
  ImageSizeValueType n = size[0] * size[1] * size[2];
  // M=size(A,2);
  ImageSizeValueType M = A.cols();

  unsigned int numGradients = A.rows();
  // the reshape takes image of size [ x y z g ] and turns it into
  // [n g ] where n = x * y * z.  Matlab reshapes columnwise, so you
  // end up with G columns of n rows.
  MatrixType _S;
  this->Reshape(inputImage,n,numGradients,_S);
  // S=S(id,:);
  MatrixType S(id.size(),numGradients);
  for(unsigned int row = 0; row < id.size(); ++row)
    {
    S.row(row) = _S.row(id[row]);
    }
  // x=zeros(size(S,1),M);
  MatrixType x; x = MatrixType::Zero(S.rows(),M);
  //
  // parfor i=1:size(S,1),
  for(unsigned i = 0; i < S.rows(); ++i)
    {
    //     x(i,:) = BPDN_homotopy_function(A, squeeze(S(i,:))', lmd, NIT);
    MatrixType squeezeS = S.row(i).transpose();
    x.row(i) = BPDN_HOMOTOPY_function(A,squeezeS,lmd,NIT);
    // end
    }
  // %matlabpool close;
  //
  // X=zeros(n,M);
  // X(id,:)=x;
  MatrixType rval; rval = MatrixType::Zero(n,M);
  for(unsigned int grad = 0; grad < numGradients; ++grad)
    {
    for(unsigned int row = 0; row < id.size(); ++row)
      {
      rval(id[row],grad) = x(row,grad);
      }
    }
  return rval;
}
// function [X] = step1(S,A,lmd,NIT,id)
//
// [nx, ny, nz, d]=size(S);
// n=nx*ny*nz;
// M=size(A,2);
//
// S=reshape(S,[n d]);
// S=S(id,:);
// x=zeros(size(S,1),M);
//
// %parallel operation possible here
// %matlabpool(12);
//
// parfor i=1:size(S,1),
//     x(i,:) = BPDN_homotopy_function(A, squeeze(S(i,:))', lmd, NIT);
// end
// %matlabpool close;
//
// X=zeros(n,M);
// X(id,:)=x;
//
// end

