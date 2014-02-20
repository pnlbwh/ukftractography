#include "DoCSEstimate.h"
#include "BuildRidges.h"
#include "BuildSensor.h"
#include "MultiSample.h"
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
#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
template <class ImageType>
void
WriteImage(const ImageType *image,
                const std::string & filename)
{
  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename  WriterType::Pointer writer = WriterType::New();
  writer->UseCompressionOn();
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
}
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
  outMatrix.conservativeResize(rows,columns);
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
  outImage = DWIVectorImageType::New();
  outImage->CopyInformation(templateImage);
  outImage->SetNumberOfComponentsPerPixel(numComponents);
  outImage->Allocate();
  DWIVectorImageType::SizeType size = templateImage->GetLargestPossibleRegion().GetSize();
  DWIVectorImageType::IndexType index;
  for(unsigned int z = 0; z < size[2]; ++z)
    {
    index[2] = z;
    for(unsigned int y = 0; y < size[1]; ++y)
      {
      index[1] = y;
      for(unsigned int x = 0; x < size[0]; ++x)
        {
        index[0] = x;
        DWIVectorImageType::PixelType curPixel;
        unsigned int row = (z * index[2]) + (x * index[0]) + y;
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
  this->m_AverageB0 = B0AvgImageType::New();
  this->m_AverageB0->SetRegions(region);
  this->m_AverageB0->Allocate();
  this->m_AverageB0->FillBuffer(0.0);

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
  this->m_IntensityData = DWIVectorImageType::New();
  this->m_IntensityData->CopyInformation(originalNrrd);
  this->m_IntensityData->SetRegions(region);
  this->m_IntensityData->SetNumberOfComponentsPerPixel(gradientIndices.size());
  this->m_IntensityData->Allocate();

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
  DWIVectorImageType::Pointer u = this->step2(this->m_IntensityData,myu,tol);
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
    typedef itk::MultiplyImageFilter<DWIVectorImageType,DWIVectorImageType,DWIVectorImageType> MultiplyFilterType;
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
    const unsigned int pixSize(Ac->GetNumberOfComponentsPerPixel());
    DWIVectorImageType::PixelType gamaVec;
    DWIVectorImageType::PixelType invGamaVec;
    gamaVec.SetSize(pixSize);
    invGamaVec.SetSize(pixSize);
    for(unsigned int i = 0; i < pixSize; ++i)
      {
      gamaVec[i] = gama;
      invGamaVec[i] = 1.0/(1.0 + gama);
      }
    AddImageFilterType::Pointer addFilter1 = AddImageFilterType::New();
    AddImageFilterType::Pointer addFilter2 = AddImageFilterType::New();
    MultiplyFilterType::Pointer multFilter1 = MultiplyFilterType::New();
    MultiplyFilterType::Pointer multFilter2 = MultiplyFilterType::New();
    addFilter1->SetInput1(Ac);
    addFilter1->SetInput1(P);
    multFilter1->SetInput(addFilter1->GetOutput());
    multFilter1->SetConstant(gamaVec);
    addFilter2->SetInput1(multFilter1->GetOutput());
    addFilter2->SetInput2(this->m_IntensityData);
    multFilter2->SetInput(addFilter2->GetOutput());
    multFilter2->SetConstant(invGamaVec);
    multFilter2->Update();
    t = multFilter2->GetOutput();
    }
    // u=step2(t,myu/(1+gama),tol);
    u = this->step2(t,myu/(1+gama),tol);
    //
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
    itk::ImageRegionConstIterator<B0AvgImageType> b0It(this->m_AverageB0,this->m_AverageB0->GetLargestPossibleRegion());
    itk::ImageRegionIterator<DWIVectorImageType> estIt(estimatedSignal,estimatedSignal->GetLargestPossibleRegion());
    for(b0It.GoToBegin(),estIt.GoToBegin(); !b0It.IsAtEnd() && !estIt.IsAtEnd(); ++b0It, ++estIt)
      {
      estIt.Set(estIt.Get() * b0It.Get());
      }
    }
    }
  // %% Insert the B0 back in
  // estimatedSignal = cat(4,averagedB0,estimatedSignal); %add B0 to the data
  // estimatedGradients=[0 0 0;new_gradients];
  DWIVectorImageType::Pointer newImage = DWIVectorImageType::New();
  newImage->CopyInformation(estimatedSignal);
  const unsigned newGradCount(estimatedSignal->GetNumberOfComponentsPerPixel() + 1);
  newImage->SetNumberOfComponentsPerPixel(newGradCount);
  newImage->Allocate();
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
  MatrixType newGrad(estimatedGradients.rows(),estimatedGradients.cols()+1);
  newGrad(0,0) = newGrad(1,0) = newGrad(2,0) = 0.0;
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
  typename ImageType::Pointer rval = ImageType::New();
  rval->CopyInformation(input);
  rval->SetRegions(input->GetLargestPossibleRegion());
  rval->Allocate();
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
  B0AvgImageType::Pointer f = B0AvgImageType::New();
  f->CopyInformation(inputImage);
  f->SetRegions(inputImage->GetLargestPossibleRegion());
  f->Allocate();

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

namespace
{
template <typename TImage>
TImage *SubtractImage(const TImage *t1, const TImage *t2)
{
  typedef itk::SubtractImageFilter<TImage,TImage,TImage>  SubtractFilterType;
  typename SubtractFilterType::Pointer sub1 = SubtractFilterType::New();
  sub1->SetInput1(t1);
  sub1->SetInput2(t2);
  sub1->Update();
  typename TImage::Pointer rval = sub1->GetOutput();
  rval.GetPointer()->Register();
  return rval.GetPointer();
}

template <typename TImage>
TImage *AddImage(const TImage *t1, const TImage *t2)
{
  typedef itk::AddImageFilter<TImage,TImage,TImage>  AddFilterType;
  typename AddFilterType::Pointer sub1 = AddFilterType::New();
  sub1->SetInput1(t1);
  sub1->SetInput2(t2);
  sub1->Update();
  typename TImage::Pointer rval = sub1->GetOutput();
  rval.GetPointer()->Register();
  return rval.GetPointer();
}

template <typename TImage>
TImage *MultiplyImage(const TImage *t1, const TImage *t2)
{
  typedef itk::MultiplyImageFilter<TImage,TImage,TImage>  MultiplyFilterType;
  typename MultiplyFilterType::Pointer sub1 = MultiplyFilterType::New();
  sub1->SetInput1(t1);
  sub1->SetInput2(t2);
  sub1->Update();
  typename TImage::Pointer rval = sub1->GetOutput();
  rval.GetPointer()->Register();
  return rval.GetPointer();
}

template <typename TImage>
TImage *MultiplyImage(const TImage *t1, typename TImage::PixelType constant)
{
  typedef itk::MultiplyImageFilter<TImage,TImage,TImage>  MultiplyFilterType;
  typename MultiplyFilterType::Pointer sub1 = MultiplyFilterType::New();
  sub1->SetInput1(t1);
  sub1->SetConstant(constant);
  sub1->Update();
  typename TImage::Pointer rval = sub1->GetOutput();
  rval.GetPointer()->Register();
  return rval.GetPointer();
}

template <typename TImage>
TImage *AllocImage(const itk::ImageBase<TImage::ImageDimension> *templateImage)
{
  typename TImage::Pointer p1 = TImage::New();
  p1->CopyInformation(templateImage);
  p1->SetRegions(templateImage->GetLargestPossibleRegion());
  p1->Allocate();
  p1.GetPointer()->Register();
  return p1.GetPointer();
}
template <typename TImage>
TImage *AllocImage(const itk::ImageBase<TImage::ImageDimension> *templateImage, typename TImage::PixelType initval)
{
  TImage *rval = AllocImage<TImage>(templateImage);
  rval->FillBuffer(initval);
  return rval;
}

template <typename TImage>
double PNormImageDiff(const TImage *im1, const TImage *im2)
{
  itk::ImageRegionConstIterator<TImage> im1It(im1,im1->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<TImage> im2It(im2,im2->GetLargestPossibleRegion());
  typename TImage::PixelType maxPixel = 0;
  for(im1It.GoToBegin(), im2It.GoToBegin(); !im1It.IsAtEnd(); ++im1It,++im2It)
    {
    const typename TImage::PixelType current = std::abs(im1It.Get() - im2It.Get());
    if(current > maxPixel)
      {
      maxPixel = current;
      }
    }
  return maxPixel;
}

}

template<typename TImage>
TImage *SubtractShift(const TImage *inImage,unsigned axis,std::vector<unsigned int> &shift)
{
  TImage *rval = AllocImage<TImage>(inImage);
  itk::ImageRegionIterator<TImage> rvalIt(rval,rval->GetLargestPossibleRegion());
  itk::ImageRegionConstIteratorWithIndex<TImage> inIt(inImage,inImage->GetLargestPossibleRegion());

  for(rvalIt.GoToBegin(), inIt.GoToBegin(); !rvalIt.IsAtEnd(); ++rvalIt, ++inIt)
    {
    typename TImage::IndexType index = inIt.GetIndex();
    index[axis] = shift[index[axis]];
    rvalIt.Set(inIt.Get() - inImage->GetPixel(index));
    }
  return rval;
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
  typedef itk::SubtractImageFilter<TImage,TImage,TImage>  SubtractFilterType;
  typedef itk::AbsImageFilter<TImage,TImage>                      AbsImageFilterType;
  typedef itk::StatisticsImageFilter<TImage>                              StatsFilterType;
  typedef itk::MultiplyImageFilter<TImage,TImage,TImage>  MultiplyImageFilterType;

  TImage::Pointer f = this->FromVecImage(inputImage,gradientIndex);

  MultiplyImageFilterType::Pointer mult = MultiplyImageFilterType::New();
  mult->SetInput1(f);
  mult->SetConstant(lambda);
  mult->Update();
  TImage::Pointer fLambda = mult->GetOutput();

  // dt = 0.25/2;
  const double dt = 0.125;

  // N = size(f);
  const ImageSizeType N = inputImage->GetLargestPossibleRegion().GetSize();

  // id = [2:N(1),N(1)];
  std::vector<unsigned int> id;
  for(unsigned int i = 1; i <= N[0] - 1; ++i)
    {
    id.push_back(i);
    }
  id.push_back(N[0] - 1);

  // iu = [1,1:N(1)-1];
  std::vector<unsigned int> iu;
  iu.push_back(0);
  for(unsigned int i = 0; i < N[0] - 1; ++i)
    {
    iu.push_back(i);
    }

  // ir = [2:N(2),N(2)];
  std::vector<unsigned int> ir;
  for(unsigned int i = 1; i < N[1]; ++i)
    {
    ir.push_back(i);
    }
  ir.push_back(N[1] - 1);

  // il = [1,1:N(2)-1];
  std::vector<unsigned int> il;
  il.push_back(0);
  for(unsigned int i = 0; i < N[1] - 1; ++i)
    {
    il.push_back(i);
    }

  // ib = [2:N(3),N(3)];
  std::vector<unsigned int> ib;
  for(unsigned int i = 1; i < N[2]; ++i)
    {
    ib.push_back(i);
    }
  ib.push_back(N[2] - 1);

  // ifr = [1,1:N(3)-1];
  std::vector<unsigned int> ifr;
  ifr.push_back(0);
  for(unsigned int i = 0; i < N[2] - 1; ++i)
    {
    ifr.push_back(i);
    }

  //
  // p1 = zeros(size(f));
  TImage::Pointer p1 = AllocImage<TImage>(f,0.0);
  // p2 = zeros(size(f));
  TImage::Pointer p2 = AllocImage<TImage>(f,0.0);
  // p3 = zeros(size(f));
  TImage::Pointer p3 = AllocImage<TImage>(f,0.0);
  //
  // divp = zeros(size(f));
  TImage::Pointer divp = AllocImage<TImage>(f,0.0);

  // lastdivp = ones(size(f));
  TImage::Pointer lastdivp;

  double pnorm(1);
//
// if (length(N) == 3), which it always will be.
  while(pnorm > tol)
    {
//     while (norm(divp(:) - lastdivp(:),inf) > Tol),
//         lastdivp = divp;
    // to avoid needless allocation, swap divp with lastdivp
    lastdivp = divp;

    // note: f doesn't change in loop so f*lambda done above and
    //         assigned to fLambda
    //         z = divp - f*lambda;
    SubtractFilterType::Pointer sub1 = SubtractFilterType::New();
    sub1->SetInput1(divp);
    sub1->SetInput2(fLambda);
    sub1->Update();
    TImage::Pointer z = sub1->GetOutput();
    // z1 = z(:,ir,:) - z;
    TImage::Pointer z1 = shiftImage<TImage>(z,1,ir);
    sub1->SetInput1(z1); sub1->SetInput2(z); sub1->Update();
    z1 = sub1->GetOutput();
    // z2 = z(id,:,:) - z;
    TImage::Pointer z2 = shiftImage<TImage>(z,0,id);
    sub1->SetInput1(z2); sub1->SetInput2(z); sub1->Update();
    z2 = sub1->GetOutput();
    // z3 = z(:,:,ib) - z;
    TImage::Pointer z3 = shiftImage<TImage>(z,2,ib);
    sub1->SetInput1(z3); sub1->SetInput2(z); sub1->Update();
    z3 = sub1->GetOutput();

    // denom = 1 + dt*sqrt(z1.^2 + z2.^2 + z3.^2);
    TImage::Pointer denom = AllocImage<TImage>(f);
    itk::ImageRegionIterator<TImage> denomIt(denom,denom->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<TImage> z1It(z1,z1->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<TImage> z2It(z2,z2->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<TImage> z3It(z3,z3->GetLargestPossibleRegion());

    for(denomIt.GoToBegin(), z1It.GoToBegin(),z2It.GoToBegin(),z3It.GoToBegin();
        !denomIt.IsAtEnd() && !z1It.IsAtEnd() && !z2It.IsAtEnd() && !z3It.IsAtEnd();
        ++denomIt, ++z1It, ++z2It, ++z3It)
      {
      TImage::PixelType z1Val = z1It.Get();
      TImage::PixelType z2Val = z2It.Get();
      TImage::PixelType z3Val = z3It.Get();
      denomIt.Set(1 + dt * std::sqrt((z1Val * z1Val) + (z2Val * z2Val)+ (z3Val * z3Val)));
      }
    itk::ImageRegionIterator<TImage> p1It(p1,p1->GetLargestPossibleRegion());
    itk::ImageRegionIterator<TImage> p2It(p2,p2->GetLargestPossibleRegion());
    itk::ImageRegionIterator<TImage> p3It(p3,p3->GetLargestPossibleRegion());

     // p1 = (p1 + dt*z1)./denom;
     // p2 = (p2 + dt*z2)./denom;
     // p3 = (p3 + dt*z3)./denom;
    for(denomIt.GoToBegin(), p1It.GoToBegin(),p2It.GoToBegin(),p3It.GoToBegin(),
          z1It.GoToBegin(), z2It.GoToBegin(), z3It.GoToBegin();
        !denomIt.IsAtEnd() && !p1It.IsAtEnd() && !p2It.IsAtEnd() && !p3It.IsAtEnd();
        ++denomIt, ++p1It, ++p2It, ++p3It, ++z1It, ++z2It, ++z3It)
      {
      double curDenom = denomIt.Get();
      double p1val = (p1It.Get() + dt * z1It.Get())/curDenom;
      double p2val = (p2It.Get() + dt * z2It.Get())/curDenom;
      double p3val = (p3It.Get() + dt * z3It.Get())/curDenom;
      p1It.Set(p1val);
      p2It.Set(p2val);
      p3It.Set(p3val);
      }
    // divp = p1 - p1(:,il,:) + p2 - p2(iu,:,:) + p3 - p3(:,:,ifr);
    divp = AllocImage<TImage>(f,0.0);
    itk::ImageRegionIteratorWithIndex<TImage>
      divpIt(divp,divp->GetLargestPossibleRegion());
    for(divpIt.GoToBegin(); !divpIt.IsAtEnd(); ++divpIt)
      {
      TImage::IndexType curIndex = divpIt.GetIndex();
      TImage::IndexType pIndex;
      double curP1 = p1->GetPixel(curIndex);
      pIndex = curIndex;
      pIndex[1] = il[curIndex[1]];
      curP1 -= p1->GetPixel(pIndex);
      double curP2 = p2->GetPixel(curIndex);
      pIndex = curIndex;
      pIndex[0] = iu[curIndex[0]];
      curP2 -= p2->GetPixel(pIndex);
      double curP3 = p3->GetPixel(curIndex);
      pIndex = curIndex;
      pIndex[2] = ifr[curIndex[2]];
      curP3 -= p3->GetPixel(pIndex);
      divp->SetPixel(curIndex, curP1 + curP2 + curP3);
      }
    /* the test moved to the bottom of the list, since the
     * first time through the loop, the test will always succeed.
     * norm(divp(:) - lastdivp(:),inf) > Tol
     * now according to matlab docs norm(X,p) is the p-norm of X.  And
     * if p = inf, p-norm(X,inf) is the unform norm which is the
     * maximum absolute value.
     * I suppose I could devise a more expensive loop test, but I'd
     * have to really work at it.
     */
    pnorm = PNormImageDiff<TImage>(divp,lastdivp);
    }
  while(pnorm > tol);
// end
//
// u = f - divp/lambda;
  SubtractFilterType::Pointer subFilter = SubtractFilterType::New();
  MultiplyImageFilterType::Pointer multFilter = MultiplyImageFilterType::New();
  multFilter->SetInput(divp);
  multFilter->SetConstant(1.0/lambda);
  subFilter->SetInput1(f);
  subFilter->SetInput2(multFilter->GetOutput());
  subFilter->Update();
  this->ToVecImage(subFilter->GetOutput(),gradientIndex,target);
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

DoCSEstimate::DWIVectorImageType *
DoCSEstimate
::step2(DWIVectorImageType::Pointer &inputImage,double myu, double tol)
{
  unsigned int numGradients = inputImage->GetNumberOfComponentsPerPixel();
  DWIVectorImageType::Pointer rval = DWIVectorImageType::New();
  rval->CopyInformation(inputImage);
  rval->SetRegions(inputImage->GetLargestPossibleRegion());
  rval->SetNumberOfComponentsPerPixel(numGradients);
  rval->Allocate();
  typedef DWIVectorImageType::SizeValueType sizetype;
  for(unsigned int k = 0; k < numGradients; ++k)
    {
    this->tvdenoise3(inputImage,k,1/myu,tol,rval);
    }
  rval.GetPointer()->Register();
  return rval.GetPointer();
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
  // [nx, ny, nz, d]=size(S);
  ImageSizeType size = inputImage->GetLargestPossibleRegion().GetSize();
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
  for(unsigned int grad = 0; grad < numGradients; ++grad)
    {
    for(unsigned int row = 0; row < id.size(); ++row)
      {
      S(row,grad) = _S(id[row],grad);
      }
    }
  // x=zeros(size(S,1),M);
  MatrixType x; x = MatrixType::Zero(S.rows(),M);
  //
  // parfor i=1:size(S,1),
  for(unsigned i = 0; i < S.rows(); ++i)
    {
    //     x(i,:) = BPDN_homotopy_function(A, squeeze(S(i,:))', lmd, NIT);
    MatrixType squeezeS = S.row(i).transpose();
    x.row(i) = this->BPDN_HOMOTOPY_function(A,squeezeS,lmd,NIT);
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

MatrixType
DoCSEstimate
::BPDN_HOMOTOPY_function(MatrixType &A,MatrixType &y,double tau, unsigned int maxiter)
{
// function [x_out, gamma_x, total_iter, total_time] = BPDN_homotopy_function(A, y, tau, maxiter)
//
// t0 = cputime;
//
// N = size(A,2);
  const unsigned int N(A.cols());
// K = size(A,1);
  const unsigned int K(A.rows());
//
// % Initialization of primal and dual sign and support
// z_x = zeros(N,1);
  MatrixType z_x; z_x = MatrixType::Zero(N,1);
// gamma_x = [];       % Primal support
  std::vector<unsigned long> gamma_x;
//
// % Initial step
// Primal_constrk = -A'*y;
  MatrixType Primal_constrk = -A.transpose() * y;
// [c i] = max(abs(Primal_constrk));
  unsigned long irow,ijunk;
  double c = Primal_constrk.cwiseAbs().maxCoeff(&irow,&ijunk);
//
// gamma_xk = i;
  std::vector<unsigned long> gamma_xk(1,irow);
  gamma_xk[0] = irow;
//
// epsilon = c;
  double epsilon(c);
// xk_1 = zeros(N,1);
  MatrixType xk_1; xk_1 = MatrixType::Zero(N,1);
//
// z_x(gamma_xk) = -sign(Primal_constrk(gamma_xk));
  if(Primal_constrk(gamma_xk[0]) > 0.0)
    {
    z_x(gamma_xk[0]) =  -1.0;
    }
  else
    {
    z_x(gamma_xk[0]) =  1.0;
    }

  // Primal_constrk(gamma_xk) = sign(Primal_constrk(gamma_xk))*epsilon;
  Primal_constrk(gamma_xk[0]) = epsilon;
  if(Primal_constrk(gamma_xk[0]) < 0.0)
    {
    Primal_constrk(gamma_xk[0]) = -epsilon;
    }
//
// z_xk = z_x;
  MatrixType z_xk = z_x;
//
// % loop parameters
// done = 0;
  bool done(false);
// iter = 0;
  unsigned int iter(0);
// data_precision = eps;   % floating point precision
  double data_precision(Eps);
//
// old_delta = 0;
  double old_delta(0);
// out_x = [];
  std::vector<unsigned long> out_x;
// count_delta_stop = 0;
  unsigned int count_delta_stop(0);
//
// constraint_plots = 1;
  unsigned constraint_plots(1);
//
// AtgxAgx = A(:,gamma_xk)'*A(:,gamma_xk);
  MatrixType AtgxAgx(gamma_xk.size(),1);
  for(unsigned int i = 0; i < gamma_xk.size(); ++i)
    {
    AtgxAgx(i) = A.col(gamma_xk[i]).transpose() * A.col(gamma_xk[i]);
    }
// iAtgxAgx = inv(A(:,gamma_xk)'*A(:,gamma_xk));
  MatrixType iAtgxAgx = AtgxAgx.inverse();
//
// while iter < maxiter
  while(iter < maxiter)
    {
//     iter = iter+1;
    ++iter;
//     % warning('off','MATLAB:divideByZero')
//
//     gamma_x = gamma_xk;
    gamma_x = gamma_xk;
//     z_x = z_xk;
    z_x = z_xk;
//     x_k = xk_1;
    MatrixType x_k = xk_1;
//
//     %%%%%%%%%%%%%%%%%%%%%
//     %%%% update on x %%%%
//     %%%%%%%%%%%%%%%%%%%%%
//
//     % Update direction
//     %del_x = inv(A(:,gamma_x)'*A(:,gamma_x))*z_x(gamma_x);
//     del_x = iAtgxAgx*z_x(gamma_x);
    MatrixType del_x;
    {
    MatrixType tmp_zx(gamma_x.size(),1);
    for(unsigned int i = 0; i < gamma_x.size(); ++i)
      {
      tmp_zx(i) = z_x(gamma_x[i]);
      }
    del_x = iAtgxAgx * tmp_zx;
    }
//     del_x_vec = zeros(N,1);
    MatrixType del_x_vec;
    del_x_vec = MatrixType::Zero(N,1);
//     del_x_vec(gamma_x) = del_x;
    for(unsigned i = 0; i < gamma_x.size(); ++i)
      {
      del_x_vec(gamma_x[i]) = del_x(i);
      }
//
//     pk = Primal_constrk;
    MatrixType pk = Primal_constrk;
//     %dk = A'*(A*del_x_vec);
//     Agdelx = A(:,gamma_x)*del_x;
    MatrixType Agdelx;
    {
    MatrixType tmpA(A.rows(),gamma_x.size());
    for(unsigned int i = 0; i < gamma_x.size(); ++i)
      {
      tmpA.col(i) = A.col(gamma_x[i]);
      }
    Agdelx = tmpA * del_x;
    }
//     dk = A'*Agdelx;
    MatrixType dk = A.transpose() * Agdelx;
//
//     %%% CONTROL THE MACHINE PRECISION ERROR AT EVERY OPERATION: LIKE BELOW.
//     pk_temp = Primal_constrk;
    MatrixType pk_temp = Primal_constrk;
//     gammaL_temp = find(abs(abs(Primal_constrk)-epsilon)<min(epsilon,2*eps));
    std::vector<unsigned long> gammaL_temp;
    for(unsigned i = 0; i < Primal_constrk.rows(); ++i)
      {
      if(std::fabs(std::fabs(Primal_constrk(i)) - epsilon) < std::min(epsilon,2*Eps))
        {
        gammaL_temp.push_back(i);
        }
      }
//     pk_temp(gammaL_temp) = sign(Primal_constrk(gammaL_temp))*epsilon;
    for(unsigned int i = 0; i < gammaL_temp.size(); ++i)
      {
      pk_temp(gammaL_temp[i]) = epsilon;
      if(Primal_constrk(gammaL_temp[i]) < 0.0)
        {
        pk_temp(gammaL_temp[i]) = -epsilon;
        }
      }
//
//     xk_temp = x_k;
    MatrixType xk_temp = x_k;
//     gammaX_temp = find(abs(x_k)<2*eps);
    std::vector<unsigned long> gammaX_temp;
    for(unsigned i = 0; i < x_k.rows(); ++i)
      {
      if(std::fabs(x_k(i)) < 2.0*Eps)
        {
        gammaX_temp.push_back(i);
        }
      }
//     xk_temp(gammaX_temp) = 0;
    for(unsigned i = 0; i < gammaX_temp.size(); ++i)
      {
      xk_temp(gammaX_temp[i]) = 0.0;
      }
//     %%%---
//
//     % Compute the step size
//     [i_delta, out_x, delta, chk_x] = update_primal(gamma_x,
//     gamma_x, z_x,  xk_temp, del_x_vec, pk_temp, dk, epsilon,  out_x);
    unsigned long i_delta;
    double        delta;
    unsigned      chk_x;
    this->update_primal(gamma_x,
                        gamma_x,
                        z_x,
                        xk_temp,
                        del_x_vec,
                        pk_temp,
                        dk,
                        epsilon,
                        out_x,
                        i_delta,
                        out_x,
                        delta,
                        chk_x);
//
//     if old_delta < 4*eps && delta < 4*eps
    if(old_delta < 4.0 * Eps && delta < 4.0 * Eps)
      {
//         count_delta_stop = count_delta_stop + 1;
      count_delta_stop++;
      }
//     else
    else
      {
//         count_delta_stop = 0;
      count_delta_stop = 0;
//     end
      }
//     if count_delta_stop >= 500
    if(count_delta_stop >= 500)
      {
//         %disp('stuck in some corner');
//         break;
      // TODO: print a message? throw an exception?
      break;
//     end
      }
//     old_delta = delta;
    old_delta = delta;
//
//     xk_1 = x_k+delta*del_x_vec;
    xk_1 = x_k + delta * del_x_vec;
//     Primal_constrk = pk+delta*dk;
    Primal_constrk = pk + delta * dk;
//     epsilon_old = epsilon;
    double epsilon_old = epsilon;
//     epsilon = epsilon-delta;
    epsilon -= delta;
//
//     if epsilon <= tau;
    if(epsilon <= tau)
      {
//         xk_1 = x_k + (epsilon_old-tau)*del_x_vec;
      xk_1 = x_k + (epsilon_old - tau) * del_x_vec;
//         total_time= cputime-t0;
//         break;
//     end
      }
//
//     if chk_x == 1
    if(chk_x == 1)
      {
//         % If an element is removed from gamma_x
//         gx_old = gamma_x;
      std::vector<unsigned long> gx_old = gamma_x;
//         len_gamma = length(gamma_x);
      unsigned long len_gamma = gamma_x.size();
//
//         outx_index = find(gamma_x==out_x);
      unsigned long outx_index(0);
      {
      unsigned long limit = std::max<unsigned long>(gamma_x.size(),out_x.size());
      for(unsigned i = 0; i < limit; ++i)
        {
        if(i < gamma_x.size() && i < out_x.size() && gamma_x[i] == out_x[i])
          {
          outx_index = i;
          break;
          }
        }
      }
//         gamma_x(outx_index) = gamma_x(len_gamma);
      gamma_x[outx_index] = gamma_x[len_gamma - 1];
//         gamma_x(len_gamma) = out_x;
      gamma_x[len_gamma - 1] = out_x[0];
//         gamma_x = gamma_x(1:len_gamma-1);
      gamma_x.pop_back();
//         gamma_xk = gamma_x;
      gamma_xk = gamma_x;
//
//         rowi = outx_index; % ith row of A is swapped with last row (out_x)
      unsigned long rowi = outx_index;
//         colj = outx_index; % jth column of A is swapped with last column (out_lambda)
      unsigned long colj = outx_index;
//         AtgxAgx_ij = AtgxAgx;
      MatrixType AtgxAgx_ij = AtgxAgx;
//         temp_row = AtgxAgx_ij(rowi,:);
      MatrixType temp_row = AtgxAgx_ij.row(rowi);
//         AtgxAgx_ij(rowi,:) = AtgxAgx_ij(len_gamma,:);
      AtgxAgx_ij.row(rowi) = AtgxAgx_ij.row(len_gamma - 1);
//         AtgxAgx_ij(len_gamma,:) = temp_row;
      AtgxAgx_ij.row(len_gamma-1) = temp_row;
//         temp_col = AtgxAgx_ij(:,colj);
      MatrixType temp_col = AtgxAgx_ij.col(colj);
//         AtgxAgx_ij(:,colj) = AtgxAgx_ij(:,len_gamma);
      AtgxAgx_ij.col(colj) = AtgxAgx_ij.col(len_gamma - 1);
//         AtgxAgx_ij(:,len_gamma) = temp_col;
      AtgxAgx_ij.col(len_gamma-1) = temp_col;
//         iAtgxAgx_ij = iAtgxAgx;
//         temp_row = iAtgxAgx_ij(colj,:);
//         iAtgxAgx_ij(colj,:) = iAtgxAgx_ij(len_gamma,:);
//         iAtgxAgx_ij(len_gamma,:) = temp_row;
//         temp_col = iAtgxAgx_ij(:,rowi);
//         iAtgxAgx_ij(:,rowi) = iAtgxAgx_ij(:,len_gamma);
//         iAtgxAgx_ij(:,len_gamma) = temp_col;
//
//         AtgxAgx = AtgxAgx_ij(1:len_gamma-1,1:len_gamma-1);
      AtgxAgx.conservativeResize(len_gamma - 1, len_gamma - 1);
//         iAtgxAgx = update_inverse(AtgxAgx_ij, iAtgxAgx_ij,2);
      iAtgxAgx = AtgxAgx.inverse();
//         xk_1(out_x) = 0;
      for(unsigned int i = 0; i < out_x.size(); ++i)
        {
        xk_1(out_x[i]) = 0.0;
        }
      }
//     else
    else
      {
//         % If an element is added to gamma_x
//         gamma_xk = [gamma_x; i_delta];
      gamma_xk = gamma_x;
      gamma_xk.push_back(i_delta);
//         new_x = i_delta;
      unsigned long new_x = i_delta;
//
//         AtgxAnx = A(:,gamma_x)'*A(:,new_x);
      MatrixType AtgxAnx;
      {
      MatrixType ATmp(A.rows(),gamma_x.size());
      for(unsigned int i = 0; i < gamma_x.size(); ++i)
        {
        ATmp.col(i) = A.col(gamma_x[i]);
        }
      AtgxAnx = ATmp.transpose() * A.col(new_x);
      }
//         AtgxAgx_mod = [AtgxAgx AtgxAnx; AtgxAnx' A(:,new_x)'*A(:,i_delta)];
      MatrixType AtgxAgx_mod(AtgxAgx.rows() + 1, AtgxAgx.cols() + 1);
      //
      // makes a new matrix of this form
      // | AtgxAgx  AtgxAnx |
      // | AtgxAnx' X       | where X = A(:,new_x)'*A(:,i_delta)
      // and this works great if AtgxAgx.rows() == AtgxAnx rows,
      // and AtgxAgx is square.
      {
      unsigned rows = AtgxAgx.rows(),
        cols = AtgxAgx.cols();
      AtgxAgx_mod.block(0,0,rows,cols) = AtgxAgx;
      AtgxAgx_mod.block(0,cols,AtgxAnx.rows(),AtgxAnx.cols()) = AtgxAnx;
      AtgxAgx_mod.block(rows,0,AtgxAnx.cols(),AtgxAnx.rows()) = AtgxAnx.transpose();
      AtgxAgx_mod(rows,cols) = A.col(new_x).transpose() * A.col(i_delta);
      }
//
//         AtgxAgx = AtgxAgx_mod;
      AtgxAgx = AtgxAgx_mod;
//         iAtgxAgx = update_inverse(AtgxAgx, iAtgxAgx,1);
      iAtgxAgx = AtgxAgx.inverse();
//         xk_1(i_delta) = 0;
      xk_1(i_delta) = 0.0;
//     end
//
//     z_xk = zeros(N,1);
      z_xk = MatrixType::Zero(N,1);
//     z_xk(gamma_xk) = -sign(Primal_constrk(gamma_xk));
      for(unsigned int i = 0; i < gamma_xk.size(); ++i)
        {
        z_xk(gamma_xk[i]) = 1.0;
        if(Primal_constrk(gamma_xk[i]) > 0.0)
          {
          z_xk(gamma_xk[i]) = -1.0;
          }
        }
//     Primal_constrk([gamma_x]) = sign(Primal_constrk([gamma_x]))*epsilon;
      for(unsigned int i = 0; i < gamma_x.size(); ++i)
        {
        Primal_constrk(gamma_x[i]) = epsilon;
        if(Primal_constrk(gamma_x[i]) < 0.0)
          {
          Primal_constrk(gamma_x[i]) = -epsilon;
          }
        }
// end
      }
    }
// total_iter = iter;
// x_out = xk_1;
  return xk_1;
}
// % BPDN_homotopy_function.m
// %
// % Solves the following basis pursuit denoising (BPDN) problem
// % min_x  \tau ||x||_1 + 1/2*||y-Ax||_2^2
// %
// % Inputs:
// % A - m x n measurement matrix
// % y - measurement vector
// % tau - final value of regularization parameter
// % maxiter - maximum number of homotopy iterations
// %
// % Outputs:
// % x_out - output for BPDN
// % gamma_x - support of the solution
// % total_iter - number of homotopy iterations taken by the solver
// % total_time - time taken by the solver
// %
// % Written by: Salman Asif, Georgia Tech
// % Email: sasif@ece.gatech.edu
// %
// %-------------------------------------------+
// % Copyright (c) 2007.  Muhammad Salman Asif
// %-------------------------------------------+
//
// function [x_out, gamma_x, total_iter, total_time] = BPDN_homotopy_function(A, y, tau, maxiter)
//
// t0 = cputime;
//
// N = size(A,2);
// K = size(A,1);
//
// % Initialization of primal and dual sign and support
// z_x = zeros(N,1);
// gamma_x = [];       % Primal support
//
// % Initial step
// Primal_constrk = -A'*y;
// [c i] = max(abs(Primal_constrk));
//
// gamma_xk = i;
//
// epsilon = c;
// xk_1 = zeros(N,1);
//
// z_x(gamma_xk) = -sign(Primal_constrk(gamma_xk));
// Primal_constrk(gamma_xk) = sign(Primal_constrk(gamma_xk))*epsilon;
//
// z_xk = z_x;
//
// % loop parameters
// done = 0;
// iter = 0;
// data_precision = eps;   % floating point precision
//
// old_delta = 0;
// out_x = [];
// count_delta_stop = 0;
//
// constraint_plots = 1;
//
// AtgxAgx = A(:,gamma_xk)'*A(:,gamma_xk);
// iAtgxAgx = inv(A(:,gamma_xk)'*A(:,gamma_xk));
//
// while iter < maxiter
//     iter = iter+1;
//     % warning('off','MATLAB:divideByZero')
//
//     gamma_x = gamma_xk;
//     z_x = z_xk;
//     x_k = xk_1;
//
//     %%%%%%%%%%%%%%%%%%%%%
//     %%%% update on x %%%%
//     %%%%%%%%%%%%%%%%%%%%%
//
//     % Update direction
//     %del_x = inv(A(:,gamma_x)'*A(:,gamma_x))*z_x(gamma_x);
//     del_x = iAtgxAgx*z_x(gamma_x);
//     del_x_vec = zeros(N,1);
//     del_x_vec(gamma_x) = del_x;
//
//     pk = Primal_constrk;
//     %dk = A'*(A*del_x_vec);
//     Agdelx = A(:,gamma_x)*del_x;
//     dk = A'*Agdelx;
//
//     %%% CONTROL THE MACHINE PRECISION ERROR AT EVERY OPERATION: LIKE BELOW.
//     pk_temp = Primal_constrk;
//     gammaL_temp = find(abs(abs(Primal_constrk)-epsilon)<min(epsilon,2*eps));
//     pk_temp(gammaL_temp) = sign(Primal_constrk(gammaL_temp))*epsilon;
//
//     xk_temp = x_k;
//     gammaX_temp = find(abs(x_k)<2*eps);
//     xk_temp(gammaX_temp) = 0;
//     %%%---
//
//     % Compute the step size
//     [i_delta, out_x, delta, chk_x] = update_primal(gamma_x,
//     gamma_x, z_x,  xk_temp, del_x_vec, pk_temp, dk, epsilon,//  out_x);
//
//     if old_delta < 4*eps && delta < 4*eps
//         count_delta_stop = count_delta_stop + 1;
//     else
//         count_delta_stop = 0;
//     end
//     if count_delta_stop >= 500
//         %disp('stuck in some corner');
//         break;
//     end
//     old_delta = delta;
//
//     xk_1 = x_k+delta*del_x_vec;
//     Primal_constrk = pk+delta*dk;
//     epsilon_old = epsilon;
//     epsilon = epsilon-delta;
//
//     if epsilon <= tau;
//         xk_1 = x_k + (epsilon_old-tau)*del_x_vec;
//         total_time= cputime-t0;
//         break;
//     end
//
//     if chk_x == 1
//         % If an element is removed from gamma_x
//         gx_old = gamma_x;
//         len_gamma = length(gamma_x);
//
//         outx_index = find(gamma_x==out_x);
//         gamma_x(outx_index) = gamma_x(len_gamma);
//         gamma_x(len_gamma) = out_x;
//         gamma_x = gamma_x(1:len_gamma-1);
//         gamma_xk = gamma_x;
//
//         rowi = outx_index; % ith row of A is swapped with last row (out_x)
//         colj = outx_index; % jth column of A is swapped with last column (out_lambda)
//         AtgxAgx_ij = AtgxAgx;
//         temp_row = AtgxAgx_ij(rowi,:);
//         AtgxAgx_ij(rowi,:) = AtgxAgx_ij(len_gamma,:);
//         AtgxAgx_ij(len_gamma,:) = temp_row;
//         temp_col = AtgxAgx_ij(:,colj);
//         AtgxAgx_ij(:,colj) = AtgxAgx_ij(:,len_gamma);
//         AtgxAgx_ij(:,len_gamma) = temp_col;
//         iAtgxAgx_ij = iAtgxAgx;
//         temp_row = iAtgxAgx_ij(colj,:);
//         iAtgxAgx_ij(colj,:) = iAtgxAgx_ij(len_gamma,:);
//         iAtgxAgx_ij(len_gamma,:) = temp_row;
//         temp_col = iAtgxAgx_ij(:,rowi);
//         iAtgxAgx_ij(:,rowi) = iAtgxAgx_ij(:,len_gamma);
//         iAtgxAgx_ij(:,len_gamma) = temp_col;
//
//         AtgxAgx = AtgxAgx_ij(1:len_gamma-1,1:len_gamma-1);
//         iAtgxAgx = update_inverse(AtgxAgx_ij, iAtgxAgx_ij,2);
//         xk_1(out_x) = 0;
//     else
//         % If an element is added to gamma_x
//         gamma_xk = [gamma_x; i_delta];
//         new_x = i_delta;
//
//         AtgxAnx = A(:,gamma_x)'*A(:,new_x);
//         AtgxAgx_mod = [AtgxAgx AtgxAnx; AtgxAnx' A(:,new_x)'*A(:,i_delta)];
//
//         AtgxAgx = AtgxAgx_mod;
//         iAtgxAgx = update_inverse(AtgxAgx, iAtgxAgx,1);
//         xk_1(i_delta) = 0;
//     end
//
//     z_xk = zeros(N,1);
//     z_xk(gamma_xk) = -sign(Primal_constrk(gamma_xk));
//     Primal_constrk([gamma_x]) = sign(Primal_constrk([gamma_x]))*epsilon;
// end
// total_iter = iter;
// x_out = xk_1;
//

void
DoCSEstimate
::update_primal(std::vector<unsigned long> & gamma_x, // current support of x
                std::vector<unsigned long> & gamma_lambda, // current support of
                                // lambda
                MatrixType &                 z_x, // sign sequence of x
                MatrixType &                 x_k, // sign sequence of lambda
                MatrixType &                 del_x_vec, // primal update direction
                MatrixType &                 pk,       //
                MatrixType &                 dk,       //
                double                       epsilon, // current value of epsilon
                std::vector<unsigned long> & out_lambda, // element removed from
                                                         // support of lambda in
                                                         // previous step if any
                                                         // OUTPUTS
                unsigned long &              i_delta, // index corresponding to
                                                      // newly active primal
                                                      // constraint (new_lambda)
                std::vector<unsigned long> & out_x, // element in x shrunk to zero;
                double  &                    delta, // primal step size
                unsigned &                   chk_x) // 1 an element is removed
                                                    // from support of x
                                                    // 0 a new element enters
                                                    // the support of lambd
{
// N = length(x_k);
  const unsigned long N(x_k.rows());
//
// % gamma_lc = setdiff([1:N]', [gamma_lambda; out_lambda]); WRONG
// % check out_lambda as well, that is if outgoing lambda switches sign in just one step
// temp_gamma = zeros(N,1);
  MatrixType temp_gamma;
  temp_gamma = MatrixType::Zero(N,1);
// temp_gamma(gamma_lambda) = gamma_lambda;
// gamma_lc = find([1:N]' ~= temp_gamma);
  // NOTE: the above two lines are just constructing a list
  // of array indices that are not in gamma_lambda
  std::vector<unsigned long> gamma_lc;
  for(unsigned int i = 0; i < N; ++i)
    {
    bool exclude(false);
    for(unsigned int j = 0; j < gamma_lambda.size(); ++j)
      {
      if(i == gamma_lambda[j])
        {
        exclude = true;
        break;
        }
      }
    if(!exclude)
      {
      gamma_lc.push_back(i);
      }
    }
//
// delta1_constr = (epsilon-pk(gamma_lc))./(1+dk(gamma_lc));
  MatrixType delta1_constr(gamma_lc.size(),1);
  for(unsigned int i = 0; i < gamma_lc.size(); ++i)
    {
    delta1_constr(i) = (epsilon - pk(gamma_lc[i])) / (1 + dk(gamma_lc[i]));
    }
// delta1_pos_ind = find(delta1_constr>0);
  std::vector<unsigned long> delta1_pos_ind;
  for(unsigned int i = 0; i < delta1_constr.rows(); ++i)
    {
    if(delta1_constr(i) > 0)
      {
      delta1_pos_ind.push_back(i);
      }
    }
// delta1_pos = delta1_constr(delta1_pos_ind);
  double delta1 = std::numeric_limits<double>::infinity();;
  unsigned long i_delta1(0), dummy(0);
  if(delta1_pos_ind.size() > 0)
    {
    MatrixType delta1_pos(delta1_pos_ind.size(),1);
    for(unsigned int i = 0; i < delta1_pos_ind.size(); ++i)
      {
      delta1_pos(i) = delta1_constr(delta1_pos_ind[i]);
      }
// [delta1 i_delta1] = min(delta1_pos);
    delta1 = delta1_pos.minCoeff(&i_delta1,&dummy);
    }
// if isempty(delta1)
//     delta1 = inf;
// end

// delta2_constr = (epsilon+pk(gamma_lc))./(1-dk(gamma_lc));
  MatrixType delta2_constr(gamma_lc.size(),1);
  for(unsigned int i = 0; i < gamma_lc.size(); ++i)
    {
    delta2_constr(i) = (epsilon + pk(gamma_lc[i])) / (1 - dk(gamma_lc[i]));
    }
// delta2_pos_ind = find(delta2_constr>0);
  std::vector<unsigned long> delta2_pos_ind;
  for(unsigned i = 0; i < delta2_constr.rows(); ++i)
    {
    if(delta2_constr(i) > 0.0)
      {
      delta2_pos_ind.push_back(i);
      }
    }
  double delta2 = std::numeric_limits<double>::infinity();
  unsigned long i_delta2(0);
  if(delta2_pos_ind.size() > 0)
    {
// delta2_pos = delta2_constr(delta2_pos_ind);
    MatrixType delta2_pos(delta2_pos_ind.size(),1);
    for(unsigned int i = 0; i < delta2_pos_ind.size(); ++i)
      {
      delta2_pos(i) = delta1_constr(delta2_pos_ind[i]);
      }
// [delta2 i_delta2] = min(delta2_pos);
    delta2 = delta2_pos.minCoeff(&i_delta2,&dummy);
// if isempty(delta2)
//     delta2 = inf;
// end
    }
// if delta1>delta2
  if(delta1 > delta2)
    {
//     delta = delta2;
    delta = delta2;
//     i_delta = gamma_lc(delta2_pos_ind(i_delta2));
    i_delta = gamma_lc[delta2_pos_ind[i_delta2]];
    }
// else
  else
    {
//     delta = delta1;
    delta = delta1;
//     i_delta = gamma_lc(delta1_pos_ind(i_delta1));
    i_delta = gamma_lc[delta1_pos_ind[i_delta1]];
// end
    }
//
// delta3_constr = (-x_k(gamma_x)./del_x_vec(gamma_x));
  MatrixType delta3_constr(gamma_x.size(),1);
  for(unsigned int i = 0; i < gamma_x.size(); ++i)
    {
    delta3_constr(i) = -x_k(gamma_x[i]) / del_x_vec(gamma_x[i]);
    }
// delta3_pos_index = find(delta3_constr>0);
  std::vector<unsigned long> delta3_pos_index;
  for(unsigned i = 0; i < delta3_constr.rows(); ++i)
    {
    if(delta3_constr(i) > 0)
      {
      delta3_pos_index.push_back(i);
      }
    }
// [delta3 i_delta3] = min(delta3_constr(delta3_pos_index));
  double delta3 = std::numeric_limits<double>::infinity();
  unsigned long i_delta3(0),ijunk(0);
  if(delta3_pos_index.size() > 0)
    {
    MatrixType delta3_pos(delta3_pos_index.size(),1);
    for(unsigned i = 0; i < delta3_pos_index.size(); ++i)
      {
      delta3_pos(i) = delta3_constr(delta3_pos_index[i]);
      }
    delta3 = delta3_pos.minCoeff(&i_delta3,&ijunk);
    }
// out_x_index = gamma_x(delta3_pos_index(i_delta3));
  unsigned long out_x_index = gamma_x[delta3_pos_index[i_delta3]];
//
// chk_x = 0;
  chk_x = 0;
// out_x = [];
  out_x.clear();
// if delta3 > 0 & delta3 <= delta
  if(delta3 > 0 && delta3 <= delta)
    {
//     chk_x = 1;
    chk_x = 1;
//     delta = delta3;
    delta = delta3;
//     out_x = out_x_index;
      out_x.push_back(out_x_index);
// end
    }
//
// %%% THESE ARE PROBABLY UNNECESSARY
// %%% NEED TO REMOVE THEM.
//
// % The following checks are just to deal with degenerate cases when more
// % than one elements want to enter or leave the support at any step
// % (e.g., Bernoulli matrix with small number of measurements)
//
// % This one is ONLY for those indices which are zero. And we don't know where
// % will its dx point in next steps, so after we calculate dx and its in opposite
// % direction to z_x, we will have to remove that index from the support.
// xk_1 = x_k+delta*del_x_vec;
  MatrixType xk_1 = x_k + delta * del_x_vec;
// xk_1(out_x) = 0;
  for(unsigned int i = 0; i < out_x.size(); ++i)
    {
    xk_1(out_x[i]) = 0;
    }
// wrong_sign = find(sign(xk_1(gamma_x)).*z_x(gamma_x)==-1);
  std::vector<unsigned long> wrong_sign;
  for(unsigned int i = 0; i < gamma_x.size(); ++i)
    {
    if((z_x(gamma_x[i]) * (xk_1(gamma_x[i]) > 0.0 ? 1.0 : -1.0)) == -1)
      {
      wrong_sign.push_back(i);
      }
    }
// if ~isempty(gamma_x(wrong_sign))
  if(wrong_sign.size() > 0)
    {
//     chk_x = 1;
    chk_x = 1;
//     delta = 0;
    delta = 0.0;
//     % can also choose specific element which became non-zero first but all
//     % that matters here is AtA(gx,gl) doesn't become singular.
//     % [val_wrong_x ind_wrong_x] =  sort(abs(del_x_vec(gamma_x(wrong_sign))),'descend');
//     out_x = gamma_x(wrong_sign(1));
    out_x.clear();
    for(unsigned int i = 0; i < wrong_sign.size(); ++i)
      {
      out_x.push_back(gamma_x[wrong_sign[i]]);
      }
// end
//
    }
// % If more than one primal constraints became active in previous iteration i.e.,
// % more than one elements wanted to enter the support and we added only one.
// % So here we need to check if those remaining elements are still active.
// i_delta_temp = gamma_lc(find(abs(pk(gamma_lc)+delta*dk(gamma_lc))-(epsilon-delta) >= 10*eps));
  std::vector<unsigned long> i_delta_temp;
  for(unsigned i = 0; i < gamma_lc.size(); ++i)
    {
    if(std::abs(pk(gamma_lc[i]) + delta * dk(gamma_lc[i])) - (epsilon - delta) >= 10*Eps)
      {
      i_delta_temp.push_back(gamma_lc[i]);
      }
    }
  std::vector<unsigned long> i_delta_more;
// if ~isempty(i_delta_temp)
  if(i_delta_temp.size() > 0)
    {
//     if ~isempty(out_lambda)
//         i_delta_more =
//         i_delta_temp;%(find(i_delta_temp~=out_lambda));
//     else
//         i_delta_more = i_delta_temp;
//     end
    // The above I guess used to do things depending on out_lambda.
    // now it does the same thing regardless
    i_delta_more = i_delta_temp;
    unsigned equalcount(0);
    for(unsigned i = 0; i < i_delta_temp.size() && i < i_delta; ++i)
      {
      if(i_delta_temp[i] == i_delta)
        {
        ++equalcount;
        }
      }
//     if length(i_delta_more)>=1 & ~sum((i_delta_temp==i_delta))
    if(i_delta_more.size() > 0 && equalcount == 0)
      {
//         % ideal way would be to check that incoming element doesn't make AtA
//         % singular!
//         [v_temp i_temp] = max(-pk(i_delta_more)./dk(i_delta_more));
      MatrixType tmp(i_delta_more.size(),1);
      for(unsigned int i = 0; i < i_delta_more.size(); ++i)
        {
        tmp(i) = -pk(i_delta_more[i]) / dk(i_delta_more[i]);
        }
      unsigned long i_temp, junk;
      double v_temp = tmp.maxCoeff(&i_temp,&junk);
//         i_delta = i_delta_more(i_temp);
      i_delta = i_delta_more[i_temp];
//         delta = 0;
      delta = 0;
//         chk_x = 0;
      chk_x = 0;
//         out_x = [];
      out_x.clear();
//     end
      }
// end
    }
}
// % update_primal.m
// %
// % This function computes the minimum step size in the primal update direction and
// % finds change in the primal or dual support with that step.
// %
// % Inputs:
// % gamma_x - current support of x
// % gamma_lambda - current support of lambda
// % z_x - sign sequence of x
// % z_lambda - sign sequence of lambda
// % del_x_vec - primal update direction
// % pk
// % dk
// % epsilon - current value of epsilon
// % out_lambda - element removed from support of lambda in previous step (if any)
// %
// % Outputs:
// % i_delta - index corresponding to newly active primal constraint (new_lambda)
// % out_x - element in x shrunk to zero
// % delta - primal step size
// % chk_x - 1  an element is removed from support of x
// %         0  a new element enters the support of lambda
// %
// % Written by: Salman Asif, Georgia Tech
// % Email: sasif@ece.gatech.edu
//
// function [i_delta, out_x, delta, chk_x] = update_primal(gamma_x,
// gamma_lambda, z_x, x_k, del_x_vec, pk, dk, epsilon, out_lambda)// ;
//
// N = length(x_k);
//
// % gamma_lc = setdiff([1:N]', [gamma_lambda; out_lambda]); WRONG
// % check out_lambda as well, that is if outgoing lambda switches sign in just one step
// temp_gamma = zeros(N,1);
// temp_gamma(gamma_lambda) = gamma_lambda;
// gamma_lc = find([1:N]' ~= temp_gamma);
//
// delta1_constr = (epsilon-pk(gamma_lc))./(1+dk(gamma_lc));
// delta1_pos_ind = find(delta1_constr>0);
// delta1_pos = delta1_constr(delta1_pos_ind);
// [delta1 i_delta1] = min(delta1_pos);
// if isempty(delta1)
//     delta1 = inf;
// end
// delta2_constr = (epsilon+pk(gamma_lc))./(1-dk(gamma_lc));
// delta2_pos_ind = find(delta2_constr>0);
// delta2_pos = delta2_constr(delta2_pos_ind);
// [delta2 i_delta2] = min(delta2_pos);
// if isempty(delta2)
//     delta2 = inf;
// end
//
// if delta1>delta2
//     delta = delta2;
//     i_delta = gamma_lc(delta2_pos_ind(i_delta2));
// else
//     delta = delta1;
//     i_delta = gamma_lc(delta1_pos_ind(i_delta1));
// end
//
// delta3_constr = (-x_k(gamma_x)./del_x_vec(gamma_x));
// delta3_pos_index = find(delta3_constr>0);
// [delta3 i_delta3] = min(delta3_constr(delta3_pos_index));
// out_x_index = gamma_x(delta3_pos_index(i_delta3));
//
// chk_x = 0;
// out_x = [];
// if delta3 > 0 & delta3 <= delta
//     chk_x = 1;
//     delta = delta3;
//     out_x = out_x_index;
// end
//
// %%% THESE ARE PROBABLY UNNECESSARY
// %%% NEED TO REMOVE THEM.
//
// % The following checks are just to deal with degenerate cases when more
// % than one elements want to enter or leave the support at any step
// % (e.g., Bernoulli matrix with small number of measurements)
//
// % This one is ONLY for those indices which are zero. And we don't know where
// % will its dx point in next steps, so after we calculate dx and its in opposite
// % direction to z_x, we will have to remove that index from the support.
// xk_1 = x_k+delta*del_x_vec;
// xk_1(out_x) = 0;
// wrong_sign = find(sign(xk_1(gamma_x)).*z_x(gamma_x)==-1);
// if ~isempty(gamma_x(wrong_sign))
//     chk_x = 1;
//     delta = 0;
//     % can also choose specific element which became non-zero first but all
//     % that matters here is AtA(gx,gl) doesn't become singular.
//     % [val_wrong_x ind_wrong_x] =  sort(abs(del_x_vec(gamma_x(wrong_sign))),'descend');
//     out_x = gamma_x(wrong_sign(1));
// end
//
// % If more than one primal constraints became active in previous iteration i.e.,
// % more than one elements wanted to enter the support and we added only one.
// % So here we need to check if those remaining elements are still active.
// i_delta_temp = gamma_lc(find(abs(pk(gamma_lc)+delta*dk(gamma_lc))-(epsilon-delta) >= 10*eps));
// if ~isempty(i_delta_temp)
//     if ~isempty(out_lambda)
//         i_delta_more = i_delta_temp;%(find(i_delta_temp~=out_lambda));
//     else
//         i_delta_more = i_delta_temp;
//     end
//     if length(i_delta_more)>=1 & ~sum((i_delta_temp==i_delta))
//         % ideal way would be to check that incoming element doesn't make AtA
//         % singular!
//         [v_temp i_temp] = max(-pk(i_delta_more)./dk(i_delta_more));
//         i_delta = i_delta_more(i_temp);
//         delta = 0;
//         chk_x = 0;
//         out_x = [];
//     end
// end
