#include "DoCSEstimate.h"
#include "BuildRidges.h"
#include "BuildSensor.h"
#include "MultiSample.h"
#include "itkSubtractImageFilter.h"
#include "itkAbsImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkVectorImageToImageAdaptor.h"
#include "itkCastImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

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

  for(b0It.Begin(), avgIt.Begin(); !b0It.IsAtEnd() && !avgIt.IsAtEnd(); ++b0It, ++avgIt)
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
  this->m_IntensityData->SetNumberOfComponentsPerPixel(gradientIndices.size());
  this->m_IntensityData->CopyInformation(originalNrrd);
  this->m_IntensityData->SetRegions(region);
  this->m_IntensityData->Allocate();

  const unsigned int numGradientIndices = gradientIndices.size();
  //
  // make new vectorimage with just the gradient values
  itk::ImageRegionIterator<DWIVectorImageType> newIt(this->m_IntensityData,region);
  for(b0It.Begin(), newIt.Begin(); !b0It.IsAtEnd() && !newIt.IsAtEnd(); ++b0It,++newIt)
    {
    DWIVectorImageType::PixelType newVal;
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

  unsigned numGradientDirections = this->m_IntensityData->GetNumberOfComponentsPerPixel();

  //
  // MATLAB CODE IS
  // DWIIntensityData = single(DWIIntensityData) ./ averagedB0(:,:,:,ones(1,numGradientDirections));
  // %remove negative values
  // DWIIntensityData(DWIIntensityData<0)=eps;
  itk::ImageRegionIterator<DWIVectorImageType>
    ivIt(this->m_IntensityData,this->m_IntensityData->GetLargestPossibleRegion());
  unsigned int numGrad = this->m_IntensityData->GetNumberOfComponentsPerPixel();

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
  rho = std::pow((rho * J),(J * p));
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
  //
  // matlabpool('open');
  //
  // u=step2(DWIIntensityData,myu,tol);
  DWIVectorImageType::Pointer u = this->step2(myu,tol);
  // c=step1(DWIIntensityData,A,lmd,NIT,id);   % initialization of ridgelet coefficients

  //
  // Ac=reshape(reshape(c,[numGradientVoxels M])*A',[nx ny nz numGradientDirections]);
  // p=Ac-u;
  //
  //
  // TNIT=3;                     % number of outer iterations
  // for itr=1:TNIT,
  // fprintf(1,'Iteration %d of %d\t',[itr TNIT]);
  //
  // t=u-p;
  // c=step1(t,A,lmd/gama,NIT,id);
  // Ac=reshape(reshape(c,[numGradientVoxels M])*A',[nx ny nz numGradientDirections]);
  //
  // t=(1/(1+gama))*(DWIIntensityData+gama*(Ac+p));
  // u=step2(t,myu/(1+gama),tol);
  //
  // p=p+(Ac-u);
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
  //
  // % Scale up  gradient direction data by averagedB0 data
  // estimatedSignal = estimatedSignal.* averagedB0(:,:,:,ones(1,n0)); %multiply by the B0 image
  // %% Insert the B0 back in
  // estimatedSignal = cat(4,averagedB0,estimatedSignal); %add B0 to the data
  // estimatedGradients=[0 0 0;new_gradients];


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
  MatrixType Y = SphericalHarmonics(u,22);

  // C=(2*(0:m)'+1)/(4*pi);
  MatrixType C(1,m+1);
  for(unsigned int i = 1; i < m+1; ++i)
    {
    C(0,i-1) = ((2 * (i-1)) + 1)/(4*this->m_Pi);
    }

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
  //for k=2:2:m,
  for(unsigned int k = 2; k <= m; k += 2)
    {
    // Lmd(k+1)=-Lmd(k-1)*(k-1)/k;
    Lmd(k) = -Lmd(k - 2) * (k - 1)/k;
    }

  std::vector<unsigned> ind(m+1);
  //ind=[1; ((1:m).^2)'+(2:m+1)'];
  ind[0] = 1;
  for(unsigned int i = 1; i < m+1; ++i)
    {
    ind[i] = (i * i) + ( i + 1);
    }
  // s=(Y'*Y)\(Y'*S);
  MatrixType s = Y.transpose() * Y;
  MatrixType s2 = Y.transpose() * S;
  Eigen::FullPivLU<MatrixType> lu(s);
  s = lu.solve(s2);

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
                      const std::vector<unsigned int> indices)
{
  typename ImageType::Pointer rval = ImageType::New();
  rval->CopyInformation(input);
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
::FromVecImage(unsigned gradientIndex)
{
  B0AvgImageType::Pointer f = B0AvgImageType::New();
  f->CopyInformation(this->m_IntensityData);
  f->Allocate();

  itk::ImageRegionConstIterator<DWIVectorImageType>
    dwiIt(this->m_IntensityData,this->m_IntensityData->GetLargestPossibleRegion());
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
::tvdenoise3(unsigned int gradientIndex,double lambda,double tol,DWIVectorImageType::Pointer &target)
{
  typedef itk::SubtractImageFilter<B0AvgImageType,B0AvgImageType,B0AvgImageType>  SubtractFilterType;
  typedef itk::AbsImageFilter<B0AvgImageType,B0AvgImageType>                      AbsImageFilterType;
  typedef itk::StatisticsImageFilter<B0AvgImageType>                              StatsFilterType;
  typedef itk::MultiplyImageFilter<B0AvgImageType,B0AvgImageType,B0AvgImageType>  MultiplyImageFilterType;

  B0AvgImageType::Pointer f = this->FromVecImage(gradientIndex);

  MultiplyImageFilterType::Pointer mult = MultiplyImageFilterType::New();
  mult->SetInput1(f);
  mult->SetConstant(lambda);
  mult->Update();
  B0AvgImageType::Pointer fLambda = mult->GetOutput();

  // dt = 0.25/2;
  const double dt = 0.125;

  // N = size(f);
  const ImageSizeType N = this->m_IntensityData->GetLargestPossibleRegion().GetSize();

  // id = [2:N(1),N(1)];
  std::vector<unsigned int> id;
  for(unsigned int i = 1; i < N[0] - 1; ++i)
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
  for(unsigned int i = 1; i < N[2]; ++i)
    {
    ir.push_back(i);
    }
  ir.push_back(N[2] - 1);

  // il = [1,1:N(2)-1];
  std::vector<unsigned int> il;
  il.push_back(0);
  for(unsigned int i = 0; i < N[2] - 1; ++i)
    {
    il.push_back(i);
    }

  // ib = [2:N(3),N(3)];
  std::vector<unsigned int> ib;
  for(unsigned int i = 1; i < N[3]; ++i)
    {
    ib.push_back(i);
    }
  ib.push_back(N[3] - 1);

  // ifr = [1,1:N(3)-1];
  std::vector<unsigned int> ifr;
  ifr.push_back(0);
  for(unsigned int i = 0; i < N[3] - 1; ++i)
    {
    ifr.push_back(i);
    }
  //
  // p1 = zeros(size(f));
  B0AvgImageType::Pointer p1 = B0AvgImageType::New();
  p1->SetRegions(N);
  p1->Allocate();
  p1->FillBuffer(0.0);
  // p2 = zeros(size(f));
  B0AvgImageType::Pointer p2 = B0AvgImageType::New();
  p2->SetRegions(N);
  p2->Allocate();
  p2->FillBuffer(0.0);
  // p3 = zeros(size(f));
  B0AvgImageType::Pointer p3 = B0AvgImageType::New();
  p3->SetRegions(N);
  p3->Allocate();
  p3->FillBuffer(0.0);
  //
  // divp = zeros(size(f));
  B0AvgImageType::Pointer divp = B0AvgImageType::New();
  divp->SetRegions(N);
  divp->Allocate();
  divp->FillBuffer(0.0);

  // lastdivp = ones(size(f));
  B0AvgImageType::Pointer lastdivp;
  lastdivp->SetRegions(N);
  lastdivp->Allocate();

  double pnorm;
//
// if (length(N) == 3), which it always will be.
  do
    {
//     while (norm(divp(:) - lastdivp(:),inf) > Tol),
//         lastdivp = divp;
    // to avoid needless allocation, swap divp with lastdivp
    B0AvgImageType::Pointer tmpvp = divp;
    divp = lastdivp;
    lastdivp = tmpvp;

    // note: f doesn't change in loop so f*lambda done above and
    //         assigned to fLambda
    //         z = divp - f*lambda;
    SubtractFilterType::Pointer sub1 = SubtractFilterType::New();
    sub1->SetInput1(divp);
    sub1->SetInput2(fLambda);
    B0AvgImageType::Pointer z = sub1->GetOutput();
    // z1 = z(:,ir,:) - z;
    B0AvgImageType::Pointer z1 = shiftImage<B0AvgImageType>(z,1,ir);
    sub1->SetInput1(z1); sub1->SetInput2(z); sub1->Update();
    z1 = sub1->GetOutput();
    // z2 = z(id,:,:) - z;
    B0AvgImageType::Pointer z2 = shiftImage<B0AvgImageType>(z,0,id);
    sub1->SetInput1(z2); sub1->SetInput2(z); sub1->Update();
    z2 = sub1->GetOutput();
    // z3 = z(:,:,ib) - z;
    B0AvgImageType::Pointer z3 = shiftImage<B0AvgImageType>(z,2,ib);
    sub1->SetInput1(z3); sub1->SetInput2(z); sub1->Update();
    z3 = sub1->GetOutput();

    // denom = 1 + dt*sqrt(z1.^2 + z2.^2 + z3.^2);
    B0AvgImageType::Pointer denom = B0AvgImageType::New();
    denom->CopyInformation(f);
    denom->Allocate();
    itk::ImageRegionIterator<B0AvgImageType> denomIt(denom,denom->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<B0AvgImageType> z1It(z1,z1->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<B0AvgImageType> z2It(z2,z2->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator<B0AvgImageType> z3It(z3,z3->GetLargestPossibleRegion());

    for(denomIt.GoToBegin(), z1It.GoToBegin(),z2It.GoToBegin(),z3It.GoToBegin();
        !denomIt.IsAtEnd() && !z1It.IsAtEnd() && !z2It.IsAtEnd() && z3It.IsAtEnd();
        ++denomIt, ++z1It, ++z2It, ++z3It)
      {
      B0AvgImageType::PixelType z1Val = z1It.Get();
      B0AvgImageType::PixelType z2Val = z2It.Get();
      B0AvgImageType::PixelType z3Val = z3It.Get();
      denomIt.Set(1 + dt * std::sqrt((z1Val * z1Val) + (z2Val * z2Val)+ (z3Val * z3Val)));
      }
    itk::ImageRegionIterator<B0AvgImageType> p1It(p1,p1->GetLargestPossibleRegion());
    itk::ImageRegionIterator<B0AvgImageType> p2It(p2,p2->GetLargestPossibleRegion());
    itk::ImageRegionIterator<B0AvgImageType> p3It(p3,p3->GetLargestPossibleRegion());

     // p1 = (p1 + dt*z1)./denom;
     // p2 = (p2 + dt*z2)./denom;
     // p3 = (p3 + dt*z3)./denom;
    for(denomIt.GoToBegin(), p1It.GoToBegin(),p2It.GoToBegin(),p3It.GoToBegin(),
          z1It.GoToBegin(), z2It.GoToBegin(), z3It.GoToBegin();
        !denomIt.IsAtEnd() && !p1It.IsAtEnd() && !p2It.IsAtEnd() && p3It.IsAtEnd();
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
    itk::ImageRegionIteratorWithIndex<B0AvgImageType>
      divpIt(divp,divp->GetLargestPossibleRegion());
    for(divpIt.GoToBegin(); !divpIt.IsAtEnd(); ++divpIt)
      {
      B0AvgImageType::IndexType curIndex = divpIt.GetIndex();
      B0AvgImageType::IndexType pIndex;
      double curP1 = p1->GetPixel(curIndex);
      pIndex = curIndex;
      pIndex[1] = curIndex[il[1]];
      curP1 -= p1->GetPixel(pIndex);
      double curP2 = p2->GetPixel(curIndex);
      pIndex = curIndex;
      pIndex[0] = curIndex[iu[0]];
      curP2 -= p2->GetPixel(pIndex);
      double curP3 = p3->GetPixel(curIndex);
      pIndex = curIndex;
      pIndex[2] = curIndex[ifr[2]];
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
    SubtractFilterType::Pointer subFilter = SubtractFilterType::New();
    subFilter->SetInput(divp);
    subFilter->SetInput(lastdivp);
    AbsImageFilterType::Pointer absImageFilter = AbsImageFilterType::New();
    absImageFilter->SetInput(subFilter->GetOutput());
    StatsFilterType::Pointer statsFilter = StatsFilterType::New();
    statsFilter->SetInput(absImageFilter->GetOutput());
    statsFilter->Update();
    pnorm = statsFilter->GetMaximum();
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

DoCSEstimate::DWIVectorImageType *DoCSEstimate
::step2(double myu, double tol)
{
  unsigned int numGradients = this->m_IntensityData->GetNumberOfComponentsPerPixel();
  DWIVectorImageType::Pointer rval = DWIVectorImageType::New();
  rval->CopyInformation(this->m_IntensityData);
  rval->SetNumberOfComponentsPerPixel(numGradients);
  rval->Allocate();
  typedef DWIVectorImageType::SizeValueType sizetype;
  for(unsigned int k = 0; k < numGradients; ++k)
    {
    this->tvdenoise3(k,1/myu,tol,rval);
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
::step1(MatrixType &A,double lmd,unsigned NIT,std::vector<unsigned int> &id)
{
  // [nx, ny, nz, d]=size(S);
  ImageSizeType size = this->m_IntensityData->GetLargestPossibleRegion().GetSize();
  // n=nx*ny*nz;
  ImageSizeValueType n = size[0] * size[1] * size[2];
  // M=size(A,2);
  ImageSizeValueType M = A.cols();

  unsigned int numGradients = this->m_IntensityData->GetNumberOfComponentsPerPixel();
  // the reshape takes image of size [ x y z g ] and turns it into
  // [n g ] where n = x * y * z.  Matlab reshapes columnwise, so you
  // end up with G columns of n rows.
  MatrixType _S(n,numGradients);
  // S=reshape(S,[n d]);
  DWIVectorImageType::IndexType index = {{0,0,0}};
  for(unsigned int grad = 0; grad < numGradients; ++grad)
    {
    for(unsigned int z = 0, row = 0; z < size[2]; ++z)
      {
      index[2] = z;
      for(unsigned int x = 0; x < size[0]; ++x)
        {
        index[0] = x;
        for(unsigned int y = 0; y < size[1]; ++y,++row)
          {
          index[1] = y;
          const DWIVectorImageType::PixelType &current = this->m_IntensityData->GetPixel(index);
          _S(row,grad) = current[grad];
          }
        }
      }
    }
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
::BPDN_HOMOTOPY_function(MatrixType &A,MatrixType &y,double lmd, unsigned int NIT)
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

//
// % Initial step
// Primal_constrk = -A'*y;
  MatrixType Primal_constrk = -A * y;
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
