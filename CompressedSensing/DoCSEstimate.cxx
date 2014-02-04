#include "DoCSEstimate.h"

DoCSEstimate
::DoCSEstimate(NrrdFile &nrrdFile,MaskImageType::Pointer &maskImage,MatrixType &newGradients) :
  m_NrrdFile(nrrdFile), m_MaskImage(maskImage), m_Gradients(newGradients), m_Eps(2.2204e-16),
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
    if(norm < this->m_Eps)
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
    if(std::fabs(avg) <= this->m_Eps)
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
    this->m_VoxelLatticeAlignedGradientDirections.row(i) = this->m_Gradients.row(gradientIndices[i]);
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
    if(std::fabs(curB0) > this->m_Eps)
      {
      for(unsigned int i = 0; i < numGrad; ++i)
        {
        vec[i] /= curB0;
        if(vec[i] < 0)
          {
          vec[i] = this->m_Eps;
          }
        }
      }
    else
      {
      for(unsigned int i = 0; i < numGrad; ++i)
        {
        vec[i] = this->m_Eps;
        }
      }
    ivIt.Set(vec);
    }

  MatrixType estimatedGradients = this->m_Gradients;
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
  // [v,M]=multisample(J);
  // A=buildsensor(gradientDirections,v,psi);
  // A0=buildsensor(new_gradients,v,psi);
  //
  // %parameter setup
  // lmd=0.06;                   % Lagrangian L1
  // myu=0.05;                   % Lagrangian TV
  // gama=0.5;                   % Bregman parameter
  // NIT=2000;                   % Max number of FISTA iterations
  // tol=1e-3;                   % Relative change in Chambolle
  //
  // id=find(mask~=0);
  //
  // matlabpool('open');
  //
  // u=step2(DWIIntensityData,myu,tol);
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

  // S=exp(-b*sum((u*D).*u,2));
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
  //Lmd(2:2:m+1)=0;
  //for k=2:2:m,
  // Lmd(k+1)=-Lmd(k-1)*(k-1)/k;
  // end
  MatrixType Lmd(m+1,1);
  Lmd = MatrixType::Ones(m+1,1);
  for(unsigned int i = 1; i < m+1; i += 2)
    {
    Lmd(i) = 0;
    }
  for(unsigned int k = 2; k <= m; k += 2)
    {
    // Lmd(k+1)=-Lmd(k-1)*(k-1)/k;
    Lmd(k) = -Lmd(k-2) * (k -1)/k;
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
  // phi=phi./phi(1);
  for(unsigned i = 0; i < m+1; i += 2)
    {
    phi(i,0) = h(i,0) / Lmd(i,0);
    }
  phi /= phi(0,0);
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
