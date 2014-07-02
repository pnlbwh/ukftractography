#include "computedispersion.h"
#include <Eigen/Dense>

#include <math.h>
#include <algorithm>
#include <vector>
/*
%COMPUTE_DISPERSION
%
%     Author: Peter Savadjiev
%     If using this code, please cite:
%     P. Savadjiev, Y. Rathi, S. Bouix, R. Verma and C.-F. Westin, Multi-scale Characterization of White Matter
%     Tract Geometry, in: Medical Image Computing and Computer Assisted Intervention (MICCAI), volume LNCS
%     7512, pp. 34–41, 2012
%
%     Implements the fiber dispersion algorithm of Savadjiev et al., MICCAI 2012.
%
%     This function requires three compulsory input arguments.
%
%     DDF = compute_dispersion(TRACT, SCALE, NUMDIRS) takes as input a cell array
%     TRACT, where each cell is a 3-by-N matrix. The rows of this matrix
%     correspond to the x, y and z components of the point coordinates
%     of one fiber. N is the number of points in that fiber.
%
%     SCALE is a numerical value that specifies the scale at which to compute
%     dispersion.
%
%     NUMDIRS is the number of directions along which to sample the Dispersion
%     Distribution Function (DDF) in the plane orthogonal to the local tangent
%     vector.
%
%     DDF is a (NUMDIRS+1)-by-P matrix. P is the total number of points in TRACT.
%     The columns of DDF correspond to points of the fiber tract, and their
%     ordering is the same as produced by a call to cell2mat(TRACT).
%     The rows of DDF correspond to the DDF computed at a given fiber point.
%     The first NUMDIRS rows give the DDF value for each sampling direction,
%     The (NUMDIRS+1)-th row gives the median value of the DDF,
%     which corresponds to the TD measure of Savadjiev et al. MICCAI 2012.
%
%     In addition, the function permits up to three optional input arguments.
%
%     DDF = compute_dispersion(TRACT, SCALE, NUMDIRS, SUB_TR)
%     subsamples the computation across the set of fibers given in TRACT, so that
%     computation is performed on only one fiber out of every SUB_TR fibers.
%     SUB_TR must be a positive integer, default value = 1 (no subsampling).
%
%
%     DDF = compute_dispersion(TRACT, SCALE, NUMDIRS, SUB_TR, SUB_FB)
%     subsamples the computation along each fiber, so that computation is
%     performed on only one point out of every SUB_FB points along the fiber.
%     SUB_FB must be a positive integer, default value = 1 (no subsampling).
%
%
%     DDF = compute_dispersion(TRACT, SCALE, NUMDIRS, SUB_TR, SUB_FB, OUTPUTFILENAME)
%     saves the output DDF matrix to a .mat file whose filename is provided
%     by argument OUTPUTFILENAME.
%     Default behavior: the output is not written to a file.

NOTE: For reference, original matlab code has been appended as a comment

*/

/**
 * PrintMat
 * NOTE: during debugging, you can call this function to show contents
 * of Eigen matrices:
 * (gdb) call PrintMat(matrixName)
 */
void PrintMat(Eigen::Matrix<ukfPrecisionType,Eigen::Dynamic,Eigen::Dynamic>  &mat, std::ostream &outfile = std::cerr)
{
  Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
  outfile << mat.format(CleanFmt);
}

void PrintMat(Eigen::Matrix<ukfPrecisionType,Eigen::Dynamic,1>  &mat, std::ostream &outfile = std::cerr)
{
  Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
  outfile<< mat.format(CleanFmt);
}

namespace
{

typedef Eigen::Matrix<ukfPrecisionType,Eigen::Dynamic,Eigen::Dynamic> MatrixType;

typedef std::vector<MatrixType> MatrixVector;

typedef Eigen::Matrix<ukfPrecisionType,Eigen::Dynamic,1> VectorType;

/** convert vector of matrices into single matrix.
 * MATLAB cell2mat is more general, I think -- this does what
 * is actually needed for computedispersion.
 *
 * fibers are stored as NRows,numPoints matrices; the output
 * needs to be NRows,numPoints*numFibers
 */
void
cell2mat(const MatrixVector &fibers /*in*/,MatrixType &tractMatrix /*out*/)
{
  const unsigned int fibersSize(fibers.size());
  unsigned int numCols(0);
  for(unsigned i = 0; i < fibersSize; ++i)
    {
    numCols += fibers[i].cols();
    }
  tractMatrix.conservativeResize(fibers[0].rows(),numCols);
  for(unsigned int i = 0, offset = 0; i < fibersSize; ++i)
    {
    const MatrixType &curFiber = fibers[i];
    tractMatrix.block(0,offset,curFiber.rows(),curFiber.cols()) = curFiber;
    offset += curFiber.cols();;
    }
}

/** Principle Component Axes
 * This is from a suggestion here:
 * http://forum.kde.org/viewtopic.php?f=74&t=110265
 * by Gael Guennebaud, an Eigen library developer
 */
MatrixType PCA(const MatrixType &mat)
{
  MatrixType centered = mat.rowwise() - mat.colwise().mean();
  MatrixType cov = centered.adjoint() * centered;

  Eigen::SelfAdjointEigenSolver<MatrixType> eig(cov);

  VectorType eigenValues = eig.eigenvalues();
  MatrixType eigenVectors = eig.eigenvectors();

  MatrixType rval(mat.cols(),mat.cols());
  for(unsigned int i = 0; i < mat.cols(); ++i)
    {
    rval.row(i) = eigenVectors.row((mat.cols() - 1) - i).normalized();
    }
  return rval;
}

/** Go from std::vector<std::vector<vec3_t> to std::vector<MatrixType>
 * the input fiber type is defined in fiberbundle as a vector of
 * vector of points.  for local computation we want a vector of
 * matrices instead, of dimension 3 x fiverSize.
 */
void
FiberVector2EigenFiberVector(const fiberbundle::FiberVector &fv /* in */,
                             MatrixVector &lv)
{
  const unsigned long fiberSize(fv.size());

  lv.clear();                   // possibly reusing existing vector
  lv.resize(fiberSize);         // same size as input fiberbundle

  for(unsigned long i = 0; i < fiberSize; ++i)
    {
    const stdVec_t &curFiber = fv[i].Points;
    const unsigned long curFiberSize(curFiber.size());

    MatrixType &curLocalFiber = lv[i];
    curLocalFiber.conservativeResize(3, curFiberSize);
    for(unsigned int j = 0; j < curFiberSize; ++j)
      {
      curLocalFiber.col(j) = curFiber[j];
      }
    }
}

/** FRAME  Compute a frame field along a curve
 *        Author: Peter Savadjiev
 *        If using this code, please cite:
 *        P. Savadjiev, Y. Rathi, S. Bouix, R. Verma and C.-F. Westin, Multi-scale Characterization of White Matter
 *        Tract Geometry, in: Medical Image Computing and Computer Assisted Intervention (MICCAI), volume LNCS
 *        7512, pp. 340-341, 2012
 *    [T, N, B] = frame(X, Y, Z, V)
 *    X, Y and Z are each 1-by-M vectors. They must respectively contain the
 *    x, y and z point coordinates of the curve points. M is the number
 *    of points along the curve.
 *    V is a 3-by-1 vector.
 *    T, N and B are 3-by-M matrices.
 *    This function returns the tangent, normal and binormal vector fields
 *    along the curve in output vectors T, N and B, respectively.
 *    Note the normal and binormal vectors are not necessarily those defined
 *    by Frenet's equations. Since the Frenet frame is quite unstable in
 *    practice - Frenet's binormal vector field is discontinuous at inflection
 *    points of the curve and vanishes on straight line segments - we use the
 *    input vector V to constrain the binormal vector to always lie in the
 *    plane defined by the tangent vector and the constraint vector V.
 *    This results in a much smoother behavior of the binormal (and hence also
 *    the normal) vector fields over the curve.
 * ORIGINAL MATLAB CODE BELOW
 */
class Frame
{
public:
  typedef Eigen::Matrix<ukfPrecisionType,3,1> FixedVectorType;
  Frame(const VectorType &xCoords,const VectorType &yCoords,const VectorType &zCoords,
        const VectorType &constraintVector) : m_Xcoords(xCoords),
                                              m_Ycoords(yCoords),
                                              m_Zcoords(zCoords),
                                              m_ConstraintVector(constraintVector)
    {
    }
  void compute();
  const MatrixType &GetTangent()  const { return m_Tangent; }
  const MatrixType &GetNormal()   const { return m_Normal; }
  const MatrixType &GetBinormal() const { return m_Binormal; }
private:
  VectorType       m_Xcoords;
  VectorType       m_Ycoords;
  VectorType       m_Zcoords;
  FixedVectorType  m_ConstraintVector;
  MatrixType       m_Tangent;
  MatrixType       m_Normal;
  MatrixType       m_Binormal;
};

void
Frame
::compute()
{
  unsigned int curveLength = this->m_Xcoords.rows();
  if(curveLength == 1)
    {
    this->m_Xcoords = this->m_Xcoords.transpose();
    this->m_Ycoords = this->m_Ycoords.transpose();
    this->m_Zcoords = this->m_Zcoords.transpose();
    curveLength = this->m_Xcoords.rows();
    }

  if(curveLength == 0)
    {
    this->m_Tangent(0) = 1.0;
    this->m_Tangent(1) = 0.0;
    this->m_Tangent(2) = 0.0;

    this->m_Normal(0) = 0.0;
    this->m_Normal(1) = 1.0;
    this->m_Normal(2) = 0.0;

    this->m_Binormal(0) = 0.0;
    this->m_Binormal(1) = 0.0;
    this->m_Binormal(2) = 1.0;
    return;
    }
  this->m_Tangent.resize(curveLength,3);
  this->m_Normal.resize(curveLength,3);
  this->m_Binormal.resize(curveLength,3);
  this->m_Tangent.setZero();
  this->m_Normal.setZero();
  this->m_Binormal.setZero();

  MatrixType points(curveLength,3);
  points.col(0) = this->m_Xcoords;
  points.col(1) = this->m_Ycoords;
  points.col(2) = this->m_Zcoords;

  for(unsigned i = 1; i < curveLength - 1; ++i)
    {
    this->m_Tangent.row(i) = (points.row(i+1) - points.row(i-1))/2;
    double tangentNorm = this->m_Tangent.row(i).norm();
    if(tangentNorm > 0.00001)
      {
      this->m_Tangent.row(i).normalize();
      }
    else
      {
      this->m_Tangent.row(i) = this->m_Tangent.row(i - 1);
      }
    }


  this->m_Tangent.row(0) = points.row(1) - points.row(0);
  double tangentNorm = this->m_Tangent.row(0).norm();
  if(tangentNorm > 0.00001)
    {
    this->m_Tangent.row(0).normalize();
    }
  else
    {
    this->m_Tangent.row(0) = this->m_Tangent.row(1);
    }

  this->m_Tangent.row(curveLength - 1) = points.row(curveLength - 1) -
    points.row(curveLength - 2);
  tangentNorm = this->m_Tangent.row(curveLength - 1).norm();
  if(tangentNorm > 0.00001)
    {
    this->m_Tangent.row(curveLength -1).normalize();
    }
  else
    {
    this->m_Tangent.row(curveLength -1) = this->m_Tangent.row(curveLength - 2);
    }



  for(unsigned i = 0; i < curveLength; ++i)
    {
    FixedVectorType curTanRow = this->m_Tangent.row(i);
    FixedVectorType orthogonalVector1 = this->m_ConstraintVector.cross(curTanRow);
    if(orthogonalVector1.norm() < 0.00001)
      {
      this->m_ConstraintVector *= 0.00002;
      this->m_ConstraintVector.normalize();
      orthogonalVector1 = this->m_ConstraintVector.cross(curTanRow);
      }
    FixedVectorType orthogonalVector2 = curTanRow.cross(orthogonalVector1).normalized();
    if(orthogonalVector2.dot(this->m_ConstraintVector) < 0.00001)
      {
      this->m_Binormal.row(i) = orthogonalVector2 * -1.0;
      }
    else
      {
      this->m_Binormal.row(i) = orthogonalVector2;
      }
    }

  for(unsigned int i = 0; i < curveLength; ++i)
    {
    FixedVectorType curTanRow = this->m_Tangent.row(i);
    FixedVectorType curBinorm = this->m_Binormal.row(i);
    this->m_Normal.row(i) = curBinorm.cross(curTanRow);
    }
}

/** RotateField
 *  Rotates the tract by rotation matrix R. As a result, the entire fiber tract
 *  is represented in the local coordinate frame at point currentPosition, such that the
 *  local fiber tangent vector is the x axis.
 *  Returns the rotated point coordinates in rotatedPointCoordinates, and their corresponding
 *  rotated tangent vectors in rotatedTangentVectorField.
 */
class RotateField
{
public:
  RotateField(const MatrixType &tractMatrix,
              const MatrixType &rotationMatrix,
              const MatrixType &currentPosition) :
    m_TractMatrix(tractMatrix),
    m_RotationMatrix(rotationMatrix),
    m_CurrentPosition(currentPosition)
    {
    }
  void compute();
  const MatrixType &GetRotatedPointCoordinates() const
    { return m_RotatedPointCoordinates; }
  const MatrixType &GetRotatedTangentVectorField() const
    { return m_RotatedTangentVectorField; }

private:
  MatrixType m_RotatedPointCoordinates;
  MatrixType m_RotatedTangentVectorField;
  const MatrixType &m_TractMatrix;
  const MatrixType &m_RotationMatrix;
  const MatrixType &m_CurrentPosition;
};

void
RotateField
::compute()
{
  const unsigned long tractCols(this->m_TractMatrix.cols());
  MatrixType currentPositionMatrix(3,tractCols);
  for(unsigned int i = 0; i < this->m_TractMatrix.cols(); ++i)
    {
    currentPositionMatrix.col(i) = this->m_CurrentPosition;
    }
  MatrixType pointCoordinates = m_TractMatrix.block(0,0,3,m_TractMatrix.cols()) - currentPositionMatrix;
  this->m_RotatedPointCoordinates = this->m_RotationMatrix * pointCoordinates;
  this->m_RotatedPointCoordinates += currentPositionMatrix;
  this->m_RotatedTangentVectorField =
    this->m_RotationMatrix * this->m_TractMatrix.block(3,0,3,tractCols);
}



/** computeMeanVector
 *  Computes the average vector within a disc-shaped neighborhood with
 *  radius scale, centered at location currentPosition. See Savadjiev et al. MICCAI
 *  2012, equation 2, and the related discussion in Section 3.1.
 */
bool computeMeanVector(const MatrixType &xCoordinates,
                      const MatrixType &yCoordinates,
                      const MatrixType &zCoordinates,
                      const MatrixType &tangentVectorField,
                      const MatrixType &currentPosition,
                      double scale,
                      Eigen::Vector3d &meanVector /* out */)
{
  const double eps(2.2204e-16); // that's what matlab thinks eps is.


  std::vector<bool> indexPointsInPlane(xCoordinates.cols(),false);
  unsigned int inPlanePointCount = 0;
  for(unsigned int i = 0; i < xCoordinates.cols(); ++i)
    {
    if(fabs(xCoordinates(i) - currentPosition(0)) < 0.5)
      {
      indexPointsInPlane[i] = true;
      ++inPlanePointCount;
      }
    }

  //
  // this implements the matlab craziness: A(b) where b
  // is a logical vector is a vector of all the points Z where b[Z] is true
  MatrixType pointCoordinates(3,inPlanePointCount);
  for(unsigned int i = 0, curPoint = 0; i < xCoordinates.cols(); ++i)
    {
    if(indexPointsInPlane[i])
      {
      pointCoordinates(0,curPoint) = xCoordinates(i) - currentPosition(0);
      pointCoordinates(1,curPoint) = yCoordinates(i) - currentPosition(1);
      pointCoordinates(2,curPoint) = zCoordinates(i) - currentPosition(2);
      ++curPoint;
      }
    }

  MatrixType pointDistance(1,inPlanePointCount);

  // indexPointsInDisk = pointDistance < scale*scale;
  std::vector<bool> indexPointsInDisk(inPlanePointCount);

  unsigned int inDiskPointCount = 0;
  for(unsigned int i = 0; i < inPlanePointCount; ++i)
    {
    double dotProd = pointCoordinates.col(i).dot(pointCoordinates.col(i));
    pointDistance(i) = dotProd;

    if(pointDistance(i) < scale * scale)
      {
      indexPointsInDisk[i] = true;
      ++inDiskPointCount;
      }
    else
      {
      indexPointsInDisk[i] = false;
      }
    }

  double n = 0.0;
  // this is a rearrangement of the matlab code;
  // if there are no points in the disk, there's no point
  // to computing the mean.
  if(inDiskPointCount > 10)
    {
    // vectorFieldInDisk = vectorFieldInDisk(:,indexPointsInDisk);
    // use second variable rather than overwrite in place
    // vectorFieldInDisk = tangentVectorField(:,indexPointsInPlane);
    MatrixType vectorFieldInDisk(3,inPlanePointCount);
    for(unsigned int i = 0, curPoint = 0; i < tangentVectorField.cols(); ++i)
      {
      if(indexPointsInPlane[i])
        {
        vectorFieldInDisk.col(curPoint) = tangentVectorField.col(i);
        curPoint++;
        }
      }
    MatrixType vectorFieldInDisk2(3,inDiskPointCount);
    for(unsigned int i = 0, curPoint = 0; i < inPlanePointCount; ++i)
      {
      if(indexPointsInDisk[i])
        {
        vectorFieldInDisk2.col(curPoint) = vectorFieldInDisk.col(i);
        ++curPoint;
        }
      }
    // indexNegativeOrientations = vectorFieldInDisk(1,:) < -eps;
    // vectorFieldInDisk(:,indexNegativeOrientations) = -vectorFieldInDisk(:,indexNegativeOrientations);
    // in other words, any  negative points multiply by scalar -1
    for(unsigned int i = 0; i < inDiskPointCount; ++i)
      {
      if(vectorFieldInDisk2(0,i) < -eps)
        {
        vectorFieldInDisk2.col(i) *= -1.0;
        }
      }
    // meanVector = mean(vectorFieldInDisk,2);
    meanVector(0) = vectorFieldInDisk2.row(0).mean();
    meanVector(1) = vectorFieldInDisk2.row(1).mean();
    meanVector(2) = vectorFieldInDisk2.row(2).mean();
    n = meanVector.norm();
    }
  //
  // the MATLAB code would take the mean of an empty matrix, which
  // resulted in a mean full of NaNs; my solution is to not try to
  // compute the mean at all, as Eigen doesn't like operations on
  // empty matrices.
  // if(~any(isnan(meanVector)) && (length(indexPointsInDisk) > 10) )
  //     n = norm(meanVector);
  //     if(n>0.1)
  if(inDiskPointCount > 10 && n > 0.1)
    {
    meanVector /= n;
    return true;
    }
  // if there were not enough points in the plane or no points in the
  // disk, then just return status of zero to indicate that.
  meanVector = Eigen::Vector3d::Zero();
  return false;
}

double median(std::vector<double> &vec)
{
  double rval;
  size_t size = vec.size();
  std::sort(vec.begin(),vec.end());
  if(size % 2 == 0)
    {
    rval = (vec[size / 2 - 1] + vec[size / 2]) / 2;
    }
  else
    {
    rval = vec[size/2];
    }
  return rval;
}

} // local namespace

int
computedispersion(fiberbundle &bundle, double scale,
                  unsigned int numberOfSamplingDirections,
                  const std::string &outputFilename,
                  unsigned int  tractSubSampling,
                  unsigned int fiberPointSubSampling)
{
  // bundle.Print();
   // according to stackoverflow, this is how to get PI
  const double Pi(std::atan(static_cast<double>(1.0))*4);

  fiberbundle::FiberVector &fibers = bundle.GetFibers();

  // vector of native eigen matrices
  MatrixVector localFibers;

  FiberVector2EigenFiberVector(fibers,localFibers);

  MatrixType tractMatrix; // 3 x x matrix accumulating all fibers tracts
  cell2mat(localFibers,tractMatrix);


  MatrixType princomp = PCA(tractMatrix.transpose());
  VectorType constraintVector = princomp.col(2);
  // compute dispersion for each point along each fiber
  // NOTE: if tractSubSampling or fiberPointSubSampling is > 1, then
  // it breaks the assumption of one dispersion per point
  //
  // NOTE: this follows matlab which stores the tan/norm/binorm after
  // the points in the curPoints matrix.  They could as easily be kept
  // separate, and it might be clearer code.
  for(unsigned i = 0; i < fibers.size(); ++i)
    {
    MatrixType &curPoints = localFibers[i];

    const unsigned int curPointsSize(curPoints.cols());

    curPoints.conservativeResize(12,curPointsSize);

    VectorType x = curPoints.row(0);
    VectorType y = curPoints.row(1);
    VectorType z = curPoints.row(2);

    Frame frame(x,y,z,constraintVector);
    frame.compute();

    curPoints.block(3,0,3,curPointsSize) = frame.GetTangent().transpose();
    curPoints.block(6,0,3,curPointsSize) = frame.GetNormal().transpose();
    curPoints.block(9,0,3,curPointsSize) = frame.GetBinormal().transpose();
    }

  cell2mat(localFibers,tractMatrix);


  // Subsampling the tract by keeping every *tractSubsampling*-th fiber %%%
  MatrixVector subSampledTract;
  for(unsigned int i = 0; i < localFibers.size(); i += tractSubSampling)
    {
    subSampledTract.push_back(localFibers[i]);
    }

  MatrixType subSampledTractMatrix;
  cell2mat(subSampledTract,subSampledTractMatrix);

  //  Initializing the DDF matrix and computing the sampling directions in a
  //  standard coordinate frame such that the sampling directions lie in the
  //  YZ plane.

  // initialize DDF matrix with -1s -- sentinal value for when there's
  // not enough points in the plane/disk to give meaninful result.
  MatrixType dispersionDistributionValues =
    MatrixType::Ones(numberOfSamplingDirections + 1, subSampledTractMatrix.cols()) * -1;

  MatrixType samplingDirections = MatrixType::Zero(3,numberOfSamplingDirections);

  // compute the sampling directions
  double theta = (2.0 * Pi)/static_cast<double>(numberOfSamplingDirections);
  for(unsigned int j = 1; j <= numberOfSamplingDirections; ++j)
    {
    samplingDirections(1,j-1) = cos(j*theta);
    samplingDirections(2,j-1) = sin(j*theta);
    }
  // Computing the DDF at every *fiberPointSubsampling*-th point along each fiber of
  // the (possibly subsampled) tract.
  for(unsigned i = 0;  i < subSampledTractMatrix.cols(); i += fiberPointSubSampling)
    {
    MatrixType rotationMatrix(3,3);
    rotationMatrix.row(0) = subSampledTractMatrix.block(3,i,3,1).transpose();
    rotationMatrix.row(1) = subSampledTractMatrix.block(6,i,3,1).transpose();
    rotationMatrix.row(2) = subSampledTractMatrix.block(9,i,3,1).transpose();

    MatrixType currentPosition = subSampledTractMatrix.block(0,i,3,1);

    RotateField rotateField(tractMatrix,rotationMatrix,currentPosition);
    rotateField.compute();
    const MatrixType &rotatedPointCoordinates =
      rotateField.GetRotatedPointCoordinates();
    const MatrixType &rotatedTangentVectorField =
      rotateField.GetRotatedTangentVectorField();

    MatrixType xCoordinates = rotatedPointCoordinates.row(0);
    MatrixType yCoordinates = rotatedPointCoordinates.row(1);
    MatrixType zCoordinates = rotatedPointCoordinates.row(2);

    Eigen::Vector3d referenceMeanVector;

    if(computeMeanVector(xCoordinates,
                         yCoordinates,
                         zCoordinates,
                         rotatedTangentVectorField,
                         currentPosition,
                         scale,
                         referenceMeanVector))
      {
      MatrixType samplingPosition(3,samplingDirections.cols());
      for(unsigned j = 0; j < samplingDirections.cols(); ++j)
        {
        samplingPosition.col(j) = currentPosition + (scale * samplingDirections.col(j));
        }
      for(unsigned int j = 0; j < samplingDirections.cols(); ++j)
        {
        Eigen::Vector3d meanVector;
        if(computeMeanVector(xCoordinates,
                             yCoordinates,
                             zCoordinates,
                             rotatedTangentVectorField,
                             samplingPosition.col(j),
                             scale,
                             meanVector))
          {
          double dot = meanVector.dot(referenceMeanVector);
          double acosDot = acos(dot);
          if(acosDot < 0.0)
            {
            acosDot *= -1.0;
            }
          dispersionDistributionValues(j,i) = acosDot;
          }
        }
      // take the median of the computed means
      MatrixType pointDDF = dispersionDistributionValues.block(0,i,numberOfSamplingDirections,1);
      std::vector<double> nonNegDDF;
      for(unsigned int j = 0; j < pointDDF.rows(); ++j)
        {
        if(pointDDF(j,0) != -1)
          {
          nonNegDDF.push_back(pointDDF(j,0));
          }
        }
      if(nonNegDDF.size() > 0)
        {
        dispersionDistributionValues(numberOfSamplingDirections,i) = median(nonNegDDF);
        }
      }
    }

  MatrixType DDFOutput = dispersionDistributionValues.row(numberOfSamplingDirections);
  //
  // So the 'punt' to handle sub-sampling is that any fibers skipped 
  for(unsigned int i = 0, curPoint = 0; i < fibers.size(); ++i)
    {
    Fiber &curFiber = fibers[i];
    std::vector<float> curDDF;
    // if this fiber was sampled, print out DDF at points for this fiber
    if(i % tractSubSampling == 0)
      {
      for(unsigned int j = 0; j < curFiber.Points.size(); ++j)
        {
        curDDF.push_back(DDFOutput(curPoint));
        ++curPoint;
        }
      }
    // if this fiber was skipped, then output -1s for each point.
    else
      {
      for(unsigned j = 0; j < curFiber.Points.size(); ++j)
        {
        curDDF.push_back(-1.0);
        }
      }
    curFiber.Fields["DDF"] = curDDF;
    }

  if(outputFilename != "")
    {
    std::ofstream outfile(outputFilename.c_str());
    PrintMat(DDFOutput, outfile);
    outfile.close();
    }
  return 0; // success
}

// ORIGINAL MATLAB CODE: COMPUTE DISPERSION
// function dispersionDistributionValues = compute_dispersion(tract, scale, numberOfSamplingDirections, varargin)
//
//     %%% Checking the optional input arguments; if unspecified,  setting
//     %%% their default values
//
//     tractSubsampling = 1;
//     fiberPointSubsampling = 1;
//     outputFilename = 0;
//
//     optargin = size(varargin,2);
//
//     if optargin == 1
//         tractSubsampling = varargin{1};
//     elseif optargin == 2
//         tractSubsampling = varargin{1};
//         fiberPointSubsampling = varargin{2};
//     elseif optargin == 3
//         tractSubsampling = varargin{1};
//         fiberPointSubsampling = varargin{2};
//         outputFilename = varargin{3};
//     else
//         fprintf('Invalid number of input arguments. Exiting.');
//         return
//     end
//
//     %%% Setting up a frame field along each fiber in the tract %%%
//     tractMatrix = cell2mat(tract);
//     % tractMatrix flattens the cell structure to a simple rectangular
//     % matrix that is 3x242 (x,y,z) X number of points in all tracts
//
//     % TODO: Kent Need to get princomp from eigenvalues
//     % The notes on this page: http://pastebin.com/2nU44TQp may provide a
//     % solution to the matlab function princomp
//     % PETER:  How do go from svd to pca?
//     % [varargout{1:nargout}]=pca(varargin{1},'Algorithm','svd','Economy',fEconomy);
//     % re-write the convience princomp to using the elmental eigen
//     % equivalents.
//
//     principalComponents = princomp([tractMatrix(1,:)' tractMatrix(2,:)' tractMatrix(3,:)']);
//
//     constraintVector = principalComponents(:,3);
//
//     for i = 1:length(tract)
//
//         x = tract{i}(1,:);
//         y = tract{i}(2,:);
//         z = tract{i}(3,:);
//
//
//         % TODO: KENT: frame is in frame.m
//         [tangent,normal,binormal] = frame(x,y,z,constraintVector);
//
//         tract{i}(4:6,:) =   tangent';
//         tract{i}(7:9,:) =   normal';
//         tract{i}(10:12,:) = binormal';
//
//     end
//
//     tractMatrix = cell2mat(tract);
//
//     %%% Subsampling the tract by keeping every *tractSubsampling*-th fiber %%%
//
//     subsampledTract = tract;
//
//     b = 1:length(tract);
//     a = 1:tractSubsampling:length(tract);
//     c = b(~ismember(b,a));
//     for i=c
//         subsampledTract{i}=[];
//     end
//     emptyCells = cellfun(@isempty,subsampledTract);
//     subsampledTract(emptyCells) = [];
//
//     subsampledTractMatrix = cell2mat(subsampledTract);
//
//
//     %%% Initializing the DDF matrix and computing the sampling directions in a
//     %%% standard coordinate frame such that the sampling directions lie in the
//     %%% YZ plane.
//
//     dispersionDistributionValues = -1*ones(numberOfSamplingDirections + 1,size(subsampledTractMatrix, 2));
//
//     samplingDirections = zeros(3,numberOfSamplingDirections);
//     for j=1:numberOfSamplingDirections
//         theta = 2*pi/numberOfSamplingDirections;
//         dir = [cos(j*theta), sin(j*theta)];
//         samplingDirections(2,j) = dir(1);
//         samplingDirections(3,j) = dir(2);
//     end
//
//     %%% Computing the DDF at every *fiberPointSubsampling*-th point along each fiber of
//     %%% the (possibly subsampled) tract.
//
//
//     for i = 1:fiberPointSubsampling:size(subsampledTractMatrix,2)
//
//         rotationMatrix = [subsampledTractMatrix(4,i) subsampledTractMatrix(5,i) subsampledTractMatrix(6,i);
//                           subsampledTractMatrix(7,i) subsampledTractMatrix(8,i) subsampledTractMatrix(9,i);
//                           subsampledTractMatrix(10,i) subsampledTractMatrix(11,i) subsampledTractMatrix(12,i);];
//
//         currentPosition = [subsampledTractMatrix(1,i);
//                            subsampledTractMatrix(2,i);
//                            subsampledTractMatrix(3,i)];
//
//         [rotatedPointCoordinates,rotatedTangentVectorField] = rotateField(tractMatrix, rotationMatrix, currentPosition);
//
//         xCoordinates = rotatedPointCoordinates(1,:);
//         yCoordinates = rotatedPointCoordinates(2,:);
//         zCoordinates = rotatedPointCoordinates(3,:);
//
//         [referenceMeanVector, statusFlag1] = computeMeanVector(xCoordinates,yCoordinates,zCoordinates,rotatedTangentVectorField,currentPosition,scale);
//
//         if(statusFlag1)
//
//             samplingPosition = currentPosition(:, ones(1,size(samplingDirections,2))) + scale * samplingDirections;
//
//             for j = 1:size(samplingDirections,2)
//                 [meanVector, statusFlag2] = computeMeanVector(xCoordinates,yCoordinates,zCoordinates,rotatedTangentVectorField,samplingPosition(:,j),scale);
//                 if(statusFlag2)
//                     dispersionDistributionValues(j,i) = abs(acos(dot(meanVector,referenceMeanVector)));
//                 end
//             end
//
//             pointDDF = dispersionDistributionValues(1:numberOfSamplingDirections, i);
//             pointDDF = pointDDF(pointDDF ~=  -1);
//
//             if(~isempty(pointDDF))
//                 dispersionDistributionValues(numberOfSamplingDirections + 1, i) = median(pointDDF);
//             end
//         end
//
//     end
//
//     dispersionDistributionValues = dispersionDistributionValues(numberOfSamplingDirections + 1, :);
//     if(outputFilename)
//         save(outputFilename, 'dispersionDistributionValues', '-v7.3');
//     end
//
// end
//
//
//
//
// function [rotatedPointCoordinates, rotatedTangentVectorField] = rotateField(tract,rotationMatrix,currentPosition)
//
// %%% Rotates the tract by rotation matrix R. As a result, the entire fiber tract
// %%% is represented in the local coordinate frame at point currentPosition, such that the
// %%% local fiber tangent vector is the x axis.
// %%% Returns the rotated point coordinates in rotatedPointCoordinates, and their corresponding
// %%% rotated tangent vectors in rotatedTangentVectorField.
//
// currentPositionMatrix =     currentPosition(:,ones(1 , size(tract,2)));
// pointCoordinates =          tract(1:3,:) - currentPositionMatrix;
// rotatedPointCoordinates =   rotationMatrix * pointCoordinates;
// rotatedPointCoordinates =   rotatedPointCoordinates + currentPositionMatrix;
// rotatedTangentVectorField = rotationMatrix * tract(4:6,:);
//
// end
//
//
// function [meanVector, statusFlag] = computeMeanVector(xCoordinates,yCoordinates,zCoordinates,tangentVectorField,currentPosition,scale)
//
// %%% Computes the average vector within a disc-shaped neighborhood with
// %%% radius scale, centered at location currentPosition. See Savadjiev et al. MICCAI
// %%% 2012, equation 2, and the related discussion in Section 3.1.
//
// indexPointsInPlane = abs(xCoordinates-currentPosition(1)) < 0.5;
// pointCoordinates = [xCoordinates(indexPointsInPlane) - currentPosition(1) ;
//                     yCoordinates(indexPointsInPlane) - currentPosition(2) ;
//                     zCoordinates(indexPointsInPlane) - currentPosition(3)];
// pointDistance = dot(pointCoordinates,pointCoordinates);
// indexPointsInDisk = pointDistance < scale*scale;
//
// vectorFieldInDisk = tangentVectorField(:,indexPointsInPlane);
// vectorFieldInDisk = vectorFieldInDisk(:,indexPointsInDisk);
// indexNegativeOrientations = vectorFieldInDisk(1,:) < -eps;
// vectorFieldInDisk(:,indexNegativeOrientations) = -vectorFieldInDisk(:,indexNegativeOrientations);
//
//
// meanVector = mean(vectorFieldInDisk,2);
//
// if(~any(isnan(meanVector)) && (length(indexPointsInDisk) > 10) )
//     n = norm(meanVector);
//     if(n>0.1)
//         meanVector = meanVector./n;
//         statusFlag=1;
//     else
//         meanVector = 0;
//         statusFlag = 0;
//     end
// else
//     meanVector = 0;
//     statusFlag = 0;
// end
//
// end
//
// ORIGINAL MATLAB CODE FRAME
//
// %FRAME  Compute a frame field along a curve
// function [tangent,normal,binormal]=frame(xCoordinates,yCoordinates,zCoordinates,constraintVector)
//
// curveLength = size(xCoordinates,1);
// if (curveLength==1)
//   xCoordinates = xCoordinates';
//   yCoordinates = yCoordinates';
//   zCoordinates = zCoordinates';
//   curveLength = size(xCoordinates,1);
// end
//
// if(curveLength ==1)
//     tangent  = [1 0 0];
//     normal   = [0 1 0];
//     binormal = [0 0 1];
//     return
// end
//
// tangent  = zeros(curveLength,3);
// normal   = zeros(curveLength,3);
// binormal = zeros(curveLength,3);
//
// points = [xCoordinates yCoordinates zCoordinates];
//
// for i = 2:(curveLength-1)
//   tangent(i,:) = (points(i+1,:)-points(i-1,:))/2;
//   tangentNorm = norm(tangent(i,:));
//   if (tangentNorm > eps)
//     tangent(i,:) = tangent(i,:)/tangentNorm;
//   else
//     tangent(i,:) = tangent(i-1,:);
//   end
//
// end
//
//
// tangent(1,:) = points(2,:)-points(1,:);
// tangentNorm = norm(tangent(1,:));
// if(tangentNorm > eps)
//     tangent(1,:) = tangent(1,:)/tangentNorm;
// else
//     tangent(1,:) = tangent(2,:);
// end
//
// tangent(curveLength,:) = points(curveLength,:) - points(curveLength-1,:);
// tangentNorm = norm(tangent(curveLength,:));
// if(tangentNorm > eps)
//     tangent(curveLength,:) = tangent(curveLength,:)/tangentNorm;
// else
//     tangent(curveLength,:) = tangent(curveLength-1,:);
// end
//
//
//
// for i = 1:curveLength
//
//     orthogonalVector1 = cross(constraintVector, tangent(i,:));
//     if(norm(orthogonalVector1) < eps)
//         constraintVector = (constraintVector + 2*eps)/norm(constraintVector + 2*eps);
//         orthogonalVector1 = cross(constraintVector, tangent(i,:));
//     end
//     orthogonalVector2 = cross(tangent(i,:),orthogonalVector1);
//     orthogonalVector2 = orthogonalVector2/norm(orthogonalVector2);
//     if(dot(orthogonalVector2,constraintVector) < -eps)
//         binormal(i,:) = -1*orthogonalVector2;
//     else
//         binormal(i,:) = orthogonalVector2;
//     end
//
// end
//
//
// for i = 1:curveLength
//
//         normal(i,:) = cross(binormal(i,:),tangent(i,:));
//
// end
//
// end



