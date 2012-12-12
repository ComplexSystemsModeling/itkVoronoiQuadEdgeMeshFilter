#ifndef __itkInCircleTestFunction_h
#define __itkInCircleTestFunction_h

#include <itkFunctionBase.h>

#include <itkMatrix.h>
#include "vnl/vnl_det.h"
#include "vnl/vnl_matrix_fixed.h"

namespace itk
{

//-------------------------------------------------------------------
// In Circle test wrapper for ITK
//
template
  <
  class MeshType,
  class OutputType =MeshType
  >
class ITK_EXPORT InCircleTestFunction :
      public FunctionBase< std::vector< typename MeshType::PointType>, bool >
{

public:
  
  /** Standart class typedefs **/
  typedef InCircleTestFunction                                             Self;
  typedef FunctionBase< std::vector< typename MeshType::PointType>, bool > Superclass;
  typedef SmartPointer< Self >                                             Pointer;
  typedef SmartPointer< const Self >                                       ConstPointer;
    
  /** Run-time type information (and related methods) **/
  itkNewMacro( Self );
  itkTypeMacro( InCircleTestFunction, FunctionBase );

  typedef MeshType Mesh;
  typedef typename Mesh::PointType PointType;
  typedef typename MeshType::PixelType PixelType;
 
  typedef std::vector< PointType > PointList;
  
  bool Evaluate( const PointList & pointList ) const
    {
    double det;
    // orientation test - determination of the sign of ad-bc
    double a = pointList[0][0] - pointList[2][0]; //ax-cx;
    double b = pointList[0][1] - pointList[2][1]; //ay-cy;
    double c = pointList[1][0] - pointList[2][0]; //bx-cx;
    double d = pointList[1][1] - pointList[2][1]; //by-cy;
    double orientation = a*d-b*c;
    
    typedef itk::Matrix< double, 4,4 > MatrixType;
    MatrixType M;
    
    M(0,0) = (double) pointList[0][0];
    M(0,1) = (double) pointList[0][1];
    M(0,2) = (double) pointList[0][0]*pointList[0][0]+pointList[0][1]*pointList[0][1];
    M(0,3) = (double) 1.0;
    M(1,0) = (double) pointList[1][0];
    M(1,1) = (double) pointList[1][1];
    M(1,2) = (double) pointList[1][0]*pointList[1][0]+pointList[1][1]*pointList[1][1];
    M(1,3) = (double) 1.0;
    M(2,0) = (double) pointList[2][0];
    M(2,1) = (double) pointList[2][1];
    M(2,2) = (double) pointList[2][0]*pointList[2][0]+pointList[2][1]*pointList[2][1];
    M(2,3) = (double) 1.0;
    M(3,0) = (double) pointList[3][0];
    M(3,1) = (double) pointList[3][1];
    M(3,2) = (double) pointList[3][0]*pointList[3][0]+pointList[3][1]*pointList[3][1];
    M(3,3) = (double) 1.0;
    
    //std::cout << "M: " << std::endl;
    //std::cout << M << std::endl;
    // determinant computation - the result is multiplied by 'orientation'
    det = vnl_det( M.GetVnlMatrix() ) * orientation;
    //std::cout << "Det(M): " << det << std::endl;
      
    // NOTE STEPH: eratum - 0 means OUT, 1 means IN ?
    return ( det < 0 ? 1 : 0 );
    }

protected:

  InCircleTestFunction() {}
  ~InCircleTestFunction() {}

private:

  InCircleTestFunction( const Self & ); // purposely not implemented
  void operator=( const Self & );       // purposely not implemented

};

} // namespace ITK
#endif //__itkInCircleTestFunction_h
