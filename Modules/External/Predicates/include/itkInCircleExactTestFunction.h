#ifndef __itkInCircleExactTestFunction_h
#define __itkInCircleExactTestFunction_h

#include <itkFunctionBase.h>
#include "itkOrientation2DFunction.h"

//--------------------------------------------------------------------
// Skewchuk code
//
//--------------------------------------------------------------------

namespace itk
{

extern "C"
  {
  double incircle( double* pa, double* pb, double* pc, double* pd );
  }

//-------------------------------------------------------------------
// In Circle test wrapper for ITK
//
template
  <
  class MeshType,
  class OutputType =MeshType
  >
class ITK_EXPORT InCircleExactTestFunction :
      public FunctionBase< std::vector< typename MeshType::PointType>, bool >
{

public:
  
  /** Standart class typedefs **/
  typedef InCircleExactTestFunction                                             Self;
  typedef FunctionBase< std::vector< typename MeshType::PointType>, bool > Superclass;
  typedef SmartPointer< Self >                                             Pointer;
  typedef SmartPointer< const Self >                                       ConstPointer;
    
  /** Run-time type information (and related methods) **/
  itkNewMacro( Self );
  itkTypeMacro( InCircleExactTestFunction, FunctionBase );

  typedef MeshType Mesh;
  typedef typename Mesh::PointType PointType;
  typedef typename MeshType::PixelType PixelType;
  typedef Orientation2DFunction< Mesh, PixelType > OrientationTest;
 
  typedef std::vector< PointType > PointList;
  
  bool Evaluate( const PointList & pointList ) const
    {
    // NOTE STEPH: number of points should be tested
    double *  pa = new double[2];
    double *  pb = new double[2];
    double *  pc = new double[2];
    double *  pd = new double[2];
    
    pa[0] = pointList[0][0];
    pa[1] = pointList[0][1];
    pb[0] = pointList[1][0];
    pb[1] = pointList[1][1];
    pc[0] = pointList[2][0];
    pc[1] = pointList[2][1];
    pd[0] = pointList[3][0];
    pd[1] = pointList[3][1];
    
    // orientation test
    typename OrientationTest::Pointer orientationTest = OrientationTest::New();
    double orientation = orientationTest->Evaluate( pointList );
    
    // incircle test - the result is multipled by the orientation test result
    double det = incircle( pa, pb, pc, pd ) * orientation;

    // zero, which means the point is ON the circle is considered IN
    // NOTE STEPH: eratum - 0 means OUT, 1 means IN ?
    return ( det < 0 ? 1 : 0 );
    }
  
protected:

  InCircleExactTestFunction() {}
  ~InCircleExactTestFunction() {}

private:

  InCircleExactTestFunction( const Self & ); // purposely not implemented
  void operator=( const Self & );       // purposely not implemented

};

} // namespace ITK
#endif //__itkInCircleExactTestFunction_h
