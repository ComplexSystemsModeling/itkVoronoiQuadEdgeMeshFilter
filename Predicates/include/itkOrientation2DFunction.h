#ifndef __itkOrientation2DFunction_h
#define __itkOrientation2DFunction_h

#include <itkFunctionBase.h>


namespace itk
{

//-----------------------------------------------------------------
// Skewchuk code
//-----------------------------------------------------------------
extern "C"
  {
  double orient2d( double* pa, double* pb, double* pc );
  }

//----------------------------------------------------------------
// Orientation test wrapper for ITK
//-----------------------------------------------------------------
template
  < 
  class MeshType,
  class OutputType = double
  >
class ITK_EXPORT Orientation2DFunction :
  public FunctionBase< std::vector< typename MeshType::PointType>, double >
{

public:
 
  /** Standart class typedefs **/
  typedef Orientation2DFunction                                                   Self;
  typedef FunctionBase< std::vector< typename MeshType::PointType>, double > Superclass;
  typedef SmartPointer< Self >                                                    Pointer;
  typedef SmartPointer< const Self >                                              ConstPointer;

  /** Run-time type information (and related methods) **/
  itkNewMacro( Self );
  itkTypeMacro( Orientation2DFunction, FunctionBase );
  
  typedef MeshType Mesh;
  typedef typename Mesh::PointType PointType;
  typedef std::vector< PointType > PointList;
 
  double Evaluate( const PointList & pointList ) const
    {
    // NOTE STEPH: number of points should be tested
    double * pa = new double[2];
    double * pb = new double[2];
    double * pc = new double[2];

    pa[0] = (double) pointList[0][0];
    pa[1] = (double) pointList[0][1];
    pb[0] = (double) pointList[1][0];
    pb[1] = (double) pointList[1][1];
    pc[0] = (double) pointList[2][0];
    pc[1] = (double) pointList[2][1];

    double orientation = orient2d( pa, pb, pc ) ;
    //std::cout << "orientation 1 : " << orientation << std::endl;
    return orientation;
    }

protected:

  Orientation2DFunction() {}
  ~Orientation2DFunction() {} 

private:

  Orientation2DFunction( const Self & ); // purposely not implemented
  void operator=( const Self & );        // purposely not implemented

};

} // namespace ITK
#endif // __itkOrientation2DFunction_h
