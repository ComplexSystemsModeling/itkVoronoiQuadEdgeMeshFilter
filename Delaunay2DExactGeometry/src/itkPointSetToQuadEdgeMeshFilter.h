
#ifndef __itkPointSetToQuadEdgeMeshFilter_h
#define __itkPointSetToQuadEdgeMeshFilter_h

// ITK includes

#include "itkQuadEdgeMeshSource.h"

namespace itk {
  
template<class TPointSet, class TQuadEdgeMesh>
class ITK_EXPORT PointSetToQuadEdgeMeshFilter : 
  public QuadEdgeMeshSource<TQuadEdgeMesh>
{
public:
  
  /** Standard class typedefs. */
  typedef PointSetToQuadEdgeMeshFilter      Self;
  typedef QuadEdgeMeshSource<TQuadEdgeMesh> Superclass;
  typedef SmartPointer<Self>                Pointer;
  typedef SmartPointer<const Self>          ConstPointer;
  
  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );
  
  /** Run-time type information (and related methods).    */
  itkTypeMacro( PointSetToQuadEdgeMeshFilter, QuadEdgeMeshSource );
  
  /** Some convenient typedefs. */
  typedef TPointSet                          PointSetType;
  typedef typename PointSetType::Pointer     PointSetPointer;
  typedef TQuadEdgeMesh                      QuadEdgeMeshType;
  typedef typename QuadEdgeMeshType::Pointer QuadEdgeMeshPointer;
  
  /** Set the pointset input of this process object.  */
  using Superclass::SetInput;
  void SetInput(const PointSetType *input);
  
  /** Get the pointset input of this process object.  */
  const PointSetType * GetInput(void) const;
  
  const PointSetType * GetInput(unsigned int idx) const;
  
protected:
  
  /** Generic Methods */
  PointSetToQuadEdgeMeshFilter();
  
  virtual ~PointSetToQuadEdgeMeshFilter() {}
  
  void CopyInputMeshToOutputMeshPoints();
  
  void CopyInputMeshToOutputMeshPointData();
  
private:
  
  PointSetToQuadEdgeMeshFilter( const Self & ); // purposely not implemented
  void operator=( const Self & );               // purposely not implemented
  
}; // PointSetToQuadEdgeMeshFilter

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPointSetToQuadEdgeMeshFilter.hxx"
#endif

#endif // __itkPointSetToQuadEdgeMeshFilter_h