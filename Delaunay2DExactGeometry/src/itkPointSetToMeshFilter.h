
#ifndef __itkPointSetToMeshFilter_h
#define __itkPointSetToMeshFilter_h

// ITK includes

#include "itkMeshSource.h"

namespace itk {
  
template< class TPointSet, class TMesh >
class ITK_EXPORT PointSetToMeshFilter : public MeshSource< TMesh >
{
public:
  
  /** Standard class typedefs. */
  typedef PointSetToMeshFilter       Self;
  typedef MeshSource<TMesh>          Superclass;
  typedef SmartPointer<Self>         Pointer;
  typedef SmartPointer<const Self>   ConstPointer;
  
  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );
  
  /** Run-time type information (and related methods).    */
  itkTypeMacro( PointSetToMeshFilter, MeshSource );
  
  /** Some convenient typedefs. */
  typedef TPointSet                       PointSetType;
  typedef typename PointSetType::Pointer  PointSetPointer;
  typedef TMesh                           MeshType;
  typedef typename MeshType::Pointer      MeshPointer;
  
  /** Set the pointset input of this process object.  */
  using Superclass::SetInput;
  void SetInput(const PointSetType *input);
  
  /** Get the pointset input of this process object.  */
  const PointSetType * GetInput(void) const;
  
  const PointSetType * GetInput(unsigned int idx) const;
  
protected:
  
  /** Generic Methods */
  PointSetToMeshFilter();
  ~PointSetToMeshFilter() {}
    
  void CopyInputMeshToOutputMeshPoints();
  
  void CopyInputMeshToOutputMeshPointData();
  
private:
  
  PointSetToMeshFilter( const Self & ); // purposely not implemented
  void operator=( const Self & );       // purposely not implemented
  
}; // PointSetToMeshFilter

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPointSetToMeshFilter.hxx"
#endif

#endif // __itkPointSetToMeshFilter_h
