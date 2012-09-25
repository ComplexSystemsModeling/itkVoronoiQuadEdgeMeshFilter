
#ifndef __itkPointSetToMeshFilter_hxx
#define __itkPointSetToMeshFilter_hxx

// ITK includes
#include "itkPointSetToMeshFilter.h"

namespace itk {

/**
*
*/
template<class TPointSet, class TMesh>
PointSetToMeshFilter< TPointSet, TMesh >::
PointSetToMeshFilter()
{
  // Modify superclass default values, can be overridden by subclasses
  this->SetNumberOfRequiredInputs(1);
}

/**
*
*/
template<class TPointSet, class TMesh>
void
PointSetToMeshFilter< TPointSet, TMesh >::
SetInput(const TPointSet *input)
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( 0,
                                 const_cast< TPointSet * >( input ) );
}

/**
*
*/
template<class TPointSet, class TMesh>
const typename PointSetToMeshFilter< TPointSet, TMesh >::PointSetType *
PointSetToMeshFilter< TPointSet, TMesh >::
GetInput() const
{
  return static_cast< const TPointSet * >( this->GetPrimaryInput() );
}

/**
*
*/
template<class TPointSet, class TMesh>
const typename PointSetToMeshFilter< TPointSet, TMesh >::PointSetType *
PointSetToMeshFilter< TPointSet, TMesh >::
GetInput(unsigned int idx) const
{
  return dynamic_cast< const TPointSet * > ( this->ProcessObject::GetInput(idx) );
}

/**
 *
 */
template<class TPointSet, class TMesh>
void
PointSetToMeshFilter< TPointSet, TMesh >::
CopyInputMeshToOutputMeshPoints(void)
{
  const PointSetType  *pointSet    =  this->GetInput();
  MeshPointer          outputMesh  =  this->GetOutput();

  typedef typename TMesh::PointsContainer      OutputPointsContainer;
  typedef typename TPointSet::PointsContainer  InputPointsContainer;

  typename OutputPointsContainer::Pointer outputPoints = OutputPointsContainer::New();
  const    InputPointsContainer           *inputPoints = pointSet->GetPoints();

  if ( inputPoints )
    {
    outputPoints->Reserve( inputPoints->Size() );
    
    typename InputPointsContainer::ConstIterator inputItr = inputPoints->Begin();
    typename InputPointsContainer::ConstIterator inputEnd = inputPoints->End();
    
    typename OutputPointsContainer::Iterator outputItr = outputPoints->Begin();
    
    while ( inputItr != inputEnd )
      {
      outputItr.Value() = inputItr.Value();
      ++inputItr;
      ++outputItr;
      }
    
    outputMesh->SetPoints(outputPoints);
    }
}

/**
 *
 */
template<class TPointSet, class TMesh>
void
PointSetToMeshFilter< TPointSet, TMesh >::
CopyInputMeshToOutputMeshPointData(void)
{
  const PointSetType  *pointSet    =  this->GetInput();
  MeshPointer          outputMesh  =  this->GetOutput();

  typedef typename TMesh::PointDataContainer      OutputPointDataContainer;
  typedef typename TPointSet::PointDataContainer  InputPointDataContainer;

  typename OutputPointDataContainer::Pointer outputPointData = OutputPointDataContainer::New();
  const    InputPointDataContainer           *inputPointData = pointSet->GetPointData();

  if ( inputPointData )
    { 
    outputPointData->Reserve( inputPointData->Size() );
    
    typename InputPointDataContainer::ConstIterator inputItr = inputPointData->Begin();
    typename InputPointDataContainer::ConstIterator inputEnd = inputPointData->End();
    
    typename OutputPointDataContainer::Iterator outputItr = outputPointData->Begin();
    
    while ( inputItr != inputEnd )
      {
      outputItr.Value() = inputItr.Value();
      ++inputItr;
      ++outputItr;
      }
    
    outputMesh->SetPointData(outputPointData);
  }
}

} // namespace itk

#endif // __itkPointSetToMeshFilter_hxx
