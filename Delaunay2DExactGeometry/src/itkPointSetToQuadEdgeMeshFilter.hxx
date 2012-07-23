
#ifndef __itkPointSetToQuadEdgeMeshFilter_hxx
#define __itkPointSetToQuadEdgeMeshFilter_hxx

// ITK includes
#include "itkPointSetToQuadEdgeMeshFilter.h"

namespace itk {

/**
*
*/
template<class TPointSet, class TQuadEdgeMesh>
PointSetToQuadEdgeMeshFilter< TPointSet, TQuadEdgeMesh >::
PointSetToQuadEdgeMeshFilter()
{
  // Modify superclass default values, can be overridden by subclasses
  this->SetNumberOfRequiredInputs(1);
}

/**
*
*/
template<class TPointSet, class TQuadEdgeMesh>
void
PointSetToQuadEdgeMeshFilter< TPointSet, TQuadEdgeMesh >::
SetInput(const TPointSet *input)
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( 0,
                                 const_cast< TPointSet * >( input ) );
}

/**
*
*/
template<class TPointSet, class TQuadEdgeMesh>
const typename PointSetToQuadEdgeMeshFilter< TPointSet, TQuadEdgeMesh >::PointSetType *
PointSetToQuadEdgeMeshFilter< TPointSet, TQuadEdgeMesh >::
GetInput() const
{
  return static_cast< const TPointSet * >( this->GetPrimaryInput() );
}

/**
*
*/
template<class TPointSet, class TQuadEdgeMesh>
  const typename PointSetToQuadEdgeMeshFilter< TPointSet, TQuadEdgeMesh >::PointSetType *
PointSetToQuadEdgeMeshFilter< TPointSet, TQuadEdgeMesh >::
GetInput(unsigned int idx) const
{
  return dynamic_cast< const TPointSet * > ( this->ProcessObject::GetInput(idx) );
}

template<class TPointSet, class TQuadEdgeMesh>
void
PointSetToQuadEdgeMeshFilter< TPointSet, TQuadEdgeMesh >::
CopyInputMeshToOutputMeshPoints(void)
{
  const PointSetType  *pointSet    =  this->GetInput();
  QuadEdgeMeshPointer  outputMesh  =  this->GetOutput();

  typedef typename TQuadEdgeMesh::PointsContainer OutputPointsContainer;
  typedef typename TPointSet::PointsContainer     InputPointsContainer;

  typename OutputPointsContainer::Pointer outputPoints = OutputPointsContainer::New();
  const InputPointsContainer *inputPoints = pointSet->GetPoints();

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

template<class TPointSet, class TQuadEdgeMesh>
void
PointSetToQuadEdgeMeshFilter< TPointSet, TQuadEdgeMesh >::
CopyInputMeshToOutputMeshPointData(void)
{
  const PointSetType  *pointSet    =  this->GetInput();
  QuadEdgeMeshPointer  outputMesh  =  this->GetOutput();

  typedef typename TQuadEdgeMesh::PointDataContainer OutputPointDataContainer;
  typedef typename TPointSet::PointDataContainer     InputPointDataContainer;

  typename OutputPointDataContainer::Pointer outputPointData = OutputPointDataContainer::New();
  const InputPointDataContainer *inputPointData = pointSet->GetPointData();

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

}

#endif // __itkPointSetToQuadEdgeMeshFilter_hxx