
/*=========================================================================
 *
 *  Portions of this file are subject to the VTK Toolkit Version 3 copyright.
 *
 *  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 *
 *  For complete copyright, license and disclaimer of warranty information
 *  please refer to the NOTICE file at the top of the ITK source tree.
 *
 *=========================================================================*/

#ifndef __itkPointSetToPointSetFilter_hxx
#define __itkPointSetToPointSetFilter_hxx

#include "itkPointSetToPointSetFilter.h"

namespace itk
{
/**
 *
 */
template< class TInputPointSet, class TOutputPointSet >
PointSetToPointSetFilter< TInputPointSet, TOutputPointSet >
::PointSetToPointSetFilter()
{
  // Modify superclass default values, can be overridden by subclasses
  this->SetNumberOfRequiredInputs( 1 );
}

/**
 *
 */
template< class TInputPointSet, class TOutputPointSet >
void
PointSetToPointSetFilter< TInputPointSet, TOutputPointSet >
::SetInput( const TInputPointSet *input )
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( 0,
                                    const_cast< TInputPointSet * >( input ) );
}

/**
 *
 */
template< class TInputPointSet, class TOutputPointSet >
const typename PointSetToPointSetFilter< TInputPointSet, TOutputPointSet >::InputPointSetType *
PointSetToPointSetFilter< TInputPointSet, TOutputPointSet >
::GetInput() const
{
  return static_cast< const TInputPointSet * >( this->GetPrimaryInput() );
}

/**
 *
 */
template< class TInputPointSet, class TOutputPointSet >
const typename PointSetToPointSetFilter< TInputPointSet, TOutputPointSet >::InputPointSetType *
PointSetToPointSetFilter< TInputPointSet, TOutputPointSet >
::GetInput( unsigned int idx ) const
{
  return dynamic_cast< const TInputPointSet * >
         ( this->ProcessObject::GetInput( idx ) );
}

template< class TInputPointSet, class TOutputPointSet >
void
PointSetToPointSetFilter< TInputPointSet, TOutputPointSet >
::CopyInputPointSetToOutputPointSetPoints( void )
{
  const InputPointSetType *inputPointSet   =  this->GetInput();
  OutputPointSetPointer    outputPointSet   =  this->GetOutput();

  typedef typename TOutputPointSet::PointsContainer OutputPointsContainer;
  typedef typename TInputPointSet::PointsContainer  InputPointsContainer;

  typename OutputPointsContainer::Pointer outputPoints = OutputPointsContainer::New();
  const InputPointsContainer *inputPoints = inputPointSet->GetPoints();

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

    outputPointSet->SetPoints(outputPoints);
    }
}

template< class TInputPointSet, class TOutputPointSet >
void
PointSetToPointSetFilter< TInputPointSet, TOutputPointSet >
::CopyInputPointSetToOutputPointSetPointData( void )
{
  const InputPointSetType *inputPointSet   =  this->GetInput();
  OutputPointSetPointer    outputPointSet   =  this->GetOutput();

  typedef typename TOutputPointSet::PointDataContainer OutputPointDataContainer;
  typedef typename TInputPointSet::PointDataContainer  InputPointDataContainer;

  typename OutputPointDataContainer::Pointer outputPointData = OutputPointDataContainer::New();
  const InputPointDataContainer *inputPointData = inputPointSet->GetPointData();

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

    outputPointSet->SetPointData(outputPointData);
    }
}

} // end namespace itk

#endif
