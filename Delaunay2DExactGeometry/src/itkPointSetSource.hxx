
#ifndef __itkPointSetSource_hxx
#define __itkPointSetSource_hxx

#include "itkPointSetSource.h"

namespace itk
{
/**
 *
 */
template< class TOutputPointSet >
PointSetSource< TOutputPointSet >
::PointSetSource()
{
  // Create the output. We use static_cast<> here because we know the default
  // output must be of type TOutputPointSet
  OutputPointSetPointer output =
    static_cast< TOutputPointSet * >( this->MakeOutput(0).GetPointer() );

  this->ProcessObject::SetNumberOfRequiredOutputs(1);
  this->ProcessObject::SetNthOutput( 0, output.GetPointer() );

  m_GenerateDataRegion = 0;
  m_GenerateDataNumberOfRegions = 0;
}

/**
 *
 */
template< class TOutputPointSet >
typename PointSetSource< TOutputPointSet >::DataObjectPointer
PointSetSource< TOutputPointSet >
::MakeOutput( DataObjectPointerArraySizeType )
{
  return static_cast< DataObject * >( TOutputPointSet::New().GetPointer() );
}

/**
 *
 */
template< class TOutputPointSet >
typename PointSetSource< TOutputPointSet >::OutputPointSetType *
PointSetSource< TOutputPointSet >
::GetOutput( void )
{
  return static_cast< TOutputPointSet * >( this->GetPrimaryOutput() );
}

/**
 *
 */
template< class TOutputPointSet >
typename PointSetSource< TOutputPointSet >::OutputPointSetType *
PointSetSource< TOutputPointSet >
::GetOutput( unsigned int idx )
{
  return static_cast< TOutputPointSet * >
         ( this->ProcessObject::GetOutput( idx ) );
}

/**
 *
 */
template< class TOutputPointSet >
void
PointSetSource< TOutputPointSet >
::SetOutput( OutputPointSetType *output )
{
  itkWarningMacro(
    <<
    "SetOutput(): This method is slated to be removed from ITK.  Please use GraftOutput() in possible combination with DisconnectPipeline() instead.");
  this->ProcessObject::SetNthOutput( 0, output );
}

/**
 *
 */
template< class TOutputPointSet >
void
PointSetSource< TOutputPointSet >
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();
}

/**
 *
 */
template< class TOutputPointSet >
void
PointSetSource< TOutputPointSet >
::GraftOutput( DataObject *graft )
{
  this->GraftNthOutput( 0, graft );
}

/**
 *
 */
template< class TOutputPointSet >
void
PointSetSource< TOutputPointSet >
::GraftOutput( const DataObjectIdentifierType & key, DataObject *graft )
{
  if ( !graft )
    {
    itkExceptionMacro(<< "Requested to graft output that is a NULL pointer");
    }

  // we use the process object method since all out output may not be
  // of the same type
  DataObject *output = this->ProcessObject::GetOutput( key );

  // Call GraftImage to copy meta-information, regions, and the pixel container
  output->Graft( graft );
}

/**
 *
 */
template< class TOutputPointSet >
void
PointSetSource< TOutputPointSet >
::GraftNthOutput( unsigned int idx, DataObject *graft )
{
  if ( idx >= this->GetNumberOfIndexedOutputs() )
    {
    itkExceptionMacro(<< "Requested to graft output " << idx
                      << " but this filter only has " << this->GetNumberOfIndexedOutputs() << " indexed Outputs.");
    }
  this->GraftOutput( this->MakeNameFromIndex(idx), graft );
}

/**
 *
 */
template< class TOutputPointSet >
void
PointSetSource< TOutputPointSet >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
}
} // end namespace itk

#endif
