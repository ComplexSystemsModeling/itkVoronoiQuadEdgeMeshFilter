
#ifndef __itkQuadEdgeMeshSource_hxx
#define __itkQuadEdgeMeshSource_hxx

#include "itkQuadEdgeMeshSource.h"

namespace itk
{
/**
 *
 */
template< class TOutputQuadEdgeMesh >
QuadEdgeMeshSource< TOutputQuadEdgeMesh >
::QuadEdgeMeshSource()
{
  // Create the output. We use static_cast<> here because we know the default
  // output must be of type TOutputQuadEdgeMesh
  OutputQuadEdgeMeshPointer output =
    static_cast< TOutputQuadEdgeMesh * >( this->MakeOutput(0).GetPointer() );

  this->ProcessObject::SetNumberOfRequiredOutputs(1);
  this->ProcessObject::SetNthOutput( 0, output.GetPointer() );

  m_GenerateDataRegion = 0;
  m_GenerateDataNumberOfRegions = 0;
}

/**
 *
 */
template< class TOutputQuadEdgeMesh >
typename QuadEdgeMeshSource< TOutputQuadEdgeMesh >::DataObjectPointer
QuadEdgeMeshSource< TOutputQuadEdgeMesh >
::MakeOutput( DataObjectPointerArraySizeType )
{
  return static_cast< DataObject * >( TOutputQuadEdgeMesh::New().GetPointer() );
}

/**
 *
 */
template< class TOutputQuadEdgeMesh >
typename QuadEdgeMeshSource< TOutputQuadEdgeMesh >::OutputQuadEdgeMeshType *
QuadEdgeMeshSource< TOutputQuadEdgeMesh >
::GetOutput( void )
{
  return static_cast< TOutputQuadEdgeMesh * >( this->GetPrimaryOutput() );
}

/**
 *
 */
template< class TOutputQuadEdgeMesh >
typename QuadEdgeMeshSource< TOutputQuadEdgeMesh >::OutputQuadEdgeMeshType *
QuadEdgeMeshSource< TOutputQuadEdgeMesh >
::GetOutput( unsigned int idx )
{
  return static_cast< TOutputQuadEdgeMesh * >
         ( this->ProcessObject::GetOutput( idx ) );
}

/**
 *
 */
template< class TOutputQuadEdgeMesh >
void
QuadEdgeMeshSource< TOutputQuadEdgeMesh >
::SetOutput( OutputQuadEdgeMeshType *output )
{
  itkWarningMacro(
    <<
    "SetOutput(): This method is slated to be removed from ITK.  Please use GraftOutput() in possible combination with DisconnectPipeline() instead.");
  this->ProcessObject::SetNthOutput( 0, output );
}

/**
 *
 */
template< class TOutputQuadEdgeMesh >
void
QuadEdgeMeshSource< TOutputQuadEdgeMesh >
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();
}

/**
 *
 */
template< class TOutputQuadEdgeMesh >
void
QuadEdgeMeshSource< TOutputQuadEdgeMesh >
::GraftOutput( DataObject *graft )
{
  this->GraftNthOutput( 0, graft );
}

/**
 *
 */
template< class TOutputQuadEdgeMesh >
void
QuadEdgeMeshSource< TOutputQuadEdgeMesh >
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
template< class TOutputQuadEdgeMesh >
void
QuadEdgeMeshSource< TOutputQuadEdgeMesh >
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
template< class TOutputQuadEdgeMesh >
void
QuadEdgeMeshSource< TOutputQuadEdgeMesh >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
}
} // end namespace itk

#endif
