
#ifndef __itkQuadEdgeQuadEdgeMeshSource_hxx
#define __itkQuadEdgeQuadEdgeMeshSource_hxx

// ITK includes
#include "itkQuadEdgeMeshSource.h"

namespace itk {

/**
 *
 */
template< class TOutputMesh >
QuadEdgeMeshSource< TOutputMesh >
::QuadEdgeMeshSource()
{
  // Create the output. We use static_cast<> here because we know the default
  // output must be of type TOutputMesh
  OutputMeshPointer output =
  static_cast< TOutputMesh * >( this->MakeOutput(0).GetPointer() );
  
  this->ProcessObject::SetNumberOfRequiredOutputs(1);
  this->ProcessObject::SetNthOutput( 0, output.GetPointer() );
  
  m_GenerateDataRegion = 0;
  m_GenerateDataNumberOfRegions = 0;
}

/**
 *
 */
template< class TOutputMesh >
typename QuadEdgeMeshSource< TOutputMesh >::DataObjectPointer
QuadEdgeMeshSource< TOutputMesh >
::MakeOutput(DataObjectPointerArraySizeType)
{
  return static_cast< DataObject * >( TOutputMesh::New().GetPointer() );
}

/**
 *
 */
template< class TOutputMesh >
typename QuadEdgeMeshSource< TOutputMesh >::OutputMeshType *
QuadEdgeMeshSource< TOutputMesh >
::GetOutput(void)
{
  return static_cast< TOutputMesh * >( this->GetPrimaryOutput() );
}

/**
 *
 */
template< class TOutputMesh >
typename QuadEdgeMeshSource< TOutputMesh >::OutputMeshType *
QuadEdgeMeshSource< TOutputMesh >
::GetOutput(unsigned int idx)
{
  return static_cast< TOutputMesh * >
  ( this->ProcessObject::GetOutput(idx) );
}

/**
 *
 */
template< class TOutputMesh >
void
QuadEdgeMeshSource< TOutputMesh >
::SetOutput(OutputMeshType *output)
{
  itkWarningMacro(
                  <<
                  "SetOutput(): This method is slated to be removed from ITK.  Please use GraftOutput() in possible combination with DisconnectPipeline() instead.");
  this->ProcessObject::SetNthOutput(0, output);
}

/**
 *
 */
template< class TOutputMesh >
void
QuadEdgeMeshSource< TOutputMesh >
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();
}

/**
 *
 */
template< class TOutputMesh >
void
QuadEdgeMeshSource< TOutputMesh >
::GraftOutput(DataObject *graft)
{
  this->GraftNthOutput(0, graft);
}

/**
 *
 */
template< class TOutputMesh >
void
QuadEdgeMeshSource< TOutputMesh >
::GraftOutput(const DataObjectIdentifierType & key, DataObject *graft)
{
  if ( !graft )
  {
    itkExceptionMacro(<< "Requested to graft output that is a NULL pointer");
  }
  
  // we use the process object method since all out output may not be
  // of the same type
  DataObject *output = this->ProcessObject::GetOutput(key);
  
  // Call GraftImage to copy meta-information, regions, and the pixel container
  output->Graft(graft);
}

/**
 *
 */
template< class TOutputMesh >
void
QuadEdgeMeshSource< TOutputMesh >
::GraftNthOutput(unsigned int idx, DataObject *graft)
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
template< class TOutputMesh >
void
QuadEdgeMeshSource< TOutputMesh >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

} // namespace itk

#endif // __itkQuadEdgeQuadEdgeMeshSource_hxx