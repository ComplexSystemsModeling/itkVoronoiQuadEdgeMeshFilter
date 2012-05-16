
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

#ifndef __itkPointSetToPointSetFilter_h
#define __itkPointSetToPointSetFilter_h

#include "itkPointSetSource.h"

namespace itk
{
/** \class PointSetToPointSetFilter
 * \brief
 *
 * PointSetToPointSetFilter is the base class for all process objects that output
 * PointSet data, and require PointSet data as input. Specifically, this class
 * defines the SetInput() method for defining the input to a filter.
 *
 * \ingroup PointSetFilters
 *
 * \ingroup ITKPointSet
 */
template< class TInputPointSet, class TOutputPointSet >
class ITK_EXPORT PointSetToPointSetFilter : public PointSetSource< TOutputPointSet >
{
  
public:
  
  /** Standard class typedefs. */
  typedef PointSetToPointSetFilter           Self;
  typedef PointSetSource< TOutputPointSet >  Superclass;
  typedef SmartPointer< Self >               Pointer;
  typedef SmartPointer< const Self >         ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( PointSetToPointSetFilter, PointSetSource );

  /** Some convenient typedefs. */
  typedef TInputPointSet                       InputPointSetType;
  typedef typename InputPointSetType::Pointer  InputPointSetPointer;
  typedef TOutputPointSet                      OutputPointSetType;
  typedef typename OutputPointSetType::Pointer OutputPointSetPointer;

  /** Set the PointSet input of this process object.  */
  using Superclass::SetInput;
  void SetInput( const InputPointSetType *input );

  /** Get the PointSet input of this process object.  */
  const InputPointSetType * GetInput( void ) const;

  const InputPointSetType * GetInput( unsigned int idx ) const;

protected:
  
  PointSetToPointSetFilter();
  ~PointSetToPointSetFilter() {}

  void CopyInputPointSetToOutputPointSetPoints();

  void CopyInputPointSetToOutputPointSetPointData();

private:
  
  PointSetToPointSetFilter( const Self & ); //purposely not implemented
  void operator=( const Self & );   //purposely not implemented
  
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPointSetToPointSetFilter.hxx"
#endif

#endif
