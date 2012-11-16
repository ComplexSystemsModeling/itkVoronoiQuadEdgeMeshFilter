/*==================================================================
 *                                                                 *
 * Barycentre Dual Point Functor                                   *
 *                                                                 *
 *   Implementation for ITK by Stéphane U. Rigaud, Humayun Irshad  *
 *   and Alexandre Gouaillard                                      *
 *                                                                 *
 *=================================================================*/

#ifndef __BarycentreDualPointFunctor_h
#define __BarycentreDualPointFunctor_h

#include "itkQuadEdgeMesh.h"

namespace itk 
{
namespace Functor 
{
  template< class TInputMesh, class TOutputMesh=TInputMesh >
  class BarycentreDualPointFunctor
  {

  public:

    typedef typename TInputMesh::PixelType PixelType;
    typedef typename TInputMesh::PointType PointType;
    typedef typename TInputMesh::PointIdentifier PointIdentifier;
    typedef typename TInputMesh::CellType CellType;
    typedef typename TInputMesh::CellIdentifier CellIdentifier;
    typedef typename TInputMesh::CellsContainer CellsContainer;  
    typedef typename CellsContainer::ConstIterator CellIterator;
    typedef typename CellType::PointIdConstIterator  PointIdConstIterator;
    
    typedef typename TOutputMesh::PixelType OutputPixelType;
    typedef typename TOutputMesh::PointType OutputPointType;    

    itkStaticConstMacro( PointDimension, unsigned int, TInputMesh::Traits::PointDimension );

    BarycentreDualPointFunctor() {};
    virtual ~BarycentreDualPointFunctor() {};

    inline OutputPointType operator() ( const TInputMesh* primalMesh, const CellIterator cellIterator ) const
    {
      OutputPointType d_point;
      PointIdConstIterator current = cellIterator.Value()->PointIdsBegin();
      PointIdConstIterator end     = cellIterator.Value()->PointIdsEnd();
      for( unsigned int i = 0; i < PointDimension; i++ )
        {
        d_point[i] = 0.0;
        }
      while( current != end )
        {
        PointType point = primalMesh->GetPoint( *current );
        for( unsigned int i = 0; i < PointDimension; i++ )
          {
          d_point[i] += static_cast<OutputPixelType>( point[i] );
          }
        current++;
        }
      for( unsigned int i =0; i < PointDimension; i++ )
        {
        d_point[i] /= PointDimension;
        }
      return d_point;
    }
  };
  
} // namespace Functor
} // namespace itk

#endif // __BarycentreDualPointFunctor_h
