/*==================================================================
 *                                                                 *
 * Circumcentre Dual Point Functor                                 *
 *                                                                 *
 *   Implementation for ITK by Stéphane U. Rigaud, Humayun Irshad  *
 *   and Alexandre Gouaillard                                      *
 *                                                                 *
 *=================================================================*/

#ifndef __CircumcentreDualPointFunctor_h
#define __CircumcentreDualPointFunctor_h

#include "itkQuadEdgeMesh.h"
#include "itkPointInCircleGeometricalPredicateFunctor.h"

namespace itk 
{
namespace Functor 
{
  template< class TInputMesh, class TOutputMesh=TInputMesh >
  class CircumcentreDualPointFunctor
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

    CircumcentreDualPointFunctor() {};
    virtual ~CircumcentreDualPointFunctor() {};

    inline OutputPointType operator() ( const TInputMesh* primalMesh, const CellIterator cellIterator ) const
    {
      OutputPointType d_point;
      PointIdConstIterator current= cellIterator.Value()->PointIdsBegin();
      PointType pointTab[3];

      unsigned int i = 0;
      while( current != cellIterator.Value()->PointIdsEnd() )
        {
          if( i < 3 )
            {
            pointTab[i] = primalMesh->GetPoint( *current );
            }
          current++;
          i++;
        }
      if( i > 3 )
        {
        // NOTE STEPH: is it possible to use the itkExceptionMacro ?
        //itkExceptionMacro("Number of vertex = " << i << " , the mesh is not a triangulation.");
        std::cerr << "itkCircumcentreDualPointFunctor" << std::endl;
        std::cerr << "Exception: Number of vertex = " << i << " , the mesh is not a triangulation." << std::endl;
        throw i;
        } 
      
      double xba, yba, zba, xca, yca, zca;
      double calength, balength;
      PointType crossbc;
      double denominator;
      double xcirca, ycirca, zcirca;
      PointType ta, tb, tc;
      PointType ba, ca;
      PointType circa;
      PointType A, B, C;
      
      A = pointTab[0];
      B = pointTab[1];
      C = pointTab[2];
      
      /* Use coordinates relative to point `a' of the triangle. */
      for( unsigned int i = 0; i < PointDimension; i++ )
        {
        ba[i] = B[i] - A[i];  
        ca[i] = C[i] - A[i];
        } 
      
      /* Squares of lengths of the edges incident to `a'. */
      balength = ba[0] * ba[0] + ba[1] * ba[1] + ba[2] * ba[2];
      calength = ca[0] * ca[0] + ca[1] * ca[1] + ca[2] * ca[2];
      
      ta[0] = B[1]; ta[1] = B[2];
      tb[0] = C[1]; tb[1] = C[2];
      tc[0] = A[1]; tc[1] = A[2];
      crossbc[0] = OrientationTest(ta, tb, tc);
      
      ta[0] = C[0]; ta[1] = C[2];
      tb[0] = B[0]; tb[1] = B[2];
      tc[0] = A[0]; tc[1] = A[2];
      crossbc[1] = OrientationTest(ta, tb, tc);
      
      ta[0] = B[0]; ta[1] = B[1];
      tb[0] = C[0]; tb[1] = C[1];
      tc[0] = A[0]; tc[1] = A[1];
      crossbc[2] = OrientationTest(ta, tb, tc);
      
      denominator = 0.5 / (crossbc[0] * crossbc[0] + crossbc[1] * crossbc[1] + crossbc[2] * crossbc[2]);
      
      circa[0] = ((balength * ca[1] - calength * ba[1]) * crossbc[2] - (balength * ca[2] - calength * ba[2]) * crossbc[1]) * denominator;
      circa[1] = ((balength * ca[2] - calength * ba[2]) * crossbc[0] - (balength * ca[0] - calength * ba[0]) * crossbc[2]) * denominator;
      circa[2] = ((balength * ca[0] - calength * ba[0]) * crossbc[1] - (balength * ca[1] - calength * ba[1]) * crossbc[0]) * denominator;

      for( unsigned int i = 0; i < PointDimension; i++ )
        {
        d_point[i] = static_cast<OutputPixelType>( circa[i] + A[i] );
        }
      
      return d_point;
    }
  };
  
} // namespace Functor
} // namespace itk

#endif // __CircumcentreDualPointFunctor_h
