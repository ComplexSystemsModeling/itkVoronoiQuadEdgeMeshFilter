#ifndef __CircumcentreDualPointFunctor_h
#define __CircumcentreDualPointFunctor_h

#include "itkQuadEdgeMesh.h"


namespace itk {
namespace Functor {

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
    ~CircumcentreDualPointFunctor() {};

    inline OutputPointType operator() ( const TInputMesh* primalMesh, CellIterator cellIterator )
    {
      OutputPointType d_point;
      PointIdConstIterator current= cellIterator.Value()->PointIdsBegin();
      PointType A, B, C;
      
      A = primalMesh->GetPoint( *current );
      current++;
      B = primalMesh->GetPoint( *current );
      current++;
      C = primalMesh->GetPoint( *current );

      /*
      B[0] = B[0] - A[0]; 
      B[1] = B[1] - A[1]; 
      C[0] = C[0] - A[0];   
      C[1] = C[1] - A[1];   

      double D = 2 * (B[0]*C[1] - B[1]*C[0]);
      point[0] = ( C[1]*(B[0]*B[0]+B[1]*B[1]) - B[1]*(C[0]*C[0]+C[1]*C[1]))/D;
      point[1] = ( B[0]*(C[0]*C[0]+C[1]*C[1]) - C[0]*(B[0]*B[0]+B[1]*B[1]))/D;      
      
      d_point[0] = static_cast<OutputPixelType>( point[0] + A[0] );
      d_point[1] = static_cast<OutputPixelType>( point[1] + A[1] );
      d_point[2] = static_cast<OutputPixelType>( ( A[2] + B[2] + C[2] ) / 3 );
      */
       
      double xba, yba, zba, xca, yca, zca;
      double calength, balength;
      PointType crossbc;
      double denominator;
      double xcirca, ycirca, zcirca;
      PointType ta, tb, tc;
      PointType ba, ca;
      PointType circa;

      std::cout << " A = " << A[0] <<" "<< A[1] <<" " << A[2] << std::endl;
      std::cout << " B = " << B[0] <<" "<< B[1] <<" " << B[2] << std::endl;
      std::cout << " C = " << C[0] <<" "<< C[1] <<" " << C[2] << std::endl;
      
      /* Use coordinates relative to point `a' of the triangle. */
      for( unsigned int i = 0; i < PointDimension; i++ )
        {
        ba[i] = B[i] - A[i];  
        ca[i] = C[i] - A[i];
        std::cout << " ba et ca de " << i << " " << ba[i] << " " << ca[i] << std::endl;
        } 
      
      /* Squares of lengths of the edges incident to `a'. */
      balength = ba[0] * ba[0] + ba[1] * ba[1] + ba[2] * ba[2];
      calength = ca[0] * ca[0] + ca[1] * ca[1] + ca[2] * ca[2];
      std::cout << "alength " << balength << " " << calength << std::endl;
      
      ta[0] = B[1]; ta[1] = B[2];
      tb[0] = C[1]; tb[1] = C[2];
      tc[0] = A[1]; tc[1] = A[2];
      crossbc[0] = OrientationTest(ta, tb, tc);
      std::cout << "orientation test 1 " << crossbc[0] << std::endl;
      
      ta[0] = C[0]; ta[1] = C[2];
      tb[0] = B[0]; tb[1] = B[2];
      tc[0] = A[0]; tc[1] = A[2];
      crossbc[1] = OrientationTest(ta, tb, tc);
      std::cout << "orientation test 2 " << crossbc[1] << std::endl;
      
      ta[0] = B[0]; ta[1] = B[1];
      tb[0] = C[0]; tb[1] = C[1];
      tc[0] = A[0]; tc[1] = A[1];
      crossbc[2] = OrientationTest(ta, tb, tc);
      std::cout << "orientation test 3 " << crossbc[2] << std::endl;
      
      denominator = 0.5 / (crossbc[0] * crossbc[0] + crossbc[1] * crossbc[1] + crossbc[2] * crossbc[2]);
      std::cout << "dominator " << denominator <<std::endl;
      
      circa[0] = ((balength * ca[1] - calength * ba[1]) * crossbc[2] - (balength * ca[2] - calength * ba[2]) * crossbc[1]) * denominator;
      circa[1] = ((balength * ca[2] - calength * ba[2]) * crossbc[0] - (balength * ca[0] - calength * ba[0]) * crossbc[2]) * denominator;
      circa[2] = ((balength * ca[0] - calength * ba[0]) * crossbc[1] - (balength * ca[1] - calength * ba[1]) * crossbc[0]) * denominator;

      for( unsigned int i = 0; i < PointDimension; i++ )
        {
        d_point[i] = static_cast<OutputPixelType>( circa[i] + A[i] );
        }
      
      std::cout << "cirmcum point is " << d_point[0] << " " << d_point[1] << " " << d_point[2] << std::endl;
      
      return d_point;
    }
  };
  
} // namespace Functor
} // namespace itk

#endif // __CircumcentreDualPointFunctor_h
