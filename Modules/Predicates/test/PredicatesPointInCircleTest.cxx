
#include <iostream>
#include "itkMesh.h"
#include "itkQuadEdgeMesh.h"

#include "itkOrientation2DFunction.h"
#include "itkInCircleTestFunction.h"
#include "itkInCircleExactTestFunction.h"

#include <itkMatrix.h>
#include "vnl/vnl_det.h"
#include "vnl/vnl_matrix_fixed.h"

template< class MeshType >
bool test( double eps_x, double eps_y, bool exact )
{
  
  typedef typename MeshType::PointType PointType;
  typedef typename MeshType::PixelType PixelType;
  
  typedef typename MeshType::CellType   CellType; 
  typedef typename MeshType::CellTraits CellTraits;
  typedef typename MeshType::CellIdentifier  CellIdentifier;
  typedef typename MeshType::PointIdentifier PointIdentifier;
  typedef typename CellType::CellAutoPointer CellAutoPointer; 

  typedef itk::CellInterface< PixelType, CellTraits > CellInterfaceType;
  typedef itk::TriangleCell< CellInterfaceType >     TriangleCellType;
  typedef typename TriangleCellType::PointIdIterator PointIdIterator;
  
  typedef itk::InCircleExactTestFunction< MeshType > InCircleExactTest;
  typedef itk::InCircleTestFunction< MeshType > InCircleTest;
  
  PointType mpa, mpb, mpc, mpd;
  mpa[0] = mpc[0] = mpc[1] = mpb[1] = 1.0;
  mpa[1] = mpb[0] =                   0.0;
  mpd[0] = eps_x;
  mpd[1] = eps_y;
    
  bool result_1;
  typename InCircleTest::Pointer circle = InCircleTest::New();
  typename InCircleExactTest::Pointer circleExact = InCircleExactTest::New();
  std::vector< PointType > pointList;
  pointList.push_back(mpa);
  pointList.push_back(mpb);
  pointList.push_back(mpc);
  pointList.push_back(mpd);

  if( exact )
    {
    result_1 = circleExact->Evaluate( pointList );
    }
  else
    {
    result_1 = circle->Evaluate( pointList );
    }

  return result_1;
}

int PredicatesPointInCircleTest( int argc, char* argv[] )
{
  if( argc < 4 )
    {
    std::cout << "usage: <exe> epsilon_x epsilon_y Exact";
    std::cout << std::endl;
    return EXIT_FAILURE;
    }
  
  char value = 0;
  char* ptr;
  ptr = &value;
  char** stringPtr;
  stringPtr = &ptr;
  
  double epsilon_x = std::strtod( argv[1], stringPtr );
  double epsilon_y = std::strtod( argv[2], stringPtr );
  bool TestExact = atoi(argv[3]);
  
  std::cout << "Epsilon_x: " << epsilon_x << std::endl;
  std::cout << "Epsilon_y: " << epsilon_y << std::endl;
  std::cout << "Exact Testing? " << TestExact << std::endl;
  
  return(
         test< itk::Mesh< float,  2 > >( epsilon_x, epsilon_y, TestExact )
         | test< itk::Mesh< float,  3 > >( epsilon_x, epsilon_y, TestExact )
         | test< itk::Mesh< float,  4 > >( epsilon_x, epsilon_y, TestExact )
         | test< itk::Mesh< double, 2 > >( epsilon_x, epsilon_y, TestExact )
         | test< itk::Mesh< double, 3 > >( epsilon_x, epsilon_y, TestExact )
         | test< itk::Mesh< double, 4 > >( epsilon_x, epsilon_y, TestExact )
         | test< itk::QuadEdgeMesh< float,  2 > >( epsilon_x, epsilon_y, TestExact )
         | test< itk::QuadEdgeMesh< float,  3 > >( epsilon_x, epsilon_y, TestExact )
         | test< itk::QuadEdgeMesh< float,  4 > >( epsilon_x, epsilon_y, TestExact )
         | test< itk::QuadEdgeMesh< double, 2 > >( epsilon_x, epsilon_y, TestExact )
         | test< itk::QuadEdgeMesh< double, 3 > >( epsilon_x, epsilon_y, TestExact )
         | test< itk::QuadEdgeMesh< double, 4 > >( epsilon_x, epsilon_y, TestExact )
         );
}
