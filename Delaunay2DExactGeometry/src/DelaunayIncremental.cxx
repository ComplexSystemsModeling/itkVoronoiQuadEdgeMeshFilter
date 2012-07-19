
#include "itkPointSet.h"
#include "itkQuadEdgeMesh.h"
#include "itkDelaunayConformingQuadEdgeMeshFilter.h"
#include "itkQuadEdgeMeshEulerOperatorsTestHelper.h"
#include "itkVTKPolyDataWriter.h"
#include "itkVTKPolyDataReader.h"

#include "itkQuadEdgeMeshToDelaunayTriangulationFilter.h"

#include <iostream>


//--------------------------------------------------------------------------------
// Random coordonates generation function
//
//--------------------------------------------------------------------------------
template< class TMesh >
std::vector< typename TMesh::PointType >
GenerateRandomCoordinates( const unsigned int& iN )
{
  typedef typename TMesh::PointType         TPointType;
  typedef typename TPointType::CoordRepType TCoordRepType;
  std::vector< TPointType > oPt( iN );

  // NOTE ALEX: is this cross platform? 
  srand(time(NULL));

  for( unsigned int i = 0; i < iN; i++ )
    {
    oPt[ i ][0] = static_cast< TCoordRepType >( rand() % 9000 - 4500 );
    oPt[ i ][1] = static_cast< TCoordRepType >( rand() % 9000 - 4500 );
    oPt[ i ][2] = static_cast< TCoordRepType >( 0. );
    }
  return oPt;
}
//--------------------------------------------------------------------------------


//--------------------------------------------------------------------------------
// Circle coordonates generation function
//--------------------------------------------------------------------------------
template< class TMesh >
std::vector< typename TMesh::PointType >
GenerateCircleCoordinates( const unsigned int& r )
{
  typedef typename TMesh::PointType         TPointType;
  typedef typename TPointType::CoordRepType TCoordRepType;
  std::vector< TPointType > oPt;
  const float DEG2RAD = 3.14159 / 180;
  TPointType p;

  for( unsigned int i = 0; i < 360; i++ )
    {
    float deg = i * DEG2RAD; 
    p[0] = static_cast< TCoordRepType >( cos(deg)*r );
    p[1] = static_cast< TCoordRepType >( sin(deg)*r );
    p[2] = static_cast< TCoordRepType >( 0. );
    oPt.push_back( p );
    }
  return oPt;
}
//--------------------------------------------------------------------------------


//--------------------------------------------------------------------------------
// Concentric Circle coordonates generation function
//-------------------------------------------------------------------------------
template< class TMesh >
std::vector< typename TMesh::PointType >
GenerateConcentricCoordinates( const unsigned int& r )
{
  typedef typename TMesh::PointType         TPointType;
  typedef typename TPointType::CoordRepType TCoordRepType;
  std::vector< TPointType > oPt;
  const float DEG2RAD = 3.14159 / 180;
  TPointType p;

  for( unsigned int j = r; j > 0; j-- )
    {
    for( unsigned int i = 0; i < 360; i++ )
      {
      float deg = i * DEG2RAD; 
      p[0] = static_cast< TCoordRepType >( cos(deg)*j );
      p[1] = static_cast< TCoordRepType >( sin(deg)*j );
      p[2] = static_cast< TCoordRepType >( 0. );
      oPt.push_back( p ); 
      }
    }
  p[0] = 0.; p[1] = 0.; p[2] = 0.;
  oPt.push_back( p ); 
  return oPt;
}
//--------------------------------------------------------------------------------


//--------------------------------------------------------------------------------
// Cross coordonates generation function
//--------------------------------------------------------------------------------
template< class TMesh >
std::vector< typename TMesh::PointType >
GenerateCrossCoordinates( const unsigned int& iN )
{
  typedef typename TMesh::PointType         TPointType;
  typedef typename TPointType::CoordRepType TCoordRepType;
  std::vector< TPointType > oPt;
  TPointType p;

  for( unsigned int i = 0 ; i < 2*iN+1; i++ )
    {
    p[0] = static_cast< TCoordRepType >( i  );
    p[1] = static_cast< TCoordRepType >( iN );
    p[2] = static_cast< TCoordRepType >( 0. );
    oPt.push_back( p );
    if( i != iN )
      {
      p[0] = static_cast< TCoordRepType >( iN );
      p[1] = static_cast< TCoordRepType >( i  );
      p[2] = static_cast< TCoordRepType >( 0. );
      oPt.push_back( p );
      }
    }
  return oPt;
}
//--------------------------------------------------------------------------------


//--------------------------------------------------------------------------------
// Main function
//--------------------------------------------------------------------------------
int
main( int argc, char* argv[] )
{
  if( argc != 3 )
    {
    std::cerr << "Exit wrong arguments" << std::endl;
    std::cerr << "Usage - arg[1] [1|2] = [\"reg\"|\"rand\"] - arg[2] int = [rowsize|nbPoint]";
    std::cerr << std::endl;

    return EXIT_FAILURE;
    }

  typedef float PixelType;
  static const unsigned int Dimensions = 3;

  typedef itk::QuadEdgeMesh< PixelType, Dimensions >  MeshType;
  typedef MeshType::PointType                         PointType;

  typedef itk::DelaunayConformingQuadEdgeMeshFilter< MeshType, MeshType > ValidityTestType;
  typedef itk::VTKPolyDataWriter< MeshType >                              MeshWriterType;

  int type           = atoi( argv[1] );
  int meshSize       = atoi( argv[2] );
  int expectedNumPts = 0;

  std::vector< PointType > pts;

  MeshType::Pointer pointSet         = MeshType::New();
  MeshType::Pointer triangulatedMesh = MeshType::New();

  // -------------------------------------------------
  // Toy Point Set creation
  
  switch(type) 
    {
    //case 1 :     
    //pts = GeneratePointCoordinates< MeshType >( meshSize ); 
    //break;
    case 2 :
      pts = GenerateRandomCoordinates< MeshType >( meshSize ); 
      break;
    case 3 :
      pts = GenerateCircleCoordinates< MeshType >( meshSize );    
      break;
    case 4 :
      pts = GenerateConcentricCoordinates< MeshType >( meshSize );
      break;
    case 5 :
      pts = GenerateCrossCoordinates< MeshType >( meshSize );
      break;
    default:
      std::cerr << "Exit wrong arguments" << std::endl;
      std::cerr << "Usage - arg[1] [1|2] = [\"reg\"|\"rand\"] - arg[2] int = [rowsize|nbPoint]";
      std::cerr << std::endl;
      
      return EXIT_FAILURE;
    }

  expectedNumPts = pts.size();
  for( int i = 0; i < expectedNumPts ; i++ )
    {
    pointSet->SetPoint( i, pts[i] );
    }

  MeshWriterType::Pointer writeToyMesh = MeshWriterType::New();
  writeToyMesh->SetFileName("./InputToyMesh.vtk");
  writeToyMesh->SetInput( pointSet );
  writeToyMesh->Update();
  
  // -------------------------------------------------- 
  // Delaunay Construction

  typedef itk::QuadEdgeMeshToDelaunayTriangulationFilter< MeshType, MeshType > MyFilter;
  MyFilter::Pointer myfilter = MyFilter::New();
  
  myfilter->SetInput( pointSet );
  triangulatedMesh = myfilter->GetOutput();
  myfilter->Update();
    
  MeshWriterType::Pointer writeDelaunay = MeshWriterType::New();
  writeDelaunay->SetFileName("./OutputDelaunayMesh.vtk");
  writeDelaunay->SetInput( triangulatedMesh );
  writeDelaunay->Update();
  
  // -------------------------------------------------
  // Delaunay Validation

  ValidityTestType::Pointer validityTest = ValidityTestType::New();
  validityTest->SetInput( triangulatedMesh );
  validityTest->GraftOutput( triangulatedMesh );
  validityTest->Update();

  if( validityTest->GetNumberOfEdgeFlips() > 0 )
    {
    return EXIT_FAILURE;
    }

  // ------------------------------------------------
  // End process ... by by

  return EXIT_SUCCESS;

}
//--------------------------------------------------------------------------------
