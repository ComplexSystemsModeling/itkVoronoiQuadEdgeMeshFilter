/*================================================================================
 *                                                                               *
 * Delaunay Triangulation Test                                                   *
 *                                                                               *
 *   Implementation for ITK by Stéphane U. Rigaud and Alexandre Gouaillard       *
 *                                                                               *
 *   Inputs: "regular" + nbRows                                                  *
 *           "random"  + nbPoints                                                *
 *           "circle"  + radius                                                  *
 *           "cross"   + axeLength                                               *
 *                                                                               *
 *===============================================================================*/

#include "itkPointSet.h"
#include "itkQuadEdgeMesh.h"
#include "itkDelaunayConformingQuadEdgeMeshFilter.h"
#include "itkQuadEdgeMeshEulerOperatorsTestHelper.h"
#include "itkVTKPolyDataWriter.h"
#include "itkVTKPolyDataReader.h"

#include "itkPointSetToDelaunayTriangulationFilter.h"

#include <iostream>


//--------------------------------------------------------------------------------
// Random coordonates generation function
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
    std::cerr << "Usage - arg[1] : PointSet Generator - arg[2] : number of points ";
    std::cerr << std::endl;

    return EXIT_FAILURE;
    }

  typedef float PixelType;
  static const unsigned int Dimension = 2;

  typedef itk::QuadEdgeMesh< PixelType, Dimension >  MeshType;

  typedef itk::PointSet< PixelType, Dimension > PointSetType;
  typedef PointSetType::PointType                PointType;

  typedef itk::DelaunayConformingQuadEdgeMeshFilter< MeshType, MeshType > ValidityTestType;
  typedef itk::VTKPolyDataWriter< MeshType >  MeshWriterType;
  
  int meshSize       = atoi( argv[2] );
  int expectedNumPts = 0;

  std::vector< PointType > pts;
  PointSetType::Pointer pointSet     = PointSetType::New();
  MeshType::Pointer triangulatedMesh = MeshType::New();

  // -------------------------------------------------
  // Toy Point Set creation
    
  if( strcmp(argv[1],"regular") == 0 )
    {
    pts = GeneratePointCoordinates< PointSetType >( meshSize ); 
    }
  else if( strcmp(argv[1],"random") == 0 )
    {
    pts = GenerateRandomCoordinates< PointSetType >( meshSize ); 
    }
  else if( strcmp(argv[1],"circle") == 0 )
    {
    pts = GenerateCircleCoordinates< PointSetType >( meshSize );   
    }
  else if( strcmp(argv[1],"cross") == 0 )
    {
    pts = GenerateCrossCoordinates< PointSetType >( meshSize );
    }
  else
    {
    std::cerr << "Main: Wrong arguments for PointSet Generator" << std::endl;
    std::cerr << std::endl;
    return EXIT_FAILURE; 
    }
  
  expectedNumPts = pts.size();
  for( int i = 0; i < expectedNumPts ; i++ )
    {
    pointSet->SetPoint( i, pts[i] );
    }
  
  // -------------------------------------------------- 
  // Delaunay Construction
  
  typedef itk::PointSetToDelaunayTriangulationFilter< PointSetType > MyFilter;
  MyFilter::Pointer myfilter = MyFilter::New();
 
  myfilter->SetInput( pointSet );
  myfilter->SetDummyPoints( false );
  triangulatedMesh = myfilter->GetOutput();
  try 
    {
    myfilter->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Main: Exception thrown while processing the input." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }
  
  // -------------------------------------------------
  // Delaunay Validation

  ValidityTestType::Pointer validityTest = ValidityTestType::New();
  validityTest->SetInput( triangulatedMesh );
  validityTest->Update();

  if( validityTest->GetNumberOfEdgeFlips() > 0 )
    {
    return EXIT_FAILURE;
    }

  MeshWriterType::Pointer writer = MeshWriterType::New();
  writer->SetFileName( "DelaunayTriangulationMesh.vtk" );
  writer->SetInput( triangulatedMesh );
  try 
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Main: Exception caught while writting the output." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }
   
  // ------------------------------------------------
  // End process ... by by

  return EXIT_SUCCESS;

}
//--------------------------------------------------------------------------------
