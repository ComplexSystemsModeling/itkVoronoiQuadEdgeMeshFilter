#include "itkQuadEdgeMeshWithDual.h"
#include "itkQuadEdgeMeshToQuadEdgeMeshWithDualFilter.h"
#include "itkQuadEdgeMeshWithDualAdaptor.h"
#include "itkVTKPolyDataWriter.h"
#include "itkQuadEdgeMeshEulerOperatorsTestHelper.h"

#include "itkBarycentreDualPointFunctor.h"
#include "itkCircumcentreDualPointFunctor.h"

#include "itkPointSetToDelaunayTriangulationFilter.h"

#include <iostream>

#include "itkRegularSphereMeshSource.h"

int main( int argc, char * argv[] )
{
  
  if( argc > 2 )
    {
    std::cerr << "Exit wrong arguments" << std::endl;
    std::cerr << "Usage - please provide a number of point." << std::endl;
    std::cerr << "Default - No argument needed" << std::endl;
    
    return EXIT_FAILURE;
    }
  
  unsigned int NumberOfPoints = 0;
  if( argc == 2 )
    {
    NumberOfPoints = (unsigned int) atoi(argv[1]);  
    }
  
  //-----------------------------------------
  //  Define all the types we will work with
  //-----------------------------------------

  // two base tyes for meshes
  const unsigned int dimension = 3;
  typedef float PixelType;

  // the mesh type
  typedef itk::QuadEdgeMeshWithDual< PixelType, dimension > SimplexMeshType;
  typedef itk::QuadEdgeMesh< PixelType, dimension > MeshType;
  // test mesh sphere
  typedef itk::RegularSphereMeshSource< MeshType > SphereMesh;
  // the pointset to primal filter  
  typedef itk::PointSetToDelaunayTriangulationFilter< MeshType > PrimalFilterType;
  // the primal to primal+dual filter
  typedef itk::Functor::CircumcentreDualPointFunctor< MeshType, SimplexMeshType > DualFunctor;
  //typedef itk::Functor::BarycentreDualPointFunctor< MeshType, SimplexMeshType > DualFunctor;
  typedef itk::QuadEdgeMeshToQuadEdgeMeshWithDualFilter< MeshType, SimplexMeshType, DualFunctor > DualFilterType;

  
  // all the filters to write the result
  typedef itk::VTKPolyDataWriter<SimplexMeshType> MeshWriterType;
  typedef itk::VTKPolyDataWriter<MeshType> WriterType;

  typedef itk::QuadEdgeMeshWithDualAdaptor< SimplexMeshType >  AdaptorType;
  typedef itk::VTKPolyDataWriter< AdaptorType > DualMeshWriterType;

  //--------------------------------------------------------------
  // Create the DataStructures that will hold the data in memory
  //--------------------------------------------------------------

  MeshType::Pointer myPrimalMesh = MeshType::New();  

  //-----------------------------------------
  // Create an input mesh to toy around with
  //-----------------------------------------

  std::cout << "Main: Create a planar PointSet." << std::endl;

  // NOTE ALEX: is this cross platform? 
  srand(time(NULL));
 
  if( NumberOfPoints > 0 )
    {
    std::cout << "Main: Create a random planar PointSet." << std::endl;
    std::vector< MeshType::PointType > oPt( NumberOfPoints );
    unsigned int i = 0;
    while( i < NumberOfPoints )
      {
      oPt[i][0] = static_cast<PixelType>( rand() % 1000 - 500 );  
      oPt[i][1] = static_cast<PixelType>( rand() % 1000 - 500 );
      oPt[i][2] = static_cast<PixelType>( rand() % 100 - 50 );
      myPrimalMesh->SetPoint( i, oPt[i] );
      i++;
      }
    }
  else 
    {
    std::cout << "Main: Create a 3D pyramid PointSet." << std::endl;
    std::vector< MeshType::PointType > oPt( 5 );
    oPt[0][0] = 0;    oPt[0][1] = 0;    oPt[0][2] = 0;
    oPt[1][0] = 0.5;  oPt[1][1] = 0.5;  oPt[1][2] = 0;
    oPt[2][0] = 0.5;  oPt[2][1] = -0.5; oPt[2][2] = 0;
    oPt[3][0] = -0.5; oPt[3][1] = -0.5; oPt[3][2] = 0;
    oPt[4][0] = -0.5; oPt[4][1] = 0.5;  oPt[4][2] = 0;
    for( unsigned int i = 0; i< 5; i++ )
      {
      myPrimalMesh->SetPoint( i, oPt[i] );
      }
    }

  WriterType::Pointer test = WriterType::New();
  test->SetInput(myPrimalMesh);
  test->SetFileName( "test.vtk" );
  test->Update();
  
  std::cout << "Main: Create a Delaunay triangulation planar mesh." << std::endl;
    
  PrimalFilterType::Pointer myPrimalFilter = PrimalFilterType::New();
  myPrimalFilter->SetInput( myPrimalMesh );
  try 
    {
    myPrimalFilter->Update();
    }
  catch (int e) 
    {
    std::cerr << "Main: Error catch in PointSetToDelaunayTriangulationFilter." << std::endl;
    }
  
  SphereMesh::Pointer sphere = SphereMesh::New();
  sphere->SetResolution( 4 );
  sphere->Update();
  
  //------------
  // Do the job
  //------------

  std::cout << "Main: Apply filter." << std::endl;
  DualFilterType::Pointer myDualFilter = DualFilterType::New();
  myDualFilter->SetInput( myPrimalFilter->GetOutput() );
  try
    {
    myDualFilter->Update( );
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Main: Exception thrown while processing the input." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  //-----------------------------------------------------
  // Write the original mesh, and the computed dual mesh
  //-----------------------------------------------------

  std::cout << "Main: Write primal mesh." << std::endl;
  MeshWriterType::Pointer writer1 = MeshWriterType::New();
  writer1->SetInput( myDualFilter->GetOutput() );
  writer1->SetFileName( "PrimalOutput.vtk" );
  writer1->Write();

  AdaptorType* adaptor = new AdaptorType();
  adaptor->SetInput( myDualFilter->GetOutput() );

  std::cout << "Main: Write dual mesh." << std::endl;
  DualMeshWriterType::Pointer writer2 = DualMeshWriterType::New();
  writer2->SetInput( adaptor );
  writer2->SetFileName( "DualOutput.vtk" );
  writer2->Write();

  //-----------------------------------------------------
  // and ... we're outta here.
  //-----------------------------------------------------

  return EXIT_SUCCESS;
}
//--------------------------------------------------------------------------------------------------
