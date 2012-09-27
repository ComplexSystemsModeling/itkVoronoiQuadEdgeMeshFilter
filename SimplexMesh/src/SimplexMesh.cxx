/*=================================================================
 *                                                                *
 * Simplex Mesh Test                                              *
 *                                                                *
 *   Implementation for ITK by Humayun Irshad, Stéphane U. Rigaud *
 *   and Alexandre Gouaillard                                     *
 *                                                                *
 *================================================================*/

// ITK Includes
#include "itkVTKPolyDataWriter.h"
#include "itkVTKPolyDataReader.h"
#include "itkQuadEdgeMeshEulerOperatorsTestHelper.h"
#include "itkRegularSphereMeshSource.h"

// Our Includes
#include "itkQuadEdgeMeshWithDual.h"
#include "itkQuadEdgeMeshToQuadEdgeMeshWithDualFilter.h"
#include "itkQuadEdgeMeshWithDualAdaptor.h"
#include "itkBarycentreDualPointFunctor.h"
#include "itkCircumcentreDualPointFunctor.h"

#include "itkPointSetToDelaunayTriangulationFilter.h"

// STD Includes
#include <iostream>

int main( int argc, char * argv[] )
{
  
  if( argc != 2 )
    {
    std::cerr << "Exit wrong arguments" << std::endl;
    std::cerr << "Usage argument - Provide a type of mesh as argument ";
    std::cerr << "\"plan\" - \"sphere\" ";
    std::cerr << "or provide an input mesh file." << std::endl;
      
    return EXIT_FAILURE;
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
  
  // the dual point functor
  typedef itk::Functor::BarycentreDualPointFunctor< MeshType, SimplexMeshType > BaryDualFunctor;
  typedef itk::Functor::CircumcentreDualPointFunctor< MeshType, SimplexMeshType > CircumDualFunctor;

  // the pointset to primal filter  
  typedef itk::PointSetToDelaunayTriangulationFilter< MeshType > PrimalFilterType;
  
  // the primal to primal+dual filter
  typedef itk::QuadEdgeMeshToQuadEdgeMeshWithDualFilter< MeshType, SimplexMeshType, BaryDualFunctor > BaryDualFilterType;
  typedef itk::QuadEdgeMeshToQuadEdgeMeshWithDualFilter< MeshType, SimplexMeshType, CircumDualFunctor > CircumDualFilterType;

  typedef itk::QuadEdgeMeshWithDualAdaptor< SimplexMeshType >  AdaptorType;
  
  // all the filters to write the result
  typedef itk::VTKPolyDataWriter<SimplexMeshType> MeshWriterType;
  typedef itk::VTKPolyDataWriter<MeshType> WriterType;
  typedef itk::VTKPolyDataReader<MeshType> ReaderType;
  typedef itk::VTKPolyDataWriter< AdaptorType > DualMeshWriterType;

  //--------------------------------------------------------------
  // Create the DataStructures that will hold the data in memory
  //--------------------------------------------------------------

  MeshType::Pointer myPrimalMesh = MeshType::New();  
  ReaderType::Pointer meshReader = ReaderType::New();

  //-----------------------------------------
  // Create an input mesh to toy around with
  //-----------------------------------------

  std::cout << "Main: Create an input mesh." << std::endl;
 
  if( strcmp(argv[1],"plan") == 0 )
    {
    std::cout << "Main: Create a regular planar mesh." << std::endl;
    CreateSquareTriangularMesh< MeshType >( myPrimalMesh );
     
    std::cout << "Main: Poke a hole in the mesh to have two boundaries." << std::endl;
    myPrimalMesh->LightWeightDeleteEdge( myPrimalMesh->FindEdge( 11, 12 ) );
    }
  else if( strcmp(argv[1],"sphere") == 0 )  
    {
    std::cout << "Main: Create a regular sphere mesh." << std::endl;
    SphereMesh::Pointer sphere = SphereMesh::New();
    sphere->SetResolution( 4 );
    sphere->Update();
    myPrimalMesh = sphere->GetOutput();
    }
  else 
    {
      meshReader->SetFileName( argv[1] );
      try 
        {
        meshReader->Update();
        myPrimalMesh = meshReader->GetOutput();
        }
      catch( itk::ExceptionObject & excp ) 
        {
        std::cerr << "Main: Exception thrown while trying to read the input." << std::endl;
        std::cerr << excp << std::endl;
        return EXIT_FAILURE; 
        }
    }
  
  //------------
  // Do the job
  //------------

  std::cout << "Main: Apply the Primal to Dual filter." << std::endl;
  BaryDualFilterType::Pointer myBaryDualFilter = BaryDualFilterType::New();
  myBaryDualFilter->SetInput( myPrimalMesh );
  myBaryDualFilter->SetMakeBorders( true );
  try
    {
    myBaryDualFilter->Update( );
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Main: Exception thrown while processing the input." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }
  
  CircumDualFilterType::Pointer myCircumDualFilter = CircumDualFilterType::New();
  myCircumDualFilter->SetInput( myPrimalMesh );
  try
    {
    myCircumDualFilter->Update( );
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
  writer1->SetInput( myBaryDualFilter->GetOutput() );
  writer1->SetFileName( "PrimalMesh.vtk" );
  writer1->Write();

  AdaptorType* adaptor = new AdaptorType();
  adaptor->SetInput( myBaryDualFilter->GetOutput() );

  std::cout << "Main: Write bary dual mesh." << std::endl;
  DualMeshWriterType::Pointer writer2 = DualMeshWriterType::New();
  writer2->SetInput( adaptor );
  writer2->SetFileName( "BaryDualMesh.vtk" );
  writer2->Write();
  
  adaptor = new AdaptorType();
  adaptor->SetInput( myCircumDualFilter->GetOutput() );
  
  std::cout << "Main: Write circum dual mesh." << std::endl;
  writer2 = DualMeshWriterType::New();
  writer2->SetInput( adaptor );
  writer2->SetFileName( "CircumDualMesh.vtk" );
  writer2->Write();

  //-----------------------------------------------------
  // and ... we're outta here.
  //-----------------------------------------------------

  return EXIT_SUCCESS;
}
//--------------------------------------------------------------------------------------------------
