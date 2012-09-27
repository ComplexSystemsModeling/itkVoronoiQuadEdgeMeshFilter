/*===============================================================================
*                                                                               *
* Walking in a Triangulation Test                                               *
*                                                                               *
*   Based on:                                                                   *
*   Devillers, O. and Pion, S. and Teillaud, M. "Walking in a triangulation",   *
*   Proceedings of the seventeenth annual symposium on Computational geometry,  *
*   pages 106-114, 2001                                                         *
*                                                                               *
*   Implementation for ITK by Stéphane U. Rigaud and Alexandre Gouaillard       *
*                                                                               *
*   Input: double       X coordinate of the destination point                   *
*          double       Y coordinate of the destination point                   *
*          int          Starting cell id                                        *
*          int          path size for validation (optional)                     *
*          N int        expected size (optional)                                *
*                                                                               *
*===============================================================================*/

//--------------------------------------------------------
// ITK includes
#include "itkPointSet.h"
#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshTraits.h"
#include "itkQuadEdgeMeshPolygonCell.h"
#include "itkVectorContainer.h"
#include "itkQuadEdgeMeshEulerOperatorsTestHelper.h"
 
//--------------------------------------------------------
// Our includes
#include "itkWalkInTriangulationFunction.h"

//--------------------------------------------------------
// STD includes
#include <iostream>

int main( int argc, char * argv[] )
{

  if( argc < 3 )
    {
    if( argc < 4 )
      {
      std::cerr<<"Usage "<<std::endl;
      std::cerr<<argv[0]<<" DestinationPointX DestinationPointY InitialCellId ( nbIds [id1 id2 ... idn] )"<<std::endl;
      return EXIT_FAILURE;
      }
    else
      {
      if( argc != 5 + atoi( argv[4] ) )
        {
        std::cerr << "Usage error " << std::endl;
        std::cerr << "Yous declared " << argv[4] << " ids but only provide " << argc-5 << std::endl;
        std::cerr << "Usage " << std::endl;
        std::cerr << argv[0] << " DestinationPointX DestinationPointY InitialCellId ( nbIds [id1 id2 ... idn] )" << std::endl;
        return EXIT_FAILURE;
        }
      }
    }

  typedef itk::QuadEdgeMesh< double, 3 >  QEMeshType;
  typedef QEMeshType::PointType           PointType;
  typedef QEMeshType::FaceRefType      FaceRefType;
  
  typedef itk::VectorContainer< unsigned int, FaceRefType > VectorContainerType;
  
  // -----------------------------------------------------
  // Initialisation of toy mesh 
  // -----------------------------------------------------    
  
  QEMeshType::Pointer mesh = QEMeshType::New();
  CreateSquareTriangularMesh< QEMeshType >( mesh );

  // -----------------------------------------------------
  // WalkInTriangulation process
  // ----------------------------------------------------- 
  
  PointType   pts;
  FaceRefType cell;
  
  pts[0] = atof( argv[1] ); 
  pts[1] = atof( argv[2] ); 
  cell.first = atoi( argv[3] );
  
  VectorContainerType::Pointer resultPath = VectorContainerType::New();
  itk::WalkInTriangulationFunction< QEMeshType >::Pointer myFunction = 
                 itk::WalkInTriangulationFunction< QEMeshType >::New();

  try
    { 
    resultPath = myFunction->Evaluate( mesh, pts, &cell );
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Main: Exception thrown while processing the input." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Main: The point (" << pts[0] << ";" << pts[1] << ") is in the cell id ";
  std::cout << resultPath->GetElement( resultPath->Size()-1 ).first << std::endl;
  
  // -----------------------------------------------------
  // WalkInTriangulation Validation 
  // -----------------------------------------------------
  
  if( argc > 5 )
    {
    std::cout << "Main: Expected path : ";
    VectorContainerType::Pointer expectedPath = VectorContainerType::New();
    FaceRefType tempRef;
    for( unsigned int i = 0; i < (unsigned int)atoi( argv[4] ); i++ )
      {
      tempRef.first = atoi(argv[5+i]);
      expectedPath->InsertElement(i,tempRef);
      std::cout << argv[5+i] << " ";
      }
    std::cout << std::endl;
    std::cout << "Main: Result path : ";
    VectorContainerType::Iterator ite = resultPath->Begin();
    while( ite != resultPath->End() )
      {
      std::cout << ite.Value().first << " ";
      ++ite;
      }
    std::cout << std::endl;
    
    VectorContainerType::Iterator ite1 = resultPath->Begin();
    VectorContainerType::Iterator ite2 = expectedPath->Begin();
    if( resultPath->Size() == expectedPath->Size() )
      {
      while( ite1 != resultPath->End() && ite2 != expectedPath->End() )
        {
        if( ite1.Value() != ite2.Value() )
          {
          return EXIT_FAILURE;
          }
        ++ite1; 
        ++ite2;      
        }
      return EXIT_SUCCESS;
      }
    else
      {
      return EXIT_FAILURE;
      }
    }
  else
    {
    return EXIT_SUCCESS;
    }
}

