/*===============================================================================
*                                                                               *
* Walking in a Triangulation algorithm implementation test                      *
*                                                                               *
*   Based on:                                                                   *
*   Devillers, O. and Pion, S. and Teillaud, M. "Walking in a triangulation",   *
*   Proceedings of the seventeenth annual symposium on Computational geometry,  *
*   pages 106-114, 2001                                                         *
*                                                                               *
*   Implementation for ITK by Stéphane Ulysse Rigaud                            *
*   IPAL (Image & Pervasive Access Lab) CNRS - A*STAR                           *
*   Singapore                                                                   *
*   http://www.ipal.cnrs.fr                                                     *
*                                                                               *
*   Input: double       X coordinate of the destination point                   *
*          double       Y coordinate of the destination point                   *
*          int          Starting cell id                                        *
*          int          path size for validation (optional)                     *
*          N int        expected size (optional)                                *
*                                                                               *
*===============================================================================*/

#include "itkPointSet.h"
#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshTraits.h"
#include "itkQuadEdgeMeshPolygonCell.h"
#include "itkVectorContainer.h"
#include "itkQuadEdgeMeshEulerOperatorsTestHelper.h"
 
#include "itkWalkInTriangulationFunction.h"

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
  typedef QEMeshType::CellIdentifier      CellIdentifier;
  
  typedef itk::VectorContainer< unsigned int, int > VectorContainerType;
  
  // -----------------------------------------------------
  // Initialisation mesh
      
  QEMeshType::Pointer mesh = QEMeshType::New();
  CreateSquareTriangularMesh< QEMeshType >( mesh );

  std::cout<<"\nNo. of Cells : "<<mesh->GetNumberOfCells()
           <<"\nNo. of Edges : "<<mesh->GetNumberOfEdges()
           <<"\nNo. of Faces : "<<mesh->GetNumberOfFaces()
           <<"\nNo. of Points : "<<mesh->GetNumberOfPoints()<<"\n\n";

  // -----------------------------------------------------
  // WalkInTriangulation 
  
  PointType      pts;
  CellIdentifier cell;
  
  pts[0] = atof( argv[1] ); 
  pts[1] = atof( argv[2] ); 
  cell   = atoi( argv[3] );
  
  VectorContainerType::Pointer resultPath = VectorContainerType::New();
  itk::WalkInTriangulationFunction< QEMeshType >::Pointer myFunction = 
	  itk::WalkInTriangulationFunction< QEMeshType >::New();

  try
    { 
    resultPath = myFunction->Evaluate( mesh, pts, cell );
    }
  catch( int e )
    {
    std::cerr << "error occured : "<< e << std::endl;
    }

  std::cout << "The point (" << pts[0] << ";" << pts[1] << ") is in the cell id ";
  std::cout << resultPath->GetElement( resultPath->Size()-1 ) << std::endl;
  
  // -----------------------------------------------------
  // Validation 
  
  if( argc > 5 )
    {
    std::cout << "expected path : ";
    VectorContainerType::Pointer expectedPath = VectorContainerType::New();
    for( unsigned int i = 0; i < (unsigned int)atoi( argv[4] ); i++ )
      {
      expectedPath->InsertElement(i,atoi(argv[5+i]));
      std::cout << argv[5+i] << " ";
      }
    std::cout << "\nresult path : ";
    VectorContainerType::Iterator ite = resultPath->Begin();
    while( ite != resultPath->End() )
      {
      std::cout << ite.Value() << " ";
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

