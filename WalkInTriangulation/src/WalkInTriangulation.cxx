/*=========================================================================
*
*    Walking in a Triangulation algorithm implementation
*   Based on :
*
*   Devillers, O. and Pion, S. and Teillaud, M. "Walking in a triangulation",
*   Proceedings of the seventeenth annual symposium on Computational geometry,
*   pages 106-114, 2001
*
*   Implementation for ITK by Stéphane Ulysse Rigaud
*   IPAL (Image & Pervasive Access Lab) CNRS - A*STAR
*   Singapore
*
*   Input string  Mtk mesh name
*         double  X coordinate of the destination point in the mesh
*         double  Y coordinate of the destination point in the mesh
*         unsigned int Starting cell id
*
*   Output int I index of the cell that contain the destination point
*          if I = -1, destination point is outside the mesh
*
*   TODO => Validation test
*        => Functorise the implementation
*
*=========================================================================*/

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
  if( argc != 5 )
    {
    if( argc < 5 )
      {
      std::cerr<<"Usage "<<std::endl;
      std::cerr<<argv[0]<<" Mesh DestinationPointX DestinationPointY InitialCellId ( nbIds id1 id2 ... idn )"<<std::endl;
      return EXIT_FAILURE;
      }
    else
      {
      if( argc != 6 + atoi( argv[5] ) )
        {
        std::cerr << "Usage error " << std::endl;
        std::cerr << "Yous declared " << argv[5] << " ids but only provide " << argc-6 << std::endl;
        std::cerr << "Usage " << std::endl;
        std::cerr << argv[0] << " Mesh DestinationPointX DestinationPointY InitialCellId ( nbIds id1 id2 ... idn )" << std::endl;
        return EXIT_FAILURE;
        }
      }
    }

  typedef itk::QuadEdgeMesh< double, 3 >  QEMeshType;
  typedef QEMeshType::PointType           PointType;
  typedef QEMeshType::CellIdentifier      CellIdentifier;
  
  typedef itk::VectorContainer< unsigned int, unsigned int > VectorContainerType;
  
  // -----------------------------------------------------
  // Initialisation
      
  QEMeshType::Pointer mesh = QEMeshType::New();
  
  std::cout<<"\nNo. of Cells : "<<mesh->GetNumberOfCells()
           <<"\nNo. of Edges : "<<mesh->GetNumberOfEdges()
           <<"\nNo. of Faces : "<<mesh->GetNumberOfFaces()
           <<"\nNo. of Points : "<<mesh->GetNumberOfPoints()<<"\n\n";


  // -----------------------------------------------------
  // WalkInTriangulation Test
  
  PointType      pts;
  CellIdentifier cell;
  
  pts[0] = atof( argv[2] ); 
  pts[1] = atof( argv[3] ); 
  cell   = atoi( argv[4] );
  
  VectorContainerType::Pointer resultPath = VectorContainerType::New();
  WalkInTriangulationFunction< QEMeshType > myFunction;
  resultPath = myFunction.Evaluate( mesh, pts, cell );
  
  std::cout << "The point (" << pts[0] << ";" << pts[1] << ") is in the cell id ";
  std::cout << resultPath->GetElement( resultPath->Size()-1 ) << std::endl;
  
  // -----------------------------------------------------
  // Path Check
  
  if( argc > 6 )
    {
    std::cout << "expected path : ";
    VectorContainerType::Pointer expectedPath = VectorContainerType::New();
    for( unsigned int i = 0; i < (unsigned int)atoi( argv[5] ); i++ )
      {
      expectedPath->InsertElement(i,atoi(argv[6+i]));
      std::cout << argv[6+i] << " ";
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
    return EXIT_SUCCESS;
    }
}

