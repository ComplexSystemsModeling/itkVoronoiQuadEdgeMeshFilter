#include "itkPointSet.h"
#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshTraits.h"
#include "itkQuadEdgeMeshPolygonCell.h"
#include "itkVectorContainer.h"
#include "itkQuadEdgeMeshEulerOperatorsTestHelper.h"
#include "itkQuadEdgeMeshEulerOperatorFlipEdgeFunction.h"

//---------------
// our code
#include "itkWalkInTriangulationFunction.h"
#include "itkPointInCircleGeometricalPredicateFunctor.h"

//---------------
#include <iostream>
#include <limits>

//TEMPORARY
#include "itkVTKPolyDataWriter.h"

//--------------------------------------------------------------------------------
// O ring check function
// TO BE DELETED
//
template< class TMeshType >
void CheckONextLink( typename TMeshType::Pointer mesh )
{
  typedef          TMeshType                         MeshType;
  typedef typename MeshType::PointsContainerIterator PointsContainerIterator;
  typedef typename MeshType::QEType                  QuadEdgeType;
   
  std::cout << "Checking consistency of all O rings." << std::endl;
  PointsContainerIterator myPointIterator = mesh->GetPoints()->Begin();
  
  while( myPointIterator != mesh->GetPoints()->End() )
    {
    
    QuadEdgeType* e = myPointIterator.Value().GetEdge();
    QuadEdgeType* next = e->GetOnext();
    while( e != next )
      {
      next = next->GetOnext();
      } 
    ++myPointIterator;
    }
  std::cout << "Checking consistency of all O rings  - DONE." << std::endl;
}
//--------------------------------------------------------------------------------


//--------------------------------------------------------------------------------
// Random coordonates generation function
//
template< class TMesh >
std::vector< typename TMesh::PointType >
GenerateRandomPointCoordinates( const unsigned int& iN )
{
  typedef typename TMesh::PointType        PointType;
  typedef typename PointType::CoordRepType CoordRepType;
  std::vector< PointType > oPt( iN * iN );
 
  // NOTE ALEX: is this cross platform? 
  srand(time(NULL));

  for( unsigned int i = 0; i < iN; i++ )
    {
    oPt[ i ][0] = static_cast< CoordRepType >( rand() % 9000 - 4500 );
    oPt[ i ][1] = static_cast< CoordRepType >( rand() % 9000 - 4500 );
    oPt[ i ][2] = static_cast< CoordRepType >( 0. );
    }
  return oPt;
}
//--------------------------------------------------------------------------------


//--------------------------------------------------------------------------------
// Dummy Mesh creation function
//
// Create a dummy square mesh for triangulation initialisation
// limit argument is to determine the size of the square
//
template< class TMeshType >
void
CreateDummyMesh( typename TMeshType::Pointer mesh, typename TMeshType::PixelType limit )
{
  // Generate a square triangulated mesh
  //
  //   2---------3
  //   |       / |
  //   |  1  /   |
  //   |   /     |
  //   | /    0  |
  //   0---------1
  //
  // Anti-ClockWise orientation
  // Cell 0 -> 013
  // Cell 1 -> 032
  
  typedef          TMeshType           MeshType;
  typedef typename MeshType::CellType  CellType;
  typedef typename MeshType::PointType PointType;
  typedef typename MeshType::PixelType PixelType;
  
  typedef typename CellType::CellAutoPointer CellAutoPointer;
  
  typedef itk::QuadEdgeMeshPolygonCell< CellType > QEPolygonCellType;
  
  if( mesh->GetNumberOfPoints() )
    {
    mesh->Clear();
    mesh->ClearFreePointAndCellIndexesLists();
    }
  
  int expectedNumPts = 4;
  int expectedNumCells = 2;
  int simpleTriangleCells[6] = { 0, 1, 3, 0, 3, 2 };
  
  PixelType min = -limit;  
  PixelType max =  limit;
  
  std::vector< PointType > pts( 4 );
  int i( 0 );
  
  pts[i][0] = min; pts[i][1] = min; pts[i++][2] = 0.;
  pts[i][0] = max; pts[i][1] = min; pts[i++][2] = 0.; 
  pts[i][0] = min; pts[i][1] = max; pts[i++][2] = 0.; 
  pts[i][0] = max; pts[i][1] = max; pts[i++][2] = 0.;
  
  for( i = 0; i < expectedNumPts; i++ )
    {
    mesh->SetPoint( i, pts[i] );
    }
  
  CellAutoPointer cellpointer;
  QEPolygonCellType *poly;
  
  for( i = 0; i < expectedNumCells; i++ )
    {
    poly = new QEPolygonCellType( 3 );
    cellpointer.TakeOwnership( poly );
    cellpointer->SetPointId( 0, simpleTriangleCells[i*3] );
    cellpointer->SetPointId( 1, simpleTriangleCells[i*3+1] );
    cellpointer->SetPointId( 2, simpleTriangleCells[i*3+2] );
    mesh->SetCell( i, cellpointer );
    }  
  
  CheckONextLink<MeshType>( mesh );
}
//--------------------------------------------------------------------------------


//--------------------------------------------------------------------------------
template<
class TMeshType
>
typename TMeshType::Pointer
RecursiveFlipEdgeTest( typename TMeshType::Pointer mesh, typename TMeshType::PointIdentifier point, typename TMeshType::CellIdentifier cell )
{
  typedef          TMeshType                              MeshType;
  typedef typename MeshType::Pointer                      MeshTypePointer;
  typedef typename MeshType::PointType                    MeshPointType;
  typedef typename MeshType::CellType                     MeshCellType;
  typedef typename MeshType::PointIdentifier              MeshPointIdentifier;
  typedef typename MeshType::CellIdentifier               MeshCellIdentifier;
  typedef typename MeshType::PointsContainer              MeshPointsContainer;
  typedef typename MeshType::CellsContainer               MeshCellsContainer;
  typedef typename MeshType::PointIdList                  MeshPointIdList;
  typedef typename MeshType::QEType                       MeshQuadEdgeType;
  typedef typename MeshType::QEPrimal                     MeshQuadEdgePrimal;

  typedef typename MeshType::CellsContainerIterator       MeshCellsContainerIteratorType;
  typedef typename MeshType::PointsContainerConstIterator MeshPointsContainerConstIterator;
  
  typedef typename MeshCellType::PointIdConstIterator  MeshCellPointIdConstIterator;
  typedef typename MeshCellType::PointIdIterator       MeshCellPointIdIterator;
  typedef typename MeshCellType::CellAutoPointer       MeshCellCellAutoPointer;

  typedef typename MeshQuadEdgeType::DualOriginRefType DualOriginRefType;
  
  typedef typename itk::QuadEdgeMeshPolygonCell< MeshCellType > QEPolygonCellType;
  typedef typename itk::QuadEdgeMeshEulerOperatorFlipEdgeFunction< MeshType, MeshQuadEdgeType > FlipEdgeFunctionType;

  typedef typename itk::VTKPolyDataWriter< MeshType > MeshWriterType;

  MeshCellCellAutoPointer cellpointer;
  MeshPointType pointCoord;
  
  MeshPointIdentifier q[3]; // points index of T
  MeshPointIdentifier p[3]; // points index of aT
  int r(0), k(0); // current point and oposite point position  
  
  if( mesh->GetCell( cell, cellpointer ) )
    {
    if( mesh->GetPoint( point, &pointCoord ) )
      {
      
      mesh->GetCell( cell, cellpointer );
      mesh->GetPoint( point, &pointCoord );        
      
      MeshCellPointIdConstIterator cellPointsIterator = cellpointer->PointIdsBegin();
      int i(0);
      while( cellPointsIterator != cellpointer->PointIdsEnd() )
        {
        p[i] = *cellPointsIterator;
        if( p[i] == point)
          {
          r=i;
          }
        ++cellPointsIterator;
        i++;
        }
      
      std::cout << std::endl;
      std::cout <<"we look at the cell " << cell << " made of the points (" << p[0] << "," << p[1] << "," << p[2] << ")";
      std::cout << std::endl;
      std::cout << "our current point is the "<<r+1<<"rd point of the triangle, which is " << p[r] << std::endl;
      std::cout << "the edge PiPj is "<< p[(r+1)%3] << "-" << p[(r+2)%3] << std::endl;
      
      MeshQuadEdgeType* e = mesh->FindEdge(p[(r+1)%3],p[(r+2)%3]);
      if( !e )
        {
        std::cerr << "ERROR: ";
        std::cerr << "The edge " << p[(r+1)%3] << "-" << p[(r+2)%3];
        std::cerr << " was not found which is impossible. Use this as a break holder for debugging.";
        std::cerr << std::endl;

        // NOTE ALEX: temporary hack: try to find the edge with reverse point order
        e = mesh->FindEdge(p[(r+2)%3],p[(r+1)%3]);
        if( !e )
          {
          std::cerr << "ERROR: ";
          std::cerr << "This edge is really not there, OMG, what shall we do now?!?!" << std::endl;
          return mesh;
          }
        }

      if( !e->IsAtBorder() )
        {
        DualOriginRefType adjCell = e->GetRight();
        if( adjCell == cell )
          {
          adjCell = e->GetLeft();
          }
        if( mesh->GetCell( adjCell, cellpointer ) )
          {

          mesh->GetCell( adjCell, cellpointer );
        
          MeshCellPointIdConstIterator adjCellPointsIterator = cellpointer->PointIdsBegin();
          int j(0);
          while( adjCellPointsIterator != cellpointer->PointIdsEnd() )
            {
            q[j] = *adjCellPointsIterator;
            if( q[j] != p[(r+1)%3] && q[j] != p[(r+2)%3] )
              {
              k = j;
              }
            ++adjCellPointsIterator;
            j++;
            }
          std::cout << "the neighbour cell is " << adjCell;
	  std::cout << " made of the points (" << q[0] <<","<< q[1] <<","<< q[2] << ")" << std::endl;
          std::cout << "the oposite point of our current point "<<p[r]<<" is the point "<< q[k] << std::endl;
            
          if( TestPointInTriangleInMesh< MeshType >( mesh, adjCell, pointCoord, true ) ) 
            {
                
	    typename MeshWriterType::Pointer writer = MeshWriterType::New();
	    writer->SetFileName( "beforeflip.vtk" );
	    writer->SetInput( mesh );
	    writer->Update();
		    
            // NOTE ALEX: use itkQuadEdgeMeshEulerOperatorFlipEdge

	    std::cout << std::endl;
	    std::cout << "I am going to flip the edge" << std::endl;
	    FlipEdgeFunctionType* flipedge = FlipEdgeFunctionType::New();
	    flipedge->SetInput( mesh );
	    e = flipedge->Evaluate( e );
	    std::cout << "I fliped the edge and I liked it" << std::endl;
            
            if( !e )
              {
	      std::cerr << "ERROR - It was just a dream" << std::endl;
              }
	   
	    std::cout << "we write a mesh"<< std::endl; 
            writer = MeshWriterType::New();
	    writer->SetFileName( "afterflip.vtk" );
	    writer->SetInput( mesh );
	    writer->Update();
	    std::cout << "mesh is writen"<< std::endl;
	
            CheckONextLink<MeshType>( mesh );
            
	    // NOTE STEF: need to retrive new cells ids for recursivity
	    //MeshCellIdentifier newCellIds[2];	 
            //mesh = RecursiveFlipEdgeTest< MeshType >( mesh, point, newCellIds[0] );
            //mesh = RecursiveFlipEdgeTest< MeshType >( mesh, point, newCellIds[1] );
              
            }
          else
            {
            std::cout<< "criterion respected, no need to flip\n\n";
            CheckONextLink<MeshType>( mesh );
            return mesh;
            }
          }
        else
          {
          std::cerr << "Error - Could not find the adjacent cell\n";
          return mesh; 
          }
        }
      else 
        {
        std::cout << "Border edge, can not be flip\n\n";
        CheckONextLink<MeshType>( mesh );
        return mesh;
        }
      }
    else
      {
      std::cerr << "Error - Could not find the point given in parameter\n";
      return mesh;
      }
    }
  else // end of if .... 
    {
    std::cerr << "Error - Could not find the cell given in parameter\n";
    return mesh;
    }

  CheckONextLink<MeshType>( mesh );
  return mesh;
}
//--------------------------------------------------------------------------------


//--------------------------------------------------------------------------------
// Delaunay triangulation process loop
//
template<
class TPointSetType,
class TMeshType
>
typename TMeshType::Pointer
DelaunayTriangulation( typename TPointSetType::Pointer myPointSet )
{
  typedef          TPointSetType                               PointSetType;
  typedef typename PointSetType::PointType                     PointSetPointType;
  typedef typename PointSetType::PointIdentifier               PointSetPointIdentifier;
  typedef typename PointSetType::PointsContainer               PointSetPointsContainer;
  typedef typename PointSetType::PointsContainerConstIterator  PointSetPointsContainerConstIterator;
  
  typedef          TMeshType                              MeshType;
  typedef typename MeshType::Pointer                      MeshTypePointer;
  typedef typename MeshType::PointType                    MeshPointType;
  typedef typename MeshType::CellType                     MeshCellType;
  typedef typename MeshType::PointIdentifier              MeshPointIdentifier;
  typedef typename MeshType::CellIdentifier               MeshCellIdentifier;
  typedef typename MeshType::PointsContainer              MeshPointsContainer;
  typedef typename MeshType::CellsContainer               MeshCellsContainer;
  typedef typename MeshType::PointIdList                  MeshPointIdList;
  typedef typename MeshType::QEType                       MeshQuadEdgeType;
  typedef typename MeshType::CellsContainerIterator       MeshCellsContainerIteratorType;
  typedef typename MeshType::PointsContainerConstIterator MeshPointsContainerConstIterator;
  
  typedef itk::VTKPolyDataWriter< MeshType > MeshWriterType;

  typedef typename MeshCellType::PointIdConstIterator  MeshCellPointIdConstIterator;
  typedef typename MeshCellType::PointIdIterator       MeshCellPointIdIterator;
  typedef typename MeshCellType::CellAutoPointer       MeshCellCellAutoPointer;
  
  typedef typename itk::WalkInTriangulationFunction< MeshType > WalkInTriangulation;
  typedef typename WalkInTriangulation::Pointer                 WalkInTriangulationPointer;
  typedef          itk::VectorContainer< unsigned int, int >    MeshCellIdVectorContainer;
  typedef typename itk::QuadEdgeMeshPolygonCell< MeshCellType > QEPolygonCellType;

  // Initialisation
  
  MeshTypePointer myMesh = MeshType::New();
  CreateDummyMesh< MeshType >( myMesh, 5000 );
  
  std::cout << std::endl << "No. of Cells  : " << myMesh->GetNumberOfCells();
  std::cout << std::endl << "No. of Edges  : " << myMesh->GetNumberOfEdges();
  std::cout << std::endl << "No. of Faces  : " << myMesh->GetNumberOfFaces();
  std::cout << std::endl << "No. of Points : " << myMesh->GetNumberOfPoints();
  std::cout << std::endl << std::endl;
  
  CheckONextLink<MeshType>( myMesh );
  
  PointSetPointsContainer              *myPoints           = myPointSet->GetPoints();
  PointSetPointsContainerConstIterator pointIterator       = myPoints->Begin();
  MeshCellIdentifier                   myStartingCellIndex = 0;
  
  int o(1);
  while( pointIterator != myPoints->End() ) 
    { 
    std::cout << "/---------------------------------------------------/" << std::endl;
    std::cout << "Delaunay iteration : " << o << "/" << myPoints->Size();
    std::cout << " - pts " << o+3 << " (" << pointIterator.Value()[0] << "," << pointIterator.Value()[1] << ")" << std::endl;
      
    PointSetPointType       myTempPoint;
    MeshPointType           myPoint;
    MeshPointIdentifier     myPointIndex;
    MeshCellIdentifier      myCellIndex;
    MeshCellCellAutoPointer myCellPointer;
    
    // Get the current point
    myTempPoint  = pointIterator.Value();
    myPoint[0]   = myTempPoint[0];
    myPoint[1]   = myTempPoint[1];
    myPoint[2]   = myTempPoint[2];
    
    // Find the triangle that include the point
    std::cout << "I walk the mesh -  starting cell " << myStartingCellIndex << std::endl;
    WalkInTriangulationPointer         letsWalk       = WalkInTriangulation::New();
    MeshCellIdVectorContainer::Pointer walkCellIdList = MeshCellIdVectorContainer::New();
    try 
      {
      walkCellIdList = letsWalk->Evaluate( myMesh, myPoint, myStartingCellIndex );
      }
    catch( int e ) 
      {
      if( e == -1 )
        { 
        std::cerr << "Error Walk - Point is out of the mesh" << std::endl;
        throw -1;
        }
      else if( e == -2 )
        {
        std::cerr << "Error Walk - Starting cell does not exist" << std::endl;
        throw -1;
        }
      else
        {
        std::cerr << "Error Walk - Unknown error" << std::endl;
        throw -2;
        }    
      break;
      }
    myCellIndex = ( --walkCellIdList->End() )->Value();
    std::cout << "I walked the mesh: " << myCellIndex << std::endl;
    
    if( myMesh->GetCell( myCellIndex, myCellPointer ) )
      { 
      MeshCellPointIdIterator pointIdIterator;
      std::vector< MeshCellIdentifier > cellPointsIds( 3 );
      std::vector< MeshCellIdentifier > newCellIds( 3 );
      MeshCellCellAutoPointer cellpointer;
      QEPolygonCellType *poly;
      
      myMesh->GetCell( myCellIndex, myCellPointer );
      
      pointIdIterator = myCellPointer->PointIdsBegin();
      cellPointsIds[0] = *pointIdIterator;
      pointIdIterator++;
      cellPointsIds[1] = *pointIdIterator;
      pointIdIterator++;
      cellPointsIds[2] = *pointIdIterator;
        
      myMesh->DeleteFace( myCellIndex );     

      myPointIndex = myMesh->FindFirstUnusedPointIndex();
      myMesh->SetPoint( myPointIndex, myPoint );

      // Create 3 new anticlockwise oriented triangles
      for( unsigned int i = 0; i < 3; i++ )
        {
        newCellIds[i] = myMesh->FindFirstUnusedCellIndex();
        poly = new QEPolygonCellType( 3 );
        cellpointer.TakeOwnership( poly );
        cellpointer->SetPointId( 0, cellPointsIds[ (i)   % 3 ] );
        cellpointer->SetPointId( 1, cellPointsIds[ (i+1) % 3 ] );
        cellpointer->SetPointId( 2, myPointIndex );
        myMesh->SetCell( newCellIds[i], cellpointer );
        } 
        
      std::cout << "Cell id " << myCellIndex <<" ( "<< cellPointsIds[0] << "," << cellPointsIds[1] << "," << cellPointsIds[2]; 
      std::cout << ") replaced by" << std::endl;
      std::cout << "\t" << newCellIds[0] << " (" << cellPointsIds[0] << "," << cellPointsIds[1] << "," << myPointIndex << ")" << std::endl;
      std::cout << "\t" << newCellIds[1] << " (" << cellPointsIds[1] << "," << cellPointsIds[2] << "," << myPointIndex << ")" << std::endl;
      std::cout << "\t" << newCellIds[2] << " (" << cellPointsIds[2] << "," << cellPointsIds[0] << "," << myPointIndex << ")" << std::endl; 
        
      myMesh = RecursiveFlipEdgeTest< MeshType >( myMesh, myPointIndex, newCellIds[0] );
      myMesh = RecursiveFlipEdgeTest< MeshType >( myMesh, myPointIndex, newCellIds[1] );
      myMesh = RecursiveFlipEdgeTest< MeshType >( myMesh, myPointIndex, newCellIds[2] );
      }
      
    std::cout << "\nNo. of Cells : " << myMesh->GetNumberOfCells()
              << "\nNo. of Edges : " << myMesh->GetNumberOfEdges()
              << "\nNo. of Faces : " << myMesh->GetNumberOfFaces()
              << "\nNo. of Points : " << myMesh->GetNumberOfPoints() <<"\n\n";
      
    myPoint = myMesh->GetPoint( myPointIndex );
    MeshQuadEdgeType* myEdge = myPoint.GetEdge();
    myStartingCellIndex = myEdge->GetLeft(); 

    if( myMesh->GetCell( myEdge->GetLeft(), myCellPointer ) )
      {
      myStartingCellIndex = myEdge->GetLeft(); 
      }
    else if( myMesh->GetCell( myEdge->GetRight(), myCellPointer ) )
      {
      myStartingCellIndex = myEdge->GetLeft();  
      }
    else 
      {
      std::cerr << "Error - Could not find a correct starting cell\n";
      break;
      }

    // check consistency    
    CheckONextLink< MeshType >( myMesh );

    // write down mesh
    std::string tempname = "./tempMesh";
    std::string ext = ".vtk";
    std::stringstream ss;
    ss << o;
    tempname = tempname + ss.str() + ext;

    typename MeshWriterType::Pointer write = MeshWriterType::New();
    write->SetFileName( tempname );
    write->SetInput( myMesh );
    write->Update();

    o++;
    ++pointIterator;
    }    
  
  CheckONextLink<MeshType>( myMesh );
  
  return myMesh;
}
//--------------------------------------------------------------------------------


//--------------------------------------------------------------------------------
// Main function
//
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
  
  const unsigned int Dimension = 3;
  typedef double     PixelType;

  typedef itk::QuadEdgeMesh< PixelType, Dimension > MeshType; 
  typedef itk::PointSet< PixelType, Dimension >     PointSetType;
  typedef PointSetType::PointType                   PointType;

  // -------------------------------------------------- //
  // Dummy QuadEdgeMesh test
  MeshType::Pointer myDummyMesh = MeshType::New();
  CreateDummyMesh< MeshType >( myDummyMesh, 5000 );

  // -------------------------------------------------- //
  // Delaunay Loop Test

  int type = atoi( argv[1] );
  int meshSize = atoi( argv[2] );
  int expectedNumPts = 0;
  std::vector<PointType> pts;

  switch(type) 
    {
    case 1 :     
      pts = GeneratePointCoordinates< PointSetType >( meshSize ); 
      expectedNumPts = meshSize*meshSize;
      break;
    case 2 :
      pts = GenerateRandomPointCoordinates< PointSetType >( meshSize ); 
      expectedNumPts = meshSize;
      break;
    default:
      std::cerr << "Exit wrong arguments\n";
      std::cerr << "Usage - arg[1] [1|2] = [\"reg\"|\"rand\"] - arg[2] int = [rowsize|nbPoint]";
      std::cerr << std::endl;
      return EXIT_FAILURE;
    }
  
  PointSetType::Pointer myPointSet = PointSetType::New();
  for( int i = 0; i < expectedNumPts ; i++ )
    {
    myPointSet->SetPoint( i, pts[i] );
    }
  
  MeshType::Pointer myTriangulatedMesh = MeshType::New();
  try
    {
    myTriangulatedMesh = DelaunayTriangulation< PointSetType, MeshType >( myPointSet );
    }
  catch( int e ) 
    {
    std::cerr << "Exception WiT Caught" << std::endl;
    return EXIT_FAILURE;
    }

  CheckONextLink< MeshType >( myTriangulatedMesh );

  typedef itk::VTKPolyDataWriter< MeshType > MeshWriter;
  MeshWriter::Pointer write = MeshWriter::New();
  write->SetFileName("./tempMesh.vtk");
  write->SetInput( myTriangulatedMesh );
  write->Update();

  return EXIT_SUCCESS;
}
//--------------------------------------------------------------------------------
