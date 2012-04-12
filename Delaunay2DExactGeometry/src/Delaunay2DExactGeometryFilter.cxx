
#include "itkPointSet.h"
#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshTraits.h"
#include "itkQuadEdgeMeshPolygonCell.h"
#include "itkVectorContainer.h"

#include "itkWalkInTriangulationFunction.h"
#include "itkPointInCircleGeometricalPredicateFunctor.h"
#include "itkQuadEdgeMeshEulerOperatorsTestHelper.h"

#include <iostream>
#include <limits>

//TEMPORARY
#include "itkVTKPolyDataWriter.h"

template< class TMeshType >
void CheckONextLink( typename TMeshType::Pointer mesh )
{
  typedef          TMeshType                         MeshType;
  typedef typename MeshType::PointType               PointType;
  typedef typename MeshType::PointIdentifier         PointIndex;
  typedef typename MeshType::PointsContainer         PointsContainer;
  typedef typename MeshType::PointsContainerIterator PointsContainerIterator;
  typedef typename MeshType::QEType                  QuadEdgeType;

  typedef typename itk::VTKPolyDataWriter< MeshType > MeshWriterType;
  typedef typename MeshWriterType::Pointer                MeshWriterPointer;

  PointsContainer*        myPoints;
  PointsContainerIterator myPointIterator;

  myPoints = mesh->GetPoints();
  myPointIterator = myPoints->Begin();

  while( myPointIterator != myPoints->End() )
    {
    PointType p = myPointIterator.Value();
    QuadEdgeType* e = p.GetEdge();
    QuadEdgeType* next = e->GetOnext();
    while( e != next )
      {
      next = next->GetOnext();
      } 
    ++myPointIterator;
    }
  
  std::cout<<"Check O ring done\n";
}

// Create a dummy mesh for triangulation initialisation
// Mesh points should be at inifinity
//
template< class TMeshType >
void CreateDummyMesh( typename TMeshType::Pointer mesh, typename TMeshType::PixelType limit )
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
 
  PixelType min = -limit; // Change to infinity if possible 
  PixelType max =  limit; // Change to infinity if possible
 
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
}

template< class TMeshType >
bool
RecursiveFlipEdgeEvaluation( typename TMeshType::Pointer mesh, 
                             typename TMeshType::PointIdentifier myPointIndex, 
                             typename TMeshType::CellIdentifier myCellIndex 
                           )
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
  typedef typename MeshType::CellsContainerIterator       MeshCellsContainerIteratorType;
  typedef typename MeshType::PointsContainerConstIterator MeshPointsContainerConstIterator;

  typedef typename MeshCellType::PointIdConstIterator  MeshCellPointIdConstIterator;
  typedef typename MeshCellType::PointIdIterator       MeshCellPointIdIterator;
  typedef typename MeshCellType::CellAutoPointer       MeshCellCellAutoPointer;

  typedef typename MeshQuadEdgeType::DualOriginRefType MeshDualOriginRefType; 

  typedef typename itk::QuadEdgeMeshPolygonCell< MeshCellType > QEPolygonCellType;

  MeshPointIdentifier     pA,pB,pC;
  MeshPointIdentifier     p[3];
  MeshPointType           myPointCoord;
  MeshDualOriginRefType   adjacentCellIndex;
  MeshCellCellAutoPointer myCellPointer;
  MeshCellPointIdIterator pointIdIterator;
  bool flip = false;

  if( mesh->GetCell( myCellIndex, myCellPointer ) ) // Check cellindex argument
    {
    //mesh->GetCell( myCellIndex, myCellPointer );
    pointIdIterator = myCellPointer->PointIdsBegin();
    p[0] = *(pointIdIterator);
    p[1] = *(++pointIdIterator);
    p[2] = *(++pointIdIterator);
    for( int i = 0; i < 3; i++ )
      {
      if( p[i%3] == myPointIndex )
        {
        pA = p[(i+1)%3];
        pB = p[(i+2)%3];
        }
      }
    std::cout << "Triangle T : "<< myCellIndex << " (" << myPointIndex << "," << pA << "," << pB <<")\n";
    
    if( mesh->FindEdge( pA, pB ) ) // Check the selected edge consistancy
      {
      if(!mesh->FindEdge( pA, pB )->IsAtBorder() ) // Check if flipable
        {
        
        MeshDualOriginRefType adjT1, adjT2;
        
        adjT1 = mesh->FindEdge( pA, pB )->GetRight();
        adjT2 = mesh->FindEdge( pA, pB )->GetLeft();        
        if( !mesh->GetCell( adjT1, myCellPointer) && !mesh->GetCell( adjT2, myCellPointer) )
          {
          std::cerr << "Right = "<<adjT1<<"\t";
          std::cerr << "Left = "<<adjT2<<"\n";
          std::cerr << "HUSTON WE HAVE A PROBLEM\n";
          return EXIT_FAILURE;
          }
        
        adjacentCellIndex = mesh->FindEdge( pA, pB )->GetRight();
        std::cout << "Right adjacency : "<< adjacentCellIndex <<"\n";
        if( adjacentCellIndex == myCellIndex )
          {
          std::cout << "It is our current Cell : "<< adjacentCellIndex <<" = "<< myCellIndex << "\n";
          adjacentCellIndex = mesh->FindEdge( pA, pB )->GetLeft();
          std::cout << "Left adjacency : "<< adjacentCellIndex <<"\n";
          }
        myPointCoord = mesh->GetPoint( myPointIndex );
        
        if( mesh->GetCell( adjacentCellIndex, myCellPointer) ) // Check the neighbour cell consistancy
          { 
          if( TestPointInTriangleInMesh< MeshType >( mesh, adjacentCellIndex, myPointCoord, true ) ) // Check Delaunay criterion
            {
            
            //mesh->GetCell( adjacentCellIndex, myCellPointer );
            pointIdIterator = myCellPointer->PointIdsBegin();
            while( pointIdIterator != myCellPointer->PointIdsEnd() )
              {
              if( pA != *pointIdIterator && pB != *pointIdIterator )
                {
                pC = *pointIdIterator;
                }
              ++pointIdIterator;
              }
            
            std::cout << "Triangle aT : "<< adjacentCellIndex << " (" << pA << "," << pC << "," << pB <<")\n";
            
            mesh->DeleteFace( myCellIndex );
            mesh->DeleteFace( adjacentCellIndex );
            
            if( mesh->GetCell( adjacentCellIndex, myCellPointer ) ) // Check cell suppression
              {
              std::cerr << "WTF is adjacentCell is still existing O_o!!\n";
              return flip;
              }
            if( mesh->GetCell( myCellIndex, myCellPointer ) ) // Check cell suppression
              {
              std::cerr << "WTF is myCell is still existing O_o!!\n";
              return flip;
              }
            
            int simpleTriangleCells[6] = 
            { myPointIndex, pA, pC, 
              myPointIndex, pC, pB };
            
            QEPolygonCellType *poly;
            MeshCellCellAutoPointer cellpointer;
            MeshCellIdentifier      cellIndexTab[2];
            for( unsigned int i = 0; i<2; i++ )
              {
              cellIndexTab[i] = mesh->FindFirstUnusedCellIndex();
              poly = new QEPolygonCellType( 3 );
              cellpointer.TakeOwnership( poly );
              cellpointer->SetPointId( 0, simpleTriangleCells[i*3] );
              cellpointer->SetPointId( 1, simpleTriangleCells[i*3+1] );
              cellpointer->SetPointId( 2, simpleTriangleCells[i*3+2] );
              mesh->SetCell( cellIndexTab[i], cellpointer ); 
              std::cout << "we create a new cell " << cellIndexTab[i] << " ("; 
              std::cout << simpleTriangleCells[i*3] << ","; 
              std::cout << simpleTriangleCells[i*3+1] << ","; 
              std::cout << simpleTriangleCells[i*3+2] << ")\n";
              }  
            
            // Check new cell edge consistancy
            if( mesh->FindEdge(pC, myPointIndex) )
              std::cout << "find edge test " << pC << "-" << myPointIndex <<" OK\n";
            if( mesh->FindEdge(myPointIndex, pC) )
              std::cout << "find edge test " << myPointIndex << "-" << pC <<" OK\n";
            if( mesh->FindEdge(pA, pC) )
              std::cout << "find edge test " << pA << "-" << pC <<" OK\n";
            if( mesh->FindEdge(pC, pA) )
              std::cout << "find edge test " << pC << "-" << pA <<" OK\n";
            if( mesh->FindEdge(myPointIndex, pA) )
              std::cout << "find edge test " << myPointIndex << "-" << pA <<" OK\n";
            if( mesh->FindEdge(pA, myPointIndex) )          
              std::cout << "find edge test " << pA << "-" << myPointIndex <<" OK\n";
            if( mesh->FindEdge(myPointIndex, pC) )
              std::cout << "find edge test " << myPointIndex << "-" << pC <<" OK\n";
            if( mesh->FindEdge(pC, myPointIndex) )
              std::cout << "find edge test " << pC << "-" << myPointIndex <<" OK\n";
            if( mesh->FindEdge(pC, pB) )
              std::cout << "find edge test " << pC << "-" << pB <<" OK\n";
            if( mesh->FindEdge(pB, pC) )
              std::cout << "find edge test " << pB << "-" << pC <<" OK\n";
            if( mesh->FindEdge(myPointIndex, pC) )
              std::cout << "find edge test " << myPointIndex << "-" << pC <<" OK\n";
            if( mesh->FindEdge(pC, myPointIndex) )
              std::cout << "find edge test " << pC << "-" << myPointIndex <<" OK\n";
            
            std::cout << "\tedge (" << pA << "-" << pB << ") is flip, we delete cells " << myCellIndex << " and " << adjacentCellIndex << "\n"
            << "\t\t new cell : " << cellIndexTab[0] << " (" << myPointIndex << "," << pA << "," << pC << ")\n"
            << "\t\t new cell : " << cellIndexTab[1] << " (" << myPointIndex << "," << pC << "," << pB << ")\n";
            
            // VERIFICATION & DEBUG  /////////////////////////////////////////////
            CheckONextLink< MeshType >( mesh );
            // VERIFICATION & DEBUG /////////////////////////////////////////////
            
            //RecursiveFlipEdgeEvaluation< MeshType >( mesh, myPointIndex, cellIndexTab[0] ); 
            //RecursiveFlipEdgeEvaluation< MeshType >( mesh, myPointIndex, cellIndexTab[1] );
            
            // VERIFICATION & DEBUG  /////////////////////////////////////////////
            CheckONextLink< MeshType >( mesh ); 
            // VERIFICATION & DEBUG /////////////////////////////////////////////
            
            flip = true;
            }
          else 
            {
            std::cout << "Not inside the circle, no need to flip the edge\n";
            return flip;
            }
          }
        else 
          {
          std::cerr << "Error mesh consistance - Adjacent cell does not exist while it should\n";
          std::cerr << "Triangle aT : "<< adjacentCellIndex << " (through the edge " << pA << "-" << pB <<")\n";
          return flip;
          }
        }
      else 
        {
        std::cout << "It is a border edge, we can not flip it\n";
        return flip;
        }
      }
    else
      {
      std::cerr << "Error mesh consistance - Edge does not exist while it should exist\n";
      std::cerr << "Edge " << pA << "-" << pB <<" has a problem\n";
      return flip;
      }
    }    
  else 
    {
    std::cerr << "Error - Cell given in parameter does not exist\n";
    return flip; 
    }

  return flip;
}

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
 
  std::cout << "\nNo. of Cells : " << myMesh->GetNumberOfCells()
            << "\nNo. of Edges : " << myMesh->GetNumberOfEdges()
            << "\nNo. of Faces : " << myMesh->GetNumberOfFaces()
            << "\nNo. of Points : " << myMesh->GetNumberOfPoints() <<"\n\n";
 
  PointSetPointsContainer              *myPoints           = myPointSet->GetPoints();
  PointSetPointsContainerConstIterator pointIterator       = myPoints->Begin();
  MeshCellIdentifier                   myStartingCellIndex = 0;
 
  // Process Loop
  int o(1);
  while( pointIterator != myPoints->End() ) 
   {
   std::cout<<"/---------------------------------------------------/\n";
   std::cout<<"Delaunay iteration : "<<o<<"/"<<myPoints->Size();
   std::cout<<" - pts "<<o+3<<" ("<<pointIterator.Value()[0]<<","<<pointIterator.Value()[1]<<")"<<std::endl;
  
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
   std::cout<<"I walk the mesh -  starting cell "<<myStartingCellIndex<<"\n";
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
       std::cerr<< "Error Walk - Point is out of the mesh\n";
       throw -1;
       }
     else if( e == -2 )
       {
       std::cerr<< "Error Walk - Starting cell does not exist\n";
       throw -1;
       }
     else
       {
       std::cerr<< "Error Walk - Unknown error\n";
       throw -2;
       }   
     break;
     }
   myCellIndex = ( --walkCellIdList->End() )->Value();
   std::cout<<"I walked the mesh\n";
   
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
   
     // Delete the triangle
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
     std::cout << ") replaced by \n";
     std::cout << "\t" << newCellIds[0] << " (" << cellPointsIds[0] << "," << cellPointsIds[1] << "," << myPointIndex << ")\n";
     std::cout << "\t" << newCellIds[1] << " (" << cellPointsIds[1] << "," << cellPointsIds[2] << "," << myPointIndex << ")\n";
     std::cout << "\t" << newCellIds[2] << " (" << cellPointsIds[2] << "," << cellPointsIds[0] << "," << myPointIndex << ")\n";
   
     // Delaunay criterion recursive test       
     RecursiveFlipEdgeEvaluation< MeshType >( myMesh, myPointIndex, newCellIds[0] );
     RecursiveFlipEdgeEvaluation< MeshType >( myMesh, myPointIndex, newCellIds[1] );
     RecursiveFlipEdgeEvaluation< MeshType >( myMesh, myPointIndex, newCellIds[2] );
   
     // VERIFICATION & DEBUG  /////////////////////////////////////////////
     CheckONextLink< MeshType >( myMesh );
     
     typedef itk::QuadEdgeMesh< double, 3 > TempMeshType;
     typedef typename itk::VTKPolyDataWriter< TempMeshType > MeshWriter;
     MeshWriter::Pointer write = MeshWriter::New();
     write->SetFileName("./tempMesh.vtk");
     write->SetInput( myMesh );
     write->Update();
     // VERIFICATION & DEBUG /////////////////////////////////////////////
       
       
     //myStartingCellIndex = newCellIds[ rand()%3 ]; 
     
     MeshQuadEdgeType* e = myMesh->FindEdge( myPointIndex );
     myStartingCellIndex = e->GetRight();
     std::cout << "We take right starting cell : " << myStartingCellIndex << "\n";
     if( !myMesh->GetCell( myStartingCellIndex, cellpointer ) )
       {
       std::cout << "\tRight cell does not exist, we take left one\n";
       myStartingCellIndex = e->GetLeft();
       std::cout << "We take left starting cell : " << myStartingCellIndex << "\n";
       }
     else if( !myMesh->GetCell( myStartingCellIndex, cellpointer ) )
       {
       std::cerr<<"Error DelaunayLoop - Could not find any starting cell\n";
       return myMesh;
       } 
     }
   ++pointIterator;
   o++;
   std::cout << "\nNo. of Cells : " << myMesh->GetNumberOfCells()
             << "\nNo. of Edges : " << myMesh->GetNumberOfEdges()
             << "\nNo. of Faces : " << myMesh->GetNumberOfFaces()
             << "\nNo. of Points : " << myMesh->GetNumberOfPoints() <<"\n\n";
  
   // Detele dummy points    
   } 
 return myMesh;
}

template< class TMesh >
std::vector< typename TMesh::PointType >
GenerateRandomPointCoordinates( const unsigned int& iN )
{
  typedef typename TMesh::PointType        PointType;
  typedef typename PointType::CoordRepType CoordRepType;
  std::vector< PointType > oPt( iN * iN );
  
  srand(time(NULL));
  for( unsigned int i = 0; i < iN; i++ )
  {
   oPt[ i ][0] = static_cast< CoordRepType >( rand() % 10000 - 5000 );
   oPt[ i ][1] = static_cast< CoordRepType >( rand() % 10000 - 5000 );
   oPt[ i ][2] = static_cast< CoordRepType >( 0. );
  }
  return oPt;
}


// Main function
//
int
main( int argc, char* argv[] )
{
  const unsigned int Dimension = 3;
  typedef double     PixelType;
 
  typedef itk::QuadEdgeMesh< PixelType, Dimension > MeshType; 
  typedef itk::PointSet< PixelType, Dimension >     PointSetType;
  typedef PointSetType::PointType                   PointType;
 
  if( argc != 3 )
    {
    std::cout << "Usage - exe reg rowsize | exe rand nbPoint\n";
    return EXIT_FAILURE;
    }

  // -------------------------------------------------- //
  // Dummy QuadEdgeMesh test
  MeshType::Pointer myDummyMesh = MeshType::New();
  CreateDummyMesh< MeshType >( myDummyMesh, 5000 );

  // -------------------------------------------------- //
  // Delaunay Loop Test

  std::string type = argv[1];

  if( type == "reg" )
    {
    int meshSize = atoi( argv[2] );
    int expectedNumPts = meshSize*meshSize;
   
    PointSetType::Pointer myPointSet = PointSetType::New();
    std::vector<PointType> pts = GeneratePointCoordinates< PointSetType >( meshSize ); 
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
      std::cerr << "Exception WiT Caught\n";
      return EXIT_FAILURE;
      }
    }
  else if( type == "rand" )
    {
    int nbPoint = atoi( argv[2] );
    PointSetType::Pointer myPointSet = PointSetType::New();
    std::vector<PointType> pts = GenerateRandomPointCoordinates< PointSetType >( nbPoint ); 
    for( int i = 0; i < nbPoint ; i++ )
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
      std::cerr << "Exception WiT Caught\n";
      return EXIT_FAILURE;
      }
    }
  else
    {
    std::cout << "usage - exe reg rowsize | exe rand nbPoint\n";
    return EXIT_FAILURE;
    }


  return EXIT_SUCCESS;
}
