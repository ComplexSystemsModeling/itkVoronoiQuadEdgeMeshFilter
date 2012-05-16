
#ifndef __itkPointSetToDelaunayTriangulationFilter_hxx
#define __itkPointSetToDelaunayTriangulationFilter_hxx

#include "itkPointSetToDelaunayTriangulationFilter.h"

namespace itk 
{

//--------------------------------------------------------------------------------
// Dummy point deletion function
//--------------------------------------------------------------------------------
template< class TInMesh, class TOutMesh = TInMesh >
void
PointSetToDelaunayTriangulationFilter::
DeleteDummyPoints( MeshType::Pointer mesh )
{
  for( unsigned int i = 0; i < 4; i++ )
    {
    QEType* e;
    while( e =  mesh->GetPoint( i ).GetEdge() )
      {
      mesh->LightWeightDeleteEdge( e );
      }
    mesh->DeletePoint( i );
    }  
}
//--------------------------------------------------------------------------------

//--------------------------------------------------------------------------------
// Dummy Mesh creation function
//--------------------------------------------------------------------------------
template< class TInMesh, class TOutMesh = TInMesh >
void
PointSetToDelaunayTriangulationFilter::
CreateDummyMesh( MeshType::Pointer mesh, 
                 PixelType         limit )
{
  if( mesh->GetNumberOfPoints() )
    {
    mesh->Clear();
    mesh->ClearFreePointAndCellIndexesLists();
    }
  
  int expectedNumPts = 4;
  int expectedNumCells = 2;
  int simpleTriangleCells[6] = { 0, 1, 3, 0, 3, 2 };
  int i( 0 );  
  std::vector< PointType > pts( 4 );
  
  PixelType min = -limit;  
  PixelType max =  limit;
  
  pts[i][0] = min; pts[i][1] = min; pts[i++][2] = 0.;
  pts[i][0] = max; pts[i][1] = min; pts[i++][2] = 0.; 
  pts[i][0] = min; pts[i][1] = max; pts[i++][2] = 0.; 
  pts[i][0] = max; pts[i][1] = max; pts[i++][2] = 0.;
  
  for( i = 0; i < expectedNumPts; i++ )
    {
    mesh->SetPoint( i, pts[i] );
    }
  
  for( i = 0; i < expectedNumCells; i++ )
    {
    mesh->AddFaceTriangle( simpleTriangleCells[ i*3     ], 
                           simpleTriangleCells[ i*3 + 1 ], 
                           simpleTriangleCells[ i*3 + 2 ] );
    } 
}
//--------------------------------------------------------------------------------

//--------------------------------------------------------------------------------
// FlipEdge criterion function
//--------------------------------------------------------------------------------
template< class TInMesh, class TOutMesh = TInMesh >
bool
PointSetToDelaunayTriangulationFilter::
RecursiveFlipEdgeTest( MeshType::Pointer mesh, 
                       PointIdentifier   point, 
                       FaceRefType       cell )
{
  CellAutoPointer      cellpointer;
  PointType            pointCoord;
  PointIdConstIterator cellPointsIterator;     
  PointIdentifier      q[3];                // points index of T
  PointIdentifier      p[3];                // points index of aT
  int r(0), k(0);                           // current point and oposite point position  
  
  if( !mesh->GetCell( cell.first, cellpointer ) )
    {
    std::cerr << "ERROR - Could not find the cell given in parameter." << std::endl;
    throw -1;
    }
  if( !mesh->GetPoint( point, &pointCoord ) )
    {
    std::cerr << "ERROR - Could not find the point given in parameter." << std::endl;
    throw -1;
    }
  
  mesh->GetCell( cell.first, cellpointer );
  mesh->GetPoint( point, &pointCoord );        
  
  cellPointsIterator = cellpointer->PointIdsBegin();
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
  
  // Border edge test
  QEType* e = mesh->FindEdge( p[(r+1)%3], p[(r+2)%3] );
  if( e->IsAtBorder() )
    {
    return false;
    }
  
  FaceRefType adjCell = e->GetRight();
  // NOTE STEF: Temporary check
  if( adjCell == cell )
    {
    adjCell = e->GetLeft();
    }
  mesh->GetCell( adjCell.first, cellpointer );       
  
  cellPointsIterator = cellpointer->PointIdsBegin();
  int j(0);
  while( cellPointsIterator != cellpointer->PointIdsEnd() )
    {
    q[j] = *cellPointsIterator;
    if( q[j] != p[(r+1)%3] && q[j] != p[(r+2)%3] )
      {
      k = j;
      }
    ++cellPointsIterator;
    j++;
    }
  
  // PointInCircle test
  // If positive, flip edge and recursively check the new created cell
  //
  // NOTE STEF: the boolean parameter of the test should be remove in a near future  
  if( TestPointInTriangleInMesh< MeshType >( mesh, adjCell.first, pointCoord, true ) ) 
    {
    FlipEdgeFunction::Pointer flipedge = FlipEdgeFunction::New();
    flipedge->SetInput( mesh );
    
    e = flipedge->Evaluate( mesh->FindEdge( p[(r+1)%3], p[(r+2)%3] ) );
    
    RecursiveFlipEdgeTest( mesh, point, e->GetLeft() );
    RecursiveFlipEdgeTest( mesh, point, e->GetRight() );
    }
  
  return true;
  
}
//--------------------------------------------------------------------------------

//--------------------------------------------------------------------------------
// Add a point to the triangulation function
//--------------------------------------------------------------------------------
template< class TInMesh, class TOutMesh = TInMesh >
PointIdentifier
PointSetToDelaunayTriangulationFilter::
AddPoint( MeshType::Pointer mesh, 
          PointType         point, 
          FaceRefType       startingCell )
{
  // Walk the mesh to find the cell containing the new point to add	
  WalkInTriangulationFunction::Pointer walk       = WalkInTriangulationFunction::New();
  CellIdVectorContainerType::Pointer   cellIdList = CellIdVectorContainerType::New();
  try 
    {
    cellIdList = walk->Evaluate( mesh, point, &startingCell );
    }
  catch( int e ) 
    {
    std::cerr << "ERROR - Exception caught in the WalkInTriangulation process" << std::endl;
    throw -1;
    }
  // NOTE STEF: Use a cleaner way to get the last value of a VectorContainer
  FaceRefType cellIndex = ( --cellIdList->End() )->Value();
  
  CellAutoPointer                cellPointer;
  PointIdIterator                pointIdIterator;
  std::vector< PointIdentifier > cellPointsIds( 3 );
  std::vector< FaceRefType >     newCellIds( 3 );
  PointIdentifier                pointIndex;
  
  // Split the cell in 3 new cell 
  if( mesh->GetCell( cellIndex.first, cellPointer ) )
    {      
    pointIdIterator  = cellPointer->PointIdsBegin();
    cellPointsIds[0] = *pointIdIterator;
    pointIdIterator++;
    cellPointsIds[1] = *pointIdIterator;
    pointIdIterator++;
    cellPointsIds[2] = *pointIdIterator;
    
    mesh->DeleteFace( cellIndex );     
    
    pointIndex = mesh->FindFirstUnusedPointIndex();
    mesh->SetPoint( pointIndex, point );
    
    MeshType::QEPrimal* p;  
    p = mesh->AddFaceTriangle( cellPointsIds[0], cellPointsIds[1], pointIndex );
    newCellIds[0] = p->GetLeft();
    p = mesh->AddFaceTriangle( cellPointsIds[1], cellPointsIds[2], pointIndex );  
    newCellIds[1] = p->GetLeft();
    p = mesh->AddFaceTriangle( cellPointsIds[2], cellPointsIds[0], pointIndex );  
    newCellIds[2] = p->GetLeft();
    
    // Check the delaunay criterion 
    RecursiveFlipEdgeTest( mesh, pointIndex, newCellIds[0] );
    RecursiveFlipEdgeTest( mesh, pointIndex, newCellIds[1] );
    RecursiveFlipEdgeTest( mesh, pointIndex, newCellIds[2] );
    }
  
  return pointIndex;
}
//--------------------------------------------------------------------------------

//--------------------------------------------------------------------------------
// Delaunay triangulation process loop
//--------------------------------------------------------------------------------
template< class TInMesh, class TOutMesh = TInMesh >
MeshType::Pointer
PointSetToDelaunayTriangulationFilter::
DelaunayTriangulation( PointSetType::Pointer pointSet )
{    
  MeshType::Pointer mesh = MeshType::New();
  
  // NOTE STEF: Replace by true infinity
  PixelType infinity = pow( 10, 10 );
  CreateDummyMesh( mesh, infinity );
  
  InPointsContainer              *points       = pointSet->GetPoints();
  InPointsContainerConstIterator pointIterator = points->Begin();
  FaceRefType startingCellIndex;
  
  while( pointIterator != points->End() ) 
    {       
    InPointType     temporaryPoint;
    PointType       currentPoint;
    PointIdentifier currentPointIndex;
    CellAutoPointer currentCellPointer;
    
    temporaryPoint  = pointIterator.Value();
    currentPoint[0] = temporaryPoint[0];
    currentPoint[1] = temporaryPoint[1];
    currentPoint[2] = temporaryPoint[2];
    
    currentPointIndex = AddPoint( mesh, currentPoint, startingCellIndex );
    
    QEType* edge =  mesh->GetPoint( currentPointIndex ).GetEdge();
    startingCellIndex = edge->GetLeft(); 
    if( !mesh->GetCell( edge->GetLeft().first, currentCellPointer ) )
      {
      startingCellIndex = edge->GetRight(); 
      }
    
    ++pointIterator;
    }    
  
  DeleteDummyPoints( mesh );
  
  return mesh;
}
//--------------------------------------------------------------------------------

} // namespace itk

#endif // __itkPointSetToDelaunayTriangulationFilter_hxx