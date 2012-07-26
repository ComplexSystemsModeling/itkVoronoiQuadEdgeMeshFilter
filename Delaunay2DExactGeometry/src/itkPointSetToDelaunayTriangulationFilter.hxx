/*================================================================================
 *                                                                               *
 * Delaunay Triangulation Incremental Algorithm                                  *
 *                                                                               *
 *                                                                               *
 *   Implementation for ITK by Stéphane Ulysse Rigaud                            *
 *   IPAL (Image & Pervasive Access Lab) CNRS - A*STAR                           *
 *   Singapore                                                                   *
 *   http://www.ipal.cnrs.fr                                                     * 
 *                                                                               *
 *   Input: Set of points itk::PointSet                                          *
 *                                                                               *
 *   Output: Delaunay triangulated mesh itk::QuadEdgeMesh                        *
 *                                                                               *
 *===============================================================================*/

#ifndef __itkPointSetToDelaunayTriangulationFilter_hxx
#define __itkPointSetToDelaunayTriangulationFilter_hxx

#include "itkPointInCircleGeometricalPredicateFunctor.h"

namespace itk {

//--------------------------------------------------------------------------------
// Constructor & Destructor
//--------------------------------------------------------------------------------
template< class TInMesh >
PointSetToDelaunayTriangulationFilter< TInMesh >::
PointSetToDelaunayTriangulationFilter()
{
  this->Superclass::SetNumberOfRequiredInputs(1);
  this->Superclass::SetNumberOfRequiredOutputs(1);
  
  this->Superclass::SetNthOutput( 0, MeshType::New() );
}
//--------------------------------------------------------------------------------
  
  
//--------------------------------------------------------------------------------
// Generate Date Function
//--------------------------------------------------------------------------------
template< class TInMesh >
void
PointSetToDelaunayTriangulationFilter< TInMesh >::
GenerateData()
{
  this->CopyInputMeshToOutputMeshPoints();
  this->CopyInputMeshToOutputMeshPointData();
  
  this->DelaunayTriangulation();
}
//--------------------------------------------------------------------------------
  

//--------------------------------------------------------------------------------
// Delete Dummy Points
//--------------------------------------------------------------------------------
template< class TInMesh >
void
PointSetToDelaunayTriangulationFilter< TInMesh >::
DeleteDummyPoints( std::vector<PointIdentifier> pts )
{
  MeshPointer mesh = this->GetOutput();
  for( unsigned int i = 0; i < pts.size(); i++ )
    {
    QEType* e;
    while( e =  mesh->GetPoint( pts[i] ).GetEdge() )
      {
      mesh->LightWeightDeleteEdge( e );
      }
    mesh->DeletePoint( pts[i] );
    } 
}
//--------------------------------------------------------------------------------

  
//--------------------------------------------------------------------------------
// Create Dummy Points
//--------------------------------------------------------------------------------  
template< class TInMesh >
std::vector< typename PointSetToDelaunayTriangulationFilter< TInMesh >::PointIdentifier >
PointSetToDelaunayTriangulationFilter< TInMesh >::
CreateDummyPoints( PixelType   limit )
{
  MeshPointer mesh = this->GetOutput();

  int expectedNumPts = 4;
  int expectedNumCells = 2;
  int simpleTriangleCells[6] = { 0, 1, 3, 0, 3, 2 };
  int i( 0 );  
  std::vector< PointType > pts( 4 );
  std::vector<PointIdentifier> idx(4);  
  
  PixelType min = -limit;  
  PixelType max =  limit;

  pts[i][0] = min; pts[i][1] = min; pts[i++][2] = 0.;
  pts[i][0] = max; pts[i][1] = min; pts[i++][2] = 0.; 
  pts[i][0] = min; pts[i][1] = max; pts[i++][2] = 0.; 
  pts[i][0] = max; pts[i][1] = max; pts[i++][2] = 0.;

  for( i = 0; i < expectedNumPts; i++ )
    {
    idx[i] = mesh->AddPoint( pts[i] );
    }
  for( i = 0; i < expectedNumCells; i++ )
    {
    mesh->AddFaceTriangle( idx[ simpleTriangleCells[ i*3     ] ], 
                           idx[ simpleTriangleCells[ i*3 + 1 ] ], 
                           idx[ simpleTriangleCells[ i*3 + 2 ] ]);
    } 
  return idx;
}
//--------------------------------------------------------------------------------

  
//--------------------------------------------------------------------------------
// Recursive Delaunay Criterion Test and Flip Edge Function
//--------------------------------------------------------------------------------
template< class TInMesh >
bool
PointSetToDelaunayTriangulationFilter< TInMesh >::
RecursiveFlipEdgeTest( PointIdentifier pointIndex, 
                         FaceRefType     cell )
{ 
  MeshPointer mesh = this->GetOutput();

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
  if( !mesh->GetPoint( pointIndex, &pointCoord ) )
    {
    std::cerr << "ERROR - Could not find the point given in parameter." << std::endl;
    throw -1;
    }

  mesh->GetCell( cell.first, cellpointer );
  mesh->GetPoint( pointIndex, &pointCoord );        

  cellPointsIterator = cellpointer->PointIdsBegin();
  int i(0);
  while( cellPointsIterator != cellpointer->PointIdsEnd() )
    {
    p[i] = *cellPointsIterator;
    if( p[i] == pointIndex)
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

  // NOTE STEF: the boolean parameter of the test should be remove in a near future  
  if( TestPointInTriangleInMesh< MeshType >( mesh, adjCell.first, pointCoord, true ) ) 
    {
    FlipEdgeFunctionPointer flipedge = FlipEdgeFunction::New();
    flipedge->SetInput( mesh );
    
    e = flipedge->Evaluate( mesh->FindEdge( p[(r+1)%3], p[(r+2)%3] ) );
    
    this->RecursiveFlipEdgeTest( pointIndex, e->GetLeft() );
    this->RecursiveFlipEdgeTest( pointIndex, e->GetRight() );
    }
  return true;
}
//--------------------------------------------------------------------------------

  
//--------------------------------------------------------------------------------
// Add a Point to the Underconstruction Delaunay Mesh
//--------------------------------------------------------------------------------  
template< class TInMesh >
typename PointSetToDelaunayTriangulationFilter< TInMesh >::PointIdentifier
PointSetToDelaunayTriangulationFilter< TInMesh >::
AddPoint( PointIdentifier pointIndex, 
          FaceRefType     startingCell )
{
  MeshPointer mesh = this->GetOutput();

  // Walk the mesh to find the cell containing the new point to add	
  WalkInTriangulationFunctionPointer walk       = WalkInTriangulationFunction::New();
  CellIdVectorContainerTypePointer   cellIdList = CellIdVectorContainerType::New();
  try 
    {
    cellIdList = walk->Evaluate( mesh, mesh->GetPoint(pointIndex), &startingCell );
    }
  catch( int e ) 
    {
    std::cerr << "ERROR - Exception caught in the WalkInTriangulation process" << std::endl;
    throw -1;
    }
  // NOTE STEF: Use a cleaner way to get the last value of a VectorContainer
  FaceRefType cellIndex = ( --cellIdList->End() )->Value();

  CellAutoPointer              cellPointer;
  PointIdIterator              pointIdIterator;
  std::vector<PointIdentifier> cellPointsIds( 3 );
  std::vector<FaceRefType>     newCellIds( 3 );

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
    
    typename MeshType::QEPrimal* p;  
    p = mesh->AddFaceTriangle( cellPointsIds[0], cellPointsIds[1], pointIndex );
    newCellIds[0] = p->GetLeft();
    p = mesh->AddFaceTriangle( cellPointsIds[1], cellPointsIds[2], pointIndex );  
    newCellIds[1] = p->GetLeft();
    p = mesh->AddFaceTriangle( cellPointsIds[2], cellPointsIds[0], pointIndex );  
    newCellIds[2] = p->GetLeft();
    
    // Check the delaunay criterion 
    this->RecursiveFlipEdgeTest( pointIndex, newCellIds[0] ); 
    this->RecursiveFlipEdgeTest( pointIndex, newCellIds[1] ); 
    this->RecursiveFlipEdgeTest( pointIndex, newCellIds[2] ); 
    }
  return pointIndex;
}
//--------------------------------------------------------------------------------

  
//--------------------------------------------------------------------------------
// Main Loop of the Delaunay Triangulation Algorithm
//--------------------------------------------------------------------------------  
template< class TInMesh >
bool
PointSetToDelaunayTriangulationFilter< TInMesh >::
DelaunayTriangulation()
{ 
  MeshPointer mesh = this->GetOutput();
      
  std::vector<PointIdentifier> pts;
  PixelType limit = pow( 10, 10 );
  pts = this->CreateDummyPoints( limit ); 
      
  PointsContainer         *points       = mesh->GetPoints();
  PointsContainerIterator pointIterator = points->Begin();
  FaceRefType             FaceIndex;
  
  FaceIndex.first = 0;
  while( pointIterator != points->End() ) 
    { 
    PointType       currentPoint;
    PointIdentifier currentPointIndex;
    CellAutoPointer currentCellPointer;
    
    currentPoint  = pointIterator.Value();
    currentPointIndex = pointIterator.Index();
    
    if( !mesh->GetPoint( currentPointIndex ).GetEdge() )
      {        
      currentPointIndex = this->AddPoint( currentPointIndex, FaceIndex );
    
      QEType* edge =  mesh->GetPoint( currentPointIndex ).GetEdge();
      FaceIndex = edge->GetLeft(); 
      if( !mesh->GetCell( FaceIndex.first, currentCellPointer ) )
        {
        FaceIndex = edge->GetRight(); 
        }
      }
    ++pointIterator;
    }
  this->DeleteDummyPoints( pts );

  return true;
}
//--------------------------------------------------------------------------------

} // namespace itk

#endif // __itkPointSetToDelaunayTriangulationFilter_hxx
