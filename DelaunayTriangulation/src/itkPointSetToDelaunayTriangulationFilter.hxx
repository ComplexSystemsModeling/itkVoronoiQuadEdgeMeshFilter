/*==========================================================================
 *                                                                         *
 * Delaunay Triangulation Incremental Algorithm                            *
 *                                                                         *
 *   Implementation for ITK by Stéphane U. Rigaud and Alexandre Gouaillard *
 *                                                                         *
 *   Input: Set of points itk::PointSet                                    *
 *                                                                         *
 *   Output: Delaunay triangulated mesh itk::QuadEdgeMesh                  *
 *                                                                         *
 *=========================================================================*/

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
  this->Superclass::SetNumberOfRequiredInputs( 1 );
  this->Superclass::SetNumberOfRequiredOutputs( 1 );
  
  this->Superclass::SetNthOutput( 0, MeshType::New() );

  m_DummyPoints = false;
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
DeleteDummyPoints( PointIdVectorType pts )
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
CreateDummyPoints()
{
  MeshPointer mesh = this->GetOutput();

  PixelType xmin(  pow(10,10) );
  PixelType ymin(  pow(10,10) );
  PixelType xmax( -pow(10,10) );
  PixelType ymax( -pow(10,10) );
  PointsContainerIterator ite = mesh->GetPoints()->Begin();
  while( ite != mesh->GetPoints()->End() )
    {
    if( ite.Value()[0] < xmin )
      {
      xmin = ite.Value()[0];
      }
    if( ite.Value()[0] > xmax )
      {
      xmax = ite.Value()[0];
      }
    if( ite.Value()[1] < ymin )
      {
      ymin = ite.Value()[1];
      }
    if( ite.Value()[1] > ymax )
      {
      ymax = ite.Value()[1];
      }
    ite++;
    }

  PixelType marge = ( ( abs(ymin) + abs(xmin) + xmax + ymax ) / 4 ) * 10 ;

  int expectedNumPts = 4;
  int expectedNumCells = 2;
  int simpleTriangleCells[6] = { 0, 1, 3, 0, 3, 2 };
  int i( 0 );  
  std::vector< PointType > pts( 4 );
  PointIdVectorType idx(4);  
  
  pts[i][0] = xmin - marge; pts[i][1] = ymin - marge; pts[i++][2] = 0.;
  pts[i][0] = xmax + marge; pts[i][1] = ymin - marge; pts[i++][2] = 0.; 
  pts[i][0] = xmin - marge; pts[i][1] = ymax + marge; pts[i++][2] = 0.; 
  pts[i][0] = xmax + marge; pts[i][1] = ymax + marge; pts[i++][2] = 0.;

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
    itkExceptionMacro( "itkPointSetToDelaunayTriangulationFilter::EdgeFlipTest" <<
                        " - Could not find the cell given in argument" );
    }
  if( !mesh->GetPoint( pointIndex, &pointCoord ) )
    {
    itkExceptionMacro( "itkPointSetToDelaunayTriangulationFilter::EdgeFlipTest" <<
                        " - Could not find the point given in argument" );
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

  // NOTE STEPH: the boolean parameter of the test should be remove in a near future  
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
  catch( itk::ExceptionObject & excp ) 
    {
    std::cerr << excp << std::endl;
    itkExceptionMacro( "itkPointSetToDelaunayTriangulationFilter::AddPoint" <<
                       " - Exception caught while walking on the mesh" );
    }
  // NOTE STEPH: Use a cleaner way to get the last value of a VectorContainer
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
    try 
      {
      this->RecursiveFlipEdgeTest( pointIndex, newCellIds[0] ); 
      }
    catch( itk::ExceptionObject & excp ) 
      {
      std::cerr << excp << std::endl;
      itkExceptionMacro( "itkPointSetToDelaunayTriangulationFilter::AddPoint" <<
                        " - Exception caught during the edgeflip test" );
      }
    try 
      {
      this->RecursiveFlipEdgeTest( pointIndex, newCellIds[1] ); 
      }
    catch( itk::ExceptionObject & excp ) 
      {
      std::cerr << excp << std::endl;
      itkExceptionMacro( "itkPointSetToDelaunayTriangulationFilter::AddPoint" <<
                        " - Exception caught during the edgeflip test" );
      }
    try 
      {
      this->RecursiveFlipEdgeTest( pointIndex, newCellIds[2] ); 
      }
    catch( itk::ExceptionObject & excp ) 
      {
      std::cerr << excp << std::endl;
      itkExceptionMacro( "itkPointSetToDelaunayTriangulationFilter::AddPoint" <<
                        " - Exception caught during the edgeflip test" );
      }
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
      
  PointIdVectorType pts;
  pts = this->CreateDummyPoints(); 
      
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
      try
        {
        currentPointIndex = this->AddPoint( currentPointIndex, FaceIndex );
        }
      catch( itk::ExceptionObject & excp ) 
        {
        std::cerr << excp << std::endl;
        itkExceptionMacro( "itkPointSetToDelaunayTriangulationFilter::GenerateData" << 
                           " - Exception caught while adding a point" );
        }
    
      QEType* edge =  mesh->GetPoint( currentPointIndex ).GetEdge();
      FaceIndex = edge->GetLeft(); 
      if( !mesh->GetCell( FaceIndex.first, currentCellPointer ) )
        {
        FaceIndex = edge->GetRight(); 
        }
      }
    ++pointIterator;
    }
  if( !m_DummyPoints )
    {
    this->DeleteDummyPoints( pts );
    }
  
  return true;
}
//--------------------------------------------------------------------------------

} // namespace itk

#endif // __itkPointSetToDelaunayTriangulationFilter_hxx
