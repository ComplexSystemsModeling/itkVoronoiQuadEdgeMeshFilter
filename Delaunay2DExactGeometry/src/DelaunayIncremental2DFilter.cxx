
//--------------
// itk code
#include "itkPointSet.h"
#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshTraits.h"
#include "itkQuadEdgeMeshPolygonCell.h"
#include "itkVectorContainer.h"
#include "itkQuadEdgeMeshEulerOperatorsTestHelper.h"
#include "itkQuadEdgeMeshEulerOperatorFlipEdgeFunction.h"
#include "itkDelaunayConformingQuadEdgeMeshFilter.h"

//--------------
// itk SimplexMesh
#include "itkQuadEdgeMeshToQuadEdgeMeshWithDualFilter.h"
#include "itkQuadEdgeMeshWithDualAdaptor.h"

//---------------
// our code
#include "itkWalkInTriangulationFunction.h"
#include "itkPointInCircleGeometricalPredicateFunctor.h"

//---------------
// std
#include <iostream>
#include <limits>

//--------------
// temporary
#include "itkVTKPolyDataWriter.h"
#include "itkVTKPolyDataReader.h"

//---------------------------------
// to be incapsulate into itk class
namespace itk
{
// TODO
// - PointSetSource class
// - pointsetToPointSetFilter class (inherits from PointSetSource)
// - QuadEdgeMeshSource class (does it exists already?)
// - rewire QuadEdgeMeshToQuadEdgeMeshFilter to use the class above.
//
// - Write GenerateData()
// - use superclass GenerateData() to copy input to output then work with getOutput()

template< class TInMesh, class TOutMesh=TInMesh >
class PointSetToDelaunayTriangulationFilter
 : public QuadEdgeMeshToQuadEdgeMeshFilter< TInMesh, TOutMesh >
{
public:
  typedef PointSetToDelaunayTriangulationFilter                       Self;
  typedef SmartPointer< Self >                                        Pointer;
  typedef SmartPointer< const Self >                                  ConstPointer;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter< TInMesh, TOutMesh > Superclass;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(PointSetToDelaunayTriangulationFilter, QuadEdgeMeshToQuadEdgeMeshFilter);

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro(Self);

  typedef float PixelType;
  static const unsigned int Dimensions = 3;

  typedef itk::PointSet< PixelType, Dimensions >      PointSetType;

  typedef PointSetType::PointType                     InputPointType;
  typedef PointSetType::PointIdentifier               InputPointIdentifier;
  typedef PointSetType::PointsContainer               InputPointsContainer;
  typedef PointSetType::PointsContainerConstIterator  InputPointsContainerConstIterator;

  typedef itk::QuadEdgeMesh< PixelType, Dimensions >  MeshType;

  typedef MeshType::PointType                         PointType;
  typedef MeshType::CellType                          CellType;
  typedef MeshType::PointIdentifier                   PointIdentifier;
  typedef MeshType::CellIdentifier                    CellIdentifier;
  typedef MeshType::PointsContainer                   PointsContainer;
  typedef MeshType::CellsContainer                    CellsContainer;
  typedef MeshType::PointIdList                       PointIdList;
  typedef MeshType::QEType                            QEType;
  typedef MeshType::CellsContainerIterator            CellsContainerIteratorType;
  typedef MeshType::PointsContainerConstIterator      PointsContainerConstIterator;

  typedef CellType::PointIdConstIterator              PointIdConstIterator;
  typedef CellType::PointIdIterator                   PointIdIterator;
  typedef CellType::CellAutoPointer                   CellAutoPointer;

  typedef QEType::DualOriginRefType                   DualOriginRefType;

  typedef itk::WalkInTriangulationFunction< MeshType >                       WalkInTriangulationFunction;
  typedef itk::VectorContainer< unsigned int, int >                          CellIdVectorContainerType;
  typedef itk::QuadEdgeMeshPolygonCell< CellType >                           QEPolygonCellType;
  typedef itk::QuadEdgeMeshEulerOperatorFlipEdgeFunction< MeshType, QEType > FlipEdgeFunction;
  //typedef itk::QuadEdgeMeshToQuadEdgeMeshWithDualFilter< MeshType >          FillDualFilterType;
  //typedef itk::QuadEdgeMeshWithDualAdaptor< MeshType >                       DualAdaptorType;


//--------------------------------------------------------------------------------
// Dummy point deletion
//
void
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
//
// Create a dummy square mesh for triangulation initialisation
// limit argument is to determine the size of the square
//
void
CreateDummyMesh( MeshType::Pointer mesh, 
                 PixelType         limit )
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
MeshType::Pointer
RecursiveFlipEdgeTest( MeshType::Pointer mesh, 
                       PointIdentifier   point, 
                       CellIdentifier    cell )
{
  CellAutoPointer      cellpointer;
  PointType            pointCoord;
  PointIdConstIterator cellPointsIterator;     
  PointIdentifier      q[3];                // points index of T
  PointIdentifier      p[3];                // points index of aT
  int r(0), k(0);                           // current point and oposite point position  

  if( !mesh->GetCell( cell, cellpointer ) )
    {
    std::cerr << "ERROR - Could not find the cell given in parameter." << std::endl;
    throw -1;
    }
  if( !mesh->GetPoint( point, &pointCoord ) )
    {
    std::cerr << "ERROR - Could not find the point given in parameter." << std::endl;
    throw -1;
    }
  
  mesh->GetCell(   cell, cellpointer );
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
    return mesh;
    }
  
  DualOriginRefType adjCell = e->GetRight();
  // NOTE STEF: Temporary check
  if( adjCell == cell )
    {
    adjCell = e->GetLeft();
    }
  mesh->GetCell( adjCell, cellpointer );       
  
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
  if( TestPointInTriangleInMesh< MeshType >( mesh, adjCell, pointCoord, true ) ) 
    {
    FlipEdgeFunction::Pointer flipedge = FlipEdgeFunction::New();
    flipedge->SetInput( mesh );
    
    e = flipedge->Evaluate( mesh->FindEdge( p[(r+1)%3], p[(r+2)%3] ) );
    
    RecursiveFlipEdgeTest( mesh, point, e->GetLeft() );
    RecursiveFlipEdgeTest( mesh, point, e->GetRight() );
    }
  
  return mesh;
}
//--------------------------------------------------------------------------------


//--------------------------------------------------------------------------------
PointIdentifier
AddPoint( MeshType::Pointer mesh, 
          PointType         point, 
          CellIdentifier    startingCell )
{
  
  // Walk the mesh to find the cell containing the new point to add	
  WalkInTriangulationFunction::Pointer walk       = WalkInTriangulationFunction::New();
  CellIdVectorContainerType::Pointer   cellIdList = CellIdVectorContainerType::New();
  try 
    {
    cellIdList = walk->Evaluate( mesh, point, startingCell );
    }
  catch( int e ) 
    {
    std::cerr << "ERROR - Exception caught in the WalkInTriangulation process" << std::endl;
    throw -1;
    }
  // NOTE STEF: Use a cleaner way to get the last value of a VectorContainer
  CellIdentifier cellIndex = ( --cellIdList->End() )->Value();

  CellAutoPointer               cellPointer;
  PointIdIterator               pointIdIterator;
  std::vector< CellIdentifier > cellPointsIds( 3 );
  std::vector< CellIdentifier > newCellIds( 3 );
  PointIdentifier               pointIndex;
   
  // Split the cell in 3 new cell 
  if( mesh->GetCell( cellIndex, cellPointer ) )
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
//
MeshType::Pointer
DelaunayTriangulation( PointSetType* pointSet )
{

  MeshType::Pointer mesh = MeshType::New();

  // NOTE STEF: Replace by true infinity
  PixelType infinity = pow( 10, 10 );
  CreateDummyMesh( mesh, infinity );

  InputPointsContainer              *points           = pointSet->GetPoints();
  InputPointsContainerConstIterator pointIterator     = points->Begin();
  CellIdentifier                    startingCellIndex = 0;

  while( pointIterator != points->End() ) 
    {       
    InputPointType  temporaryPoint;
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
    if( !mesh->GetCell( edge->GetLeft(), currentCellPointer ) )
      {
      startingCellIndex = edge->GetRight(); 
      }
    
    ++pointIterator;
    }    

  DeleteDummyPoints( mesh );

  return mesh;
}
//--------------------------------------------------------------------------------
};
} // endof namespace


//--------------------------------------------------------------------------------
// Random coordonates generation function
//
template< class TMesh >
std::vector< typename TMesh::PointType >
GenerateRandomCoordinates( const unsigned int& iN )
{
  typedef typename TMesh::PointType         TPointType;
  typedef typename TPointType::CoordRepType TCoordRepType;
  std::vector< TPointType > oPt( iN );
  
  // NOTE ALEX: is this cross platform? 
  srand(time(NULL));
  
  for( unsigned int i = 0; i < iN; i++ )
    {
    oPt[ i ][0] = static_cast< TCoordRepType >( rand() % 9000 - 4500 );
    oPt[ i ][1] = static_cast< TCoordRepType >( rand() % 9000 - 4500 );
    oPt[ i ][2] = static_cast< TCoordRepType >( 0. );
    }
  return oPt;
}
//--------------------------------------------------------------------------------


//--------------------------------------------------------------------------------
// Circle coordonates generation function
//
template< class TMesh >
std::vector< typename TMesh::PointType >
GenerateCircleCoordinates( const unsigned int& r )
{
  typedef typename TMesh::PointType         TPointType;
  typedef typename TPointType::CoordRepType TCoordRepType;
  std::vector< TPointType > oPt;
  const float DEG2RAD = 3.14159 / 180;
  TPointType p;
  
  for( unsigned int i = 0; i < 360; i++ )
    {
    float deg = i * DEG2RAD; 
    p[0] = static_cast< TCoordRepType >( cos(deg)*r );
    p[1] = static_cast< TCoordRepType >( sin(deg)*r );
    p[2] = static_cast< TCoordRepType >( 0. );
    oPt.push_back( p );
    }
  return oPt;
}
//--------------------------------------------------------------------------------


//--------------------------------------------------------------------------------
// Concentric Circle coordonates generation function
//
template< class TMesh >
std::vector< typename TMesh::PointType >
GenerateConcentricCoordinates( const unsigned int& r )
{
  typedef typename TMesh::PointType         TPointType;
  typedef typename TPointType::CoordRepType TCoordRepType;
  std::vector< TPointType > oPt;
  const float DEG2RAD = 3.14159 / 180;
  TPointType p;

  for( unsigned int j = r; j > 0; j-- )
    {
    for( unsigned int i = 0; i < 360; i++ )
      {
      float deg = i * DEG2RAD; 
      p[0] = static_cast< TCoordRepType >( cos(deg)*j );
      p[1] = static_cast< TCoordRepType >( sin(deg)*j );
      p[2] = static_cast< TCoordRepType >( 0. );
      oPt.push_back( p ); 
      }
    }
  p[0] = 0.; p[1] = 0.; p[2] = 0.;
  oPt.push_back( p ); 
  return oPt;
}
//--------------------------------------------------------------------------------


//--------------------------------------------------------------------------------
// Cross coordonates generation function
//
template< class TMesh >
std::vector< typename TMesh::PointType >
GenerateCrossCoordinates( const unsigned int& iN )
{
  typedef typename TMesh::PointType         TPointType;
  typedef typename TPointType::CoordRepType TCoordRepType;
  std::vector< TPointType > oPt;
  TPointType p;

  for( unsigned int i = 0 ; i < 2*iN+1; i++ )
    {
    p[0] = static_cast< TCoordRepType >( i  );
    p[1] = static_cast< TCoordRepType >( iN );
    p[2] = static_cast< TCoordRepType >( 0. );
    oPt.push_back( p );
    if( i != iN )
      {
      p[0] = static_cast< TCoordRepType >( iN );
      p[1] = static_cast< TCoordRepType >( i  );
      p[2] = static_cast< TCoordRepType >( 0. );
      oPt.push_back( p );
      }
    }
  return oPt;
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

  typedef float PixelType;
  static const unsigned int Dimensions = 3;
  typedef itk::QuadEdgeMesh< PixelType, Dimensions >      MeshType;
  typedef MeshType::PointType                             InputPointType;

  typedef itk::DelaunayConformingQuadEdgeMeshFilter< MeshType, MeshType > ValidityTestType;
  typedef itk::VTKPolyDataWriter< MeshType >                              MeshWriterType;

  int type = atoi( argv[1] );
  int meshSize = atoi( argv[2] );
  int expectedNumPts = 0;
  std::vector< InputPointType > pts;
  
  typedef itk::PointSet< PixelType, Dimensions >      PointSetType;
  PointSetType::Pointer pointSet = PointSetType::New();
  MeshType::Pointer triangulatedMesh = MeshType::New();
      
  // -------------------------------------------------
  // Toy Point Set creation

  switch(type) 
    {
    case 1 :     
      pts = GeneratePointCoordinates< MeshType >( meshSize ); 
      break;
    case 2 :
      pts = GenerateRandomCoordinates< MeshType >( meshSize ); 
      break;
    case 3 :
      pts = GenerateCircleCoordinates< MeshType >( meshSize );    
      break;
    case 4 :
      pts = GenerateConcentricCoordinates< MeshType >( meshSize );
      break;
    case 5 :
      pts = GenerateCrossCoordinates< MeshType >( meshSize );
      break;
    default:
      std::cerr << "Exit wrong arguments" << std::endl;
      std::cerr << "Usage - arg[1] [1|2] = [\"reg\"|\"rand\"] - arg[2] int = [rowsize|nbPoint]";
      std::cerr << std::endl;
      
      return EXIT_FAILURE;
    }

  expectedNumPts = pts.size();
  for( int i = 0; i < expectedNumPts ; i++ )
    {
    pointSet->SetPoint( i, pts[i] );
    }

  // -------------------------------------------------- 
  // Delaunay Construction
  typedef itk::PointSetToDelaunayTriangulationFilter< MeshType > myFilterType;
  myFilterType::Pointer myFilter = myFilterType::New();
  try
    {
    triangulatedMesh = myFilter->DelaunayTriangulation( pointSet );
    }
  catch( int e ) 
    {
    std::cerr << "Exception WiT Caught" << std::endl;
    return EXIT_FAILURE;
    }
  
  MeshWriterType::Pointer write = MeshWriterType::New();
  write->SetFileName("./OutputDelaunayMesh.vtk");
  write->SetInput( triangulatedMesh );
  write->Update();
 
  // -------------------------------------------------
  // Delaunay Validation

  ValidityTestType::Pointer test = ValidityTestType::New();
  test->SetInput( triangulatedMesh );
  test->GraftOutput( triangulatedMesh );
  test->Update();
  
  if( test->GetNumberOfEdgeFlips() > 0 )
    {
    return EXIT_FAILURE;
    }

  // -------------------------------------------------
  // Delaunay Dual => Voronoi

  FillDualFilterType::Pointer fillDual = FillDualFilterType::New();
  fillDual->SetInput( triangulatedMesh );
  try
    {
    fillDual->Update();
    }
  catch( ... )
    {
    return EXIT_FAILURE;
    }

  DualAdaptorType* adaptor = new DualAdaptorType();
  adaptor->SetInput( fillDual->GetOutput() );

  DualMeshWriterType::Pointer dualwrite = DualMeshWriterType::New();
  dualwrite->SetFileName("./OutputVoronoiTesselation.vtk");
  dualwrite->SetInput( adaptor );
  dualwrite->Update();

  // ------------------------------------------------
  // End process ... by by

  return EXIT_SUCCESS;

}
//--------------------------------------------------------------------------------
