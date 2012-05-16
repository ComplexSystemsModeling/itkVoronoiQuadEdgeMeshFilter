#ifndef __itkPointSetToDelaunayTriangulationFilter_h
#define __itkPointSetToDelaunayTriangulationFilter_h

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
//#include "itkQuadEdgeMeshToQuadEdgeMeshWithDualFilter.h"
//#include "itkQuadEdgeMeshWithDualAdaptor.h"

//---------------
// our code
#include "itkWalkInTriangulationFunction.h"
#include "itkPointInCircleGeometricalPredicateFunctor.h"
#include "itkPointSetSource.h"
#include "itkPointSetToPointSetFilter.h"
#include "itkQuadEdgeMeshSource.h"


//---------------
// std
#include <iostream>
#include <limits>

namespace itk 
{
// TOCHECK
// - PointSetSource class
// - pointsetToPointSetFilter class
// - QuadEdgeMeshSource class
	
// TODO
// - rewire QuadEdgeMeshToQuadEdgeMeshFilter to use the class above.
//
// - Write GenerateData()
// - use superclass GenerateData() to copy input to output then work with getOutput()

template< 
  class TInMesh, 
  class TOutMesh = TInMesh 
  >
class ITK_EXPORT PointSetToDelaunayTriangulationFilter
  : public QuadEdgeMeshToQuadEdgeMeshFilter< TInMesh, TOutMesh >
{

public:

  /** Basic type     */
  typedef PointSetToDelaunayTriangulationFilter                 Self;
  typedef SmartPointer< Self >                                  Pointer;
  typedef SmartPointer< const Self >                            ConstPointer;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter< TInMesh, TOutMesh > Superclass;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(PointSetToDelaunayTriangulationFilter, QuadEdgeMeshToQuadEdgeMeshFilter);

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro(Self);

  /** Our type      */
  typedef float PixelType;
  static const unsigned int Dimensions = 3;

  typedef itk::PointSet< PixelType, Dimensions >      PointSetType;
  
  typedef PointSetType::PointType                     InPointType;
  typedef PointSetType::PointsContainer               InPointsContainer;
  typedef PointSetType::PointsContainerConstIterator  InPointsContainerConstIterator;
  
  typedef itk::QuadEdgeMesh< PixelType, Dimensions >  MeshType;

  typedef MeshType::PointType                         PointType;
  typedef MeshType::PointIdentifier                   PointIdentifier;
  typedef MeshType::CellType                          CellType;
  typedef MeshType::CellIdentifier                    CellIdentifier;
  typedef MeshType::PointsContainer                   PointsContainer;
  typedef MeshType::CellsContainer                    CellsContainer;
  typedef MeshType::PointIdList                       PointIdList;
  typedef MeshType::QEType                            QEType;
  typedef MeshType::CellsContainerIterator            CellsContainerIteratorType;
  typedef MeshType::PointsContainerConstIterator      PointsContainerConstIterator;
  
  typedef MeshType::FaceRefType                       FaceRefType;
  typedef MeshType::VertexRefType                     VertexRefType;

  typedef CellType::PointIdConstIterator              PointIdConstIterator;
  typedef CellType::PointIdIterator                   PointIdIterator;
  typedef CellType::CellAutoPointer                   CellAutoPointer;

  typedef QEType::DualOriginRefType                   DualOriginRefType;
  typedef QEType::OriginRefType                       OriginRefType;

  typedef itk::WalkInTriangulationFunction< MeshType >                       WalkInTriangulationFunction;
  typedef itk::VectorContainer< unsigned int, FaceRefType >                  CellIdVectorContainerType;
  typedef itk::QuadEdgeMeshPolygonCell< CellType >                           QEPolygonCellType;
  typedef itk::QuadEdgeMeshEulerOperatorFlipEdgeFunction< MeshType, QEType > FlipEdgeFunction;

  /** Our methods     */
  void DeleteDummyPoints( MeshType::Pointer mesh );

  void CreateDummyMesh( MeshType::Pointer mesh, 
                   PixelType         limit );
  
  bool RecursiveFlipEdgeTest( MeshType::Pointer mesh, 
                              PointIdentifier   point, 
                              FaceRefType       cell );
  
  PointIdentifier AddPoint( MeshType::Pointer mesh, 
                            PointType         point, 
                            FaceRefType       startingCell );

  MeshType::Pointer DelaunayTriangulation( PointSetType::Pointer pointSet );

}; // class PointSetToDelaunayTriangulationFilter

} // namespace itk

#endif // __itkPointSetToDelaunayTriangulationFilter_h