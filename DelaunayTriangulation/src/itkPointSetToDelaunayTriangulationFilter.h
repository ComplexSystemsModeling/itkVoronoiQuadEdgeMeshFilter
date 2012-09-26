/*================================================================================
 *                                                                               *
 * Delaunay Triangulation Incremental Algorithm                                  *
 *                                                                               *
 *                                                                               *
 *   Implementation for ITK by Stéphane Ulysse Rigaud and Alexandre Gouaillard   *
 *   IPAL (Image & Pervasive Access Lab) CNRS - A*STAR                           *
 *   Singapore                                                                   *
 *   http://www.ipal.cnrs.fr                                                     * 
 *                                                                               *
 *   Input: Set of points itk::PointSet                                          *
 *                                                                               *
 *   Output: Delaunay triangulated mesh itk::QuadEdgeMesh                        *
 *                                                                               *
 *===============================================================================*/

#ifndef __itkPointSetToDelaunayTriangulationFilter_h
#define __itkPointSetToDelaunayTriangulationFilter_h

// ITK includes
#include "itkQuadEdgeMeshEulerOperatorFlipEdgeFunction.h"
#include "itkWalkInTriangulationFunction.h"
#include "itkPointInCircleGeometricalPredicateFunctor.h"

#include "itkQuadEdgeMesh.h"
#include "itkPointSetToMeshFilter.h"

// TODO
// - Modify dummy points coordinates to be at infinity 
//
// BUG
// - segmentation fault with some particular configuration of points

namespace itk {
  
template< class TInMesh >
class ITK_EXPORT PointSetToDelaunayTriangulationFilter : 
  public PointSetToMeshFilter< TInMesh, QuadEdgeMesh< typename TInMesh::PixelType, GetPointSetDimension<TInMesh>::PointDimension > >
{
public:
  
  /** Force OutputType */  
  typedef QuadEdgeMesh< typename TInMesh::PixelType, GetPointSetDimension<TInMesh>::PointDimension > MeshType;
  
  /** Standard class typedefs. */
  typedef PointSetToDelaunayTriangulationFilter      Self;
  typedef PointSetToMeshFilter< TInMesh, MeshType >  Superclass;
  typedef SmartPointer< Self >                       Pointer;
  typedef SmartPointer< const Self >                 ConstPointer;
  
  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );
  
  /** Run-time type information (and related methods).    */
  itkTypeMacro( PointSetToDelaunayTriangulationFilter, PointSetToMeshFilter );
  
  /** Convenient macro   */
  itkStaticConstMacro( PointDimension, unsigned int, MeshType::Traits::PointDimension );
  
  /** Some convenient typedefs.    */
  typedef typename MeshType::Pointer                 MeshPointer;
  typedef typename MeshType::PixelType               PixelType;
  typedef typename MeshType::PointType               PointType;
  typedef typename MeshType::PointIdentifier         PointIdentifier;
  typedef typename MeshType::PointIdIterator         PointIdIterator;
  typedef typename MeshType::PointsContainer         PointsContainer;
  typedef typename MeshType::PointsContainerIterator PointsContainerIterator;
  typedef typename MeshType::CellAutoPointer         CellAutoPointer;  
  typedef typename MeshType::FaceRefType             FaceRefType;
  typedef typename MeshType::QEType                  QEType;
  typedef typename MeshType::CellType                CellType;
  typedef typename CellType::PointIdConstIterator    PointIdConstIterator;
  
  typedef          QuadEdgeMeshEulerOperatorFlipEdgeFunction< MeshType, QEType > FlipEdgeFunction;
  typedef typename FlipEdgeFunction::Pointer                                     FlipEdgeFunctionPointer;
  typedef          WalkInTriangulationFunction< MeshType >                       WalkInTriangulationFunction;
  typedef typename WalkInTriangulationFunction::Pointer                          WalkInTriangulationFunctionPointer;
  typedef          VectorContainer< unsigned int, FaceRefType >                  CellIdVectorContainerType;
  typedef typename CellIdVectorContainerType::Pointer                            CellIdVectorContainerTypePointer;

  typedef std::vector< PointIdentifier >    PointIdVectorType;
  
protected:
  
  /** Generic Methods */
  PointSetToDelaunayTriangulationFilter();
  ~PointSetToDelaunayTriangulationFilter() {}
  
  void GenerateData();
    
  /** Process Methods */
  void DeleteDummyPoints( PointIdVectorType pts );
  
  PointIdVectorType CreateDummyPoints( PixelType limit );
  
  bool RecursiveFlipEdgeTest( PointIdentifier pointIndex, FaceRefType cell );
  
  PointIdentifier AddPoint( PointIdentifier pointIndex, FaceRefType startingCell );

  bool DelaunayTriangulation();
  
private:
  
  PointSetToDelaunayTriangulationFilter( const Self & ); // purposely not implemented
  void operator=( const Self & );                        // purposely not implemented
  
}; // QuadEdgeMeshToDelaunayTriangulationFilter

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPointSetToDelaunayTriangulationFilter.hxx"
#endif

#endif // __itkPointSetToDelaunayTriangulationFilter_h
