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

#ifndef __itkPointSetToDelaunayTriangulationFilter_h
#define __itkPointSetToDelaunayTriangulationFilter_h

// ITK includes
#include "itkQuadEdgeMeshEulerOperatorFlipEdgeFunction.h"
#include "itkWalkInTriangulationFunction.h"
#include "itkPointInCircleGeometricalPredicateFunctor.h"

#include "itkQuadEdgeMeshSource.h"
#include "itkPointSetToQuadEdgeMeshFilter.h"

// TODO
// - Modify dummy points coordinates to be at infinity
//
// BUG
// - Random point generation with more than 2000 points can do segmentation fault

namespace itk {
  
template<class TInMesh, class TOutMesh>
class ITK_EXPORT PointSetToDelaunayTriangulationFilter : 
  public PointSetToQuadEdgeMeshFilter<TInMesh, TOutMesh>
{
public:
  
  /** Standard class typedefs. */
  typedef PointSetToDelaunayTriangulationFilter               Self;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<TInMesh, TOutMesh> Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;
  
  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );
  
  /** Run-time type information (and related methods).    */
  itkTypeMacro( PointSetToDelaunayTriangulationFilter, QuadEdgeMeshToQuadEdgeMeshFilter );
  
  typedef          TOutMesh                          MeshType;
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
  
  typedef          QuadEdgeMeshEulerOperatorFlipEdgeFunction<MeshType, QEType> FlipEdgeFunction;
  typedef typename FlipEdgeFunction::Pointer                                   FlipEdgeFunctionPointer;
  typedef          WalkInTriangulationFunction<MeshType>                       WalkInTriangulationFunction;
  typedef typename WalkInTriangulationFunction::Pointer                        WalkInTriangulationFunctionPointer;
  typedef          VectorContainer<unsigned int, FaceRefType>                  CellIdVectorContainerType;
  typedef typename CellIdVectorContainerType::Pointer                          CellIdVectorContainerTypePointer;
  
protected:
  
  /** Generic Methods */
  PointSetToDelaunayTriangulationFilter();
  
  virtual ~PointSetToDelaunayTriangulationFilter();
  
  void GenerateData();
  
  void PrintSelf( std::ostream & os, Indent indent ) const;
  
  /** Process Methods */
  void DeleteDummyPoints( std::vector<PointIdentifier> pts );
  
  std::vector<PointIdentifier> CreateDummyPoints( PixelType limit );
  
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