#include "itkQuadEdgeMesh.h"
#include "itkDefaultStaticMeshTraits.h"

#include <iostream>

namespace itk
{
/**
 * \class QuadEdgeMeshWithDual
 *
 * \brief Mesh class for 2D manifolds embedded in ND space and their dual.
 *
 * \author Alexandre Gouaillard, Humayun
 *
 * This implementation was contributed as a paper to the Insight Journal
 *
 * \ingroup ITK-QuadEdgeMesh
 */
template< typename TPixel, unsigned int VDimension,
          typename TTraits = QuadEdgeMeshTraits< TPixel, VDimension, bool, bool > >
class ITK_EXPORT QuadEdgeMeshWithDual : public QuadEdgeMesh< TPixel, VDimension, TTraits >
{
public:
  /** Input template parameters. */
  typedef TTraits Traits;
  typedef TPixel  PixelType;

  /** Standard typedefs. */
  typedef QuadEdgeMeshWithDual                       Self;
  typedef QuadEdgeMesh< TPixel, VDimension, Traits > Superclass;
  typedef SmartPointer< Self >                       Pointer;
  typedef SmartPointer< const Self >                 ConstPointer;

  /** Convenient constants obtained from MeshTraits. */
  itkStaticConstMacro(PointDimension, unsigned int,
                      Superclass::PointDimension);
  itkStaticConstMacro(MaxTopologicalDimension, unsigned int,
                      Superclass::MaxTopologicalDimension);

  /** Basic Object interface. */
  itkNewMacro(Self);
  itkTypeMacro(QuadEdgeMeshWithDual, QuadEdgeMesh);

  typedef typename Superclass::CellPixelType           CellPixelType;
  typedef typename Superclass::CoordRepType            CoordRepType;
  typedef typename Superclass::PointHashType           PointHashType;
  typedef typename Superclass::CellTraits              CellTraits;
  typedef typename Superclass::PointIdInternalIterator PointIdInternalIterator;
  typedef typename Superclass::PointIdIterator         PointIdIterator;

  typedef typename Superclass::PointIdentifier               PointIdentifier;
  typedef typename Superclass::PointType                     PointType    ;
  typedef typename Superclass::PointIdList                   PointIdList;
  typedef typename Superclass::PointsContainer               PointsContainer;
  typedef typename Superclass::PointsContainerPointer        PointsContainerPointer;
  typedef typename Superclass::PointsContainerConstIterator  PointsContainerConstIterator;
  typedef typename Superclass::PointsContainerIterator       PointsContainerIterator;
  typedef typename Superclass::PointDataContainer            PointDataContainer;
  typedef typename Superclass::PointDataContainerPointer     PointDataContainerPointer;
  typedef typename Superclass::PointDataContainerIterator    PointDataContainerIterator;

  typedef typename Superclass::CellIdentifier CellIdentifier;
  typedef typename Superclass::CellType      CellType;
  typedef typename Superclass::CellAutoPointer CellAutoPointer;
  typedef typename Superclass::CellMultiVisitorType CellMultiVisitorType;
  typedef typename Superclass::CellsContainer      CellsContainer;
  typedef typename Superclass::CellsContainerPointer CellsContainerPointer;
  typedef typename Superclass::CellsContainerConstIterator CellsConatinerConstIterator;
  typedef typename Superclass::CellsContainerIterator     CellsContainerIterator;
  typedef typename Superclass::CellDataContainer          CellDataContainer;
  typedef typename Superclass::CellDataContainerPointer   CellDataContainerPointer;
  typedef typename Superclass::CellDataContainerIterator  CellDataContainerIterator;

  typedef typename Superclass::PolygonCellType PolygonCellType;
  typedef typename Superclass::EdgeCellType    EdgeCellType;

  /** Accessors */
  const CellsContainerPointer  GetDualCells()     const { return m_DualCellsContainer;     };
  const CellsContainerPointer  GetDualEdgeCells() const { return m_DualEdgeCellsContainer; };
  const PointsContainerPointer GetDualPoints()    const { return m_DualPointsContainer;    };

  /** Add cells / edges / points */

  void SetDualPoint( PointIdentifier id, PointType p )
    {
    return m_DualPointsContainer->InsertElement( id, p );
    }
  PointIdentifier AddDualPoint( PointType p )
    {
    PointIdentifier pid = m_DualPointsContainer->size();
    this->SetDualPoint( pid, p );
    return( pid );
    }

  // NOTE ALEX: this does not create underlying QE layer for now
  CellIdentifier AddDualFace( const PointIdList & points )
    {
    // Create the cell without underlying QE layer and add it to the container
    PolygonCellType *faceCell = new PolygonCellType( points.size() );

    CellIdentifier fid = m_DualCellsContainer->size();
    faceCell->SetIdent( fid );
    CellAutoPointer face;
    face.TakeOwnership( faceCell );
    for( unsigned int i = 0; i < points.size(); i++ )
      face->SetPointId( i, points[i] );

    m_DualCellsContainer->InsertElement( fid, face.ReleaseOwnership() );

    return fid;
    }

  void AddDualEdge( const PointIdentifier& pid1, const PointIdentifier& pid2 )
    {
    // Create the cell without unterlying QE layer and add it to the container
    // NOTE ALEX: is it actually possible to have an edgecell type without QE?
    EdgeCellType *edgeCell = new EdgeCellType( );

    CellIdentifier   eid = m_DualEdgeCellsContainer->size();
    edgeCell->SetIdent( eid );
    CellAutoPointer edge;
    edge.TakeOwnership( edgeCell );
    edge->SetPointId( 0, pid1 );
    edge->SetPointId( 1, pid2 );

    m_DualEdgeCellsContainer->InsertElement( eid, edge.ReleaseOwnership() );
    }

  PointIdentifier GetNumberOfDualPoints() { return m_DualPointsContainer->size(); };

private:
  CellsContainerPointer  m_DualCellsContainer;
  CellsContainerPointer  m_DualEdgeCellsContainer;
  PointsContainerPointer m_DualPointsContainer;

protected:
  /** Constructor and Destructor. */
  QuadEdgeMeshWithDual()
    {
    m_DualCellsContainer     = CellsContainer::New();
    m_DualEdgeCellsContainer = CellsContainer::New();
    m_DualPointsContainer    = PointsContainer::New();
    };

  virtual ~QuadEdgeMeshWithDual() { };

};
}
