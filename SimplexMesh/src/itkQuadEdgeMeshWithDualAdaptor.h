#include <iostream>

namespace itk
{

template< class TMesh >
class QuadEdgeMeshWithDualAdaptor
{
public:
  /** Standard typedefs. */
  typedef QuadEdgeMeshWithDualAdaptor  Self;
  typedef Self*                        Pointer;
  typedef const Self*                  ConstPointer;

  /** to foul the itkSuperclassTraitMacro */
  typedef TMesh Superclass;

  /** Convenient constants obtained from MeshTraits. */
  itkStaticConstMacro(PointDimension, unsigned int,
                      Superclass::PointDimension);
  itkStaticConstMacro(MaxTopologicalDimension, unsigned int,
                      Superclass::MaxTopologicalDimension);

  /** those types need to be defined for the itkVTKPolydataWriter */
  typedef typename Superclass::PixelType               PixelType;
  typedef typename Superclass::PointType               PointType;
  typedef typename Superclass::CellType                CellType;
  typedef typename Superclass::PointIdentifier         PointIdentifier;
  typedef typename Superclass::PointIdIterator         PointIdIterator;
  typedef typename Superclass::CellTraits              CellTraits;
  typedef typename Superclass::PointsContainer         PointsContainer;
  typedef typename Superclass::CellsContainer          CellsContainer;
  typedef typename Superclass::PointsContainerPointer  PointsContainerPointer;
  typedef typename Superclass::CellsContainerPointer   CellsContainerPointer;

  /** API that will be used by itkVTKPolydataWriter */
  const PointIdentifier  GetNumberOfPoints() const { return GetPoints()->size(); };
  const CellsContainerPointer  GetCells()    const { return m_Mesh->GetDualCells();  };
  const PointsContainerPointer GetPoints()   const { return m_Mesh->GetDualPoints(); };
  void SetInput( TMesh* mesh ) { m_Mesh = mesh; };

private:
  TMesh* m_Mesh;

};

}
