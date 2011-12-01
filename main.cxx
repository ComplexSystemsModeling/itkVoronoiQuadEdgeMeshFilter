
#include <iostream>


#include "itkMesh.h"
#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshTraits.h"
#include "itkPointSet.h"
#include "itkVTKPolyDataWriter.h"
#include "itkVTKPolyDataReader.h"


template< class TPointType >
double orientation ( TPointType p , TPointType q , TPointType r )
{
	
	typedef typename TPointType::VectorType VectorType;
	
	VectorType pq;
	pq[0] = q[0] - p[0] ;
	pq[1] = q[1] - p[1] ;

	VectorType pr;
	pr[0] = r[0] - p[0] ;
	pr[1] = r[1] - p[1] ;
	
	double scalar = pq[0] * pr[0] - pq[1] * pr[1] ; 
	
	return scalar;
}


int main(int argc, char * argv[] )
{
	
	//-----------------------------------
	// Typedef
	//-----------------------------------
	
	typedef double			 PixelType;
	const		unsigned int Dimension = 3;
	
	typedef itk::QuadEdgeMeshTraits< PixelType, Dimension, PixelType, PixelType, PixelType, PixelType > QEMTraits;
	typedef itk::QuadEdgeMesh< PixelType, Dimension, QEMTraits >																				QEMeshType;

	typedef QEMeshType::PointType					PointType;
	typedef QEMeshType::PointIdentifier		PointId;
	typedef QEMeshType::VectorType				VectorType;
	typedef QEMeshType::PointIdList				PointIdList;
	
	typedef QEMeshType::CellType												CellType;					// abstract
	typedef QEMeshType::CellIdentifier									CellId;
  typedef CellType::CellAutoPointer										CellAutoPointer;	// abstract	
	typedef QEMeshType::CellsContainerIterator					CellsContainerIteratorType;
  
	typedef QEMeshType::PixelType												RealType;
	typedef QEMeshType::CellTraits											CellTraits;
	typedef itk::CellInterface< RealType, CellTraits >	CellInterfaceType;
  typedef itk::TriangleCell< CellInterfaceType >			TriangleCellType;
	typedef TriangleCellType::PointIdIterator						PointIdIterator;
	
	typedef itk::VTKPolyDataWriter< QEMeshType >				MeshWriterType;
	typedef itk::VTKPolyDataReader< QEMeshType >				MeshReaderType;

	
	//-----------------------------------
	// Test Mesh 
	//-----------------------------------
	
	QEMeshType::Pointer myMesh = QEMeshType::New();
	MeshWriterType::Pointer writer = MeshWriterType::New();
	MeshReaderType::Pointer reader = MeshReaderType::New();
	
	PointType localPoint;
	PointId pId;
	PointIdList IdList (3);

	// SetMeshPoint
	for (unsigned int i = 0; i<9; i++) {
		pId = i;
		
		localPoint[0] = (i % 3) * 10 ;
		localPoint[1] = floor( i / 3 ) * 10;
		
		myMesh->SetPoint( pId, localPoint );
	}
	
	// SetMesh Edge & Face
	myMesh->AddFaceTriangle(0, 1, 4);
	myMesh->AddFaceTriangle(4, 3, 0);
	myMesh->AddFaceTriangle(1, 2, 5);
	myMesh->AddFaceTriangle(5, 4, 1);
	myMesh->AddFaceTriangle(3, 4, 7);
	myMesh->AddFaceTriangle(7, 6, 3);
	myMesh->AddFaceTriangle(4, 5, 8);
	myMesh->AddFaceTriangle(8, 7, 4);
	
	writer->SetFileName("./myMesh.vtk");
	writer->SetInput(myMesh);
	writer->Update();
	
	reader->SetFileName("./myMesh.vtk");
	reader->Update();
	myMesh = reader->GetOutput();
	
	// Test point
	PointType seekedPoint;
	seekedPoint[0] = 15;
	seekedPoint[1] = 18;
	
	//-----------------------------------
	// Walk a Triangulation
	// algorithm
	// 
	//	Get the first cell ID in T
	//	While the seeked point is not inside T
	//	do
	//		Go through all the vertices of T
	//		Calculate the barycentre P
	//		Go through all the edge ER of the cell
	//		Until
	//			orientation(PQR) && orientation(PQE) >=0
	//		Get the ID of the cell that share the edge ER
	//		Save the ID in T
	//	end while
	//		
	//-----------------------------------
	
	// Initialisation
	bool win			 = false, 
	     edgeFound = false;
	double scalarA, 
	       scalarB;
	
	VectorType directionVector, edgeVector;
	
	PointType barycentre;
	PointType mpa, mpb,	mpc, epa, epb;
	PointId		e1,	e2;
	
	CellAutoPointer  cellIterator;  // Get the first CellId of the mesh
	CellsContainerIteratorType CellIterator = myMesh->GetCells()->Begin();
	
	while( !win )  {
		
		// Bertrand's InCirclePredicateFunctor
		// ERR can be in circle but not in triangle...
		// win = IsInside( TriangleId t, PointId p ) 
		
		if (!win ) {
	 	 
		
			if( myMesh->GetCell( CellIterator.Index(), cellIterator ) )   // Test if the selected Cell id exist
			{
				
				// TODO => Automation regardless the Cell type
				// WARNING => Type cell is Polygone and not Triangle ??
				if( cellIterator->GetType() == 2 )													// Test type of Cell (triangle, squarre, polygone ... )
				{ 
					
					// Calculate the barycentre of the current cell
					// TODO => Upgrade for any type of cell
					PointIdIterator pointIdIterator = CellIterator.Value()->PointIdsBegin();
					myMesh->GetPoint( *pointIdIterator, &mpa );
					pointIdIterator++;
					myMesh->GetPoint( *pointIdIterator, &mpb );
					pointIdIterator++;
					myMesh->GetPoint( *pointIdIterator, &mpc );
					
					barycentre[0] = ( mpa[0] + mpb[0] + mpc[0] ) / 3;
					barycentre[1] = ( mpa[1] + mpb[1] + mpc[1] ) / 3;

				}
				else
				{
					std::cout<<" err - this is not a triangle \n";
				}
			}
			else
			{
				std::cout<<" err - the cell ID was not found in the container \n";
			}
			
			// Determined the Cell edge that cross the direction we need to go
			// Get the id of the first point in the cell
			PointIdIterator pointIdIterator = CellIterator.Value()->PointIdsBegin();
			myMesh->GetPoint( *pointIdIterator, &epa );
			pointIdIterator++;

			while ( !edgeFound )  {

				// Get the id of the next point in the cell
				myMesh->GetPoint( *pointIdIterator, &epb );
				
				// Calculate the orientation test
				// if the 2 scalar product of the direction vector given by the barycentre and the seekedPoint
				// and the 2 vectors given by the barycentre and 2 vertices of an edge of the cell
				// are positive, then the edge given by the 2 current vertices is the edge we need to cross
				scalarA = orientation(barycentre, seekedPoint, epa) ;
				scalarB = orientation(barycentre, seekedPoint, epb) ;
				
				if (scalarA >= 0 && scalarB >= 0) {
					// We stock the 2 vertices id
					e2 = *pointIdIterator;
					e1 = *--pointIdIterator;
					
					// We set the flag to true
					edgeFound = true;
				}
				else 
				{
					// if the edge is not found
					// we go to the next edge
					epa = epb;
					if (pointIdIterator != CellIterator.Value()->PointIdsEnd() ) { 
						pointIdIterator++;
					}
					else {
						// if we arrived at the end of the vertices list, we go back at begin
						pointIdIterator = CellIterator.Value()->PointIdsBegin();
					}
				}
									
			}
	 
			// using the saved id of the selected vertices => e1 and e2
			// and the current cell id => CellIterator.Index() & CellIterator.Value()
			// we go to the cell that contain this edge and that is not the current cell we are.
			// Use the CellIterator.GoTo() ?
			
		}
	}

	
  return EXIT_SUCCESS;

}

