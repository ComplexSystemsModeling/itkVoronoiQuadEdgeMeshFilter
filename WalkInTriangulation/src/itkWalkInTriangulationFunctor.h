/*=========================================================================
 *
 *	 Walking in a Triangulation algorithm implementation
 *   Based on :
 *
 *   Devillers, O. and Pion, S. and Teillaud, M. "Walking in a triangulation",
 *   Proceedings of the seventeenth annual symposium on Computational geometry,
 *   pages 106-114, 2001
 *
 *   Implementation for ITK by Stéphane Ulysse Rigaud
 *   IPAL (Image & Pervasive Access Lab) CNRS - A*STAR
 *   Singapore
 *
 *   Input double  X coordinate of the destination point in the mesh
 *         double  Y coordinate of the destination point in the mesh
 *         unsigned int I starting cell id
 *         bool debug display
 *
 *   Output int I index of the cell that contain the destination point
 *          if I = -1, destination point is outside the mesh
 *          if I = -2, errors was encounter
 *
 *   TODO => Validation test
 *        => Functorise the implementation
 *
 *=========================================================================*/

#ifndef __itkWalkInTriangulation_h__
#define __itkWalkInTriangulation_h__

#include "itkPointSet.h"
#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshTraits.h"
#include "itkQuadEdgeMeshPolygonCell.h"
#include "itkVectorContainer.h"
#include "itkVTKPolyDataWriter.h"

#include "vnl/vnl_det.h" 
#include "vnl/vnl_matrix_fixed.h" 

#include <iostream>

//------------------------------------------------------------------------------
// Skewchuck code
//
extern "C"
{
  double orient2d(double* pa, double* pb, double* pc);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// The test functor to map skewchuck code to ITK's API
//
template< typename TPointType >
double
orientation2d( const TPointType& TrianglePoint1,
							 const TPointType& TrianglePoint2,
							 const TPointType& PointToTest   ) 
{	
  double * pa = new double[2];
  double * pb = new double[2];
  double * pc = new double[2];
  pa[0] = TrianglePoint1[0];
  pa[1] = TrianglePoint1[1];
  pb[0] = TrianglePoint2[0];
  pb[1] = TrianglePoint2[1];
  pc[0] = PointToTest[0];
  pc[1] = PointToTest[1];
	
  // orientation test
  double orientation = orient2d( pa, pb, pc );
	return( orientation >= 0 ? 1 : -1 );
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// WalkInTriangulation Functor
//
template<
typename TMeshType, 
typename TMeshPointerType, 
typename TCellIdentifier, 
typename TPointType >
itk::VectorContainer<unsigned int,unsigned int>::Pointer
WalkInTriangulation( TMeshPointerType myMesh, 
										 const TPointType& myPts, 
										 const TCellIdentifier& myCell = -1,
										 bool Debug = false                   )
{
	
	// -----------------------------------------------------
	// Typedef
	
	typedef typename TMeshType::PointType              PointType;
	typedef typename TMeshType::CellType               CellType;
	typedef typename TMeshType::PointIdentifier        PointIdentifier;
	typedef typename TMeshType::CellIdentifier         CellIdentifier;
	typedef typename TMeshType::PointsContainer        PointsContainer;
	typedef typename TMeshType::CellsContainer         CellsContainer;
	typedef typename TMeshType::PointIdList            PointIdList;
	typedef typename TMeshType::QEType                 QuadEdgeType;
	typedef typename TMeshType::CellsContainerIterator CellsContainerIteratorType;
	
	typedef typename CellType::PointIdConstIterator     PointIdConstIterator;
	typedef typename CellType::PointIdIterator          PointIdIterator;
	typedef typename CellType::CellAutoPointer          CellAutoPointer;
	
	typedef itk::VectorContainer<unsigned int,unsigned int> VectorContainerType;
	
	// -----------------------------------------------------
	// Initialisation
	
	CellAutoPointer             myCellPointer;  
	CellsContainer              *myCellsContainer = myMesh->GetCells();
	CellsContainerIteratorType  myCellIterator;
	CellIdentifier              myCellIndex;
	CellIdentifier              myOldCellIndex;
	
	PointType destination = myPts;	
	unsigned int triangleVisitedCompter = 0;
	unsigned int orientationTestCompter = 0;
	
	VectorContainerType::Pointer path = VectorContainerType::New();
	
	if ( myCell >= 0 ) { 
		myCellIndex = myCell ;
	}
	else { // Default start
		myCellIterator = myCellsContainer->Begin();
		myCellIndex = myCellIterator.Index();
	}
	
	// -----------------------------------------------------
	// WalkInTriangulation Algorithm
	
	if (myMesh->GetCell( myCellIndex, myCellPointer)) { 
		
		PointType	pointQ, pointA, pointB, pointC;
		PointIdIterator pointIdIterator = myCellPointer->PointIdsBegin();
		
		PointIdentifier pointIdQ = *pointIdIterator;
		myMesh->GetPoint( *pointIdIterator, &pointQ );
		pointIdIterator++;
		PointIdentifier pointIdB = *pointIdIterator;
		myMesh->GetPoint( *pointIdIterator, &pointB );
		pointIdIterator++;
		PointIdentifier pointIdC = *pointIdIterator;
		myMesh->GetPoint( *pointIdIterator, &pointC );
		PointIdentifier pointIdA;
		
		path->InsertElement(triangleVisitedCompter,myCellIndex);
		triangleVisitedCompter+=1;
		
		if ( orientation2d( pointB, pointQ, destination ) < 0 )  {
			orientationTestCompter += 1;
			
			while ( orientation2d( pointC, pointQ, destination ) < 0 ) {
				orientationTestCompter += 1;
				
				// r = l
				pointIdB = pointIdC;
				pointB = pointC;
				
				// t = neighbour( t through ql )
				if (myMesh->FindEdge( pointIdQ, pointIdC )->IsAtBorder()) {
					//std::cout<<"This is a border edge, the point is out of the mesh\n";
					myCellIndex = -1; 
					
					break;
				}
				else {
					myOldCellIndex = myCellIndex;
					typename QuadEdgeType::DualOriginRefType rightCell = myMesh->FindEdge( pointIdQ, pointIdC )->GetRight();
					typename QuadEdgeType::DualOriginRefType leftCell = myMesh->FindEdge( pointIdQ, pointIdC )->GetLeft();
					if (leftCell == myCellIndex) {
						myCellIndex = rightCell;
					}
					else {
						myCellIndex = leftCell;
					}		
					
					//std::cout<<"We go from cell "<<myOldCellIndex<<" to cell "<<myCellIndex<<" "
					//<<"trough the edge "<<pointIdQ<<" - "<<pointIdC<<"\n";
				}
				path->InsertElement(triangleVisitedCompter,myCellIndex);
				triangleVisitedCompter+=1;
				
				// l = vertex of t, l!=q, l!=r
				if( myMesh->GetCell( myCellIndex, myCellPointer) ) {
					PointIdIterator pointIdIterator = myCellPointer->PointIdsBegin();
					pointIdC = *pointIdIterator;
					myMesh->GetPoint( *pointIdIterator, &pointC );
					pointIdIterator++;
					while (pointIdC == pointIdB || pointIdC == pointIdQ) {
						if (pointIdIterator != myCellPointer->PointIdsEnd()) {
							pointIdC = *pointIdIterator;
							myMesh->GetPoint( *pointIdIterator, &pointC );
							pointIdIterator++;
						}
					}
				}
			}
		}
		else {
			orientationTestCompter += 1;
			
			do {
				orientationTestCompter += 1;
				
				// l = r
				pointIdC = pointIdB;
				pointC = pointB;
				
				
				// t = neighbour( t through qr )
				if (myMesh->FindEdge( pointIdQ, pointIdB )->IsAtBorder()) {
					//std::cout<<"This is a border edge, the point is out of the mesh\n";
					myCellIndex = -1; 
					break;
				}
				else {
					myOldCellIndex = myCellIndex;
					typename QuadEdgeType::DualOriginRefType rightCell = myMesh->FindEdge( pointIdQ, pointIdB )->GetRight();
					typename QuadEdgeType::DualOriginRefType leftCell = myMesh->FindEdge( pointIdQ, pointIdB )->GetLeft();
					if (leftCell == myCellIndex) {
						myCellIndex = rightCell;
					}
					else {
						myCellIndex = leftCell;
					}		
					
					//std::cout<<"We go from cell "<<myOldCellIndex<<" to cell "<<myCellIndex<<" "
					//<<"trough the edge "<<pointIdQ<<" - "<<pointIdB<<"\n";
				}
				path->InsertElement(triangleVisitedCompter,myCellIndex);
				triangleVisitedCompter+=1;
				
				// r = vertex of t, r!=q, r!=l
				if( myMesh->GetCell( myCellIndex, myCellPointer) ) {
					PointIdIterator pointIdIterator = myCellPointer->PointIdsBegin();
					pointIdB = *pointIdIterator;
					myMesh->GetPoint( *pointIdIterator, &pointB );
					pointIdIterator++;
					while (pointIdB == pointIdC || pointIdB == pointIdQ) {
						if (pointIdIterator != myCellPointer->PointIdsEnd()) {
							pointIdB = *pointIdIterator;
							myMesh->GetPoint( *pointIdIterator, &pointB );
							pointIdIterator++;
						}
					}
				}
			}
			while ( orientation2d( pointB, pointQ, destination ) > 0 );			
		}
		
		// End of initialisation step
		// Q-destination vector has B on its right and C on its left
		while ( orientation2d( destination, pointB, pointC ) < 0 ) {
			orientationTestCompter += 1;
			
			// t = neighbour( t through rl )
			if (myMesh->FindEdge( pointIdB, pointIdC )->IsAtBorder()) {
				//std::cout<<"This is a border edge, the point is out of the mesh\n";
				myCellIndex = -1; 
				break;
			}
			else {
				myOldCellIndex = myCellIndex;
				typename QuadEdgeType::DualOriginRefType rightCell = myMesh->FindEdge( pointIdB, pointIdC )->GetRight();
				typename QuadEdgeType::DualOriginRefType leftCell = myMesh->FindEdge( pointIdB, pointIdC )->GetLeft();
				if (leftCell == myCellIndex) {
					myCellIndex = rightCell;
				}
				else {
					myCellIndex = leftCell;
				}
				
				//std::cout<<"We go from cell "<<myOldCellIndex<<" to cell "<<myCellIndex<<" "
				//<<"trough the edge "<<pointIdB<<" - "<<pointIdC<<"\n";
			}
			path->InsertElement(triangleVisitedCompter,myCellIndex);
			triangleVisitedCompter+=1;
			
			// s = vertex of t, s!=r, s!=l
			if( myMesh->GetCell( myCellIndex, myCellPointer) ) {
				PointIdIterator pointIdIterator = myCellPointer->PointIdsBegin();
				pointIdA = *pointIdIterator;
				myMesh->GetPoint( *pointIdIterator, &pointA );
				pointIdIterator++;
				while (pointIdA == pointIdB || pointIdA == pointIdC) {
					if (pointIdIterator != myCellPointer->PointIdsEnd()) {
						pointIdA = *pointIdIterator;
						myMesh->GetPoint( *pointIdIterator, &pointA );
						pointIdIterator++;
					}
				}
			}
			if ( orientation2d( pointA, pointQ, destination ) < 0 ) {
				orientationTestCompter += 1;
				
				// r = s
				pointIdB = pointIdA;
				pointB = pointA;
			}
			else {
				orientationTestCompter += 1;
				
				// l = s
				pointIdC = pointIdA;
				pointC = pointA;
			}
		}
		
		// destination reached
		//std::cout<<"We arrived at destination : "<<myCellIndex<<"\n";
		//std::cout<< orientationTestCompter <<" orientation test were made \n";
		//std::cout<< triangleVisitedCompter <<" triangle were visited \n\n";
	}
	else {
		std::cerr<<"Wrong initial cell id - the id specified does not exist"<<std::endl;
		myCellIndex = -2;
		
		path->InsertElement(triangleVisitedCompter,myCellIndex);
		triangleVisitedCompter+=1;
	}
	
	return path;
}

#endif // __itkWalkInTriangulation_h__
