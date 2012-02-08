/*================================================================================
 *                                                                               *
 * Walking in a Triangulation algorithm implementation                           *
 *                                                                               *
 *   Based on :                                                                  *
 *   Devillers, O. and Pion, S. and Teillaud, M. "Walking in a triangulation",   *
 *   Proceedings of the seventeenth annual symposium on Computational geometry,  *
 *   pages 106-114, 2001                                                         *
 *                                                                               *
 *   Implementation for ITK by Stéphane Ulysse Rigaud                            *
 *   IPAL (Image & Pervasive Access Lab) CNRS - A*STAR                           *
 *   Singapore                                                                   *
 *                                                                               *
 *   Input: itk::Mesh*          input mesh                                       *
 *          Mesh::PointType&    point coordinate                                 *
 *          unsigned int        starting cell id                                 *
 *                                                                               *
 *   Output: itk::vectorContainer  list of visited triangle index                *
 *                                                                               *
 *   NOTE ALEX: no -1 or -2: throw an exception                                  *
 *                                                                               *
 *===============================================================================*/

#ifndef __itkWalkInTriangulationFunction_h__
#define __itkWalkInTriangulationFunction_h__

#include "itkPointSet.h"
#include "itkVectorContainer.h"

#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshTraits.h"
#include "itkQuadEdgeMeshPolygonCell.h"
#include "itkQuadEdgeMeshFunctionBase.h"

#include "itkPointInCircleGeometricalPredicateFunctor.h"

#include <iostream>

//
// WalkInTriangulation Function
//

template<
  class    TMeshType, 
  typename TOutputType = itk::VectorContainer< unsigned int, unsigned int >::Pointer
  >
class WalkInTriangulationFunction: public itk::QuadEdgeMeshFunctionBase< TMeshType, TOutputType >
{
 
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

  typedef typename QuadEdgeType::DualOriginRefType DualOriginRefType;

  typedef itk::VectorContainer< unsigned int, unsigned int > VectorContainerType;

public:
 
  TOutputType Evaluate( 
    TMeshType* myMesh,
    const PointType& myPts, 
    const CellIdentifier& myCell = -1
    )
  {
  
	// 
  // Initialisation
  //
	 
  CellAutoPointer             myCellPointer;  
  CellsContainer              *myCellsContainer = myMesh->GetCells();
  CellsContainerIteratorType  myCellIterator;
  CellIdentifier              myCellIndex;
  CellIdentifier              myOldCellIndex;

  PointType destination = myPts;
  unsigned int triangleVisitedCompter = 0;
  unsigned int orientationTestCompter = 0;

  VectorContainerType::Pointer path = VectorContainerType::New();

  if( myCell >= 0 )
    { 
    myCellIndex = myCell ;
    }
  else // Default start
    { 
    myCellIterator = myCellsContainer->Begin();
    myCellIndex = myCellIterator.Index();
    }
  
	// 
  // WalkInTriangulation Algorithm
  //
	 
  if( myMesh->GetCell( myCellIndex, myCellPointer ) )
    { 

    PointType pointQ, pointA, pointB, pointC;
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

    path->InsertElement( triangleVisitedCompter, myCellIndex );
    triangleVisitedCompter += 1;

    if( orient2d( pointB, pointQ, destination ) < 0 )
      {
      orientationTestCompter += 1;

      while( orient2d( pointC, pointQ, destination ) < 0 )
        {
        orientationTestCompter += 1;

        pointIdB = pointIdC;
        pointB = pointC;

        if( myMesh->FindEdge( pointIdQ, pointIdC )->IsAtBorder() )
          {
          throw "Point is outside the mesh structure.";
          myCellIndex = -1; 
          break;
          }
        else
          {
          myOldCellIndex = myCellIndex;
          DualOriginRefType rightCell = myMesh->FindEdge( pointIdQ, pointIdC )->GetRight();
          DualOriginRefType leftCell = myMesh->FindEdge( pointIdQ, pointIdC )->GetLeft();
          if( leftCell == myCellIndex )
            {
            myCellIndex = rightCell;
            }
          else
            {
            myCellIndex = leftCell;
            }
          }

        path->InsertElement( triangleVisitedCompter, myCellIndex );
        triangleVisitedCompter += 1;

        if( myMesh->GetCell( myCellIndex, myCellPointer) )
          {
          PointIdIterator pointIdIterator = myCellPointer->PointIdsBegin();
          pointIdC = *pointIdIterator;
          myMesh->GetPoint( *pointIdIterator, &pointC );
          pointIdIterator++;
          while( pointIdC == pointIdB || pointIdC == pointIdQ )
            {
            if( pointIdIterator != myCellPointer->PointIdsEnd() )
              {
              pointIdC = *pointIdIterator;
              myMesh->GetPoint( *pointIdIterator, &pointC );
              pointIdIterator++;
              }
            }
          }
        }
      }
    else
      {
      orientationTestCompter += 1;
      do
        {
        orientationTestCompter += 1;
 
        pointIdC = pointIdB;
        pointC = pointB;
 
        if( myMesh->FindEdge( pointIdQ, pointIdB )->IsAtBorder() )
          {
					throw "Point is outside of the mesh structure";
          myCellIndex = -1; 
          break;
          }
        else
          {
          myOldCellIndex = myCellIndex;
          DualOriginRefType rightCell = myMesh->FindEdge( pointIdQ, pointIdB )->GetRight();
          DualOriginRefType leftCell = myMesh->FindEdge( pointIdQ, pointIdB )->GetLeft();
          if( leftCell == myCellIndex )
            {
            myCellIndex = rightCell;
            }
          else
            {
            myCellIndex = leftCell;
            }
          }
        path->InsertElement( triangleVisitedCompter, myCellIndex );
        triangleVisitedCompter += 1;
 
        if( myMesh->GetCell( myCellIndex, myCellPointer) )
          {
          PointIdIterator pointIdIterator = myCellPointer->PointIdsBegin();
          pointIdB = *pointIdIterator;
          myMesh->GetPoint( *pointIdIterator, &pointB );
          pointIdIterator++;
          while( pointIdB == pointIdC || pointIdB == pointIdQ )
            {
            if( pointIdIterator != myCellPointer->PointIdsEnd() )
              {
              pointIdB = *pointIdIterator;
              myMesh->GetPoint( *pointIdIterator, &pointB );
              pointIdIterator++;
              }
            }
          }
        }
      while ( orient2d( pointB, pointQ, destination ) > 0 );
      }

    while( orient2d( destination, pointB, pointC ) < 0 )
      {
      orientationTestCompter += 1;

      if( myMesh->FindEdge( pointIdB, pointIdC )->IsAtBorder() )
        {
				throw "Point is outside of the mesh structure";
        myCellIndex = -1; 
        break;
        }
      else
        {
        myOldCellIndex = myCellIndex;
        DualOriginRefType rightCell = myMesh->FindEdge( pointIdB, pointIdC )->GetRight();
        DualOriginRefType leftCell = myMesh->FindEdge( pointIdB, pointIdC )->GetLeft();
        if( leftCell == myCellIndex )
          {
          myCellIndex = rightCell;
          }
        else
          {
          myCellIndex = leftCell;
          }
        }
      path->InsertElement( triangleVisitedCompter, myCellIndex );
      triangleVisitedCompter += 1;

      if( myMesh->GetCell( myCellIndex, myCellPointer) )
        {
        PointIdIterator pointIdIterator = myCellPointer->PointIdsBegin();
        pointIdA = *pointIdIterator;
        myMesh->GetPoint( *pointIdIterator, &pointA );
        pointIdIterator++;
        while( pointIdA == pointIdB || pointIdA == pointIdC )
          {
          if( pointIdIterator != myCellPointer->PointIdsEnd() )
            {
            pointIdA = *pointIdIterator;
            myMesh->GetPoint( *pointIdIterator, &pointA );
            pointIdIterator++;
            }
          }
        }
      if( orient2d( pointA, pointQ, destination ) < 0 )
        {
        orientationTestCompter += 1;
        pointIdB = pointIdA;
        pointB = pointA;
        }
      else
        {
        orientationTestCompter += 1;
        pointIdC = pointIdA;
        pointC = pointA;
        }
      }
    }
  else
    {
    throw "The starting triangle does not exist";
    myCellIndex = -2;
    path->InsertElement( triangleVisitedCompter, myCellIndex );
    triangleVisitedCompter += 1;
    }

  return path;
  }
 
};

#endif // __itkWalkInTriangulationFunction_h__
