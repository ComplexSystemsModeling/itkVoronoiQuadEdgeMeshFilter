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
 *   http://www.ipal.cnrs.fr                                                     * 
 *                                                                               *
 *   Input: itk::Mesh*          input mesh                                       *
 *          Mesh::PointType&    point coordinate                                 *
 *          unsigned int        starting cell id                                 *
 *                                                                               *
 *   Output: itk::vectorContainer  list of visited triangle index                *
 *                                                                               *
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

namespace itk
{

//
// WalkInTriangulation Function
//

template<
  class    TMeshType, 
  typename TOutputType = itk::VectorContainer< unsigned int, int >::Pointer
  >
class WalkInTriangulationFunction: public QuadEdgeMeshFunctionBase< TMeshType, TOutputType >
{

public:

  typedef WalkInTriangulationFunction                         Self;
  typedef QuadEdgeMeshFunctionBase< TMeshType, TOutputType >  Superclass;
  typedef SmartPointer< Self >                                Pointer;
  typedef SmartPointer< const Self >                          ConstPointer;

  itkNewMacro( WalkInTriangulationFunction );
  itkTypeMacro( WalkInTriangulationFunction, QuadEdgeMeshFunctionBase );

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

  typedef VectorContainer< unsigned int, int > VectorContainerType;

  TOutputType Evaluate( 
    TMeshType* myMesh,
    const PointType& myPts, 
    const CellIdentifier& myCell = 0
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

  if( myCell > 0 )
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

    if( OrientationTest( pointB, pointQ, destination ) < 0 )
      {
      orientationTestCompter += 1;

      while( OrientationTest( pointC, pointQ, destination ) < 0 )
        {
        orientationTestCompter += 1;

        pointIdB = pointIdC;
        pointB = pointC;

        if( myMesh->FindEdge( pointIdQ, pointIdC )->IsAtBorder() )
          {
          throw -1;
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
	  throw -1;
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
      while ( OrientationTest( pointB, pointQ, destination ) > 0 );
      }

    while( OrientationTest( destination, pointB, pointC ) < 0 )
      {
      orientationTestCompter += 1;

      if( myMesh->FindEdge( pointIdB, pointIdC )->IsAtBorder() )
        {
	throw -1;
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
      if( OrientationTest( pointA, pointQ, destination ) < 0 )
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
    throw -2;
    myCellIndex = -2;
    path->InsertElement( triangleVisitedCompter, myCellIndex );
    triangleVisitedCompter += 1;
    }

  return path;
  }

 protected:

  WalkInTriangulationFunction()
  {}

  ~WalkInTriangulationFunction()
  {}

 private:

  WalkInTriangulationFunction( const Self & ); // purposely not implemented
  void operator=( const Self & );              // purposely not implemented

};

} // namespace itk

#endif // __itkWalkInTriangulationFunction_h__
