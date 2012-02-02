/*=========================================================================
 *
 * Walking in a Triangulation algorithm implementation
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
 *
 *   Output int I index of the cell that contain the destination point
 *          if I = -1, destination point is outside the mesh
 *          if I = -2, errors was encounter
 *   NOTE ALEX: no -1 or -2: throw an exception*
 *
 *=========================================================================*/

#ifndef __itkWalkInTriangulationFunction_h__
#define __itkWalkInTriangulationFunction_h__

#include "itkPointSet.h"
#include "itkVectorContainer.h"

#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshTraits.h"
#include "itkQuadEdgeMeshPolygonCell.h"
#include "itkQuadEdgeMeshFunctionBase.h"

#include "itkPointInCircleGeometricalPredicateFunctor.h"

//#include "vnl/vnl_det.h" 
//#include "vnl/vnl_matrix_fixed.h" 

#include <iostream>

//------------------------------------------------------------------------------
// WalkInTriangulation Functor
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

  typedef itk::VectorContainer<unsigned int,unsigned int> VectorContainerType;

public:
  TOutputType Evaluate( 
    TMeshType* myMesh,
    const PointType& myPts, 
    const CellIdentifier& myCell = -1
    )
  {

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

  if( myCell >= 0 )
    { 
    myCellIndex = myCell ;
    }
  else
    { // Default start
    myCellIterator = myCellsContainer->Begin();
    myCellIndex = myCellIterator.Index();
    }

  // -----------------------------------------------------
  // WalkInTriangulation Algorithm

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

    path->InsertElement(triangleVisitedCompter,myCellIndex);
    triangleVisitedCompter+=1;

    if( orient2d( pointB, pointQ, destination ) < 0 )
      {
      orientationTestCompter += 1;

      while( orient2d( pointC, pointQ, destination ) < 0 )
        {
        orientationTestCompter += 1;

        // r = l
        pointIdB = pointIdC;
        pointB = pointC;

        // t = neighbour( t through ql )
        if( myMesh->FindEdge( pointIdQ, pointIdC )->IsAtBorder() )
          {
          myCellIndex = -1; 
          break;
          }
        else
          {
          myOldCellIndex = myCellIndex;
          typename QuadEdgeType::DualOriginRefType rightCell = myMesh->FindEdge( pointIdQ, pointIdC )->GetRight();
          typename QuadEdgeType::DualOriginRefType leftCell = myMesh->FindEdge( pointIdQ, pointIdC )->GetLeft();
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

        // l = vertex of t, l!=q, l!=r
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

       // l = r
       pointIdC = pointIdB;
       pointC = pointB;

       // t = neighbour( t through qr )
       if( myMesh->FindEdge( pointIdQ, pointIdB )->IsAtBorder() )
         {
         myCellIndex = -1; 
         break;
         }
       else
         {
         myOldCellIndex = myCellIndex;
         typename QuadEdgeType::DualOriginRefType rightCell = myMesh->FindEdge( pointIdQ, pointIdB )->GetRight();
         typename QuadEdgeType::DualOriginRefType leftCell = myMesh->FindEdge( pointIdQ, pointIdB )->GetLeft();
         if( leftCell == myCellIndex )
           {
           myCellIndex = rightCell;
           }
         else
           {
           myCellIndex = leftCell;
           }
         }
       path->InsertElement(triangleVisitedCompter,myCellIndex);
       triangleVisitedCompter+=1;

       // r = vertex of t, r!=q, r!=l
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

  // End of initialisation step
  // Q-destination vector has B on its right and C on its left
  while( orient2d( destination, pointB, pointC ) < 0 )
    {
    orientationTestCompter += 1;

    // t = neighbour( t through rl )
    if( myMesh->FindEdge( pointIdB, pointIdC )->IsAtBorder() )
      {
      myCellIndex = -1; 
      break;
      }
    else
      {
      myOldCellIndex = myCellIndex;
      typename QuadEdgeType::DualOriginRefType rightCell = myMesh->FindEdge( pointIdB, pointIdC )->GetRight();
      typename QuadEdgeType::DualOriginRefType leftCell = myMesh->FindEdge( pointIdB, pointIdC )->GetLeft();
      if( leftCell == myCellIndex )
        {
        myCellIndex = rightCell;
        }
      else
        {
        myCellIndex = leftCell;
        }
      }
    path->InsertElement(triangleVisitedCompter,myCellIndex);
    triangleVisitedCompter+=1;

    // s = vertex of t, s!=r, s!=l
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

      // r = s
      pointIdB = pointIdA;
      pointB = pointA;
      }
    else
      {
      orientationTestCompter += 1;

      // l = s
      pointIdC = pointIdA;
      pointC = pointA;
      }
    }
  }
else
  {
  myCellIndex = -2;

  path->InsertElement(triangleVisitedCompter,myCellIndex);
  triangleVisitedCompter+=1;
  }

return path;
}
};
#endif // __itkWalkInTriangulation_h__
