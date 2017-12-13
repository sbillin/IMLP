#ifndef TEST_CLOSESTPOINTSOLVER_H
#define TEST_CLOSESTPOINTSOLVER_H

#include <cisstVector.h>
#include <cisstCommon.h>

#include "triangleClosestPointSolver.h";
#include "cisstTriangleClosestPointSolver.h";

void test_ClosestPointSolver()
{
  int numRandTrials = 10;

  vctDynamicVector<vct3> testPoints;
  testPoints.SetSize(20);

  testPoints[0].Assign(0.1, 0.5, 0.0);
  testPoints[1].Assign(1, 0.5, 0.0);
  testPoints[2].Assign(0.5, 0.2, 0.0);
  testPoints[3].Assign(-0.1, 0.5, 0.0);
  testPoints[4].Assign(-1.1, -2.0, 0.0);
  testPoints[5].Assign(-0.5, -2.0, 0.0);
  testPoints[6].Assign(0.5, -2.0, 0.0);
  testPoints[7].Assign(-0.2, 1.1, 0.0);
  testPoints[8].Assign(-0.1, 1.2, 0.0);
  testPoints[9].Assign(-0.1, 1.1, 0.0);
  testPoints[10].Assign(0.0, 2.0, 0.0);
  testPoints[11].Assign(2.0, 3.1, 0.0);
  testPoints[12].Assign(2.1, 1.0, 0.0);
  testPoints[13].Assign(0.0, 0.0, 0.0);
  testPoints[14].Assign(0.0, 1.0, 0.0);
  testPoints[15].Assign(1.0, 2.0, 0.0);
  testPoints[16].Assign(1.0, -1.0, 0.0);
  testPoints[17].Assign(2.0, -2.1, 0.0);
  testPoints[18].Assign(2.0, -3.1, 0.0);
  testPoints[19].Assign(2.0, 0.5, 0.0);


  vct3 P1, P2;
  P1.Assign(0.0, 0.0, 0.0);
  P2.Assign(0.0, 1.0, 0.0);

  vctDynamicVector<vct3> vctP3;
  vctP3.SetSize(3);
  vctP3[0].Assign(1.0, 2.0, 0.0);
  vctP3[1].Assign(1.0, -1.0, 0.0);
  vctP3[2].Assign(1.0, 0.4, 0.0);

  vctDynamicVector<vct3> vertices(3);
  vertices[0] = P1;
  vertices[1] = P2;

  vctDynamicVector<vctInt3> triangles(1);
  triangles[0][0] = 0;
  triangles[0][1] = 1;
  triangles[0][2] = 2;


  // Process

  for (unsigned int indexP3 = 0; indexP3 < vctP3.size(); indexP3++)
  {
    std::cout << "indexP3: " << indexP3 << std::endl;

    vertices[2] = vctP3[indexP3];

    for (unsigned int trial = 0; trial < numRandTrials; trial++)
    {
      // compute a random transformation of the vertices and test points
      vctDynamicVector<double> randVect(7);
      vctRandom(randVect, -1.0, 1.0);

      vct3 axis(randVect[0], randVect[1], randVect[2]);
      axis = axis.Normalized() * randVect[3] * cmnPI;
      vctRodRot3 rodriguesRotation(axis[0], axis[1], axis[2]);

      vct3 translation(randVect[4], randVect[5], randVect[6]);
      translation = translation * 10.0;

      vctFrm3 xfm(vctRot3(rodriguesRotation), translation);

      vctDynamicVector<vct3> verticesXfmd(3);
      verticesXfmd[0] = xfm * vertices[0];
      verticesXfmd[1] = xfm * vertices[1];
      verticesXfmd[2] = xfm * vertices[2];

      vct3 closestPoint1, closestPoint2;
      triangleClosestPointSolver solver1( verticesXfmd, triangles );
      cisstTriangleClosestPointSolver solver2;

      for (unsigned int indexPt = 0; indexPt < testPoints.size(); indexPt++)
      {
        vct3 testPointXfmd = xfm * testPoints[indexPt];

        solver1.FindClosestPointOnTriangle( testPointXfmd, 0, closestPoint1 );

        solver2.FindClosestPointOnTriangle(
            testPointXfmd,
            verticesXfmd[0], verticesXfmd[1], verticesXfmd[2],
            -1,
            closestPoint2 );

        if ((closestPoint1 - closestPoint2).Norm() > 1.0e-6)
        {
          std::cout << "ERROR: closest point computed for indexP3 " << indexP3 << " trial " << trial
              << " point index " << indexPt << " do not match" << std::endl
              << " closestPoint1: " << closestPoint1 << std::endl
              << " closestPoint2: " << closestPoint2 << std::endl
              << " dist1: " << (closestPoint1 - testPointXfmd).Norm() << std::endl
              << " dist2: " << (closestPoint2 - testPointXfmd).Norm() << std::endl
              << " testPointXfmd: " << testPointXfmd << std::endl
              << " verticesXfmd: " << std::endl << verticesXfmd << std::endl;

          vctDynamicVector<vct3> verticesLocal(3);
          verticesLocal[0] = solver1.triXfm[0] * verticesXfmd[0];
          verticesLocal[1] = solver1.triXfm[0] * verticesXfmd[1];
          verticesLocal[2] = solver1.triXfm[0] * verticesXfmd[2];

          std::cout << "In Local Coordinates: " << std::endl
              << " point: " << solver1.triXfm[0] * testPointXfmd << std::endl
              << " closestPoint1: " << solver1.triXfm[0] * closestPoint1 << std::endl
              << " closestPoint2: " << solver1.triXfm[0] * closestPoint2 << std::endl
              << " vertices: " << std::endl << verticesLocal
              << " P2: " << solver1.P2[0] << std::endl
              << " P2P3: " << solver1.P2P3[0] << std::endl
              << " E13: " << solver1.E13[0] << std::endl
              << " E23: " << solver1.E23[0] << std::endl
              << std::endl;
        }
      }
    }
  }
}

#endif // TEST_CLOSESTPOINTSOLVER_H
