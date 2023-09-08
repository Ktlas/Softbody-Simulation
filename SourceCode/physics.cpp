/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"

void computeHooksForce(point pointAPos, point pointBPos, double k, double restLength, point& resultForce)
{
    point Lvect;
    double Llen;
    pDIFFERENCE(pointAPos, pointBPos, Lvect);
    //printf("%f, %f, %f \n", pointAPos.x, pointAPos.y, pointAPos.z);
    //printf("%f, %f, %f \n", pointBPos.x, pointBPos.y, pointBPos.z);
    //printf("%f, %f, %f \n", Lvect.x, Lvect.y, Lvect.z);
    //printf("Rest Length: %f\n", restLength);
    Llen = pLength(Lvect);
    //printf("Current Length: %f\n", Llen);
    // hook's law
    double springLengthDiff = Llen - restLength;
    //printf("Length Diff: %f\n", springLengthDiff);
    resultForce = pMultiply(Lvect, -k * springLengthDiff / Llen);
    //printf("%f, %f, %f \n", resultForce.x, resultForce.y, resultForce.z);
    //system("pause");
}

void computeDampingForce(point pointAPos, point pointBPos, point pointAVelocity, point pointBVelocity, double d, point& resultForce)
{
    point Lvect;
    double Llen;
    pDIFFERENCE(pointBPos, pointAPos, Lvect);
    Llen = pLength(Lvect);
    // damping
    point velocityDiff;
    pDIFFERENCE(pointAVelocity, pointBVelocity, velocityDiff);
    double tempResult = pDot(velocityDiff, Lvect);
    resultForce = pMultiply(Lvect, -d * tempResult / Llen / Llen);
    //printf("%f, %f, %f \n", resultForce.x, resultForce.y, resultForce.z);
}

void computeIndividualSpringForce(int srcx, int srcy, int srcz, int destx, int desty, int destz, double restLength, struct world* jello, point*** f)
{
    point posA = jello->p[srcx][srcy][srcz];
    point posB = jello->p[destx][desty][destz];
    point vA = jello->v[srcx][srcy][srcz];
    point vB = jello->v[destx][desty][destz];

    point hookForce, dampingForce, srcPos, destPos, srcV, destV;
    computeHooksForce(posA, posB, jello->kElastic, restLength, hookForce);
    computeDampingForce(posA, posB, vA, vB, jello->dElastic, dampingForce);

    point totalForce, negativeTotalForce;
    totalForce = pSum(hookForce, dampingForce);
    negativeTotalForce = pNegate(totalForce);
    
    f[srcx][srcy][srcz] = pSum(f[srcx][srcy][srcz], totalForce);
    f[destx][desty][destz] = pSum(f[destx][desty][destz], negativeTotalForce);
}

void computeSpringForce(world* jello, point*** f) {
    double standardRestLength = 1.0 / 7.0;
    double diagRestLength = 1.0 / 7.0 * sqrt(2);
    double longDiagRestLength = 1.0 / 7.0 * sqrt(3);
    double doubleDiagRestLength = diagRestLength * 2;
    int counter = 0;
    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            for (int k = 0; k < 8; k++)
            {
                // structural spring
                // six directions
                if (i > 0) {
                    computeIndividualSpringForce(i, j, k, i - 1, j, k, standardRestLength, jello, f);
                }

                /*if (i < 7) {
                    computeIndividualSpringForce(i, j, k, i + 1, j, k, standardRestLength, jello, f);
                }*/

                if (j > 0) {
                    computeIndividualSpringForce(i, j, k, i, j - 1, k, standardRestLength, jello, f);
                }

                /*if (j < 7) {
                    computeIndividualSpringForce(i, j, k, i, j + 1, k, standardRestLength, jello, f);
                }*/

                if (k > 0) {
                    computeIndividualSpringForce(i, j, k, i, j, k - 1, standardRestLength, jello, f);
                }

                /*if (k < 7) {
                    computeIndividualSpringForce(i, j, k, i, j, k + 1, standardRestLength, jello, f);
                }*/

                // shear spring
                // 20 directions
                // 4 each on each plane, which means there are 12 on x, y, z plane.
                // 8 for diagonal direction
                // x
                if (j > 0 && k > 0) {
                    computeIndividualSpringForce(i, j, k, i, j - 1, k - 1, diagRestLength, jello, f);
                }

                if (j > 0 && k < 7) {
                    computeIndividualSpringForce(i, j, k, i, j - 1, k + 1, diagRestLength, jello, f);
                }

                /*if (j < 7 && k > 0) {
                    computeIndividualSpringForce(i, j, k, i, j + 1, k - 1, diagRestLength, jello, f);
                }

                if (j < 7 && k < 7) {
                    computeIndividualSpringForce(i, j, k, i, j + 1, k + 1, diagRestLength, jello, f);
                }*/

                //y
                if (i > 0 && k > 0) {
                    computeIndividualSpringForce(i, j, k, i - 1, j, k - 1, diagRestLength, jello, f);
                }

                if (i > 0 && k < 7) {
                    computeIndividualSpringForce(i, j, k, i - 1, j, k + 1, diagRestLength, jello, f);
                }

                /*if (i < 7 && k > 0) {
                    computeIndividualSpringForce(i, j, k, i + 1, j, k - 1, diagRestLength, jello, f);
                }

                if (i < 7 && k < 7) {
                    computeIndividualSpringForce(i, j, k, i + 1, j, k + 1, diagRestLength, jello, f);
                }*/

                //z
                if (i > 0 && j > 0) {
                    computeIndividualSpringForce(i, j, k, i - 1, j - 1, k, diagRestLength, jello, f);
                }

                if (i > 0 && j < 7) {
                    computeIndividualSpringForce(i, j, k, i - 1, j + 1, k, diagRestLength, jello, f);
                }

                /*if (i < 7 && j > 0) {
                    computeIndividualSpringForce(i, j, k, i + 1, j - 1, k, diagRestLength, jello, f);
                }

                if (i < 7 && j < 7) {
                    computeIndividualSpringForce(i, j, k, i + 1, j + 1, k, diagRestLength, jello, f);
                }*/

                // diag
                if (i > 0 && j > 0 && k > 0) {
                    computeIndividualSpringForce(i, j, k, i - 1, j - 1, k - 1, longDiagRestLength, jello, f);
                }

                if (i > 0 && j > 0 && k < 7) {
                    computeIndividualSpringForce(i, j, k, i - 1, j - 1, k + 1, longDiagRestLength, jello, f);
                }

                if (i > 0 && j < 7 && k > 0) {
                    computeIndividualSpringForce(i, j, k, i - 1, j + 1, k - 1, longDiagRestLength, jello, f);
                }

                if (i > 0 && j < 7 && k < 7) {
                    computeIndividualSpringForce(i, j, k, i - 1, j + 1, k + 1, longDiagRestLength, jello, f);
                }

                /*if (i < 7 && j > 0 && k > 0) {
                    computeIndividualSpringForce(i, j, k, i + 1, j - 1, k - 1, longDiagRestLength, jello, f);
                }

                if (i < 7 && j > 0 && k < 7) {
                    computeIndividualSpringForce(i, j, k, i + 1, j - 1, k + 1, longDiagRestLength, jello, f);
                }

                if (i < 7 && j < 7 && k > 0) {
                    computeIndividualSpringForce(i, j, k, i + 1, j + 1, k - 1, longDiagRestLength, jello, f);
                }

                if (i < 7 && j < 7 && k < 7) {
                    computeIndividualSpringForce(i, j, k, i + 1, j + 1, k + 1, longDiagRestLength, jello, f);
                }*/

                // Bend
                // six directions
                if (i > 1) {
                    computeIndividualSpringForce(i, j, k, i - 2, j, k, doubleDiagRestLength, jello, f);
                }

                /*if (i < 6) {
                    computeIndividualSpringForce(i, j, k, i + 2, j, k, doubleDiagRestLength, jello, f);
                }*/

                if (j > 1) {
                    computeIndividualSpringForce(i, j, k, i, j - 2, k, doubleDiagRestLength, jello, f);
                }

                /*if (j < 6) {
                    computeIndividualSpringForce(i, j, k, i, j + 2, k, doubleDiagRestLength, jello, f);
                }*/

                if (k > 1) {
                    computeIndividualSpringForce(i, j, k, i, j, k - 2, doubleDiagRestLength, jello, f);
                }

                /*if (k < 6) {
                    computeIndividualSpringForce(i, j, k, i, j, k + 2, doubleDiagRestLength, jello, f);
                }*/
            }
        }
    }
}

void computeCollisionForce(world* jello, point*** f)
{
    //do collision detection with bounding box
    //printf("%f\n", jello->kCollision);
    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            for (int k = 0; k < 8; k++)
            {
                //printf("%f, %f, %f \n", jello->p[i][j][k].x, jello->p[i][j][k].y, jello->p[i][j][k].z);
                point pos = jello->p[i][j][k];
                point velocity = jello->v[i][j][k];
                point staticVelocity;
                pZeros(staticVelocity);
                // x = 2
                if (pos.x > 2) {
                    point contactPoint;
                    double panelty;
                    contactPoint.x = 2;
                    contactPoint.y = pos.y;
                    contactPoint.z = pos.z;
                    panelty = pos.x - 2;

                    point hooksForce, dampingForce;
                    computeHooksForce(pos, contactPoint, jello->kCollision * panelty, 0.0, hooksForce);
                    computeDampingForce(pos, contactPoint, velocity, staticVelocity, jello->dCollision * panelty, dampingForce);

                    f[i][j][k] = pSum(f[i][j][k], hooksForce);
                    f[i][j][k] = pSum(f[i][j][k], dampingForce);
                }
                // x = -2
                if (pos.x < -2) {
                    point contactPoint;
                    double panelty;
                    contactPoint.x = -2;
                    contactPoint.y = pos.y;
                    contactPoint.z = pos.z;
                    panelty = -(pos.x + 2);

                    point hooksForce, dampingForce;
                    computeHooksForce(pos, contactPoint, jello->kCollision * panelty, 0.0, hooksForce);
                    computeDampingForce(pos, contactPoint, velocity, staticVelocity, jello->dCollision * panelty, dampingForce);

                    f[i][j][k] = pSum(f[i][j][k], hooksForce);
                    f[i][j][k] = pSum(f[i][j][k], dampingForce);
                }
                // y = 2
                if (pos.y > 2) {
                    point contactPoint;
                    double panelty;
                    contactPoint.x = pos.x;
                    contactPoint.y = 2;
                    contactPoint.z = pos.z;
                    panelty = pos.y - 2;

                    point hooksForce, dampingForce;
                    computeHooksForce(pos, contactPoint, jello->kCollision * panelty, 0.0, hooksForce);
                    computeDampingForce(pos, contactPoint, velocity, staticVelocity, jello->dCollision * panelty, dampingForce);

                    f[i][j][k] = pSum(f[i][j][k], hooksForce);
                    f[i][j][k] = pSum(f[i][j][k], dampingForce);
                }
                // y = -2
                if (pos.y < -2) {
                    point contactPoint;
                    double panelty;
                    contactPoint.x = pos.x;
                    contactPoint.y = -2;
                    contactPoint.z = pos.z;
                    panelty = -(pos.y + 2);

                    point hooksForce, dampingForce;
                    computeHooksForce(pos, contactPoint, jello->kCollision * panelty, 0.0, hooksForce);
                    computeDampingForce(pos, contactPoint, velocity, staticVelocity, jello->dCollision * panelty, dampingForce);

                    f[i][j][k] = pSum(f[i][j][k], hooksForce);
                    f[i][j][k] = pSum(f[i][j][k], dampingForce);
                }
                // z = 2
                if (pos.z > 2) {
                    point contactPoint;
                    double panelty;
                    contactPoint.x = pos.x;
                    contactPoint.y = pos.y;
                    contactPoint.z = 2;
                    panelty = pos.z - 2;

                    point hooksForce, dampingForce;
                    computeHooksForce(pos, contactPoint, jello->kCollision * panelty, 0.0, hooksForce);
                    computeDampingForce(pos, contactPoint, velocity, staticVelocity, jello->dCollision * panelty, dampingForce);

                    f[i][j][k] = pSum(f[i][j][k], hooksForce);
                    f[i][j][k] = pSum(f[i][j][k], dampingForce);
                }
                // z = -2
                if (pos.z < -2) {
                    point contactPoint;
                    double panelty;
                    contactPoint.x = pos.x;
                    contactPoint.y = pos.y;
                    contactPoint.z = -2;
                    panelty = -(pos.z + 2);

                    point hooksForce, dampingForce;
                    computeHooksForce(pos, contactPoint, jello->kCollision * panelty, 0.0, hooksForce);
                    computeDampingForce(pos, contactPoint, velocity, staticVelocity, jello->dCollision * panelty, dampingForce);

                    f[i][j][k] = pSum(f[i][j][k], hooksForce);
                    f[i][j][k] = pSum(f[i][j][k], dampingForce);
                }
            }
        }
    }
    //do collision detection with arbitary plane
    if (jello->incPlanePresent == 1) {
        for (int i = 0; i < 8; i++)
        {
            for (int j = 0; j < 8; j++)
            {
                for (int k = 0; k < 8; k++)
                {
                    point pos = jello->p[i][j][k];
                    point velocity = jello->v[i][j][k];
                    point staticVelocity;
                    pZeros(staticVelocity);
                    double upper = jello->a * pos.x + jello->b * pos.y + jello->c * pos.z + jello->d;
                    double lower = jello->a * jello->a + jello->b * jello->b + jello->c * jello->c;
                    if (upper > 0)
                    {
                        double dist = upper / sqrt(lower);
                        point planeNormalRaw;
                        planeNormalRaw.x = jello->a;
                        planeNormalRaw.y = jello->b;
                        planeNormalRaw.z = jello->c;

                        point planeNormal = pNormalize(planeNormalRaw);
                        point temp = pMultiply(planeNormal, dist);
                        point contactPoint = pDifference(pos, temp);

                        point hooksForce, dampingForce;
                        computeHooksForce(pos, contactPoint, jello->kCollision * dist, 0.0, hooksForce);
                        computeDampingForce(pos, contactPoint, velocity, staticVelocity, jello->dCollision * dist, dampingForce);

                        f[i][j][k] = pSum(f[i][j][k], hooksForce);
                        f[i][j][k] = pSum(f[i][j][k], dampingForce);
                    }
                }
            }
        }
    }
}

void computeForceFieldForce(world* jello, point*** f)
{
    //interpolate force field force at each position
    if (jello->resolution == 0) {
        return;
    }
    //printf("%d\n", jello->resolution);
    double forceFieldUnit = 4.0 / double(jello->resolution - 1);
    // in input, force field have a fixed offset of -2, -2, -2
    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            for (int k = 0; k < 8; k++)
            {
                //find which cell current point is
                point pos = jello->p[i][j][k];
                int forceFieldCellX = int((pos.x + 2) / forceFieldUnit);
                int forceFieldCellY = int((pos.y + 2) / forceFieldUnit);
                int forceFieldCellZ = int((pos.z + 2) / forceFieldUnit);
                //clamp to force field boundary (cell's lower boundary)
                forceFieldCellX = forceFieldCellX < 0 ? 0 : forceFieldCellX;
                forceFieldCellX = forceFieldCellX >= (jello->resolution - 1) ? (jello->resolution - 2) : forceFieldCellX;
                forceFieldCellY = forceFieldCellY < 0 ? 0 : forceFieldCellY;
                forceFieldCellY = forceFieldCellY >= (jello->resolution - 1) ? (jello->resolution - 2) : forceFieldCellY;
                forceFieldCellZ = forceFieldCellZ < 0 ? 0 : forceFieldCellZ;
                forceFieldCellZ = forceFieldCellZ >= (jello->resolution - 1) ? (jello->resolution - 2) : forceFieldCellZ;
                //extract force field forces at cell's corner
                point forceFieldCorner000, forceFieldCorner001, forceFieldCorner010, forceFieldCorner100, forceFieldCorner011, forceFieldCorner101, forceFieldCorner110, forceFieldCorner111;
                pCPY(jello->forceField[forceFieldCellX * jello->resolution * jello->resolution + forceFieldCellY * jello->resolution + forceFieldCellZ], forceFieldCorner000);
                pCPY(jello->forceField[forceFieldCellX * jello->resolution * jello->resolution + forceFieldCellY * jello->resolution + (forceFieldCellZ + 1)], forceFieldCorner001);
                pCPY(jello->forceField[forceFieldCellX * jello->resolution * jello->resolution + (forceFieldCellY + 1) * jello->resolution + forceFieldCellZ], forceFieldCorner010);
                pCPY(jello->forceField[(forceFieldCellX + 1) * jello->resolution * jello->resolution + forceFieldCellY * jello->resolution + forceFieldCellZ], forceFieldCorner100);
                pCPY(jello->forceField[forceFieldCellX * jello->resolution * jello->resolution + (forceFieldCellY + 1) * jello->resolution + (forceFieldCellZ + 1)], forceFieldCorner011);
                pCPY(jello->forceField[(forceFieldCellX + 1) * jello->resolution * jello->resolution + forceFieldCellY * jello->resolution + (forceFieldCellZ + 1)], forceFieldCorner101);
                pCPY(jello->forceField[(forceFieldCellX + 1) * jello->resolution * jello->resolution + (forceFieldCellY + 1) * jello->resolution + forceFieldCellZ], forceFieldCorner110);
                pCPY(jello->forceField[(forceFieldCellX + 1) * jello->resolution * jello->resolution + (forceFieldCellY + 1) * jello->resolution + (forceFieldCellZ + 1)], forceFieldCorner111);
                //normalize current point coordinate for trilinear interpolation
                point normalizedPoint;
                normalizedPoint.x = (jello->p[i][j][k].x - (forceFieldCellX * forceFieldUnit - 2)) / forceFieldUnit;
                normalizedPoint.y = (jello->p[i][j][k].y - (forceFieldCellY * forceFieldUnit - 2)) / forceFieldUnit;
                normalizedPoint.z = (jello->p[i][j][k].z - (forceFieldCellZ * forceFieldUnit - 2)) / forceFieldUnit;
                // trilinear interpolation
                point interpolatedForce, tempForceComponent;
                pZeros(interpolatedForce);
                pMULTIPLY(forceFieldCorner000, (1 - normalizedPoint.x) * (1 - normalizedPoint.y) * (1 - normalizedPoint.z), tempForceComponent);
                pSUM(interpolatedForce, tempForceComponent, interpolatedForce);
                pMULTIPLY(forceFieldCorner001, (1 - normalizedPoint.x) * (1 - normalizedPoint.y) * (normalizedPoint.z), tempForceComponent);
                pSUM(interpolatedForce, tempForceComponent, interpolatedForce);
                pMULTIPLY(forceFieldCorner010, (1 - normalizedPoint.x) * (normalizedPoint.y) * (1 - normalizedPoint.z), tempForceComponent);
                pSUM(interpolatedForce, tempForceComponent, interpolatedForce);
                pMULTIPLY(forceFieldCorner100, (normalizedPoint.x) * (1 - normalizedPoint.y) * (1 - normalizedPoint.z), tempForceComponent);
                pSUM(interpolatedForce, tempForceComponent, interpolatedForce);
                pMULTIPLY(forceFieldCorner011, (1 - normalizedPoint.x) * (normalizedPoint.y) * (normalizedPoint.z), tempForceComponent);
                pSUM(interpolatedForce, tempForceComponent, interpolatedForce);
                pMULTIPLY(forceFieldCorner101, (normalizedPoint.x) * (1 - normalizedPoint.y) * (normalizedPoint.z), tempForceComponent);
                pSUM(interpolatedForce, tempForceComponent, interpolatedForce);
                pMULTIPLY(forceFieldCorner110, (normalizedPoint.x) * (normalizedPoint.y) * (1 - normalizedPoint.z), tempForceComponent);
                pSUM(interpolatedForce, tempForceComponent, interpolatedForce);
                pMULTIPLY(forceFieldCorner111, (normalizedPoint.x) * (normalizedPoint.y) * (normalizedPoint.z), tempForceComponent);
                pSUM(interpolatedForce, tempForceComponent, interpolatedForce);
                //printf("%f, %f, %f \n", interpolatedForce.x, interpolatedForce.y, interpolatedForce.z);
                // update force to result force array
                f[i][j][k] = pSum(f[i][j][k], interpolatedForce);
            }
        }
    }
}

void computeMouseDragForce(world* jello, point*** f)
{
    //interpolate force field force at each position
    if (jello->resolution == 0) {
        return;
    }
    //printf("%d\n", jello->resolution);
    double forceFieldUnit = 4.0 / double(jello->resolution - 1);
    // in input, force field have a fixed offset of -2, -2, -2
    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            for (int k = 0; k < 8; k++)
            {
                //find which cell current point is
                point pos = jello->p[i][j][k];
                int forceFieldCellX = int((pos.x + 2) / forceFieldUnit);
                int forceFieldCellY = int((pos.y + 2) / forceFieldUnit);
                int forceFieldCellZ = int((pos.z + 2) / forceFieldUnit);
                //clamp to force field boundary (cell's lower boundary)
                forceFieldCellX = forceFieldCellX < 0 ? 0 : forceFieldCellX;
                forceFieldCellX = forceFieldCellX >= (jello->resolution - 1) ? (jello->resolution - 2) : forceFieldCellX;
                forceFieldCellY = forceFieldCellY < 0 ? 0 : forceFieldCellY;
                forceFieldCellY = forceFieldCellY >= (jello->resolution - 1) ? (jello->resolution - 2) : forceFieldCellY;
                forceFieldCellZ = forceFieldCellZ < 0 ? 0 : forceFieldCellZ;
                forceFieldCellZ = forceFieldCellZ >= (jello->resolution - 1) ? (jello->resolution - 2) : forceFieldCellZ;
                //extract force field forces at cell's corner
                point draggedForce;
                draggedForce.x = draggedX;
                draggedForce.y = draggedY;
                draggedForce.z = 0;
                //printf("%f, %f\n", draggedX, draggedY);
                dragged = false;
                point forceFieldCorner000, forceFieldCorner001, forceFieldCorner010, forceFieldCorner100, forceFieldCorner011, forceFieldCorner101, forceFieldCorner110, forceFieldCorner111;
                pCPY(draggedForce, forceFieldCorner000);
                pCPY(draggedForce, forceFieldCorner001);
                pCPY(draggedForce, forceFieldCorner010);
                pCPY(draggedForce, forceFieldCorner100);
                pCPY(draggedForce, forceFieldCorner011);
                pCPY(draggedForce, forceFieldCorner101);
                pCPY(draggedForce, forceFieldCorner110);
                pCPY(draggedForce, forceFieldCorner111);
                //normalize current point coordinate for trilinear interpolation
                point normalizedPoint;
                normalizedPoint.x = (jello->p[i][j][k].x - (forceFieldCellX * forceFieldUnit - 2)) / forceFieldUnit;
                normalizedPoint.y = (jello->p[i][j][k].y - (forceFieldCellY * forceFieldUnit - 2)) / forceFieldUnit;
                normalizedPoint.z = (jello->p[i][j][k].z - (forceFieldCellZ * forceFieldUnit - 2)) / forceFieldUnit;
                // trilinear interpolation
                point interpolatedForce, tempForceComponent;
                pZeros(interpolatedForce);
                pMULTIPLY(forceFieldCorner000, (1 - normalizedPoint.x) * (1 - normalizedPoint.y) * (1 - normalizedPoint.z), tempForceComponent);
                pSUM(interpolatedForce, tempForceComponent, interpolatedForce);
                pMULTIPLY(forceFieldCorner001, (1 - normalizedPoint.x) * (1 - normalizedPoint.y) * (normalizedPoint.z), tempForceComponent);
                pSUM(interpolatedForce, tempForceComponent, interpolatedForce);
                pMULTIPLY(forceFieldCorner010, (1 - normalizedPoint.x) * (normalizedPoint.y) * (1 - normalizedPoint.z), tempForceComponent);
                pSUM(interpolatedForce, tempForceComponent, interpolatedForce);
                pMULTIPLY(forceFieldCorner100, (normalizedPoint.x) * (1 - normalizedPoint.y) * (1 - normalizedPoint.z), tempForceComponent);
                pSUM(interpolatedForce, tempForceComponent, interpolatedForce);
                pMULTIPLY(forceFieldCorner011, (1 - normalizedPoint.x) * (normalizedPoint.y) * (normalizedPoint.z), tempForceComponent);
                pSUM(interpolatedForce, tempForceComponent, interpolatedForce);
                pMULTIPLY(forceFieldCorner101, (normalizedPoint.x) * (1 - normalizedPoint.y) * (normalizedPoint.z), tempForceComponent);
                pSUM(interpolatedForce, tempForceComponent, interpolatedForce);
                pMULTIPLY(forceFieldCorner110, (normalizedPoint.x) * (normalizedPoint.y) * (1 - normalizedPoint.z), tempForceComponent);
                pSUM(interpolatedForce, tempForceComponent, interpolatedForce);
                pMULTIPLY(forceFieldCorner111, (normalizedPoint.x) * (normalizedPoint.y) * (normalizedPoint.z), tempForceComponent);
                pSUM(interpolatedForce, tempForceComponent, interpolatedForce);
                //printf("%f, %f, %f \n", interpolatedForce.x, interpolatedForce.y, interpolatedForce.z);
                // update force to result force array
                f[i][j][k] = pSum(f[i][j][k], interpolatedForce);
            }
        }
    }
}

/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(world * jello, point*** a)
{
    //struct point springForces[8][8][8];
    //struct point forceFieldForces[8][8][8];
    //struct point collisionForces[8][8][8];
    //struct point f[8][8][8];

    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            for (int k = 0; k < 8; k++)
            {
                //pZeros(springForces[i][j][k]);
                //pZeros(forceFieldForces[i][j][k]);
                //pZeros(collisionForces[i][j][k]);
                pZeros(a[i][j][k]);
            }
        }
    }

    world buffer = *jello;
    computeSpringForce(&buffer, a);
    computeCollisionForce(&buffer, a);
    computeForceFieldForce(&buffer, a);
    computeMouseDragForce(&buffer, a);

    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            for (int k = 0; k < 8; k++)
            {
                a[i][j][k].x = a[i][j][k].x / jello->mass;
                a[i][j][k].y = a[i][j][k].y / jello->mass;
                a[i][j][k].z = a[i][j][k].z / jello->mass;
            }
        }
    }
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point*** a = new point * *[8];

  for (int i = 0; i < 8; i++)
  {
      a[i] = new point * [8];
      for (int j = 0; j < 8; j++) {
          a[i][j] = new point[8];
      }
  }

  // assign values to the allocated memory
  for (int i = 0; i < 8; i++)
  {
      for (int j = 0; j < 8; j++)
      {
          for (int k = 0; k < 8; k++) {
              pZeros(a[i][j][k]);
          }
      }
  }

  computeAcceleration(jello, a);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

      }

  // deallocate memory
  for (int i = 0; i < 8; i++)
  {
      for (int j = 0; j < 8; j++) {
          delete[] a[i][j];
      }
      delete[] a[i];
  }

  delete[] a;
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
    point*** F1p = new point * *[8];
    point*** F2p = new point * *[8];
    point*** F3p = new point * *[8];
    point*** F4p = new point * *[8];
    point*** F1v = new point * *[8];
    point*** F2v = new point * *[8];
    point*** F3v = new point * *[8];
    point*** F4v = new point * *[8];
    point*** a = new point * *[8];

    for (int i = 0; i < 8; i++)
    {
        F1p[i] = new point * [8];
        F2p[i] = new point * [8];
        F3p[i] = new point * [8];
        F4p[i] = new point * [8];
        F1v[i] = new point * [8];
        F2v[i] = new point * [8];
        F3v[i] = new point * [8];
        F4v[i] = new point * [8];
        a[i] = new point * [8];
        for (int j = 0; j < 8; j++) {
            F1p[i][j] = new point[8];
            F2p[i][j] = new point[8];
            F3p[i][j] = new point[8];
            F4p[i][j] = new point[8];
            F1v[i][j] = new point[8];
            F2v[i][j] = new point[8];
            F3v[i][j] = new point[8];
            F4v[i][j] = new point[8];
            a[i][j] = new point[8];
        }
    }

    // assign values to the allocated memory
    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            for (int k = 0; k < 8; k++) {
                pZeros(F1p[i][j][k]);
                pZeros(F2p[i][j][k]);
                pZeros(F3p[i][j][k]);
                pZeros(F4p[i][j][k]);
                pZeros(F1v[i][j][k]);
                pZeros(F2v[i][j][k]);
                pZeros(F3v[i][j][k]);
                pZeros(F4v[i][j][k]);
                pZeros(a[i][j][k]);
            }
        }
    }
  //point F1p[8][8][8], F1v[8][8][8], 
  //      F2p[8][8][8], F2v[8][8][8],
  //      F3p[8][8][8], F3v[8][8][8],
  //      F4p[8][8][8], F4v[8][8][8];

  //point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],1.0,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],1.0,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }

  // deallocate memory
  for (int i = 0; i < 8; i++)
  {
      for (int j = 0; j < 8; j++) {
          delete[] F1p[i][j];
          delete[] F2p[i][j];
          delete[] F3p[i][j];
          delete[] F4p[i][j];
          delete[] F1v[i][j];
          delete[] F2v[i][j];
          delete[] F3v[i][j];
          delete[] F4v[i][j];
          delete[] a[i][j];
      }
      delete[] F1p[i];
      delete[] F2p[i];
      delete[] F3p[i];
      delete[] F4p[i];
      delete[] F1v[i];
      delete[] F2v[i];
      delete[] F3v[i];
      delete[] F4v[i];
      delete[] a[i];
  }

  delete[] F1p;
  delete[] F2p;
  delete[] F3p;
  delete[] F4p;
  delete[] F1v;
  delete[] F2v;
  delete[] F3v;
  delete[] F4v;
  delete[] a;

  return;  
}
