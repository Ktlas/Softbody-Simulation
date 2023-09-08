/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#ifndef _PHYSICS_H_
#define _PHYSICS_H_

void computeHooksForce(point pointAPos, point pointBPos, double k, double restLength, point& resultForce);

void computeDampingForce(point pointAPos, point pointBPos, point pointAVelocity, point pointBVelocity, double d, point& resultForce);

void computeSpringForce(world* jello, point*** f);

void computeIndividualSpringForce(int srcx, int srcy, int srcz, int destx, int desty, int destz, double restLength, const world* jello, point*** f);

void computeForceFieldForce(world* jello, point*** f);

void computeCollisionForce(world* jello, point*** f);

void computeMouseDragForce(world* jello, point*** f);

void computeAcceleration(world* jello, point*** a);

// perform one step of Euler and Runge-Kutta-4th-order integrators
// updates the jello structure accordingly
void Euler(struct world * jello);
void RK4(struct world * jello);

#endif

