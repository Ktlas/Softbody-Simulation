/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include <vector>
#ifndef _SHOWCUBE_H_
#define _SHOWCUBE_H_

struct Plane {
    point normal;
    point point;
};

struct AABB {
    point min;
    point max;
};

struct Segment {
    point start;
    point end;
};

point getIntersection(Plane& plane, Segment& segment, bool flag);

std::vector<point> getIntersectionPolygon(Plane& plane, AABB& aabb);

void showCube(struct world * jello);

void showBoundingBox();

void showInclinedPlane(struct world * jello);

#endif
