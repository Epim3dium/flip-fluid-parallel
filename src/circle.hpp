#ifndef EMP_CIRCLE_HPP
#define EMP_CIRCLE_HPP
#include "AABB.hpp"

struct Circle {
    vec2f pos;
    float radius;
    Circle(vec2f p = vec2f(0, 0), float r = 1.f) : pos(p), radius(r) {}
};
#endif
