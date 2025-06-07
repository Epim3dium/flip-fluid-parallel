#pragma once
#include <SFML/System/Vector2.hpp>
#include <cmath>
typedef sf::Vector2f vec2f;
template<class Vec>
float perp_dot(Vec a, Vec b) {
    return a.x * b.y - b.x * a.y;
}
template<class Vec>
float qlen(Vec v) {
    return v.x * v.x + v.y * v.y;
}
template<class Vec>
float angle(Vec a, Vec b) {
    return atan2(perp_dot(a, b), dot(a, b));
}
template<class Vec>
vec2f rotate(Vec vec, float angle) {
    return {
            cosf(angle) * vec.x - sinf(angle) * vec.y,
            sinf(angle) * vec.x + cosf(angle) * vec.y,
    };
}
template<class Vec>
float length(Vec v) {
    return sqrt(v.x * v.x + v.y * v.y);
}
template<class Vec>
float dot(Vec a, Vec b) {
    return a.x * b.x + a.y * b.y;
}
template<class Vec>
vec2f normal(Vec v) {
    return v / length(v);
}
template<class Vec>
vec2f proj(Vec a, Vec plane_norm) {
    return (dot(a, plane_norm) / dot(plane_norm, plane_norm)) * plane_norm;
}
