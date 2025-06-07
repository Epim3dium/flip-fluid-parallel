#ifndef EMP_AABB_HPP
#define EMP_AABB_HPP

#include <SFML/Graphics/Color.hpp>
#include <SFML/System/Vector2.hpp>
typedef sf::Color Color;
#include "vec2.hpp"
#include <concepts>
#include <limits>
#include <cmath>
#include <math.h>
#include <numeric>
#include <vector>



struct AABB {
    vec2f min;
    vec2f max;

    vec2f bl() const {
        return min;
    }
    vec2f br() const {
        return vec2f(max.x, min.y);
    }
    vec2f tl() const {
        return vec2f(min.x, max.y);
    }
    vec2f tr() const {
        return max;
    }

    vec2f center() const {
        return min + (max-min) / 2.f;
    }
    float right() const {
        return max.x;
    };
    float left() const {
        return min.x;
    };
    float bottom() const {
        return min.y;
    };
    float top() const {
        return max.y;
    };
    vec2f size() const {
        return max - min;
    }
    AABB& setCenter(vec2f c) {
        auto t = size();
        min = c - t / 2.f;
        max = c + t / 2.f;
        return *this;
    }
    AABB& move(vec2f offset) {
        min += offset;
        max += offset;
        return *this;
    }
    AABB& expandToContain(vec2f point) {
        min.x = std::fmin(min.x, point.x);
        min.y = std::fmin(min.y, point.y);
        max.x = std::fmax(max.x, point.x);
        max.y = std::fmax(max.y, point.y);
        return *this;
    }
    AABB& setSize(vec2f s) {
        auto t = center();
        min = t - s / 2.f;
        max = t + s / 2.f;
        return *this;
    }
    AABB combine(AABB val) {
        val.min.x = std::min(min.x, val.min.x);
        val.min.y = std::min(min.y, val.min.y);
        val.max.x = std::max(max.x, val.max.x);
        val.max.y = std::max(max.y, val.max.y);
        return val;
    }

    static AABB Expandable() {
        return {{INFINITY, INFINITY}, {-INFINITY, -INFINITY}};
    }
    static AABB CreateMinMax(vec2f min, vec2f max) {
        return {min, max};
    }
    static AABB CreateCenterSize(vec2f center, vec2f size) {
        AABB  result;
        result.setCenter(center);
        result.setSize(size);
        return result;
    }
    static AABB CreateMinSize(vec2f min, vec2f size) {
        AABB  result{min, min+size};
        return result;
    }
};

#endif
