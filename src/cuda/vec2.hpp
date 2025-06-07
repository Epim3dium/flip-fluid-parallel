#pragma once
#include <cmath>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
struct vec2f {
    float x, y;

    // ======= Constructors =======
    __host__ __device__ vec2f() : x(0), y(0) {}
    __host__ __device__ vec2f(float x_, float y_) : x(x_), y(y_) {}

    // ======= Unary =======
    __host__ __device__ vec2f operator-() const {
        return {-x, -y};
    }

    // ======= Binary: vec2f <op> vec2f =======
    __host__ __device__ vec2f operator+(const vec2f& rhs) const {
        return {x + rhs.x, y + rhs.y};
    }

    __host__ __device__ vec2f operator-(const vec2f& rhs) const {
        return {x - rhs.x, y - rhs.y};
    }

    __host__ __device__ vec2f operator*(const vec2f& rhs) const {
        return {x * rhs.x, y * rhs.y};
    }

    __host__ __device__ vec2f operator/(const vec2f& rhs) const {
        return {x / rhs.x, y / rhs.y};
    }

    // ======= Compound assignment: vec2f <op>= vec2f =======
    __host__ __device__ vec2f& operator+=(const vec2f& rhs) {
        x += rhs.x; y += rhs.y;
        return *this;
    }

    __host__ __device__ vec2f& operator-=(const vec2f& rhs) {
        x -= rhs.x; y -= rhs.y;
        return *this;
    }

    __host__ __device__ vec2f& operator*=(const vec2f& rhs) {
        x *= rhs.x; y *= rhs.y;
        return *this;
    }

    __host__ __device__ vec2f& operator/=(const vec2f& rhs) {
        x /= rhs.x; y /= rhs.y;
        return *this;
    }

    // ======= Binary: vec2f <op> float =======
    __host__ __device__ vec2f operator+(float scalar) const {
        return {x + scalar, y + scalar};
    }

    __host__ __device__ vec2f operator-(float scalar) const {
        return {x - scalar, y - scalar};
    }

    __host__ __device__ vec2f operator*(float scalar) const {
        return {x * scalar, y * scalar};
    }

    __host__ __device__ vec2f operator/(float scalar) const {
        return {x / scalar, y / scalar};
    }

    // ======= Compound assignment: vec2f <op>= float =======
    __host__ __device__ vec2f& operator+=(float scalar) {
        x += scalar; y += scalar;
        return *this;
    }

    __host__ __device__ vec2f& operator-=(float scalar) {
        x -= scalar; y -= scalar;
        return *this;
    }

    __host__ __device__ vec2f& operator*=(float scalar) {
        x *= scalar; y *= scalar;
        return *this;
    }

    __host__ __device__ vec2f& operator/=(float scalar) {
        x /= scalar; y /= scalar;
        return *this;
    }
        // === Comparison Operators (Component-wise) ===
    __host__ __device__ bool operator==(const vec2f& rhs) const { return x == rhs.x && y == rhs.y; }
    __host__ __device__ bool operator!=(const vec2f& rhs) const { return !(*this == rhs); }
};

// ======= Global: float <op> vec2f =======
__host__ __device__ inline vec2f operator+(float scalar, const vec2f& v) {
    return {scalar + v.x, scalar + v.y};
}

__host__ __device__ inline vec2f operator-(float scalar, const vec2f& v) {
    return {scalar - v.x, scalar - v.y};
}

__host__ __device__ inline vec2f operator*(float scalar, const vec2f& v) {
    return {scalar * v.x, scalar * v.y};
}

__host__ __device__ inline vec2f operator/(float scalar, const vec2f& v) {
    return {scalar / v.x, scalar / v.y};
}
template <class Vec>
__host__ __device__ float perp_dot(Vec a, Vec b) {
    return a.x * b.y - b.x * a.y;
}

template <class Vec>
__host__ __device__ float qlen(Vec v) {
    return v.x * v.x + v.y * v.y;
}

template <class Vec>
__host__ __device__ float dot(Vec a, Vec b) {
    return a.x * b.x + a.y * b.y;
}

template <class Vec>
__host__ __device__ float length(Vec v) {
    return sqrtf(qlen(v));
}

template <class Vec>
__host__ __device__ float angle(Vec a, Vec b) {
    return atan2f(perp_dot(a, b), dot(a, b));
}

template <class Vec>
__host__ __device__ vec2f rotate(Vec v, float angle) {
    float ca = cosf(angle);
    float sa = sinf(angle);
    return { ca * v.x - sa * v.y, sa * v.x + ca * v.y };
}

template <class Vec>
__host__ __device__ vec2f normal(Vec v) {
    float len = length(v);
    return len > 0.0f ? v / len : vec2f{0.0f, 0.0f};
}

template <class Vec>
__host__ __device__ vec2f proj(Vec a, Vec plane_norm) {
    float d = dot(plane_norm, plane_norm);
    if (d == 0.0f) return vec2f{0.0f, 0.0f}; // avoid divide by zero
    return (dot(a, plane_norm) / d) * plane_norm;
}
