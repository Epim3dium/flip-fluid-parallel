#include "geometry_func.hpp"
#include <array>

float angleAround(vec2f a, vec2f pivot, vec2f b) {
    return angle(a - pivot, b - pivot);
}
vec2f sign(vec2f x) {
    return {std::copysign(1.f, x.x), std::copysign(1.f, x.y)};
}
vec2f rotateVec(vec2f vec, float angle) {
    return vec2f(
            cos(angle) * vec.x - sin(angle) * vec.y,
            sin(angle) * vec.x + cos(angle) * vec.y
    );
}

float qlen(vec2f v) {
    return v.x * v.x + v.y * v.y;
}
float angle(vec2f a, vec2f b) {
    return atan2(perp_dot(a, b), dot(a, b));
}
vec2f rotate(vec2f vec, float angle) {
    return {
            cosf(angle) * vec.x - sinf(angle) * vec.y,
            sinf(angle) * vec.x + cosf(angle) * vec.y,
    };
}
float length(vec2f v) {
    return sqrt(v.x * v.x + v.y * v.y);
}
float dot(vec2f a, vec2f b) {
    return a.x * b.x + a.y * b.y;
}
vec2f normal(vec2f v) {
    return v / length(v);
}
vec2f proj(vec2f a, vec2f plane_norm) {
    return (dot(a, plane_norm) / dot(plane_norm, plane_norm)) * plane_norm;
}
float perp_dot(vec2f a, vec2f b) {
    return a.x * b.y - b.x * a.y;
}
#define SQR(x) ((x) * (x))
bool isOverlappingPointAABB(const vec2f& p, const AABB& r) {
    return (p.x >= r.center().x - r.size().x / 2 &&
            p.y >= r.center().y - r.size().y / 2 &&
            p.x <= r.center().x + r.size().x / 2 &&
            p.y <= r.center().y + r.size().y / 2);
}
bool isOverlappingPointCircle(const vec2f& p, const Circle& c) {
    return length(p - c.pos) <= c.radius;
}
bool isOverlappingPointPoly(const vec2f& p, const std::vector<vec2f>& points) {
    int i, j, c = 0;
    for (i = 0, j = points.size() - 1; i < points.size(); j = i++) {
        auto& vi = points[i];
        auto& vj = points[j];
        if (((vi.y > p.y) != (vj.y > p.y)) &&
            (p.x < (vj.x - vi.x) * (p.y - vi.y) / (vj.y - vi.y) + vi.x))
            c = !c;
    }
    return c;
}
bool isOverlappingAABBAABB(const AABB& r1, const AABB& r2) {
    return (r1.min.x <= r2.max.x && r1.max.x >= r2.min.x &&
            r1.min.y <= r2.max.y && r1.max.y >= r2.min.y);
}
bool AABBcontainsAABB(const AABB& r1, const AABB& r2) {
    return (r2.min.x >= r1.min.x) && (r2.max.x <= r1.max.x) &&
           (r2.min.y >= r1.min.y) && (r2.max.y <= r1.max.y);
}
IntersectionRayAABBResult intersectRayAABB(
        vec2f ray_origin, vec2f ray_dir, const AABB& target
) {
    IntersectionRayAABBResult result;
    vec2f invdir = {1.0f / ray_dir.x, 1.0f / ray_dir.y};
    vec2f t_size = target.size();
    // VVVVVVVVVVVVV
    // if((int)target.size.y % 2 == 0 && target.pos.y > ray_origin.y)
    // t_size -= vec2f(0, 1);
    //^^^^^^^^^^^^^
    vec2f t_near = (target.center() - t_size / 2.f - ray_origin);
    t_near.x *= invdir.x;
    t_near.y *= invdir.y;
    vec2f t_far = (target.center() + t_size / 2.f - ray_origin);
    t_far.x *= invdir.x;
    t_far.y *= invdir.y;

    if (std::isnan(t_far.y) || std::isnan(t_far.x))
        return {false};
    if (std::isnan(t_near.y) || std::isnan(t_near.x))
        return {false};
    if (t_near.x > t_far.x)
        std::swap(t_near.x, t_far.x);
    if (t_near.y > t_far.y)
        std::swap(t_near.y, t_far.y);

    if (t_near.x > t_far.y || t_near.y > t_far.x)
        return {false};
    float t_hit_near = std::max(t_near.x, t_near.y);
    result.time_hit_near = t_hit_near;
    float t_hit_far = std::min(t_far.x, t_far.y);
    result.time_hit_far = t_hit_far;

    if (t_hit_far < 0)
        return {false};
    result.contact_point = ray_origin + ray_dir * t_hit_near;
    if (t_near.x > t_near.y) {
        if (invdir.x < 0)
            result.contact_normal = {1, 0};
        else
            result.contact_normal = {-1, 0};
    } else if (t_near.x < t_near.y) {
        if (invdir.y < 0)
            result.contact_normal = {0, 1};
        else
            result.contact_normal = {0, -1};
    }
    result.detected = true;
    return result;
}
IntersectionRayRayResult intersectRayRay(
        vec2f ray0_origin, vec2f ray0_dir, vec2f ray1_origin, vec2f ray1_dir
) {
    if (ray0_origin == ray1_origin) {
        return {true, false, ray0_origin, 0.f, 0.f};
    }
    auto dx = ray1_origin.x - ray0_origin.x;
    auto dy = ray1_origin.y - ray0_origin.y;
    float det = ray1_dir.x * ray0_dir.y - ray1_dir.y * ray0_dir.x;
    if (det != 0) { // near parallel line will yield noisy results
        float u = (dy * ray1_dir.x - dx * ray1_dir.y) / det;
        float v = (dy * ray0_dir.x - dx * ray0_dir.y) / det;
        return {u >= 0 && v >= 0 && u <= 1.f && v <= 1.f,
                false,
                ray0_origin + ray0_dir * u,
                u,
                v};
    }
    return {false, true};
}
/*
    bool detected;
    vec2f contact_point;
    float t_hit_near0;
    float t_hit_near1;
*/
vec2f findClosestPointOnRay(vec2f ray_origin, vec2f ray_dir, vec2f point) {
    float ray_dir_len = length(ray_dir);
    vec2f ray_dir_normal = ray_dir / ray_dir_len;
    float proj = dot(point - ray_origin, ray_dir_normal);
    if (proj <= 0)
        return ray_origin;
    if (proj >= ray_dir_len)
        return ray_origin + ray_dir;
    return ray_dir_normal * proj + ray_origin;
}
vec2f findClosestPointOnEdge(vec2f point, const std::vector<vec2f>& points) {
    vec2f closest(INFINITY, INFINITY);
    float closest_dist = INFINITY;
    for (size_t i = 0; i < points.size(); i++) {
        vec2f a = points[i];
        vec2f b = points[(i + 1) % points.size()];
        vec2f adir = b - a;
        vec2f t = findClosestPointOnRay(a, adir, point);
        float dist = length(t - point);
        if (dist < closest_dist) {
            closest_dist = dist;
            closest = t;
        }
    }
    return closest;
}
IntersectionCircleCircleResult intersectCircleCircle(
        const Circle& c1, const Circle& c2
) {
    vec2f dist = c1.pos - c2.pos;
    float dist_len = length(dist);
    if (dist_len > c1.radius + c2.radius) {
        return {false};
    }
    float overlap = c1.radius + c2.radius - dist_len;
    vec2f contact_point = dist / dist_len * c2.radius + c2.pos;
    return {true, dist / dist_len, contact_point, overlap};
}
float calculateInertia(const std::vector<vec2f>& model, float mass) {
    float area = 0;
    float mmoi = 0;

    auto prev = model.back();
    for (auto next : model) {
        float area_step = abs(perp_dot(prev, next)) / 2.f;
        float mmoi_step =
                area_step *
                (dot(prev, prev) + dot(next, next) + abs(dot(prev, next))) /
                6.f;

        area += area_step;
        mmoi += mmoi_step;

        prev = next;
    }

    double density = mass / area;
    mmoi *= density;
    // mmoi -= mass * dot(center, center);
    if (std::isnan(mmoi)) {
        mmoi = 0.f;
    }
    return abs(mmoi);
}

