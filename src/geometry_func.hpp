#ifndef EMP_GEOMETRY_FUNC_HPP
#define EMP_GEOMETRY_FUNC_HPP
#include <vector>
#include "AABB.hpp"
#include "circle.hpp"


float angleAround(vec2f a, vec2f pivot, vec2f b);
vec2f sign(vec2f x);
vec2f rotateVec(vec2f vec, float angle);
float qlen(vec2f v);
float angle(vec2f a, vec2f b);
vec2f rotate(vec2f vec, float angle);
float length(vec2f v);
float dot(vec2f a, vec2f b);
vec2f normal(vec2f v);
vec2f proj(vec2f a, vec2f plane_norm);
float perp_dot(vec2f a, vec2f b);
// returns true if r1 contains the whole of r2
bool AABBcontainsAABB(const AABB& r1, const AABB& r2);
// finds the closest vector to point that lies on ray
vec2f findClosestPointOnRay(vec2f ray_origin, vec2f ray_dir, vec2f point);
// finds the closest vetor to point that lies on one of poly's edges
vec2f findClosestPointOnEdge(vec2f point, const std::vector<vec2f>& poly);
std::vector<vec2f> findContactPoints(
        const std::vector<vec2f>& p0, const std::vector<vec2f>& p1
);
// returns true if p is within aabb
bool isOverlappingPointAABB(const vec2f& p, const AABB& r);
// returns true if p is within circle
bool isOverlappingPointCircle(const vec2f& p, const Circle& c);
// returns true if p is within polygon
bool isOverlappingPointPoly(
        const vec2f& p, const std::vector<vec2f>& poly_points
);
// returns true if aabb and aabb are overlapping
bool isOverlappingAABBAABB(const AABB& r1, const AABB& r2);

/**
 * structure containing all info returned by Ray and AABB intersection
 *
 * detected - true if intersection occured
 * time_hit_near - time along ray_dir where first intersection occured
 * time_hit_far - time along ray_dir where second intersection occured
 * contact normal - normal of first collision
 * contact point - point where first intersection took place
 */
struct IntersectionRayAABBResult {
    bool detected;
    float time_hit_near;
    float time_hit_far;
    vec2f contact_normal;
    vec2f contact_point;
};
/**
 * Calculates all information connected to Ray and AABB intersection
 * @param ray_origin is the origin point of ray
 * @param ray_dir is the direction that after adding to ray_origin gives other
 * point on ray
 * @return IntersectionRayAABBResult that contains: (in order) [bool]detected,
 * [float]time_hit_near, [float]time_hit_far, [vec2f]contact_normal,
 * [vec2f]contact_point
 */
IntersectionRayAABBResult intersectRayAABB(
        vec2f ray_origin, vec2f ray_dir, const AABB& target
);

/**
 * structure containing all info returned by Ray and Ray intersection
 *
 * detected - true if intersection occured [note that even when time_hit_near is
 * larger than 1, aka it 'goes out of ray' this still returns true] contact
 * point - point where first intersection took place time_hit_near0 - time along
 * ray_dir0 where intersection occured time_hit_near1 - time along ray_dir1
 * where second intersection occured
 */
struct IntersectionRayRayResult {
    bool detected;
    bool nearParallel;
    vec2f contact_point;
    float t_hit_near0;
    float t_hit_near1;
};
/**
 * Calculates all information connected to Ray and Ray intersection
 * @return IntersectionRayRayResult that contains: (in order) [bool]detected,
 * [vec2f]contact_point, [float]t_hit_near0, [float]t_hit_near1
 */
IntersectionRayRayResult intersectRayRay(
        vec2f ray0_origin, vec2f ray0_dir, vec2f ray1_origin, vec2f ray1_dir
);

struct IntersectionCircleCircleResult {
    bool detected;
    vec2f contact_normal;
    vec2f contact_point;
    float overlap;
};
/**
 * Calculates all information connected to Polygon and Polygon intersection
 * @return IntersectionPolygonPolygonResult that contains: (in order)
 * [bool]detected, [vec2f]contact_normal, [float]overlap
 */
IntersectionCircleCircleResult intersectCircleCircle(
        const Circle& c1, const Circle& c2
);

#endif
