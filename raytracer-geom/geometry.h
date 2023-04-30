#include <vector.h>
#include <sphere.h>
#include <intersection.h>
#include <triangle.h>
#include <ray.h>

#include <optional>
#include <cmath>

inline std::optional<Intersection> GetIntersection(const Ray& ray, const Sphere& sphere) {
    double eps = 1e-8;
    Vector l = Difference(sphere.GetCenter(), ray.GetOrigin());
    if (std::fabs(Length(l) - sphere.GetRadius()) < eps) {
        if (DotProduct(l, ray.GetDirection()) < 0) {
            return std::nullopt;
        }
        double cross_prod = Length(CrossProduct(l, ray.GetDirection()));
        cross_prod /= Length(ray.GetDirection());
        double length = sqrt(Length(l) * Length(l) - cross_prod * cross_prod) +
            sqrt(sphere.GetRadius() * sphere.GetRadius() - cross_prod * cross_prod);

        double original_length = Length(ray.GetDirection());
        double coef = length / original_length;
        Vector answer_point = ray.GetDirection();
        for (size_t i = 0; i < 3; ++i) {
            answer_point[i] *= coef;
        }
        Vector result_point = Sum(ray.GetOrigin(), answer_point);
        //assert(std::fabs(Length(Difference(result_point, sphere.GetCenter())) - sphere.GetRadius()) < 1e-2);
        Vector normal = Difference(result_point, sphere.GetCenter());
        normal.Normalize();
        {
            for (size_t i = 0; i < 3; ++i) {
                normal[i] *= -1;
            }
        }
        return Intersection(result_point, normal, Length(answer_point));
    }
    if (DotProduct(l, ray.GetDirection()) < 0 && Length(l) - sphere.GetRadius() > eps) {
        return std::nullopt;
    }
    double cross_prod = Length(CrossProduct(l, ray.GetDirection()));
    cross_prod /= Length(ray.GetDirection());
    if (cross_prod > eps + sphere.GetRadius()) {
        return std::nullopt;
    }
    double length = sqrt(Length(l) * Length(l) - cross_prod * cross_prod) -
                    sqrt(sphere.GetRadius() * sphere.GetRadius() - cross_prod * cross_prod);
    if (Length(l) + eps < sphere.GetRadius()) {
        if (DotProduct(l, ray.GetDirection()) > 0) {
            length = sqrt(Length(l) * Length(l) - cross_prod * cross_prod) +
                     sqrt(sphere.GetRadius() * sphere.GetRadius() - cross_prod * cross_prod);
        } else {
            length = -sqrt(Length(l) * Length(l) - cross_prod * cross_prod) +
                sqrt(sphere.GetRadius() * sphere.GetRadius() - cross_prod * cross_prod);
        }
    }
    double original_length = Length(ray.GetDirection());
    double coef = length / original_length;
    Vector answer_point = ray.GetDirection();
    for (size_t i = 0; i < 3; ++i) {
        answer_point[i] *= coef;
    }
    Vector result_point = Sum(ray.GetOrigin(), answer_point);
    Vector normal = Difference(result_point, sphere.GetCenter());
    normal.Normalize();
    if (Length(l) + eps < sphere.GetRadius()) {
        for (size_t i = 0; i < 3; ++i) {
            normal[i] *= -1;
        }
    }
    return Intersection(result_point, normal, Length(answer_point));
}

inline bool Inside(const Vector& a, const Triangle& triangle) {
    double eps = 1e-8;
    auto a1 = Triangle({a, triangle.GetVertex(0), triangle.GetVertex(1)});
    auto a2 = Triangle({a, triangle.GetVertex(0), triangle.GetVertex(2)});
    auto a3 = Triangle({a, triangle.GetVertex(1), triangle.GetVertex(2)});

    double sum = a1.Area() + a2.Area() + a3.Area();
    double required = triangle.Area();

    return std::fabs(sum - required) < eps;
}

inline std::optional<Intersection> GetIntersection(const Ray& ray, const Triangle& triangle) {
    Vector normal = CrossProduct(Difference(triangle.GetVertex(0), triangle.GetVertex(1)),
                                 Difference(triangle.GetVertex(0), triangle.GetVertex(2)));
    if (DotProduct(normal, ray.GetDirection()) < 0) {
        for (size_t i = 0; i < 3; ++i) {
            normal[i] *= -1;
        }
    }
    normal.Normalize();
    Vector sample = Difference(triangle.GetVertex(0), ray.GetOrigin());

    double height = DotProduct(sample, normal);
    double mult_coef = height / DotProduct(ray.GetDirection(), normal);
    if (mult_coef < 0) {
        return std::nullopt;
    }
    Vector intersection_point = ray.GetDirection();
    for (size_t i = 0; i < 3; ++i) {
        intersection_point[i] *= mult_coef;
    }
    Vector intersection_point_0 = Sum(intersection_point, ray.GetOrigin());
    if (!Inside(intersection_point_0, triangle) || Length(intersection_point) < 1e-8) {
        return std::nullopt;
    }
    for (size_t i = 0; i < 3; ++i) {
        normal[i] *= -1;
    }
    return Intersection(intersection_point_0, normal, Length(intersection_point));
}

std::optional<Vector> Refract(const Vector& ray, const Vector& normal, double eta) {
    Vector normalised_normal = normal;
    normalised_normal.Normalize();
    double mult_coef = DotProduct(normalised_normal, ray);
    Vector vx = normalised_normal;
    double cos_alpha = std::fabs(mult_coef / Length(ray)) - 1e-8;
    double sin_alpha = sqrt(1.0 - (cos_alpha * cos_alpha));
    double sin_betta = eta * sin_alpha;
    if (std::fabs(sin_betta) > 1) {
        return std::nullopt;
    }
    double cos_betta = sqrt(1.0 - (sin_betta * sin_betta));
    for (size_t i = 0; i < 3; ++i) {
        vx[i] *= mult_coef;
    }
    Vector vy = Difference(ray, vx);
    vx.Normalize();
    vy.Normalize();
    vx = Mult(vx, cos_betta);
    vy = Mult(vy, sin_betta);
    Vector answer = Sum(vx, vy);
    answer.Normalize();
    return answer;
}
Vector Reflect(const Vector& ray, const Vector& normal) {
    Vector add = normal;
    add.Normalize();
    if (DotProduct(ray, add) > 0) {
        add = Mult(add, -1);
    }
    double coef = -DotProduct(ray, add);
    for (size_t i = 0; i < 3; ++i) {
        add[i] *= (2 * coef);
    }
    return Sum(ray, add);
}
Vector GetBarycentricCoords(const Triangle& triangle, const Vector& point) {
    double eps = 1e-8;
    double x00 = triangle.GetVertex(0)[0] - triangle.GetVertex(2)[0];
    double x01 = triangle.GetVertex(1)[0] - triangle.GetVertex(2)[0];
    double x10 = triangle.GetVertex(0)[1] - triangle.GetVertex(2)[1];
    double x11 = triangle.GetVertex(1)[1] - triangle.GetVertex(2)[1];
    double y0 = point[0] - triangle.GetVertex(2)[0];
    double y1 = point[1] - triangle.GetVertex(2)[1];

    if (std::fabs(x01 * x10 - x11 * x00) > eps) {
        double alpha2 = (y0 * x10 - y1 * x00) / (x01 * x10 - x11 * x00 + eps);
        double alpha1 = 1;
        if (std::fabs(x00) > eps) {
            alpha1 = (y0 - alpha2 * x01) / x00;
        } else {
            alpha1 = (y1 - alpha2 * x11) / x10;
        }
        return Vector({alpha1, alpha2, 1 - alpha1 - alpha2});
    } else {
        x00 = triangle.GetVertex(0)[0] - triangle.GetVertex(2)[0];
        x01 = triangle.GetVertex(1)[0] - triangle.GetVertex(2)[0];
        x10 = triangle.GetVertex(0)[2] - triangle.GetVertex(2)[2];
        x11 = triangle.GetVertex(1)[2] - triangle.GetVertex(2)[2];
        y0 = point[0] - triangle.GetVertex(2)[0];
        y1 = point[2] - triangle.GetVertex(2)[2];
        double alpha2 = (y0 * x10 - y1 * x00) / (x01 * x10 - x11 * x00 + eps);
        double alpha1 = 1;
        if (std::fabs(x00) > eps) {
            alpha1 = (y0 - alpha2 * x01) / (x00 + eps);
        } else {
            alpha1 = (y1 - alpha2 * x11) / (x10 + eps);
        }
        return Vector({alpha1,  1 - alpha1 - alpha2, alpha2});
    }
}
