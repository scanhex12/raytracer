#pragma once

#include <image.h>
#include <camera_options.h>
#include <render_options.h>
#include <string>
#include <scene.h>
#include <ray.h>
#include <intersection.h>
#include <geometry.h>
#include <cmath>
#include <map>
#include <vector.h>
#include <cassert>

std::vector<Vector> RotateCamera(Vector from, Vector to) {
    double eps = 1e-8;
    Vector diff = Difference(to, from);
    diff.Normalize();
    Vector basis_x = CrossProduct(Vector({0, 1, 0}), diff);
    if (Length(basis_x) < eps) {
        basis_x = Vector({-1, 0, 0});
    }
    basis_x.Normalize();
    Vector basis_y = CrossProduct(diff, basis_x);
    if (Length(basis_y) < eps) {
        basis_y = Vector({0, -1, 0});
    }
    basis_y.Normalize();
    assert(DotProduct(basis_y, diff) < eps);
    assert(DotProduct(basis_x, diff) < eps);
    assert(DotProduct(basis_y, basis_x) < eps);

    return {basis_y, basis_x, diff};
}

Image PostProcessing(std::vector<std::vector<std::vector<double>>>& image) {
    double c;
    for (size_t i = 0; i < image.size(); ++i) {
        for (size_t j = 0; j < image[0].size(); ++j) {
            for (int k = 0; k < 3; ++k) {
                c = std::max(c, image[i][j][k]);
            }
        }
    }
    for (size_t i = 0; i < image.size(); ++i) {
        for (size_t j = 0; j < image[0].size(); ++j) {
            for (int k = 0; k < 3; ++k) {
                image[i][j][k] =
                    (image[i][j][k] * (1 + (image[i][j][k] / (c * c)))) / (1 + image[i][j][k]);
                image[i][j][k] = std::pow(image[i][j][k], 1.0 / 2.2);
                image[i][j][k] *= 255;
            }
        }
    }
    Image answer = Image(image[0].size(), image.size());
    for (size_t i = 0; i < image.size(); ++i) {
        for (size_t j = 0; j < image[0].size(); ++j) {
            answer.SetPixel(RGB(image[i][j][0], image[i][j][1], image[i][j][2]), i, j);
        }
    }
    return answer;
}

Vector CalculateIntensivity(Ray ray_normal, std::vector<Object>& objects,
                            std::vector<SphereObject>& spheres, std::vector<Light>& lights,
                            size_t max_recursive_height, double n = 1.,
                            size_t recursive_height = 0) {
    if (max_recursive_height == recursive_height) {
        return Vector({0, 0, 0});
    }
    Vector intersection_point;
    Vector normal_intersection;
    Material material_intersection;
    double d = 1e9;
    ray_normal.Step(1e-6);
    Vector l_from = ray_normal.GetOrigin();
    bool inside = false;
    for (size_t k = 0; k < objects.size(); ++k) {
        auto x = GetIntersection(ray_normal, objects[k].polygon);
        if (x == std::nullopt) {
            continue;
        }
        if (x->GetDistance() < d) {
            d = x->GetDistance();
            intersection_point = x->GetPosition();
            if (!objects[k].normals[0]) {
                normal_intersection = x->GetNormal();
            } else {
                Vector coefs = GetBarycentricCoords(objects[k].polygon, x->GetPosition());
                normal_intersection = Vector({0, 0, 0});
                for (size_t ind = 0; ind < 3; ++ind) {
                    normal_intersection =
                        Sum(normal_intersection, Mult(*objects[k].normals[ind], coefs[ind]));
                }
                normal_intersection.Normalize();
                if (DotProduct(x->GetNormal(), normal_intersection) < 1e-8) {
                    normal_intersection = Mult(normal_intersection, -1);
                }
            }
            material_intersection = *objects[k].material;
        }
    }
    for (size_t k = 0; k < spheres.size(); ++k) {
        auto x = GetIntersection(ray_normal, spheres[k].sphere);
        if (x == std::nullopt) {
            continue;
        }
        if (x->GetDistance() < d) {
            inside = true;
            d = x->GetDistance();
            intersection_point = x->GetPosition();
            normal_intersection = x->GetNormal();
            material_intersection = *spheres[k].material;
        }
    }
    if (d > 9e8) {
        return Vector({0, 0, 0});
    }
    Vector i_result = Sum(material_intersection.ambient_color, material_intersection.intensity);
    {
        for (auto cur_light : lights) {
            auto direction = Difference(cur_light.position, intersection_point);
            Ray ray(intersection_point, direction);
            ray.Step(1e-5);
            double min_d = Length(direction);
            for (size_t id_obj = 0; id_obj < objects.size(); ++id_obj) {
                auto q = GetIntersection(ray, objects[id_obj].polygon);
                if (!q) {
                    continue;
                }
                min_d = std::min(min_d, q->GetDistance());
            }
            for (size_t id_obj = 0; id_obj < spheres.size(); ++id_obj) {
                auto q = GetIntersection(ray, spheres[id_obj].sphere);
                if (!q) {
                    continue;
                }
                min_d = std::min(min_d, q->GetDistance());
            }
            if (std::fabs(min_d - Length(direction)) < 5e-6) {
                direction = ray.GetDirection();
                direction.Normalize();
                normal_intersection.Normalize();
                double scalar_coef = std::max(0.0, DotProduct(normal_intersection, direction));
                Vector result_intensivity = Mult(cur_light.intensity, scalar_coef);
                for (size_t dim = 0; dim < 3; ++dim) {
                    result_intensivity[dim] *= material_intersection.diffuse_color[dim];
                }
                i_result =
                    Sum(i_result, (Mult(result_intensivity, material_intersection.albedo[0])));

                Vector reflected = Reflect(Mult(direction, -1), normal_intersection);
                reflected.Normalize();
                Vector ve = Difference(l_from, intersection_point);
                ve.Normalize();
                double cos_theta = DotProduct(reflected, ve);
                scalar_coef =
                    std::pow(std::max(0.0, cos_theta), material_intersection.specular_exponent);
                result_intensivity = Mult(cur_light.intensity, scalar_coef);
                for (size_t dim = 0; dim < 3; ++dim) {
                    result_intensivity[dim] *= material_intersection.specular_color[dim];
                }
                i_result =
                    Sum(i_result, (Mult(result_intensivity, material_intersection.albedo[0])));
            }
        }
    }
    if (n < 1 + 1e-4) {
        Vector start_dir = ray_normal.GetDirection();
        start_dir.Normalize();
        Vector reflect = Reflect(start_dir, normal_intersection);
        reflect.Normalize();
        Ray refracted_ray = Ray(intersection_point, reflect);
        Vector refracted_intensivity = CalculateIntensivity(
            refracted_ray, objects, spheres, lights, max_recursive_height, n, recursive_height + 1);

        i_result = Sum(i_result, Mult(refracted_intensivity, material_intersection.albedo[1]));
    }
    double next_refract_coef = material_intersection.refraction_index;
    if (n > 1 + 1e-4) {
        next_refract_coef = 1;
    }
    auto refract =
        Refract(ray_normal.GetDirection(), normal_intersection, n / (next_refract_coef + 0.0));
    if (refract) {
        refract->Normalize();
        Ray refracted_ray = Ray(intersection_point, *refract);
        Vector refracted_intensivity =
            CalculateIntensivity(refracted_ray, objects, spheres, lights, max_recursive_height,
                                 next_refract_coef, recursive_height + 1);

        double cur_tr = (n > 1 + 1e-4 && inside) ? 1 : material_intersection.albedo[2];
        refracted_intensivity = Mult(refracted_intensivity, cur_tr);
        i_result = Sum(i_result, refracted_intensivity);
    }

    return i_result;
}

Image Render(const std::string& filename, const CameraOptions& camera_options,
             const RenderOptions& render_options) {
    double h = camera_options.screen_height;
    double w = camera_options.screen_width;

    Image answer = Image(w, h);
    Scene scene = ReadScene(filename);
    auto objects = scene.GetObjects();
    auto spheres = scene.GetSphereObjects();
    auto lights = scene.GetLights();

    if (render_options.mode == RenderMode::kDepth) {
        Vector l_from = Vector(camera_options.look_from);
        Vector l_to = Vector(camera_options.look_to);

        auto basises = RotateCamera(l_from, l_to);

        std::vector<std::vector<double>> pixels(h, std::vector<double>(w));
        double scale = 2 * std::tan(camera_options.fov / 2) / (std::min(h, w) + 0.0);
        double max_d = 0;
        for (size_t i = 0; i < h; ++i) {
            for (size_t j = 0; j < w; ++j) {
                double d = 1e9;

                Vector direction = Vector({0, 0, 0});
                direction = Sum(direction, Mult(basises[0], (i - (h - 1.0) / 2.0) * scale));
                direction = Sum(direction, Mult(basises[1], (j - (w - 1.0) / 2.0) * scale));
                direction = Sum(direction, Mult(basises[2], 1));

                Vector start = direction;
                start.Normalize();

                Ray ray_normal = Ray(l_from, start);
                for (size_t k = 0; k < objects.size(); ++k) {
                    auto x = GetIntersection(ray_normal, objects[k].polygon);
                    if (x == std::nullopt) {
                        continue;
                    }
                    if (x->GetDistance() < d) {
                        d = x->GetDistance();
                    }
                }
                for (size_t k = 0; k < spheres.size(); ++k) {
                    auto x = GetIntersection(ray_normal, spheres[k].sphere);
                    if (x == std::nullopt) {
                        continue;
                    }
                    if (x->GetDistance() < d) {
                        d = x->GetDistance();
                    }
                }
                if (d > 9e8) {
                    continue;
                }
                max_d = std::max(max_d, d);
                pixels[i][j] = d;
            }
        }
        for (size_t i = 0; i < h; ++i) {
            for (size_t j = 0; j < w; ++j) {
                if (std::fabs(pixels[h - i - 1][w - j - 1]) < 1e-9) {
                    pixels[h - i - 1][w - j - 1] = max_d;
                }
                pixels[h - i - 1][w - j - 1] /= max_d;
                answer.SetPixel(
                    RGB(pixels[h - i - 1][w - j - 1] * 255, pixels[h - i - 1][w - j - 1] * 255,
                        pixels[h - i - 1][w - j - 1] * 255),
                    i, j);
            }
        }
        return answer;
    }

    if (render_options.mode == RenderMode::kNormal) {
        Vector l_from = Vector(camera_options.look_from);
        Vector l_to = Vector(camera_options.look_to);

        auto basises = RotateCamera(l_from, l_to);

        double scale = 2 * std::tan(camera_options.fov / 2) / (std::min(h, w) + 0.0);
        double max_d = 0;
        for (size_t i = 0; i < h; ++i) {
            for (size_t j = 0; j < w; ++j) {
                double d = 1e9;

                Vector direction = Vector({0, 0, 0});
                direction = Sum(direction, Mult(basises[0], (i - (h - 1.0) / 2.0) * scale));
                direction = Sum(direction, Mult(basises[1], (j - (w - 1.0) / 2.0) * scale));
                direction = Sum(direction, Mult(basises[2], 1));

                Vector start = direction;
                start.Normalize();

                Vector normal_intersection = Vector({0, 0, 0});

                Ray ray_normal = Ray(l_from, start);
                for (size_t k = 0; k < objects.size(); ++k) {
                    auto x = GetIntersection(ray_normal, objects[k].polygon);
                    if (x == std::nullopt) {
                        continue;
                    }
                    if (x->GetDistance() < d) {
                        d = x->GetDistance();
                        if (!objects[k].normals[0]) {
                            normal_intersection = x->GetNormal();
                            normal_intersection.Normalize();
                        } else {
                            Vector coefs =
                                GetBarycentricCoords(objects[k].polygon, x->GetPosition());
                            normal_intersection = Vector({0, 0, 0});
                            for (size_t ind = 0; ind < 3; ++ind) {
                                normal_intersection =
                                    Sum(normal_intersection,
                                        Mult(*objects[k].normals[ind], coefs[ind]));
                            }
                        }
                    }
                }
                for (size_t k = 0; k < spheres.size(); ++k) {
                    auto x = GetIntersection(ray_normal, spheres[k].sphere);
                    if (x == std::nullopt) {
                        continue;
                    }
                    if (x->GetDistance() < d) {
                        d = x->GetDistance();
                        normal_intersection = x->GetNormal();
                        normal_intersection.Normalize();
                    }
                }
                if (d > 9e8) {
                    answer.SetPixel(RGB(0.0, 0.0, 0.0), h - i - 1, w - j - 1);
                    continue;
                }
                max_d = std::max(max_d, d);
                answer.SetPixel(RGB((0.5 + 0.5 * normal_intersection[0]) * 255,
                                    (0.5 + 0.5 * normal_intersection[1]) * 255,
                                    (0.5 + 0.5 * normal_intersection[2]) * 255),
                                h - i - 1, w - j - 1);
            }
        }
        return answer;
    }
    if (render_options.mode == RenderMode::kFull) {
        Vector l_from = Vector(camera_options.look_from);
        Vector l_to = Vector(camera_options.look_to);

        auto basises = RotateCamera(l_from, l_to);

        std::vector<std::vector<std::vector<double>>> pixels(
            h, std::vector<std::vector<double>>(w, std::vector<double>(3)));
        double scale = 2 * std::tan(camera_options.fov / 2) / (std::min(h, w) + 0.0);
        for (size_t i = 0; i < h; ++i) {
            for (size_t j = 0; j < w; ++j) {
                Vector direction = Vector({0, 0, 0});
                direction = Sum(direction, Mult(basises[0], (i - (h - 1.0) / 2.0) * scale));
                direction = Sum(direction, Mult(basises[1], (j - (w - 1.0) / 2.0) * scale));
                direction = Sum(direction, Mult(basises[2], 1));

                Vector start = direction;
                start.Normalize();
                Ray ray_normal = Ray(l_from, start);
                Vector i_result = CalculateIntensivity(ray_normal, objects, spheres, lights,
                                                       render_options.depth);
                pixels[h - i - 1][w - j - 1][0] = i_result[0];
                pixels[h - i - 1][w - j - 1][1] = i_result[1];
                pixels[h - i - 1][w - j - 1][2] = i_result[2];
            }
        }
        return PostProcessing(pixels);
    }
    throw "Not implemented";
}
