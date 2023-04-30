#pragma once

#include <object.h>
#include <vector>
#include <string>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <memory>

void PreProcessString(std::string& s) {
    std::string final;
    for (auto elem : s) {
        if (elem == '\t') {
            final.push_back(' ');
        } else {
            final.push_back(elem);
        }
    }
    s = final;
}

bool Count(const std::string& x, char c) {
    for (auto elem : x) {
        if (elem == c) {
            return true;
        }
    }
    return false;
}

bool CountSubstr(const std::string& x, std::string t) {
    for (size_t i = 0; i < x.size() - t.size() + 1; ++i) {
        if (x.substr(i, t.size()) == t) {
            return true;
        }
    }
    return false;
}

std::vector<std::string> Split(const std::string& strin, const std::string& delimiter = " ") {
    std::vector<std::string> answer;
    std::string back;
    for (size_t i = 0; i < strin.size();) {
        if (i < strin.size() - delimiter.size() + 1 &&
            strin.substr(i, delimiter.size()) == delimiter) {
            if (!back.empty()) {
                answer.push_back(back);
            }
            back = "";
            i += delimiter.size();
        } else {
            back.push_back(strin[i]);
            i++;
        }
    }
    answer.push_back(back);
    return answer;
}

std::pair<bool, double> Str2Double(const std::string& s) {
    try {
        std::pair<bool, double> answer = {true, std::stod(s)};
        return answer;
    } catch (...) {
        return {false, 0};
    }
}

std::vector<Vector> GetPoints(std::string filename) {
    std::ifstream myfile(filename);
    myfile.tie(nullptr);
    std::string line;
    std::vector<Vector> answer;
    while (std::getline(myfile, line)) {
        PreProcessString(line);
        auto objects = Split(line);
        if (objects[0] != "v") {
            continue;
        }
        answer.push_back(Vector({Str2Double(objects[1]).second, Str2Double(objects[2]).second,
                                 Str2Double(objects[3]).second}));
    }
    return answer;
}

std::vector<std::shared_ptr<Vector>> GetNormals(std::string filename) {
    std::ifstream myfile(filename);
    myfile.tie(nullptr);
    std::string line;
    std::vector<std::shared_ptr<Vector>> answer;
    while (std::getline(myfile, line)) {
        PreProcessString(line);
        auto objects = Split(line);
        if (objects[0] != "vn") {
            continue;
        }
        answer.push_back(std::shared_ptr<Vector>(
            new Vector({Str2Double(objects[1]).second, Str2Double(objects[2]).second,
                        Str2Double(objects[3]).second})));
    }
    return answer;
}

std::map<std::string, Material> FileMaterials(std::string filename) {
    std::ifstream myfile(filename);
    myfile.tie(nullptr);
    std::string line;
    std::map<std::string, Material> answer;

    std::string material_name;
    Vector ka = Vector({0, 0, 0}), kd = Vector({0, 0, 0}), ks = Vector({0, 0, 0}),
           ke = Vector({0, 0, 0});
    double ns = 1, ni = 1;
    double tr = 1;
    std::array<double, 3> albedo = {0, 0, 0};
    albedo[0] = 1;

    while (std::getline(myfile, line)) {
        PreProcessString(line);
        auto objects = Split(line);
        if (objects[0] == "newmtl") {
            if (!material_name.empty()) {
                answer[material_name] = Material(material_name, ka, kd, ks, ke, ns, ni, albedo, tr);
                material_name = "";
                ka = Vector();
                kd = Vector();
                ks = Vector();
                ke = Vector();
                ns = 1, ni = 1;
                albedo[0] = 1;
                albedo[1] = 0;
                albedo[2] = 0;
                tr = 1;
            }
            material_name = objects[1];
            continue;
        }
        if (objects[0] == "Ka") {
            ka = Vector({Str2Double(objects[1]).second, Str2Double(objects[2]).second,
                         Str2Double(objects[3]).second});
            continue;
        }
        if (objects[0] == "Kd") {
            kd = Vector({Str2Double(objects[1]).second, Str2Double(objects[2]).second,
                         Str2Double(objects[3]).second});
            continue;
        }
        if (objects[0] == "Ks") {
            ks = Vector({Str2Double(objects[1]).second, Str2Double(objects[2]).second,
                         Str2Double(objects[3]).second});
            continue;
        }
        if (objects[0] == "Ke") {
            ke = Vector({Str2Double(objects[1]).second, Str2Double(objects[2]).second,
                         Str2Double(objects[3]).second});
            continue;
        }
        if (objects[0] == "Ns") {
            ns = Str2Double(objects[1]).second;
            continue;
        }
        if (objects[0] == "Ni") {
            ni = Str2Double(objects[1]).second;
            continue;
        }
        if (objects[0] == "al") {
            albedo[0] = std::stod(objects[1]);
            albedo[1] = std::stod(objects[2]);
            albedo[2] = std::stod(objects[3]);
            continue;
        }
    }
    if (!material_name.empty()) {
        answer[material_name] = Material(material_name, ka, kd, ks, ke, ns, ni, albedo, tr);
    }

    return answer;
}

std::vector<Object> FileObjects(std::string filename,
                                std::map<std::string, Material>& all_materials) {
    std::vector<Vector> all_points = GetPoints(filename);
    std::vector<std::shared_ptr<Vector>> all_norms = GetNormals(filename);
    size_t cur_poins = 0, cur_norms = 0;

    std::ifstream myfile(filename);
    std::string line;
    std::vector<Object> answer;

    std::string last_material;
    while (std::getline(myfile, line)) {
        PreProcessString(line);
        auto objects = Split(line);
        if (objects[0] == "v") {
            cur_poins++;
            continue;
        }
        if (objects[0] == "vn") {
            cur_norms++;
            continue;
        }
        if (objects[0] == "usemtl") {
            last_material = objects[1];
        }
        if (objects[0] == "f") {
            std::vector<Vector> triangle_points;
            std::vector<std::shared_ptr<Vector>> triangle_norms;
            for (size_t j = 1; j < objects.size(); ++j) {
                if (objects[j].empty()) {
                    continue;
                }
                if (Count(objects[j], '/')) {
                    if (CountSubstr(objects[j], "//")) {
                        auto point_normal = Split(objects[j], "//");
                        int point_lag = std::stoi(point_normal[0]);
                        int norm_lag = std::stoi(point_normal[1]);

                        if (point_lag > 0) {
                            --point_lag;
                        } else {
                            point_lag = cur_poins + point_lag;
                        }

                        if (norm_lag > 0) {
                            --norm_lag;
                        } else {
                            norm_lag = cur_norms + norm_lag;
                        }
                        triangle_points.push_back(all_points[point_lag]);
                        triangle_norms.push_back(all_norms[norm_lag]);
                    } else {
                        auto point_normal = Split(objects[j], "/");
                        int point_lag = std::stoi(point_normal[0]);
                        int norm_lag = std::stoi(point_normal[2]);

                        if (point_lag > 0) {
                            --point_lag;
                        } else {
                            point_lag = cur_poins + point_lag;
                        }

                        if (norm_lag > 0) {
                            --norm_lag;
                        } else {
                            norm_lag = cur_norms + norm_lag;
                        }
                        triangle_points.push_back(all_points[point_lag]);
                        triangle_norms.push_back(all_norms[norm_lag]);
                    }
                } else {
                    int point_lag = std::stoi(objects[j]);
                    if (point_lag > 0) {
                        --point_lag;
                    } else {
                        point_lag = cur_poins + point_lag;
                    }
                    triangle_points.push_back(all_points[point_lag]);
                    triangle_norms.push_back(nullptr);
                }
            }
            for (size_t j = 2; j < triangle_points.size(); ++j) {
                answer.push_back(
                    Object(std::shared_ptr<Material>(new Material(all_materials[last_material])),
                           triangle_points[0], triangle_points[j - 1], triangle_points[j],
                           triangle_norms[0], triangle_norms[j - 1], triangle_norms[j]));
            }
        }
    }
    return answer;
}

std::vector<SphereObject> FileSpheres(std::string filename,
                                      std::map<std::string, Material>& all_materials) {
    std::ifstream myfile(filename);
    myfile.tie(nullptr);
    std::string line;
    std::vector<SphereObject> answer;

    std::string last_material;
    while (std::getline(myfile, line)) {
        PreProcessString(line);
        auto objects = Split(line);
        if (objects[0] == "usemtl") {
            last_material = objects[1];
            continue;
        }
        if (objects[0] == "S") {
            double x = std::stod(objects[1]);
            double y = std::stod(objects[2]);
            double z = std::stod(objects[3]);
            double r = std::stod(objects[4]);
            answer.push_back(SphereObject(
                std::shared_ptr<Material>(new Material(all_materials[last_material])), x, y, z, r));
        }
    }
    return answer;
}

std::vector<Light> FileLights(std::string filename) {
    std::ifstream myfile(filename);
    myfile.tie(nullptr);
    std::string line;
    std::vector<Light> answer;

    std::string last_material;
    while (std::getline(myfile, line)) {
        PreProcessString(line);
        auto objects = Split(line);
        if (objects[0] == "P") {
            double x = std::stod(objects[1]);
            double y = std::stod(objects[2]);
            double z = std::stod(objects[3]);
            double r = std::stod(objects[4]);
            double g = std::stod(objects[5]);
            double b = std::stod(objects[6]);
            answer.push_back(Light(x, y, z, r, g, b));
        }
    }
    return answer;
}

std::string GetMtlFileName(std::string obj_filename) {
    std::ifstream myfile(obj_filename);
    myfile.tie(nullptr);
    std::string line;

    std::string dir = obj_filename;
    while (dir.back() != '/') {
        dir.pop_back();
    }

    while (std::getline(myfile, line)) {
        PreProcessString(line);
        auto objects = Split(line);
        if (objects[0] == "mtllib") {
            dir += objects[1];
        }
    }
    return dir;
}
