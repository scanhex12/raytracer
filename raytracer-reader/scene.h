#pragma once

#include <material.h>
#include <vector.h>
#include <object.h>
#include <light.h>
#include <reader_utils.h>

#include <vector>
#include <string>
#include <map>
#include <string>
#include <fstream>
#include <sstream>

class Scene {
public:
    Scene(std::vector<Object> objects, std::vector<SphereObject> spheres, std::vector<Light> lights,
          std::map<std::string, Material> materials)
        : objects_(objects), spheres_(spheres), lights_(lights), materials_(materials) {
    }

    const std::vector<Object>& GetObjects() const {
        return objects_;
    }
    const std::vector<SphereObject>& GetSphereObjects() const {
        return spheres_;
    }
    const std::vector<Light>& GetLights() const {
        return lights_;
    }
    const std::map<std::string, Material>& GetMaterials() const {
        return materials_;
    }

private:
    std::vector<Object> objects_;
    std::vector<SphereObject> spheres_;
    std::vector<Light> lights_;
    std::map<std::string, Material> materials_;
};

inline std::map<std::string, Material> ReadMaterials(std::string filename) {
    return FileMaterials(filename);
}
inline Scene ReadScene(std::string filename) {
    std::string mtl_filename = GetMtlFileName(filename);
    auto materials = FileMaterials(mtl_filename);
    auto objects = FileObjects(filename, materials);
    auto sphere_objects = FileSpheres(filename, materials);
    auto lights = FileLights(filename);
    return Scene(objects, sphere_objects, lights, materials);
}
