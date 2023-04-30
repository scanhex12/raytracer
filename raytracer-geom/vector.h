#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <initializer_list>
#include <algorithm>

class Vector {
public:
    Vector() {
        for (int i = 0; i < 3; ++i) {
            data_[i] = 0;
        }
    }
    Vector(std::initializer_list<double> list) {
        auto it = list.begin();
        for (int i = 0; i < 3; ++i) {
            data_[i] = *it;
            it++;
        }
    }
    Vector(std::array<double, 3> data) : data_(data) {
    }

    double& operator[](size_t ind) {
        return data_[ind];
    }
    double operator[](size_t ind) const {
        return data_[ind];
    }

    void Normalize() {
        double len = sqrt(data_[0] * data_[0] + data_[1] * data_[1] + data_[2] * data_[2]);
        data_[0] /= len;
        data_[1] /= len;
        data_[2] /= len;
    }

    Vector(const Vector& another) {
        data_[0] = another[0];
        data_[1] = another[1];
        data_[2] = another[2];
    }

    Vector& operator=(const Vector& another) {
        data_[0] = another[0];
        data_[1] = another[1];
        data_[2] = another[2];
        return *this;
    }

private:
    std::array<double, 3> data_;
};

inline double DotProduct(const Vector& lhs, const Vector& rhs) {
    double answer = 0;
    for (size_t i = 0; i < 3; ++i) {
        answer += lhs[i] * rhs[i];
    }
    return answer;
}
inline Vector CrossProduct(const Vector& a, const Vector& b) {
    return Vector(
        {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]});
}
inline double Length(const Vector& vec) {
    double ans = 0;
    for (int i = 0; i < 3; ++i) {
        ans += vec[i] * vec[i];
    }
    return sqrt(ans);
}

inline Vector Difference(const Vector& a, const Vector& b) {
    return Vector({a[0] - b[0], a[1] - b[1], a[2] - b[2]});
}

inline Vector Sum(const Vector& a, const Vector& b) {
    return Vector({a[0] + b[0], a[1] + b[1], a[2] + b[2]});
}

inline Vector Mult(const Vector& a, double x) {
    return Vector({a[0] * x, a[1] * x, a[2] * x});
}