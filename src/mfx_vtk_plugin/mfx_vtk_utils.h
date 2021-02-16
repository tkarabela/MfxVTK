/*
MfxVTK Open Mesh Effect plug-in
Copyright (c) 2020 Tomas Karabela

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#pragma once

#include <cmath>
#include <array>

static inline constexpr bool is_positive_double(double x) {
    return x >= DBL_EPSILON;
}

template <typename T>
static inline constexpr T clamp(T value, T minimum, T maximum) {
    return std::min(maximum, std::max(minimum, value));
}

template <typename T>
static inline constexpr T vec3_squared_distance(const T x[3], const T y[3]) {
    T acc = 0;
    for (int i = 0; i < 3; i++) {
        acc += (x[i]-y[i]) * (x[i]-y[i]);
    }
    return acc;
}

template <typename T>
static inline constexpr T vec3_dot(const T x[3], const T y[3]) {
    T acc = 0;
    for (int i = 0; i < 3; i++) {
        acc += x[i]*y[i];
    }
    return acc;
}

template <typename T>
static inline constexpr void vec3_min(T x[3], const T y[3]) {
    for (int i = 0; i < 3; i++) {
        x[i] = std::min(x[i], y[i]);
    }
}

template <typename T>
static inline constexpr void vec3_max(T x[3], const T y[3]) {
    for (int i = 0; i < 3; i++) {
        x[i] = std::max(x[i], y[i]);
    }
}

// TODO improve, it's a bit too much dependent
template <int N>
class AdditiveRecurrence {
public:
    AdditiveRecurrence() {
        for (int i = 0; i < N; i++) {
            s[i] = 0.0;
        }
        Next();
    }

    double GetValue(int i) const {
        return s[i];
    }

    double GetRangeValue(int i, double low, double high) const {
        return (GetValue(i) * (high - low)) + low;
    }

    void Next() {
        const double alpha[8] = {
                1.4142135623730951, // sqrt(2)
                1.7320508075688772, // sqrt(3)
                2.2360679774997898, // sqrt(5)
                2.6457513110645907, // sqrt(7)
                3.3166247903553998, // sqrt(11)
                3.6055512754639891, // sqrt(13)
                4.1231056256176606, // sqrt(17)
                4.3588989435406740, // sqrt(19)
        };

        double dummy;
        for (int i = 0; i < N; i++) {
            s[i] = modf(alpha[i] + s[i], &dummy); // (s_i * alpha_i) mod 1
        }
    }

protected:
    double s[N];
    static_assert(N > 0 && N < 8, "wrong dimension");
};
