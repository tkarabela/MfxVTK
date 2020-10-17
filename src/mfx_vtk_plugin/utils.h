#pragma once

#include <cmath>
#include <array>

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
