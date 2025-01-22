#pragma once

//Derived from DFTCXX
//https://github.com/ifilot/dftcxx.git

namespace SCESS {
    //Define Lebedev orders
    enum class LebedevOrder :int {
        LEBEDEV_6 = 6,
        LEBEDEV_14 = 14,
        LEBEDEV_26 = 26,
        LEBEDEV_38 = 38,
        LEBEDEV_50 = 50,
        LEBEDEV_74 = 74,
        LEBEDEV_86 = 86,
        LEBEDEV_110 = 110,
        LEBEDEV_146 = 146,
        LEBEDEV_170 = 170,
        LEBEDEV_194 = 194
    };

    namespace Quadrature {
        //Number of coefficients per order
        constexpr std::array<uint32_t, 11> nLebedevPoints = {6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194};

        // Total size of lebedev_coefficient matrix
        constexpr uint32_t nLebedevCoefficients = []() constexpr {
            uint32_t sum = 0;
            for (uint32_t nLebedevPoint: nLebedevPoints)
                sum += nLebedevPoint;

            return sum;
        }();

        // Define Lebedev coefficients
        extern double LEBEDEV_COEFFICIENTS[nLebedevCoefficients][4];
    };
}
