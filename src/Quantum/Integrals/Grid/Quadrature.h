#pragma once

//Derived from DFTCXX
//https://github.com/ifilot/dftcxx.git

namespace SCESS {
    //Define Lebedev orders
    enum class LebedevOrder :int {
        LEBEDEV_6 = 0,
        LEBEDEV_14,
        LEBEDEV_26,
        LEBEDEV_38,
        LEBEDEV_50,
        LEBEDEV_74,
        LEBEDEV_86,
        LEBEDEV_110,
        LEBEDEV_146,
        LEBEDEV_170,
        LEBEDEV_194,

        NUM_LEBEDEV_POINTS
    };

    class Quadrature {
    public:
        //Number of coefficients per order
        static constexpr std::array<uint32_t, 11> nLebedevPoints = {6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194};

        // Total size of lebedev_coefficient matrix
        static constexpr uint32_t nLebedevCoefficients = []() constexpr {
            uint32_t sum = 0;
            for (uint32_t nLebedevPoint: nLebedevPoints)
                sum += nLebedevPoint;

            return sum;
        }();

        // Define Lebedev coefficients
        static double LEBEDEV_COEFFICIENTS[nLebedevCoefficients][4];
    };
}
