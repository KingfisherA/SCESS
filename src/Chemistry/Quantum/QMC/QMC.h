#pragma once
#include <cstdint>
#include <string>
#include<vector>


namespace SCESS {
    class RHFResult;

    struct Particle {
        public:
        int population;
    };
class QuantumMonteCarlo {
public:
    //Full CI quantum monte carlo
    static double SimpleFCIQMC(const RHFResult& rhfWfn, double Vnn, std::stringstream* console = nullptr);

};

} // SCESS
