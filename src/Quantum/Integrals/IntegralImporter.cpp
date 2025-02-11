#include<vector>
#include<string>
#include<boost/algorithm/string.hpp>
#include<filesystem>
#include<fstream>
#include "GraphicsEngine/Tools/Strings.h"
#include<unsupported/Eigen/CXX11/Tensor>

#include "IntegralsF12.h"
#include"Integrals.h"


//
//
namespace SCESS::Integrals {
    using namespace Eigen;

    void F12::ImportF12IntegralsMRCC(const std::string &fileName, uint32_t nbf, F12Integrals& f12Ints, TEITensor& eri) {
        if (!std::filesystem::exists(fileName)) {
            throw std::runtime_error("Cannot import a F12 integral file; The file does not exist");
        }

        std::fstream file(fileName.data(), std::ios::in);

        if (!file.is_open())
            throw std::runtime_error("Cannot open file");


        array<Index, 4> dimensions{nbf, nbf, nbf, nbf};
        eri = Eigen::Tensor<double, 4>(dimensions);
        f12Ints.f12 = Eigen::Tensor<double, 4>(dimensions);
        f12Ints.f12Sq = Eigen::Tensor<double, 4>(dimensions);
        f12Ints.f12DelSq = Eigen::Tensor<double, 4>(dimensions);
        f12Ints.f12Repulsion = Eigen::Tensor<double, 4>(dimensions);

        std::string line;
        std::string intType;
        while (std::getline(file, line)) {
            if ((line.find("vint") != std::string::npos) ||
                (line.find("f12int") != std::string::npos) ||
                (line.find("f122int") != std::string::npos) ||
                (line.find("delf122int") != std::string::npos) ||
                (line.find("f12r12int") != std::string::npos)) {
                boost::trim(line);
                intType = line;
                continue;
            }

            if (!intType.empty()) {
                //Remove whitespace from the start and the end
                boost::trim(line);

                //Checking if the end of the integral block section has been found
                if (!std::isdigit(line[0]) && line[0] != '-')
                    break;

                auto tokens = GEng::Strings::Split(line, ' ');

                double integralVal = std::stod(tokens[0]);
                array<Index, 4> orbIndices = {
                    std::stoi(tokens[1]),
                    std::stoi(tokens[2]),
                    std::stoi(tokens[3]),
                    std::stoi(tokens[4])
                };

                if (intType == "vint") {
                    eri(orbIndices) = integralVal;
                } else if (intType == "f12int") {
                    f12Ints.f12(orbIndices) = integralVal;
                } else if (intType == "f122int") {
                    f12Ints.f12Sq(orbIndices) = integralVal;
                } else if (intType == "delf122int") {
                    f12Ints.f12DelSq(orbIndices) = integralVal;
                } else if (intType == "f12r12int") {
                    f12Ints.f12Repulsion(orbIndices) = integralVal;
                }
            }
        }
        file.close();

    }
}
