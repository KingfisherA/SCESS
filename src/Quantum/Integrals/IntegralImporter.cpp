#include<vector>
#include<string>
#include<boost/algorithm/string.hpp>
#include<filesystem>
#include<fstream>
#include<GraphicsEngine/Utilities/StringOperations.h>
#include<unsupported/Eigen/CXX11/Tensor>

#include "IntegralsF12.h"
#include"Integrals.h"


//
//
namespace SCESS::Integrals {
    using namespace Eigen;

    F12::F12Integrals F12::ImportF12IntegralsMRCC(const std::string &fileName, uint32_t nbf) {
        if (!std::filesystem::exists(fileName)) {
            throw std::runtime_error("Cannot import a F12 integral file; The file does not exist");
        }

        std::fstream file(fileName.data(), std::ios::in);

        if (!file.is_open())
            throw std::runtime_error("Cannot open file");

        std::vector<std::pair<Eigen::Vector4i, double> > eriInts, f12Ints, f12SqInts, f12DelSqInts, f12RepulsionInts;

        std::string line;
        std::string intType = "";
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

            if (intType != "") {
                //Remove whitespace from the start and the end
                boost::trim(line);

                //Checking if the end of the integral block section has been found
                if (!std::isdigit(line[0]) && line[0] != '-')
                    break;

                auto tokens = GEng::StringOperations::Split(line, ' ');

                double integralVal = std::stod(tokens[0]);
                Eigen::Vector4i orbIndices;
                orbIndices.x() = std::stoi(tokens[1]);
                orbIndices.y() = std::stoi(tokens[2]);
                orbIndices.z() = std::stoi(tokens[3]);
                orbIndices.w() = std::stoi(tokens[4]);

                if (intType == "vint") {
                    eriInts.emplace_back(orbIndices, integralVal);
                } else if (intType == "f12int") {
                    f12Ints.emplace_back(orbIndices, integralVal);
                } else if (intType == "f122int") {
                    f12SqInts.emplace_back(orbIndices, integralVal);
                } else if (intType == "delf122int") {
                    f12DelSqInts.emplace_back(orbIndices, integralVal);
                } else if (intType == "f12r12int") {
                    f12RepulsionInts.emplace_back(orbIndices, integralVal);
                }
            }
        }
        file.close();

        Eigen::array<Eigen::Index, 4> dimensions{nbf, nbf, nbf, nbf};
        F12Integrals returnInts;
        returnInts.eri = Eigen::Tensor<double, 4>(dimensions);
        returnInts.f12 = Eigen::Tensor<double, 4>(dimensions);
        returnInts.f12Sq = Eigen::Tensor<double, 4>(dimensions);
        returnInts.f12DelSq = Eigen::Tensor<double, 4>(dimensions);
        returnInts.f12Repulsion = Eigen::Tensor<double, 4>(dimensions);

        for (const auto& ele : eriInts) {
            const Vector4i& orbIndices = ele.first;
            returnInts.eri(orbIndices.x(), orbIndices.y(), orbIndices.z(), orbIndices.w()) = ele.second;
        }

        for (const auto& ele : f12Ints) {
            const Vector4i& orbIndices = ele.first;
            returnInts.f12(orbIndices.x(), orbIndices.y(), orbIndices.z(), orbIndices.w()) = ele.second;
        }

        for (const auto& ele : f12SqInts) {
            const Vector4i& orbIndices = ele.first;
            returnInts.f12Sq(orbIndices.x(), orbIndices.y(), orbIndices.z(), orbIndices.w()) = ele.second;
        }

        for (const auto& ele : f12DelSqInts) {
            const Vector4i& orbIndices = ele.first;
            returnInts.f12DelSq(orbIndices.x(), orbIndices.y(), orbIndices.z(), orbIndices.w()) = ele.second;
        }

        for (const auto& ele : f12RepulsionInts) {
            const Vector4i& orbIndices = ele.first;
            returnInts.f12Repulsion(orbIndices.x(), orbIndices.y(), orbIndices.z(), orbIndices.w()) = ele.second;
        }

        return returnInts;
    }
}
