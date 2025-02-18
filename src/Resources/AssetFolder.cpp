#include "AssetFolder.h"
#include<nlohmann/json.hpp>
#include<iostream>
#include<fstream>


//
//
namespace SCESS::Res {
    //
    //
    using namespace boost;
    //
    //
    AssetsFolder &AssetsFolder::Get() {
        static AssetsFolder instance;
        return instance;
    }
    //
    bool AssetsFolder::Init(const char *settingsFile) {
        using namespace nlohmann;

        if (filesystem::exists(settingsFile)) {
            std::ifstream jsonFile(settingsFile);
            try {
                json settingsJSON = json::parse(jsonFile);
                std::string tempPath;
                settingsJSON["assetsPath"].get_to(tempPath);
                mAssetsPath = tempPath;
            } catch (const std::exception &e) {
                std::cout << "Could not parse settings JSON. \n" << e.what() << '\n';
            }
        } else
            std::cout << "Could not find settings file! Defaulting to ./Assets directory\n";


        if (!filesystem::exists(mAssetsPath)) {
            std::cout << "Error! Required assets directory not found at \"" << mAssetsPath << "\"\n";
            return false;
        }

        std::cout << "Found assets directory at " << mAssetsPath << "\n";
        return true;
    }
    //
    //
}
