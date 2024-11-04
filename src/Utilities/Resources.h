#pragma once

#include <boost/filesystem.hpp>

namespace SCESS {
    class Resources {
    public:
        //Relative Assets: "./Assets"
#define ASSETS_PATH "../../../../Assets"
        [[nodiscard]] static const char *GetAssetsPath() { return ASSETS_PATH; }

        [[nodiscard]] static boost::filesystem::path GetResourcePath(const std::string &resource) {
            return boost::filesystem::path(GetAssetsPath()).append(resource);
        };

        static const char *GetFontPath() {
            return ASSETS_PATH "/Fonts/Ubuntu-Light.ttf";
        }
    };
}
