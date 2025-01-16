#pragma once

#include <boost/filesystem.hpp>

namespace SCESS::Res {
    class AssetsFolder {
    public:
        static constexpr const char *SETTINGS_FILE = "./settings.json";

        static AssetsFolder &Get();

        bool Init(const char *settingsFile = SETTINGS_FILE);
        [[nodiscard]] const boost::filesystem::path &GetAssetsPath() const { return mAssetsPath; }
        [[nodiscard]] const std::string &GetAssetsPathStr() const { return mAssetsPath.string(); }
        [[nodiscard]] const char *GetAssetsPathCStr() const { return (const char*)mAssetsPath.c_str(); }

        [[nodiscard]] boost::filesystem::path GetResourcePath(const std::string &resource) const {
            return boost::filesystem::path(GetAssetsPath()).append(resource);
        }

    private:
        boost::filesystem::path mAssetsPath = "./Assets";
    };
}
