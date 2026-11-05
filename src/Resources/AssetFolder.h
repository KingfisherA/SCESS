#pragma once
#include <boost/filesystem.hpp>


//
//
namespace SCESS::Res {
    //
    //
    class AssetsFolder {
    public:
        static constexpr const char *SETTINGS_FILE = "settings.json";

        static AssetsFolder &Get();

        // Finds the Assets folder, in order:
        //  1. "assetsPath" in a settings.json in the working directory or next to the executable
        //  2. Assets inside the macOS app bundle (Contents/Resources/Assets)
        //  3. Assets next to the executable or in the working directory
        //  4. The source tree's Assets folder (SCESS_SOURCE_ASSETS_DIR, for builds run in place)
        bool Init(const char *settingsFile = SETTINGS_FILE);
        [[nodiscard]] const boost::filesystem::path &GetAssetsPath() const { return mAssetsPath; }
        [[nodiscard]] std::string GetAssetsPathStr() const { return mAssetsPath.string(); }
        [[nodiscard]] const char *GetAssetsPathCStr() const { return (const char*)mAssetsPath.c_str(); }

        [[nodiscard]] boost::filesystem::path GetResourcePath(const std::string &resource) const {
            return boost::filesystem::path(GetAssetsPath()).append(resource);
        }

        // Directory containing the running executable (empty if it cannot be determined)
        static boost::filesystem::path GetExecutableDirectory();
        // Per-user directory for files SCESS writes (e.g. the ImGui layout); created if missing.
        // macOS: ~/Library/Application Support/SCESS, Windows: %APPDATA%\SCESS,
        // Linux: $XDG_CONFIG_HOME/SCESS or ~/.config/SCESS. Falls back to the working directory.
        static boost::filesystem::path GetUserDataDirectory();

    private:
        boost::filesystem::path mAssetsPath = "./Assets";
    };
    //
    //
}
