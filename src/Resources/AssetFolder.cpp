#include "AssetFolder.h"
#include<nlohmann/json.hpp>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<vector>

#if defined(__APPLE__)
#include <mach-o/dyld.h>
#elif defined(_WIN32)
#include <windows.h>
#else
#include <unistd.h>
#endif


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

        const filesystem::path exeDir = GetExecutableDirectory();
        std::vector<filesystem::path> candidates;

        // 1. An explicit assetsPath from settings.json (relative paths are relative to the working directory)
        for (const filesystem::path &dir: {filesystem::path("."), exeDir}) {
            const filesystem::path settingsPath = dir / settingsFile;
            system::error_code ec;
            if (dir.empty() || !filesystem::exists(settingsPath, ec))
                continue;

            std::ifstream jsonFile(settingsPath.string());
            try {
                json settingsJSON = json::parse(jsonFile);
                std::string tempPath;
                settingsJSON["assetsPath"].get_to(tempPath);
                candidates.emplace_back(tempPath);
            } catch (const std::exception &e) {
                std::cout << "Could not parse settings file " << settingsPath << ".\n" << e.what() << '\n';
            }
            break;
        }

        // 2-3. Assets inside the app bundle, next to the executable or in the working directory
        if (!exeDir.empty()) {
            candidates.push_back(exeDir / ".." / "Resources" / "Assets");
            candidates.push_back(exeDir / "Assets");
        }
        candidates.emplace_back("./Assets");

        // 4. The source tree's Assets folder
#ifdef SCESS_SOURCE_ASSETS_DIR
        candidates.emplace_back(SCESS_SOURCE_ASSETS_DIR);
#endif

        for (const filesystem::path &candidate: candidates) {
            system::error_code ec;
            if (filesystem::is_directory(candidate, ec)) {
                const filesystem::path canonical = filesystem::canonical(candidate, ec);
                mAssetsPath = ec ? candidate : canonical;
                std::cout << "Found assets directory at " << mAssetsPath << "\n";
                return true;
            }
        }

        std::cout << "Error! Required assets directory not found. Searched:\n";
        for (const filesystem::path &candidate: candidates)
            std::cout << "    " << candidate << "\n";
        return false;
    }
    //
    filesystem::path AssetsFolder::GetExecutableDirectory() {
        filesystem::path exePath;
#if defined(__APPLE__)
        uint32_t size = 0;
        _NSGetExecutablePath(nullptr, &size);
        std::string buffer(size, '\0');
        if (_NSGetExecutablePath(buffer.data(), &size) != 0)
            return {};
        exePath = buffer.c_str();
#elif defined(_WIN32)
        std::wstring buffer(MAX_PATH, L'\0');
        DWORD length = 0;
        while ((length = GetModuleFileNameW(nullptr, buffer.data(), (DWORD) buffer.size())) == buffer.size())
            buffer.resize(buffer.size() * 2);
        if (length == 0)
            return {};
        buffer.resize(length);
        exePath = buffer;
#else
        std::string buffer(4096, '\0');
        const ssize_t length = readlink("/proc/self/exe", buffer.data(), buffer.size());
        if (length <= 0)
            return {};
        buffer.resize((size_t) length);
        exePath = buffer;
#endif
        system::error_code ec;
        const filesystem::path canonical = filesystem::canonical(exePath, ec);
        return (ec ? exePath : canonical).parent_path();
    }
    //
    filesystem::path AssetsFolder::GetUserDataDirectory() {
        const auto env = [](const char *name) -> filesystem::path {
            const char *value = std::getenv(name);
            return value && *value ? filesystem::path(value) : filesystem::path();
        };

        filesystem::path dir;
#if defined(__APPLE__)
        if (const filesystem::path home = env("HOME"); !home.empty())
            dir = home / "Library" / "Application Support" / "SCESS";
#elif defined(_WIN32)
        if (const filesystem::path appData = env("APPDATA"); !appData.empty())
            dir = appData / "SCESS";
#else
        if (const filesystem::path config = env("XDG_CONFIG_HOME"); !config.empty())
            dir = config / "SCESS";
        else if (const filesystem::path home = env("HOME"); !home.empty())
            dir = home / ".config" / "SCESS";
#endif

        system::error_code ec;
        if (dir.empty() || (filesystem::create_directories(dir, ec), ec))
            return filesystem::current_path(ec);
        return dir;
    }
    //
    //
}
