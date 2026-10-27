#include "models.h"

#include <algorithm>
#include <map>

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

namespace {

std::filesystem::path executableDir() {
#ifdef __APPLE__
    char buf[4096];
    uint32_t size = sizeof(buf);
    if (_NSGetExecutablePath(buf, &size) == 0) {
        return std::filesystem::weakly_canonical(buf).parent_path();
    }
#endif
    return {};
}

void addDir(std::vector<std::filesystem::path>& dirs, const std::filesystem::path& path) {
    std::error_code error;
    if (!std::filesystem::is_directory(path, error)) return;
    const std::filesystem::path canonical = std::filesystem::weakly_canonical(path, error);
    if (error) return;
    if (std::find(dirs.begin(), dirs.end(), canonical) == dirs.end()) dirs.push_back(canonical);
}

} // namespace

std::vector<ModelFile> listModels() {
    std::vector<std::filesystem::path> dirs;
    const std::filesystem::path exe = executableDir();
    const std::filesystem::path relatives[] = {
        "input",
        "../../input",
        "../../../input",
        "../../../../assets/Models2",
        "assets/Models2",
    };
    for (const std::filesystem::path& relative : relatives) addDir(dirs, relative);
    if (!exe.empty()) {
        addDir(dirs, exe / "../../../../input");
        addDir(dirs, exe / "../../../../../../assets/Models2");
    }

    std::map<std::string, std::vector<ModelFile>> grouped;
    for (const std::filesystem::path& dir : dirs) {
        for (const std::filesystem::directory_entry& entry : std::filesystem::directory_iterator(dir)) {
            if (!entry.is_regular_file()) continue;
            if (entry.path().extension() != ".obj") continue;
            ModelFile model;
            model.label = entry.path().filename().string();
            model.path = entry.path();
            grouped[model.label].push_back(std::move(model));
        }
    }

    std::vector<ModelFile> models;
    for (auto& item : grouped) {
        if (item.second.size() == 1) {
            models.push_back(std::move(item.second[0]));
            continue;
        }
        for (ModelFile& model : item.second) {
            model.label += "  (" + model.path.parent_path().filename().string() + ")";
            models.push_back(std::move(model));
        }
    }
    std::sort(models.begin(), models.end(), [](const ModelFile& a, const ModelFile& b) {
        return a.label < b.label;
    });
    return models;
}

int modelIndex(const std::vector<ModelFile>& models, const std::filesystem::path& path) {
    std::error_code error;
    const std::filesystem::path canonical = std::filesystem::weakly_canonical(path, error);
    for (int i = 0; i < static_cast<int>(models.size()); ++i) {
        if (!error && std::filesystem::weakly_canonical(models[i].path) == canonical) return i;
        if (models[i].path.filename() == path.filename() && models[i].path == path) return i;
    }
    return -1;
}

int modelIndexByName(const std::vector<ModelFile>& models, const std::string& filename) {
    for (int i = 0; i < static_cast<int>(models.size()); ++i) {
        if (models[i].path.filename() == filename) return i;
    }
    return models.empty() ? -1 : 0;
}
