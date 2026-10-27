#pragma once

#include <filesystem>
#include <string>
#include <vector>

struct ModelFile {
    std::string label;
    std::filesystem::path path;
};

std::vector<ModelFile> listModels();
int modelIndex(const std::vector<ModelFile>& models, const std::filesystem::path& path);
int modelIndexByName(const std::vector<ModelFile>& models, const std::string& filename);
