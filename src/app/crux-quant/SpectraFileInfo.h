#pragma once

#include <string>

namespace CruxQuant{
    class SpectraFileInfo{
    public:
        const std::string fullFilePathWithExtension;
        const std::string filenameWithoutExtension;
        std::string condition;
        const int biologicalReplicate;
        const int fraction;
        const int technicalReplicate;

        SpectraFileInfo(const std::string& fullFilePathWithExtension, const std::string& condition, int biorep, int techrep, int fraction);

        // Files are considered the same if the absolute file path is the same
        bool operator==(const SpectraFileInfo& other) const;

        bool operator!=(const SpectraFileInfo& other) const;

        // Use the full file path hashcode as the hash code
        size_t GetHashCode() const;

    private:
        // Extract the filename without extension from a given file path
        static std::string GetFilenameWithoutExtension(const std::string& filePath);
    };
} // namespace CruxQuant
