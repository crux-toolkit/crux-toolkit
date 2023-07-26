#include "SpectraFileInfo.h"
#include <algorithm>

namespace CruxQuant{
    SpectraFileInfo::SpectraFileInfo(const std::string& fullFilePathWithExtension, const std::string& condition, int biorep, int techrep, int fraction)
        : fullFilePathWithExtension(fullFilePathWithExtension),
          filenameWithoutExtension(GetFilenameWithoutExtension(fullFilePathWithExtension)),
          condition(condition),
          biologicalReplicate(biorep),
          fraction(fraction),
          technicalReplicate(techrep){}

    // Files are considered the same if the absolute file path is the same
    bool SpectraFileInfo::operator==(const SpectraFileInfo& other) const{
        return fullFilePathWithExtension == other.fullFilePathWithExtension;
    }

    bool SpectraFileInfo::operator!=(const SpectraFileInfo& other) const{
        return !(*this == other);
    }

    // Use the full file path hashcode as the hash code
    size_t SpectraFileInfo::GetHashCode() const{
        return std::hash<std::string>()(fullFilePathWithExtension);
    }

    // Extract the filename without extension from a given file path
    std::string SpectraFileInfo::GetFilenameWithoutExtension(const std::string& filePath){
        size_t lastSlash = filePath.find_last_of("/\\");
        size_t lastDot = filePath.find_last_of(".");
        if (lastDot != std::string::npos && (lastSlash == std::string::npos || lastDot > lastSlash))
        {
            return filePath.substr(lastSlash + 1, lastDot - lastSlash - 1);
        }
        return filePath.substr(lastSlash + 1);
    }
} // namespace CruxQuant
