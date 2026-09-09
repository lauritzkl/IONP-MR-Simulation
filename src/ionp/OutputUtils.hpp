/**
 * @file OutputUtils.hpp
 *
 * Output numbering based on directory contents. This is a count-based naming
 * convention, not a reservation mechanism: gaps, unrelated files or concurrent
 * runs can reuse names. Use separate empty output directories for each run.
 */

/*****************************************************************************
*Utility helpers for output file naming based on directory contents
******************************************************************************/

#ifndef OUTPUT_UTILS_HPP_
#define OUTPUT_UTILS_HPP_

#include <cstddef>
#include <filesystem>

namespace output {

inline std::size_t count_regular_files(const std::filesystem::path& dir) {
    if(!std::filesystem::exists(dir)) {
        return 0;
    }

    std::size_t count = 0;
    for(const auto& entry : std::filesystem::directory_iterator(dir)) {
        if(entry.is_regular_file()) {
            ++count;
        }
    }
    return count;
}

// Count regular files, not the largest numeric filename prefix. Non-atomic.
inline std::size_t next_file_index(const std::filesystem::path& dir) {
    return count_regular_files(dir) + 1;
}

} // namespace output

#endif
