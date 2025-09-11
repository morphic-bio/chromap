#ifndef OVERFLOW_READER_H_
#define OVERFLOW_READER_H_

#include <string>
#include <cstdio>
#include <cstdint>

// Simple reader for overflow files
// Reads length-prefixed records sequentially
class OverflowReader {
public:
    explicit OverflowReader(const std::string& path);
    ~OverflowReader();

    // Read next record header and payload
    // Returns false on EOF, true on success
    bool ReadNext(uint32_t& out_rid, std::string& out_payload);

    // Check if reader is valid (file opened successfully)
    bool IsValid() const { return file_ != nullptr; }

    // Get current file path
    const std::string& GetPath() const { return path_; }

private:
    std::string path_;
    FILE* file_;
};

#endif  // OVERFLOW_READER_H_
