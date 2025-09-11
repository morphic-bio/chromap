#include "overflow_reader.h"

#include <cassert>

OverflowReader::OverflowReader(const std::string& path) 
    : path_(path), file_(nullptr) {
    file_ = fopen(path.c_str(), "rb");
}

OverflowReader::~OverflowReader() {
    if (file_) {
        fclose(file_);
        file_ = nullptr;
    }
}

bool OverflowReader::ReadNext(uint32_t& out_rid, std::string& out_payload) {
    if (!file_) {
        return false;
    }
    
    // Read header: rid (4 bytes) + byte_len (4 bytes)
    uint32_t rid, byte_len;
    
    if (fread(&rid, sizeof(uint32_t), 1, file_) != 1) {
        // EOF or error
        return false;
    }
    
    if (fread(&byte_len, sizeof(uint32_t), 1, file_) != 1) {
        // Incomplete header - this is an error
        return false;
    }
    
    // Read payload
    out_payload.resize(byte_len);
    if (byte_len > 0) {
        if (fread(&out_payload[0], 1, byte_len, file_) != byte_len) {
            // Incomplete payload - this is an error
            return false;
        }
    }
    
    out_rid = rid;
    return true;
}
