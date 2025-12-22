#ifndef OVERFLOW_WRITER_H_
#define OVERFLOW_WRITER_H_

#include <atomic>
#include <string>
#include <unordered_map>
#include <vector>
#include <cstdio>
#include <cstdint>

// Thread-safe overflow writer for mapping records
// Each thread writes to its own files, avoiding contention
class OverflowWriter {
public:
    explicit OverflowWriter(const std::string& base_dir = "", const std::string& prefix = "chromap_of");
    ~OverflowWriter();

    // Thread-safe: write a record to thread-local file for the given rid
    template<typename MappingRecord>
    void Write(uint32_t rid, const MappingRecord& rec) {
        FILE* fp = GetFileForRid(rid);
        if (!fp) {
            return;
        }

        const uint32_t byte_len = static_cast<uint32_t>(rec.SerializedSize());

        fwrite(&rid, sizeof(uint32_t), 1, fp);
        fwrite(&byte_len, sizeof(uint32_t), 1, fp);
        rec.WriteToFile(fp);
    }

    // Finalize all thread-local files and return file paths
    std::vector<std::string> Close();

private:
    std::string base_dir_;
    std::string prefix_;
    static std::atomic<uint32_t> global_counter_;
    
    // Thread-local storage for file handles
    static thread_local std::unordered_map<uint32_t, FILE*> tls_files_;
    static thread_local std::unordered_map<uint32_t, std::string> tls_file_paths_;
    static thread_local bool tls_initialized_;
    
    // Get or create thread-local file for given rid
    FILE* GetFileForRid(uint32_t rid);
    
    // Generate unique filename
    std::string GenerateFilename(uint32_t rid);
};


#endif  // OVERFLOW_WRITER_H_
