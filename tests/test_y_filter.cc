// Unit tests for Y-chromosome filtering functionality
// Tests Y contig detection and basic routing logic

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_set>
#include <unistd.h>
#include <vector>

#include "../src/y_contig_detector.h"
#include "../src/sequence_batch.h"

using namespace chromap;

// Helper: Create a temporary FASTA file with given contig names
std::string CreateTempFasta(const std::vector<std::string> &contig_names) {
  char tmpfile[] = "/tmp/test_y_filter_XXXXXX";
  int fd = mkstemp(tmpfile);
  if (fd == -1) {
    std::cerr << "Failed to create temp file\n";
    exit(1);
  }
  close(fd);
  
  std::ofstream out(tmpfile);
  for (size_t i = 0; i < contig_names.size(); ++i) {
    out << ">" << contig_names[i] << "\n";
    out << "ACGTACGTACGTACGTACGTACGTACGTACGT\n";  // Dummy sequence
  }
  out.close();
  return std::string(tmpfile);
}

// Test 1: Y contig detection with various naming conventions
void TestYContigDetection() {
  std::cout << "Test 1: Y contig detection...\n";
  
  // Test cases: (contig_name, should_match)
  std::vector<std::pair<std::string, bool>> test_cases = {
    {"chrY", true},      // Standard UCSC format
    {"Y", true},         // No prefix
    {"chr_y", true},     // Lowercase with underscore
    {"CHRY", true},      // Uppercase
    {"chrY_random", false},  // Should NOT match (exact match required)
    {"Y_alt", false},    // Should NOT match
    {"chrY_hap1", false}, // Should NOT match
    {"chr1", false},     // Different chromosome
    {"chrX", false},     // Different chromosome
    {"chrM", false},     // Mitochondria
  };
  
  for (const auto &test_case : test_cases) {
    std::vector<std::string> contigs = {test_case.first};
    std::string fasta_file = CreateTempFasta(contigs);
    
    SequenceBatch reference;
    reference.InitializeLoading(fasta_file);
    reference.LoadAllSequences();
    
    std::unordered_set<uint32_t> y_rids = BuildYContigRidMask(
        reference.GetNumSequences(), reference);
    
    bool found = (y_rids.count(0) > 0);
    if (found != test_case.second) {
      std::cerr << "FAIL: " << test_case.first 
                << " expected " << (test_case.second ? "match" : "no match")
                << " but got " << (found ? "match" : "no match") << "\n";
      unlink(fasta_file.c_str());
      exit(1);
    }
    
    reference.FinalizeLoading();
    unlink(fasta_file.c_str());
  }
  
  std::cout << "  PASS: All Y contig detection cases passed\n";
}

// Test 2: Multiple contigs including Y
void TestMultipleContigs() {
  std::cout << "Test 2: Multiple contigs with Y...\n";
  
  std::vector<std::string> contigs = {"chr1", "chr2", "chrY", "chr3", "chrX"};
  std::string fasta_file = CreateTempFasta(contigs);
  
  SequenceBatch reference;
  reference.InitializeLoading(fasta_file);
  reference.LoadAllSequences();
  
  std::unordered_set<uint32_t> y_rids = BuildYContigRidMask(
      reference.GetNumSequences(), reference);
  
  // Should find exactly one Y contig at index 2
  if (y_rids.size() != 1 || y_rids.count(2) == 0) {
    std::cerr << "FAIL: Expected Y contig at index 2, got " << y_rids.size()
              << " Y contig(s)\n";
    unlink(fasta_file.c_str());
    exit(1);
  }
  
  reference.FinalizeLoading();
  unlink(fasta_file.c_str());
  
  std::cout << "  PASS: Correctly identified Y contig among multiple contigs\n";
}

// Test 3: No Y contigs
void TestNoYContigs() {
  std::cout << "Test 3: Reference with no Y contigs...\n";
  
  std::vector<std::string> contigs = {"chr1", "chr2", "chr3", "chrX"};
  std::string fasta_file = CreateTempFasta(contigs);
  
  SequenceBatch reference;
  reference.InitializeLoading(fasta_file);
  reference.LoadAllSequences();
  
  std::unordered_set<uint32_t> y_rids = BuildYContigRidMask(
      reference.GetNumSequences(), reference);
  
  if (!y_rids.empty()) {
    std::cerr << "FAIL: Expected no Y contigs, found " << y_rids.size() << "\n";
    unlink(fasta_file.c_str());
    exit(1);
  }
  
  reference.FinalizeLoading();
  unlink(fasta_file.c_str());
  
  std::cout << "  PASS: Correctly identified no Y contigs\n";
}

// Test 4: Case-insensitive matching
void TestCaseInsensitive() {
  std::cout << "Test 4: Case-insensitive matching...\n";
  
  std::vector<std::pair<std::string, bool>> cases = {
    {"Y", true},
    {"y", true},
    {"ChrY", true},
    {"CHRY", true},
    {"chrY", true},
  };
  
  for (const auto &test_case : cases) {
    std::vector<std::string> contigs = {test_case.first};
    std::string fasta_file = CreateTempFasta(contigs);
    
    SequenceBatch reference;
    reference.InitializeLoading(fasta_file);
    reference.LoadAllSequences();
    
    std::unordered_set<uint32_t> y_rids = BuildYContigRidMask(
        reference.GetNumSequences(), reference);
    
    bool found = (y_rids.count(0) > 0);
    if (!found) {
      std::cerr << "FAIL: Case-insensitive match failed for " << test_case.first << "\n";
      unlink(fasta_file.c_str());
      exit(1);
    }
    
    reference.FinalizeLoading();
    unlink(fasta_file.c_str());
  }
  
  std::cout << "  PASS: Case-insensitive matching works correctly\n";
}

int main() {
  std::cout << "Running Y-filter unit tests...\n\n";
  
  TestYContigDetection();
  TestMultipleContigs();
  TestNoYContigs();
  TestCaseInsensitive();
  
  std::cout << "\nAll unit tests PASSED!\n";
  return 0;
}

