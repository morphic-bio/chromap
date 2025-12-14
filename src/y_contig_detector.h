#ifndef Y_CONTIG_DETECTOR_H_
#define Y_CONTIG_DETECTOR_H_

#include <algorithm>
#include <cctype>
#include <iostream>
#include <string>
#include <unordered_set>

#include "sequence_batch.h"

namespace chromap {

// Build a set of reference sequence IDs (RIDs) that represent Y chromosome contigs.
// 
// Detection rules (intentional design choices):
// - Case-insensitive matching
// - Strips "chr" or "chr_" prefix if present
// - EXACT match for "y" after stripping (no partial/substring matching)
// 
// Matches: "Y", "chrY", "CHR_Y", "chr_y", "y"
// Does NOT match: "chrY_random", "Y_alt", "chrY_hap1" (exact match required)
//
// This is intentional: we want to filter reads mapping to the primary Y chromosome,
// not to decoy/random/alt contigs which may have different filtering requirements.
inline std::unordered_set<uint32_t> BuildYContigRidMask(
    uint32_t num_reference_sequences,
    const SequenceBatch &reference) {
  
  std::unordered_set<uint32_t> y_rids;
  
  for (uint32_t rid = 0; rid < num_reference_sequences; ++rid) {
    std::string name(reference.GetSequenceNameAt(rid));
    
    // Convert to lowercase
    std::transform(name.begin(), name.end(), name.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    
    // Strip "chr" or "chr_" prefix if present
    if (name.length() > 3 && name.substr(0, 3) == "chr") {
      name = name.substr(3);
      if (!name.empty() && name[0] == '_') {
        name = name.substr(1);
      }
    }
    
    // Exact match for "y" only
    if (name == "y") {
      y_rids.insert(rid);
      std::cerr << "  Y contig detected: " << reference.GetSequenceNameAt(rid)
                << " (rid=" << rid << ")\n";
    }
  }
  
  if (y_rids.empty()) {
    std::cerr << "WARNING: No Y chromosome contigs found in reference. "
              << "Y-filtering will have no effect.\n";
  } else {
    std::cerr << "Found " << y_rids.size() << " Y chromosome contig(s).\n";
  }
  
  return y_rids;
}

}  // namespace chromap

#endif  // Y_CONTIG_DETECTOR_H_

