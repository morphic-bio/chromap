#include "mapping_writer.h"

#include <queue>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <unistd.h>

namespace chromap {

#ifndef LEGACY_OVERFLOW
// Static member definitions for thread-local and shared overflow handling
template <typename MappingRecord>
thread_local std::unique_ptr<OverflowWriter> MappingWriter<MappingRecord>::tls_overflow_writer_;

template <typename MappingRecord>
std::vector<std::string> MappingWriter<MappingRecord>::shared_overflow_file_paths_;

template <typename MappingRecord>
std::mutex MappingWriter<MappingRecord>::overflow_paths_mutex_;
#endif

// Specialization for BED format.
template <>
void MappingWriter<MappingWithBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference) {}

template <>
void MappingWriter<MappingWithBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const MappingWithBarcode &mapping) {
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BED) {
    std::string strand = mapping.IsPositiveStrand() ? "+" : "-";
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t mapping_end_position = mapping.GetEndPosition();
    this->AppendMappingOutput(std::string(reference_sequence_name) + "\t" +
                              std::to_string(mapping.GetStartPosition()) +
                              "\t" + std::to_string(mapping_end_position) +
                              "\t" +
                              barcode_translator_.Translate(
                                  mapping.cell_barcode_, cell_barcode_length_) +
                              "\t" + std::to_string(mapping.num_dups_) + "\n");
  } else {
    std::string strand = mapping.IsPositiveStrand() ? "+" : "-";
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t mapping_end_position = mapping.GetEndPosition();
    this->AppendMappingOutput(std::string(reference_sequence_name) + "\t" +
                              std::to_string(mapping.GetStartPosition()) +
                              "\t" + std::to_string(mapping_end_position) +
                              "\tN\t" + std::to_string(mapping.mapq_) + "\t" +
                              strand + "\n");
  }
}

template <>
void MappingWriter<MappingWithoutBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference) {}

template <>
void MappingWriter<MappingWithoutBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const MappingWithoutBarcode &mapping) {
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BED) {
    std::string strand = mapping.IsPositiveStrand() ? "+" : "-";
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t mapping_end_position = mapping.GetEndPosition();
    this->AppendMappingOutput(std::string(reference_sequence_name) + "\t" +
                              std::to_string(mapping.GetStartPosition()) +
                              "\t" + std::to_string(mapping_end_position) +
                              "\tN\t" + std::to_string(mapping.mapq_) + "\t" +
                              strand + "\t" + std::to_string(mapping.num_dups_) + "\n");
  } else {
    std::string strand = mapping.IsPositiveStrand() ? "+" : "-";
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t mapping_end_position = mapping.GetEndPosition();
    this->AppendMappingOutput(std::string(reference_sequence_name) + "\t" +
                              std::to_string(mapping.GetStartPosition()) +
                              "\t" + std::to_string(mapping_end_position) +
                              "\tN\t" + std::to_string(mapping.mapq_) + "\t" +
                              strand + "\t" + std::to_string(mapping.num_dups_) + "\n");
  }
}

// Specialization for BEDPE format.
template <>
void MappingWriter<PairedEndMappingWithoutBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference) {}

template <>
void MappingWriter<PairedEndMappingWithoutBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const PairedEndMappingWithoutBarcode &mapping) {
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BED) {
    std::string strand = mapping.IsPositiveStrand() ? "+" : "-";
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t mapping_end_position = mapping.GetEndPosition();
    this->AppendMappingOutput(std::string(reference_sequence_name) + "\t" +
                              std::to_string(mapping.GetStartPosition()) +
                              "\t" + std::to_string(mapping_end_position) +
                              "\tN\t" + std::to_string(mapping.mapq_) + "\t" +
                              strand + "\t" + std::to_string(mapping.num_dups_) + "\n");
  } else {
    bool positive_strand = mapping.IsPositiveStrand();
    uint32_t positive_read_end =
        mapping.fragment_start_position_ + mapping.positive_alignment_length_;
    uint32_t negative_read_end =
        mapping.fragment_start_position_ + mapping.fragment_length_;
    uint32_t negative_read_start =
        negative_read_end - mapping.negative_alignment_length_;
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    if (positive_strand) {
      this->AppendMappingOutput(
          std::string(reference_sequence_name) + "\t" +
          std::to_string(mapping.fragment_start_position_) + "\t" +
          std::to_string(positive_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t+\n" +
          std::string(reference_sequence_name) + "\t" +
          std::to_string(negative_read_start) + "\t" +
          std::to_string(negative_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t-\t" + 
          std::to_string(mapping.num_dups_) + "\n");
    } else {
      this->AppendMappingOutput(
          std::string(reference_sequence_name) + "\t" +
          std::to_string(negative_read_start) + "\t" +
          std::to_string(negative_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t-\n" +
          std::string(reference_sequence_name) + "\t" +
          std::to_string(mapping.fragment_start_position_) + "\t" +
          std::to_string(positive_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t+\t" +
          std::to_string(mapping.num_dups_) + "\n");
    }
  }
}

template <>
void MappingWriter<PairedEndMappingWithBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference) {}

template <>
void MappingWriter<PairedEndMappingWithBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const PairedEndMappingWithBarcode &mapping) {
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BED) {
    std::string strand = mapping.IsPositiveStrand() ? "+" : "-";
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t mapping_end_position = mapping.GetEndPosition();
    this->AppendMappingOutput(std::string(reference_sequence_name) + "\t" +
                              std::to_string(mapping.GetStartPosition()) +
                              "\t" + std::to_string(mapping_end_position) +
                              "\t" +
                              barcode_translator_.Translate(
                                  mapping.cell_barcode_, cell_barcode_length_) +
                              "\t" + std::to_string(mapping.num_dups_) + "\n");
  } else {
    bool positive_strand = mapping.IsPositiveStrand();
    uint32_t positive_read_end =
        mapping.fragment_start_position_ + mapping.positive_alignment_length_;
    uint32_t negative_read_end =
        mapping.fragment_start_position_ + mapping.fragment_length_;
    uint32_t negative_read_start =
        negative_read_end - mapping.negative_alignment_length_;
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    if (positive_strand) {
      this->AppendMappingOutput(
          std::string(reference_sequence_name) + "\t" +
          std::to_string(mapping.fragment_start_position_) + "\t" +
          std::to_string(positive_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t+\n" +
          std::string(reference_sequence_name) + "\t" +
          std::to_string(negative_read_start) + "\t" +
          std::to_string(negative_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t-\n");
    } else {
      this->AppendMappingOutput(
          std::string(reference_sequence_name) + "\t" +
          std::to_string(negative_read_start) + "\t" +
          std::to_string(negative_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t-\n" +
          std::string(reference_sequence_name) + "\t" +
          std::to_string(mapping.fragment_start_position_) + "\t" +
          std::to_string(positive_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t+\n");
    }
  }
}

// Specialization for PAF format.
template <>
void MappingWriter<PAFMapping>::OutputHeader(uint32_t num_reference_sequences,
                                             const SequenceBatch &reference) {}

template <>
void MappingWriter<PAFMapping>::AppendMapping(uint32_t rid,
                                              const SequenceBatch &reference,
                                              const PAFMapping &mapping) {
  const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
  uint32_t reference_sequence_length = reference.GetSequenceLengthAt(rid);
  std::string strand = mapping.IsPositiveStrand() ? "+" : "-";
  uint32_t mapping_end_position =
      mapping.fragment_start_position_ + mapping.fragment_length_;
  this->AppendMappingOutput(
      mapping.read_name_ + "\t" + std::to_string(mapping.read_length_) + "\t" +
      std::to_string(0) + "\t" + std::to_string(mapping.read_length_) + "\t" +
      strand + "\t" + std::string(reference_sequence_name) + "\t" +
      std::to_string(reference_sequence_length) + "\t" +
      std::to_string(mapping.fragment_start_position_) + "\t" +
      std::to_string(mapping_end_position) + "\t" +
      std::to_string(mapping.read_length_) + "\t" +
      std::to_string(mapping.fragment_length_) + "\t" +
      std::to_string(mapping.mapq_) + "\n");
}

template <>
void MappingWriter<PAFMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<PAFMapping> > &mappings) {
  FILE *temp_mapping_output_file =
      fopen(temp_mapping_output_file_path.c_str(), "wb");
  assert(temp_mapping_output_file != NULL);
  for (size_t ri = 0; ri < num_reference_sequences; ++ri) {
    // make sure mappings[ri] exists even if its size is 0
    size_t num_mappings = mappings[ri].size();
    fwrite(&num_mappings, sizeof(size_t), 1, temp_mapping_output_file);
    if (mappings[ri].size() > 0) {
      for (size_t mi = 0; mi < num_mappings; ++mi) {
        mappings[ri][mi].WriteToFile(temp_mapping_output_file);
      }
      // fwrite(mappings[ri].data(), sizeof(MappingRecord), mappings[ri].size(),
      // temp_mapping_output_file);
    }
  }
  fclose(temp_mapping_output_file);
}

// Specialization for PairedPAF format.
template <>
void MappingWriter<PairedPAFMapping>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference) {}

template <>
void MappingWriter<PairedPAFMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<PairedPAFMapping> > &mappings) {
  FILE *temp_mapping_output_file =
      fopen(temp_mapping_output_file_path.c_str(), "wb");
  assert(temp_mapping_output_file != NULL);
  for (size_t ri = 0; ri < num_reference_sequences; ++ri) {
    // make sure mappings[ri] exists even if its size is 0
    size_t num_mappings = mappings[ri].size();
    fwrite(&num_mappings, sizeof(size_t), 1, temp_mapping_output_file);
    if (mappings[ri].size() > 0) {
      for (size_t mi = 0; mi < num_mappings; ++mi) {
        mappings[ri][mi].WriteToFile(temp_mapping_output_file);
      }
      // fwrite(mappings[ri].data(), sizeof(MappingRecord), mappings[ri].size(),
      // temp_mapping_output_file);
    }
  }
  fclose(temp_mapping_output_file);
}

template <>
void MappingWriter<PairedPAFMapping>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const PairedPAFMapping &mapping) {
  bool positive_strand = mapping.IsPositiveStrand();
  uint32_t positive_read_end =
      mapping.fragment_start_position_ + mapping.positive_alignment_length_;
  uint32_t negative_read_end =
      mapping.fragment_start_position_ + mapping.fragment_length_;
  uint32_t negative_read_start =
      negative_read_end - mapping.negative_alignment_length_;
  const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
  uint32_t reference_sequence_length = reference.GetSequenceLengthAt(rid);
  if (positive_strand) {
    this->AppendMappingOutput(
        mapping.read1_name_ + "\t" + std::to_string(mapping.read1_length_) +
        "\t" + std::to_string(0) + "\t" +
        std::to_string(mapping.read1_length_) + "\t" + "+" + "\t" +
        std::string(reference_sequence_name) + "\t" +
        std::to_string(reference_sequence_length) + "\t" +
        std::to_string(mapping.fragment_start_position_) + "\t" +
        std::to_string(positive_read_end) + "\t" +
        std::to_string(mapping.read1_length_) + "\t" +
        std::to_string(mapping.positive_alignment_length_) + "\t" +
        std::to_string(mapping.mapq1_) + "\n");
    this->AppendMappingOutput(
        mapping.read2_name_ + "\t" + std::to_string(mapping.read2_length_) +
        "\t" + std::to_string(0) + "\t" +
        std::to_string(mapping.read2_length_) + "\t" + "-" + "\t" +
        std::string(reference_sequence_name) + "\t" +
        std::to_string(reference_sequence_length) + "\t" +
        std::to_string(negative_read_start) + "\t" +
        std::to_string(negative_read_end) + "\t" +
        std::to_string(mapping.read2_length_) + "\t" +
        std::to_string(mapping.negative_alignment_length_) + "\t" +
        std::to_string(mapping.mapq2_) + "\n");
  } else {
    this->AppendMappingOutput(
        mapping.read1_name_ + "\t" + std::to_string(mapping.read1_length_) +
        "\t" + std::to_string(0) + "\t" +
        std::to_string(mapping.read1_length_) + "\t" + "-" + "\t" +
        std::string(reference_sequence_name) + "\t" +
        std::to_string(reference_sequence_length) + "\t" +
        std::to_string(negative_read_start) + "\t" +
        std::to_string(negative_read_end) + "\t" +
        std::to_string(mapping.read1_length_) + "\t" +
        std::to_string(mapping.negative_alignment_length_) + "\t" +
        std::to_string(mapping.mapq1_) + "\n");
    this->AppendMappingOutput(
        mapping.read2_name_ + "\t" + std::to_string(mapping.read2_length_) +
        "\t" + std::to_string(0) + "\t" +
        std::to_string(mapping.read2_length_) + "\t" + "+" + "\t" +
        std::string(reference_sequence_name) + "\t" +
        std::to_string(reference_sequence_length) + "\t" +
        std::to_string(mapping.fragment_start_position_) + "\t" +
        std::to_string(positive_read_end) + "\t" +
        std::to_string(mapping.read2_length_) + "\t" +
        std::to_string(mapping.positive_alignment_length_) + "\t" +
        std::to_string(mapping.mapq2_) + "\n");
  }
}

// Specialization for SAM format.
template <>
void MappingWriter<SAMMapping>::OutputHeader(uint32_t num_reference_sequences,
                                             const SequenceBatch &reference) {
  for (uint32_t rid = 0; rid < num_reference_sequences; ++rid) {
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t reference_sequence_length = reference.GetSequenceLengthAt(rid);
    std::string header_line = "@SQ\tSN:" + std::string(reference_sequence_name) +
                              "\tLN:" + std::to_string(reference_sequence_length) + "\n";
    
    // Write to primary output
    this->AppendMappingOutput(header_line);
    
    // Mirror to secondary streams (must be open before this call)
    if (noY_output_file_) {
      fwrite(header_line.data(), 1, header_line.size(), noY_output_file_);
    }
    if (Y_output_file_) {
      fwrite(header_line.data(), 1, header_line.size(), Y_output_file_);
    }
  }
}

template <>
void MappingWriter<SAMMapping>::AppendMapping(uint32_t rid,
                                              const SequenceBatch &reference,
                                              const SAMMapping &mapping) {
  const char *reference_sequence_name =
      (mapping.flag_ & BAM_FUNMAP) > 0 ? "*" : reference.GetSequenceNameAt(rid);
  const char *mate_ref_sequence_name =
      mapping.mrid_ < 0 ? "*" : 
      ((uint32_t)mapping.mrid_ == rid ? "=" : reference.GetSequenceNameAt(mapping.mrid_));
  const uint32_t mapping_start_position = mapping.GetStartPosition();
  const uint32_t mate_mapping_start_position = mapping.mrid_ < 0 ? 0 : (mapping.mpos_ + 1);

  // Build the entire SAM record in one string to ensure atomic output
  std::string out;
  out.reserve(256 + mapping.sequence_.size() + mapping.sequence_qual_.size() + mapping.MD_.size());
  
  out.append(mapping.read_name_);
  out.push_back('\t');
  out.append(std::to_string(mapping.flag_));
  out.push_back('\t');
  out.append(reference_sequence_name);
  out.push_back('\t');
  out.append(std::to_string(mapping_start_position));
  out.push_back('\t');
  out.append(std::to_string(mapping.mapq_));
  out.push_back('\t');
  out.append(mapping.GenerateCigarString());
  out.push_back('\t');
  out.append(mate_ref_sequence_name);
  out.push_back('\t');
  out.append(std::to_string(mate_mapping_start_position));
  out.push_back('\t');
  out.append(std::to_string(mapping.tlen_));
  out.push_back('\t');
  out.append(mapping.sequence_);
  out.push_back('\t');
  out.append(mapping.sequence_qual_);
  out.push_back('\t');
  out.append(mapping.GenerateIntTagString("NM", mapping.NM_));
  out.append("\tMD:Z:");
  out.append(mapping.MD_);
  
  if (cell_barcode_length_ > 0) {
    out.append("\tCB:Z:");
    out.append(barcode_translator_.Translate(mapping.cell_barcode_, cell_barcode_length_));
  }
  
  out.push_back('\n');
  
  // Write to primary output
  this->AppendMappingOutput(out);
  
  // Route to Y-filter streams based on read ID
  if (reads_with_y_hit_ && (noY_output_file_ || Y_output_file_)) {
    bool is_y_hit = reads_with_y_hit_->count(mapping.read_id_) > 0;
    
    if (Y_output_file_ && is_y_hit) {
      fwrite(out.data(), 1, out.size(), Y_output_file_);
    }
    if (noY_output_file_ && !is_y_hit) {
      fwrite(out.data(), 1, out.size(), noY_output_file_);
    }
  }
}

template <>
void MappingWriter<SAMMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<SAMMapping> > &mappings) {
  FILE *temp_mapping_output_file =
      fopen(temp_mapping_output_file_path.c_str(), "wb");
  assert(temp_mapping_output_file != NULL);
  for (size_t ri = 0; ri < num_reference_sequences; ++ri) {
    // make sure mappings[ri] exists even if its size is 0
    size_t num_mappings = mappings[ri].size();
    fwrite(&num_mappings, sizeof(size_t), 1, temp_mapping_output_file);
    if (mappings[ri].size() > 0) {
      for (size_t mi = 0; mi < num_mappings; ++mi) {
        mappings[ri][mi].WriteToFile(temp_mapping_output_file);
      }
      // fwrite(mappings[ri].data(), sizeof(MappingRecord), mappings[ri].size(),
      // temp_mapping_output_file);
    }
  }
  fclose(temp_mapping_output_file);
}

// Specialization for pairs format.
template <>
void MappingWriter<PairsMapping>::OutputHeader(uint32_t num_reference_sequences,
                                               const SequenceBatch &reference) {
  std::vector<uint32_t> rid_order;
  rid_order.resize(num_reference_sequences);
  uint32_t i;
  for (i = 0; i < num_reference_sequences; ++i) {
    rid_order[pairs_custom_rid_rank_[i]] = i;
  }
  this->AppendMappingOutput("## pairs format v1.0.0\n#shape: upper triangle\n");
  for (i = 0; i < num_reference_sequences; ++i) {
    uint32_t rid = rid_order[i];
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t reference_sequence_length = reference.GetSequenceLengthAt(rid);
    this->AppendMappingOutput(
        "#chromsize: " + std::string(reference_sequence_name) + " " +
        std::to_string(reference_sequence_length) + "\n");
  }
  this->AppendMappingOutput(
      "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n");
}

template <>
void MappingWriter<PairsMapping>::AppendMapping(uint32_t rid,
                                                const SequenceBatch &reference,
                                                const PairsMapping &mapping) {
  const char *reference_sequence_name1 =
      reference.GetSequenceNameAt(mapping.rid1_);
  const char *reference_sequence_name2 =
      reference.GetSequenceNameAt(mapping.rid2_);
  this->AppendMappingOutput(mapping.read_name_ + "\t" +
                            std::string(reference_sequence_name1) + "\t" +
                            std::to_string(mapping.GetPosition(1)) + "\t" +
                            std::string(reference_sequence_name2) + "\t" +
                            std::to_string(mapping.GetPosition(2)) + "\t" +
                            std::string(1, mapping.GetStrand(1)) + "\t" +
                            std::string(1, mapping.GetStrand(2)) + "\tUU\n");
}

template <>
void MappingWriter<PairsMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<PairsMapping> > &mappings) {
  FILE *temp_mapping_output_file =
      fopen(temp_mapping_output_file_path.c_str(), "wb");
  assert(temp_mapping_output_file != NULL);
  for (size_t ri = 0; ri < num_reference_sequences; ++ri) {
    // make sure mappings[ri] exists even if its size is 0
    size_t num_mappings = mappings[ri].size();
    fwrite(&num_mappings, sizeof(size_t), 1, temp_mapping_output_file);
    if (mappings[ri].size() > 0) {
      for (size_t mi = 0; mi < num_mappings; ++mi) {
        mappings[ri][mi].WriteToFile(temp_mapping_output_file);
      }
      // fwrite(mappings[ri].data(), sizeof(MappingRecord), mappings[ri].size(),
      // temp_mapping_output_file);
    }
  }
  fclose(temp_mapping_output_file);
}

#ifndef LEGACY_OVERFLOW
// Overflow writer methods
template <typename MappingRecord>
void MappingWriter<MappingRecord>::OutputTempMappingsToOverflow(
    uint32_t num_reference_sequences,
    std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs) {
  
  // Initialize thread-local overflow writer if not already done
  if (!tls_overflow_writer_) {
    // Use user-specified temp directory, or let OverflowWriter choose optimal location
    std::string base_dir = mapping_parameters_.temp_directory_path;
    tls_overflow_writer_ = std::unique_ptr<OverflowWriter>(new OverflowWriter(base_dir, "chromap"));
  }
  
  // Write all mappings to overflow files using thread-local writer
  for (uint32_t rid = 0; rid < num_reference_sequences; ++rid) {
    for (const auto& mapping : mappings_on_diff_ref_seqs[rid]) {
      tls_overflow_writer_->Write(rid, mapping);
    }
    mappings_on_diff_ref_seqs[rid].clear();
  }
}

template <typename MappingRecord>
void MappingWriter<MappingRecord>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t num_mappings_in_mem, uint32_t num_reference_sequences,
    const SequenceBatch &reference,
    const khash_t(k64_seq) * barcode_whitelist_lookup_table) {
  
  // All thread-local writers should already be closed by now
  // Use the shared collection of overflow file paths
  if (shared_overflow_file_paths_.empty()) {
    return;
  }
  
  std::cerr << "Processing " << shared_overflow_file_paths_.size() << " overflow files for k-way merge\n";

  double sort_and_dedupe_start_time = GetRealTime();

  // Structure to hold current record from each file for k-way merge
  struct FileRecord {
    uint32_t rid;
    MappingRecord mapping;
    size_t file_index;
    
    bool operator<(const FileRecord& other) const {
      if (rid != other.rid) {
        return rid > other.rid;  // Min-heap: smaller rid first
      }
      return mapping < other.mapping ? false : true;  // Min-heap: smaller mapping first
    }
  };

  // First, scan files to determine which rids they contain
  // This allows us to group files by rid and process rids in ascending order
  // Note: OverflowWriter creates one file per rid, but we scan to handle edge cases
  std::unordered_map<uint32_t, std::vector<size_t>> rid_to_files;
  std::unordered_set<uint32_t> all_rids;
  
  // Scan all files to build rid -> file_index mapping
  for (size_t fi = 0; fi < shared_overflow_file_paths_.size(); ++fi) {
    OverflowReader scanner(shared_overflow_file_paths_[fi]);
    if (!scanner.IsValid()) {
      continue;
    }
    
    // Scan file to find all rids it contains
    uint32_t rid;
    std::string payload;
    std::unordered_set<uint32_t> rids_in_file;
    while (scanner.ReadNext(rid, payload)) {
      rids_in_file.insert(rid);
      all_rids.insert(rid);
    }
    
    // Add this file to the list for each rid it contains
    for (uint32_t r : rids_in_file) {
      rid_to_files[r].push_back(fi);
    }
  }
  
  // Create readers for actual merge (will be opened per-rid to limit FDs)
  std::vector<std::unique_ptr<OverflowReader>> readers;
  readers.resize(shared_overflow_file_paths_.size());

  // Merge and dedupe (replicating logic from ProcessAndOutputMappingsInLowMemory)
  uint32_t last_rid = std::numeric_limits<uint32_t>::max();
  MappingRecord last_mapping = MappingRecord();
  uint32_t num_last_mapping_dups = 0;
  uint64_t num_uni_mappings = 0;
  uint64_t num_multi_mappings = 0;
  uint64_t num_mappings_passing_filters = 0;
  uint64_t num_total_mappings = 0;
  std::vector<MappingRecord> temp_dups_for_bulk_level_dedup;
  temp_dups_for_bulk_level_dedup.reserve(255);

  const bool deduplicate_at_bulk_level_for_single_cell_data =
      mapping_parameters_.remove_pcr_duplicates &&
      !mapping_parameters_.is_bulk_data &&
      mapping_parameters_.remove_pcr_duplicates_at_bulk_level;

  // Process rids in ascending order to preserve coordinate-sorted output
  std::vector<uint32_t> rids(all_rids.begin(), all_rids.end());
  std::sort(rids.begin(), rids.end());

  // For each rid, k-way merge all files containing that rid
  for (uint32_t current_rid : rids) {
    const auto& file_indices = rid_to_files[current_rid];
    
    // Priority queue for k-way merge within this rid
    std::priority_queue<FileRecord> merge_heap;
    
    // Initialize: read first record from each file for this rid
    // Note: OverflowWriter creates one file per rid, so each file contains only records for one rid
    for (size_t fi : file_indices) {
      readers[fi].reset(new OverflowReader(shared_overflow_file_paths_[fi]));
      if (!readers[fi] || !readers[fi]->IsValid()) {
        continue;
      }
      
      uint32_t rid;
      std::string payload;
      if (readers[fi]->ReadNext(rid, payload)) {
        // Verify rid matches (should always be true, but check for safety)
        if (rid == current_rid) {
          FILE* mem_file = fmemopen(const_cast<char*>(payload.data()), payload.size(), "rb");
          if (mem_file) {
            MappingRecord mapping;
            mapping.LoadFromFile(mem_file);
            fclose(mem_file);
            
            FileRecord rec;
            rec.rid = rid;
            rec.mapping = mapping;
            rec.file_index = fi;
            merge_heap.push(rec);
          }
        }
      }
    }

    // K-way merge for this rid
    while (!merge_heap.empty()) {
      FileRecord min_rec = merge_heap.top();
      merge_heap.pop();
      
      ++num_total_mappings;
      
      const MappingRecord& current_min_mapping = min_rec.mapping;
      const uint32_t min_rid = min_rec.rid;
      
      const bool is_first_iteration = num_total_mappings == 1;
      const bool current_mapping_is_duplicated_at_cell_level =
          !is_first_iteration && current_min_mapping == last_mapping;
      const bool current_mapping_is_duplicated_at_bulk_level =
          !is_first_iteration && deduplicate_at_bulk_level_for_single_cell_data &&
          current_min_mapping.IsSamePosition(last_mapping);
      const bool current_mapping_is_duplicated =
          last_rid == min_rid && (current_mapping_is_duplicated_at_cell_level ||
                                  current_mapping_is_duplicated_at_bulk_level);
      
      if (mapping_parameters_.remove_pcr_duplicates &&
          current_mapping_is_duplicated) {
        ++num_last_mapping_dups;
        if (deduplicate_at_bulk_level_for_single_cell_data) {
          if (!temp_dups_for_bulk_level_dedup.empty() &&
              current_min_mapping == temp_dups_for_bulk_level_dedup.back()) {
            temp_dups_for_bulk_level_dedup.back() = current_min_mapping;
            temp_dups_for_bulk_level_dedup.back().num_dups_ += 1;
          } else {
            temp_dups_for_bulk_level_dedup.push_back(current_min_mapping);
            temp_dups_for_bulk_level_dedup.back().num_dups_ = 1;
          }
        }
        if (current_min_mapping.mapq_ > last_mapping.mapq_) {
          last_mapping = current_min_mapping;
        }
      } else {
        if (!is_first_iteration) {
          if (deduplicate_at_bulk_level_for_single_cell_data) {
            size_t best_mapping_index = FindBestMappingIndexFromDuplicates(
                barcode_whitelist_lookup_table, temp_dups_for_bulk_level_dedup);
            last_mapping = temp_dups_for_bulk_level_dedup[best_mapping_index];
            temp_dups_for_bulk_level_dedup.clear();
          }

          if (last_mapping.mapq_ >= mapping_parameters_.mapq_threshold) {
            last_mapping.num_dups_ =
                std::min((uint32_t)std::numeric_limits<uint8_t>::max(),
                         num_last_mapping_dups);
            if (mapping_parameters_.Tn5_shift) {
              last_mapping.Tn5Shift();
            }

            AppendMapping(last_rid, reference, last_mapping);
            ++num_mappings_passing_filters;
            if (!mapping_parameters_.summary_metadata_file_path.empty())
              summary_metadata_.UpdateCount(last_mapping.GetBarcode(), SUMMARY_METADATA_DUP,
                num_last_mapping_dups - 1);
          } else {
            if (!mapping_parameters_.summary_metadata_file_path.empty())
              summary_metadata_.UpdateCount(last_mapping.GetBarcode(), SUMMARY_METADATA_LOWMAPQ, 
                        num_last_mapping_dups);
          }
          if (!mapping_parameters_.summary_metadata_file_path.empty())
            summary_metadata_.UpdateCount(last_mapping.GetBarcode(), SUMMARY_METADATA_MAPPED, 
                  num_last_mapping_dups);

          if (last_mapping.is_unique_ == 1) {
            ++num_uni_mappings;
          } else {
            ++num_multi_mappings;
          }
        }

        last_mapping = current_min_mapping;
        last_rid = min_rid;
        num_last_mapping_dups = 1;

        if (deduplicate_at_bulk_level_for_single_cell_data) {
          temp_dups_for_bulk_level_dedup.push_back(current_min_mapping);
          temp_dups_for_bulk_level_dedup.back().num_dups_ = 1;
        }
      }

      // Read next record from the file we just processed
      // Note: Since OverflowWriter creates one file per rid, all records in this file are for current_rid
      if (min_rec.file_index < readers.size() && readers[min_rec.file_index]) {
        uint32_t rid;
        std::string payload;
        if (readers[min_rec.file_index]->ReadNext(rid, payload)) {
          // Should always be current_rid, but verify for safety
          if (rid == current_rid) {
            FILE* mem_file = fmemopen(const_cast<char*>(payload.data()), payload.size(), "rb");
            if (mem_file) {
              MappingRecord mapping;
              mapping.LoadFromFile(mem_file);
              fclose(mem_file);
              
              FileRecord rec;
              rec.rid = rid;
              rec.mapping = mapping;
              rec.file_index = min_rec.file_index;
              merge_heap.push(rec);
            }
          }
        }
      }
    }
    
    // Close readers for this rid to free file descriptors before moving to next rid
    for (size_t fi : file_indices) {
      readers[fi].reset();
    }
  }

  // Output the last mapping if any
  if (num_total_mappings > 0) {
    if (last_mapping.mapq_ >= mapping_parameters_.mapq_threshold) {
      if (deduplicate_at_bulk_level_for_single_cell_data) {
        size_t best_mapping_index = FindBestMappingIndexFromDuplicates(
            barcode_whitelist_lookup_table, temp_dups_for_bulk_level_dedup);
        last_mapping = temp_dups_for_bulk_level_dedup[best_mapping_index];
        temp_dups_for_bulk_level_dedup.clear();
      }

      last_mapping.num_dups_ = std::min(
          (uint32_t)std::numeric_limits<uint8_t>::max(), num_last_mapping_dups);
      if (mapping_parameters_.Tn5_shift) {
        last_mapping.Tn5Shift();
      }
      AppendMapping(last_rid, reference, last_mapping);
      ++num_mappings_passing_filters;
      
      if (!mapping_parameters_.summary_metadata_file_path.empty())
        summary_metadata_.UpdateCount(last_mapping.GetBarcode(), SUMMARY_METADATA_DUP,
            num_last_mapping_dups - 1);
    } else {
      if (!mapping_parameters_.summary_metadata_file_path.empty())
        summary_metadata_.UpdateCount(last_mapping.GetBarcode(), SUMMARY_METADATA_LOWMAPQ, 
            num_last_mapping_dups);
    }
    if (!mapping_parameters_.summary_metadata_file_path.empty())
      summary_metadata_.UpdateCount(last_mapping.GetBarcode(), SUMMARY_METADATA_MAPPED, 
            num_last_mapping_dups);

    if (last_mapping.is_unique_ == 1) {
      ++num_uni_mappings;
    } else {
      ++num_multi_mappings;
    }
  }

  // Clean up overflow files
  for (const std::string& file_path : shared_overflow_file_paths_) {
    unlink(file_path.c_str());
  }
  shared_overflow_file_paths_.clear();

  if (mapping_parameters_.remove_pcr_duplicates) {
    std::cerr << "Sorted, deduped and outputed mappings in "
              << GetRealTime() - sort_and_dedupe_start_time << "s.\n";
  } else {
    std::cerr << "Sorted and outputed mappings in "
              << GetRealTime() - sort_and_dedupe_start_time << "s.\n";
  }
  std::cerr << "# uni-mappings: " << num_uni_mappings
            << ", # multi-mappings: " << num_multi_mappings
            << ", total: " << num_uni_mappings + num_multi_mappings << ".\n";
  std::cerr << "Number of output mappings (passed filters): "
            << num_mappings_passing_filters << "\n";
}

// Explicit template instantiations for overflow methods
// Only instantiate for types that have WriteToFile/LoadFromFile methods
template void MappingWriter<SAMMapping>::OutputTempMappingsToOverflow(
    uint32_t, std::vector<std::vector<SAMMapping>>&);
template void MappingWriter<SAMMapping>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t, uint32_t, const SequenceBatch&, const khash_t(k64_seq)*);

template void MappingWriter<PAFMapping>::OutputTempMappingsToOverflow(
    uint32_t, std::vector<std::vector<PAFMapping>>&);
template void MappingWriter<PAFMapping>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t, uint32_t, const SequenceBatch&, const khash_t(k64_seq)*);

template void MappingWriter<PairsMapping>::OutputTempMappingsToOverflow(
    uint32_t, std::vector<std::vector<PairsMapping>>&);
template void MappingWriter<PairsMapping>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t, uint32_t, const SequenceBatch&, const khash_t(k64_seq)*);

template void MappingWriter<PairedPAFMapping>::OutputTempMappingsToOverflow(
    uint32_t, std::vector<std::vector<PairedPAFMapping>>&);
template void MappingWriter<PairedPAFMapping>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t, uint32_t, const SequenceBatch&, const khash_t(k64_seq)*);

// Helper function to close thread-local writer and collect paths
template <typename MappingRecord>
void MappingWriter<MappingRecord>::CloseThreadOverflowWriter() {
  if (!tls_overflow_writer_) {
    return;
  }
  
  // Close the thread-local writer and get its file paths
  auto file_paths = tls_overflow_writer_->Close();
  tls_overflow_writer_.reset();
  
  // Add the paths to the shared collection under lock
  {
    std::lock_guard<std::mutex> lock(overflow_paths_mutex_);
    shared_overflow_file_paths_.insert(shared_overflow_file_paths_.end(),
                                       file_paths.begin(), file_paths.end());
  }
}

// Helper function to rotate thread-local writer after each flush
// Closes current writer (collecting file paths) and resets it so next spill creates fresh files
// This ensures each overflow file contains exactly one sorted run for correct k-way merge
template <typename MappingRecord>
void MappingWriter<MappingRecord>::RotateThreadOverflowWriter() {
  if (!tls_overflow_writer_) {
    return;
  }
  
  // Close the thread-local writer and get its file paths
  auto file_paths = tls_overflow_writer_->Close();
  tls_overflow_writer_.reset();  // Next spill will create new files
  
  // Add the paths to the shared collection under lock
  {
    std::lock_guard<std::mutex> lock(overflow_paths_mutex_);
    shared_overflow_file_paths_.insert(shared_overflow_file_paths_.end(),
                                       file_paths.begin(), file_paths.end());
  }
}

// MappingWithBarcode and MappingWithoutBarcode don't have WriteToFile/LoadFromFile methods
// so we provide fallback implementations that use the legacy temp file system
template <>
void MappingWriter<MappingWithBarcode>::OutputTempMappingsToOverflow(
    uint32_t num_reference_sequences,
    std::vector<std::vector<MappingWithBarcode>> &mappings_on_diff_ref_seqs) {
  // Fallback to legacy temp file system for barcode types
  // This should not be called in practice since barcode types typically don't use low-memory mode
  // But we need the symbol to exist for linking
}

template <>
void MappingWriter<MappingWithBarcode>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t num_mappings_in_mem, uint32_t num_reference_sequences,
    const SequenceBatch &reference,
    const khash_t(k64_seq) * barcode_whitelist_lookup_table) {
  // Fallback - should not be called for barcode types
}

template <>
void MappingWriter<MappingWithoutBarcode>::OutputTempMappingsToOverflow(
    uint32_t num_reference_sequences,
    std::vector<std::vector<MappingWithoutBarcode>> &mappings_on_diff_ref_seqs) {
  // Fallback to legacy temp file system for barcode types
  // This should not be called in practice since barcode types typically don't use low-memory mode
  // But we need the symbol to exist for linking
}

template <>
void MappingWriter<MappingWithoutBarcode>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t num_mappings_in_mem, uint32_t num_reference_sequences,
    const SequenceBatch &reference,
    const khash_t(k64_seq) * barcode_whitelist_lookup_table) {
  // Fallback - should not be called for barcode types
}

template <>
void MappingWriter<PairedEndMappingWithBarcode>::OutputTempMappingsToOverflow(
    uint32_t num_reference_sequences,
    std::vector<std::vector<PairedEndMappingWithBarcode>> &mappings_on_diff_ref_seqs) {
  // Fallback - should not be called for paired-end barcode types
}

template <>
void MappingWriter<PairedEndMappingWithBarcode>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t num_mappings_in_mem, uint32_t num_reference_sequences,
    const SequenceBatch &reference,
    const khash_t(k64_seq) * barcode_whitelist_lookup_table) {
  // Fallback - should not be called for barcode types
}

template <>
void MappingWriter<PairedEndMappingWithoutBarcode>::OutputTempMappingsToOverflow(
    uint32_t num_reference_sequences,
    std::vector<std::vector<PairedEndMappingWithoutBarcode>> &mappings_on_diff_ref_seqs) {
  // Fallback - should not be called for paired-end barcode types
}

template <>
void MappingWriter<PairedEndMappingWithoutBarcode>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t num_mappings_in_mem, uint32_t num_reference_sequences,
    const SequenceBatch &reference,
    const khash_t(k64_seq) * barcode_whitelist_lookup_table) {
  // Fallback - should not be called for barcode types
}

// Add explicit instantiation for CloseThreadOverflowWriter
template void MappingWriter<SAMMapping>::CloseThreadOverflowWriter();
template void MappingWriter<PAFMapping>::CloseThreadOverflowWriter();
template void MappingWriter<PairedPAFMapping>::CloseThreadOverflowWriter();
template void MappingWriter<PairsMapping>::CloseThreadOverflowWriter();
template void MappingWriter<MappingWithBarcode>::CloseThreadOverflowWriter();
template void MappingWriter<MappingWithoutBarcode>::CloseThreadOverflowWriter();
template void MappingWriter<PairedEndMappingWithBarcode>::CloseThreadOverflowWriter();
template void MappingWriter<PairedEndMappingWithoutBarcode>::CloseThreadOverflowWriter();

// Add explicit instantiation for RotateThreadOverflowWriter
template void MappingWriter<SAMMapping>::RotateThreadOverflowWriter();
template void MappingWriter<PAFMapping>::RotateThreadOverflowWriter();
template void MappingWriter<PairedPAFMapping>::RotateThreadOverflowWriter();
template void MappingWriter<PairsMapping>::RotateThreadOverflowWriter();
template void MappingWriter<MappingWithBarcode>::RotateThreadOverflowWriter();
template void MappingWriter<MappingWithoutBarcode>::RotateThreadOverflowWriter();
template void MappingWriter<PairedEndMappingWithBarcode>::RotateThreadOverflowWriter();
template void MappingWriter<PairedEndMappingWithoutBarcode>::RotateThreadOverflowWriter();

#endif  // LEGACY_OVERFLOW

}  // namespace chromap
