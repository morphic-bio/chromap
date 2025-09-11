#include "mapping_writer.h"

namespace chromap {

#ifdef NEW_OVERFLOW
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
    this->AppendMappingOutput(
        "@SQ\tSN:" + std::string(reference_sequence_name) +
        "\tLN:" + std::to_string(reference_sequence_length) + "\n");
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
  
  this->AppendMappingOutput(out);
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

#ifdef NEW_OVERFLOW
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
  
  // Simple implementation: read all overflow files and output directly
  // This is a basic version - could be optimized with proper k-way merge
  for (const std::string& file_path : shared_overflow_file_paths_) {
    OverflowReader reader(file_path);
    if (!reader.IsValid()) {
      continue;
    }
    
    uint32_t rid;
    std::string payload;
    while (reader.ReadNext(rid, payload)) {
      // Create a temporary file from the payload to use existing LoadFromFile
      FILE* mem_file = fmemopen(const_cast<char*>(payload.data()), payload.size(), "rb");
      if (mem_file) {
        MappingRecord mapping;
        mapping.LoadFromFile(mem_file);
        fclose(mem_file);
        
        // Output the mapping using existing logic
        AppendMapping(rid, reference, mapping);
      }
    }
    
    // Clean up the overflow file
    unlink(file_path.c_str());
  }
  
  shared_overflow_file_paths_.clear();
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
void MappingWriter<PairedEndMappingWithBarcode>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t num_mappings_in_mem, uint32_t num_reference_sequences,
    const SequenceBatch &reference,
    const khash_t(k64_seq) * barcode_whitelist_lookup_table) {
  // Fallback - should not be called for barcode types
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

#endif

}  // namespace chromap
