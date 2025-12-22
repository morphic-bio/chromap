CXX=g++

# htslib submodule
HTSLIB_DIR=third_party/htslib
ifeq (,$(wildcard $(HTSLIB_DIR)/htslib/sam.h))
$(error "htslib submodule not initialized. Run: git submodule update --init --recursive")
endif

CXXFLAGS=-std=c++11 -Wall -O3 -fopenmp -msse4.1 -I$(HTSLIB_DIR)
LDFLAGS=-L$(HTSLIB_DIR) -lhts -lm -lz -lpthread -lcurl -lcrypto -lbz2 -llzma -ldeflate

cpp_source=sequence_batch.cc index.cc minimizer_generator.cc candidate_processor.cc alignment.cc feature_barcode_matrix.cc ksw.cc draft_mapping_generator.cc mapping_generator.cc mapping_writer.cc overflow_writer.cc overflow_reader.cc bam_sorter.cc chromap.cc chromap_driver.cc
src_dir=src
objs_dir=objs
objs+=$(patsubst %.cc,$(objs_dir)/%.o,$(cpp_source))

exec=chromap

ifneq ($(asan),)
	CXXFLAGS+=-fsanitize=address -g
	LDFLAGS+=-fsanitize=address -ldl -g
endif

ifneq ($(LEGACY_OVERFLOW),)
	CXXFLAGS+=-DLEGACY_OVERFLOW
endif

all: dir $(exec) 
	
dir:
	mkdir -p $(objs_dir)

$(exec): $(objs)
	$(CXX) $(CXXFLAGS) $(objs) -o $(exec) $(LDFLAGS)
	
$(objs_dir)/%.o: $(src_dir)/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean test-unit
clean:
	-rm -rf $(exec) $(objs_dir)

# Unit test for Y-filtering
test-unit: dir
	@mkdir -p tests
	$(CXX) $(CXXFLAGS) -I$(src_dir) tests/test_y_filter.cc $(src_dir)/sequence_batch.cc -o tests/test_y_filter $(LDFLAGS)
	./tests/test_y_filter
