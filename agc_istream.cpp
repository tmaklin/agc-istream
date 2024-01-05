#include <exception>

#include "zstd.h"

#include "agc_decompressor_lib.h" // CAGCDecompressorLibrary and sample_contig_data_t and contig_task_t
#include "collection_v3.h" // sample_desc_t
#include "queue.h" // CBoundedQueue and CPriorityQueue
#include "genome_io.h" // CGenomeIO
#include "defs.h" // contig_t

class AgcStreamer : public CAGCDecompressorLibrary {
private:
    ZSTD_DCtx* zstd_ctx;

    sample_desc_t sample_desc;
    std::vector<std::vector<uint8_t>> contigs;
    std::vector<std::string> contig_names;

    size_t line_length = 80;

public:
    AgcStreamer(std::string _archive_path) : CAGCDecompressorLibrary(true) {
	this->Open(_archive_path);
	this->zstd_ctx = ZSTD_createDCtx();
    }

    ~AgcStreamer() {
	ZSTD_freeDCtx(this->zstd_ctx);
    }

    std::istringstream get(const std::string &sample_name) {
	if (!this->collection_desc->get_sample_desc(sample_name, this->sample_desc)) {
	    throw std::runtime_error("There is no sample " + sample_name);
	}

	size_t n_contigs = this->sample_desc.size();
	this->contigs.resize(n_contigs);

	for (size_t i = 0; i < n_contigs; ++i) {
	    contig_task_t contig_desc(i, "", this->sample_desc[i].first, this->sample_desc[i].second);

	    this->decompress_contig(contig_desc, zstd_ctx, this->contigs[i]);
	    this->convert_to_alpha(this->contigs[i]);

	    // Add contig descriptor line
	    contig_names.emplace_back(">" + contig_desc.name_range.str());
	}

	ZSTD_DCtx_reset(this->zstd_ctx, ZSTD_reset_session_only);

	std::vector<std::string> contig_strings(n_contigs);
	for (size_t i = 0; i < n_contigs; ++i) {
	    contig_strings[i] += this->contig_names[i];
	    contig_strings[i] += '\n';
	    contig_strings[i] += std::string(this->contigs[i].begin(), this->contigs[i].end());
	    contig_strings[i] += '\n';
	}
	return std::istringstream(std::accumulate(contig_strings.begin(), contig_strings.end(), std::string("")));
    }
};

int main() {
    std::string sample_name = "ref";

    AgcStreamer streamer("ref.agc");
    std::istringstream in = streamer.get(sample_name);
    size_t line_nr = 0;
    std::string line;
    while (std::getline(in, line)) {
	std::cerr << line_nr << '\t' << line << std::endl;
	++line_nr;
    }

    sample_name = "ref2";
    in = streamer.get(sample_name);
    line_nr = 0;
    while (std::getline(in, line)) {
	std::cerr << line_nr << '\t' << line << std::endl;
	++line_nr;
    }

    return 1;
}
