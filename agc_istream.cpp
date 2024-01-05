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

    std::string sample_data;
    sample_desc_t sample_desc;

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
	// Reset the sample data in case reusing the same object
	this->sample_data = "";
	ZSTD_DCtx_reset(this->zstd_ctx, ZSTD_reset_session_only);

	if (!this->collection_desc->get_sample_desc(sample_name, this->sample_desc)) {
	    throw std::runtime_error("There is no sample " + sample_name);
	}

	for (size_t i = 0; i < this->sample_desc.size(); ++i) {
	    contig_task_t contig_desc(i, "", this->sample_desc[i].first, this->sample_desc[i].second);
	    std::vector<uint8_t> ctg;

	    this->decompress_contig(contig_desc, zstd_ctx, ctg);
	    this->convert_to_alpha(ctg);

	    // Add contig descriptor line
	    this->sample_data += ">";
	    this->sample_data += contig_desc.name_range.str();
	    this->sample_data += '\n';

	    size_t nt_count = 0;
	    for (auto nt : ctg) {
		this->sample_data += nt;
		++nt_count;
		if (nt_count%80 == 0) {
		    this->sample_data += '\n';
		}
	    }
	    if (nt_count%80 != 0) {
		this->sample_data += '\n';
	    }
	}

	return std::istringstream(this->sample_data);
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
