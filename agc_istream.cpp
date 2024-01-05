#include <exception>

#include "zstd.h"

#include "agc_decompressor_lib.h" // CAGCDecompressorLibrary and sample_contig_data_t and contig_task_t
#include "collection_v3.h" // sample_desc_t
#include "queue.h" // CBoundedQueue and CPriorityQueue
#include "genome_io.h" // CGenomeIO
#include "defs.h" // contig_t

class AgcStreamer : public CAGCDecompressorLibrary {
private:
    size_t line_length = 80;
    size_t no_threads = 1;

    std::string sample_data;
    sample_desc_t sample_desc;

public:
    AgcStreamer(std::string _archive_path) : CAGCDecompressorLibrary(true) {
	this->Open(_archive_path);
    }

    std::istringstream get(const std::string &sample_name) {
	if (!this->collection_desc->get_sample_desc(sample_name, this->sample_desc)) {
	    throw std::runtime_error("There is no sample " + sample_name);
	}

	this->q_contig_tasks = std::make_unique<CBoundedQueue<CAGCDecompressorLibrary::contig_task_t>>(1, 1);
	this->pq_contigs_to_save = std::make_unique<CPriorityQueue<CAGCDecompressorLibrary::sample_contig_data_t>>(no_threads * 1);

	std::vector<CAGCDecompressorLibrary::contig_task_t> v_tasks;
	std::vector<std::thread> v_threads;

	v_tasks.reserve(this->sample_desc.size());
	v_tasks.shrink_to_fit();

	for (size_t i = 0; i < this->sample_desc.size(); ++i)
	    v_tasks.emplace_back(i, "", this->sample_desc[i].first, this->sample_desc[i].second);

	std::sort(v_tasks.begin(), v_tasks.end(), [](auto& x, auto& y) {return x.segments.size() > y.segments.size(); });

	this->q_contig_tasks->Restart(1);

	size_t i = 0;
	bool converted_to_alpha = true;
	v_threads.emplace_back([this, i, converted_to_alpha] {

	    ZSTD_DCtx *zstd_ctx = ZSTD_createDCtx();

	    std::vector<uint8_t> ctg;
	    CAGCDecompressorLibrary::contig_task_t contig_desc;

	    while (!this->q_contig_tasks->IsCompleted())
		{
		    if (!this->q_contig_tasks->Pop(contig_desc))
			break;

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

	    this->pq_contigs_to_save->MarkCompleted();

	    ZSTD_freeDCtx(zstd_ctx);
	});

	for (auto& task : v_tasks)
	    this->q_contig_tasks->Push(task, 0);
	this->q_contig_tasks->MarkCompleted();

	this->join_threads(v_threads);

	v_threads.clear();

	this->q_contig_tasks.release();
	this->pq_contigs_to_save.release();

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
    return 1;
}
