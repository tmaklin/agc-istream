#define protected public // I don't want to make yet another derived class from CAGCDecompressorLibrary

#include <exception>

#include "zstd.h"

#include "agc_decompressor_lib.h" // CAGCDecompressorLibrary and sample_contig_data_t and contig_task_t
#include "collection_v3.h" // sample_desc_t
#include "queue.h" // CBoundedQueue and CPriorityQueue
#include "genome_io.h" // CGenomeIO
#include "defs.h" // contig_t

std::istringstream open_agc_istream(const std::string &sample_name) {
    size_t line_length = 80;
    size_t no_threads = 1;

    CAGCDecompressorLibrary agc_d(true);
    agc_d.Open("ref.agc", false);

    sample_desc_t sample_desc;

    if (!agc_d.collection_desc->get_sample_desc(sample_name, sample_desc)) {
	throw std::runtime_error("There is no sample " + sample_name);
    }

    agc_d.q_contig_tasks = std::make_unique<CBoundedQueue<CAGCDecompressorLibrary::contig_task_t>>(1, 1);
    agc_d.pq_contigs_to_save = std::make_unique<CPriorityQueue<CAGCDecompressorLibrary::sample_contig_data_t>>(no_threads * 1);

    std::vector<CAGCDecompressorLibrary::contig_task_t> v_tasks;
    std::vector<std::thread> v_threads;

    // Saving thread
    std::thread gio_thread([&] {
	CGenomeIO gio;
	CAGCDecompressorLibrary::sample_contig_data_t ctg;

	gio.Open(sample_name, true);

	while (!agc_d.pq_contigs_to_save->IsCompleted())
	    {
		if (!agc_d.pq_contigs_to_save->Pop(ctg))
		    break;

		gio.SaveContig(ctg.contig_name, ctg.contig_data, line_length);
	    }

	gio.Close();
    });

    v_threads.clear();
    v_threads.reserve(no_threads);

    std::ostringstream out;

    v_tasks.clear();
    v_tasks.shrink_to_fit();
    v_tasks.reserve(sample_desc.size());

    for (size_t i = 0; i < sample_desc.size(); ++i)
	v_tasks.emplace_back(i, "", sample_desc[i].first, sample_desc[i].second);

    std::sort(v_tasks.begin(), v_tasks.end(), [](auto& x, auto& y) {return x.segments.size() > y.segments.size(); });

    agc_d.q_contig_tasks->Restart(1);

    size_t i = 0;
    bool converted_to_alpha = true;
    v_threads.emplace_back([&agc_d, i, converted_to_alpha, &out] {

	ZSTD_DCtx *zstd_ctx = ZSTD_createDCtx();

	std::vector<uint8_t> ctg;
	CAGCDecompressorLibrary::contig_task_t contig_desc;

	while (!agc_d.q_contig_tasks->IsCompleted())
	    {
		if (!agc_d.q_contig_tasks->Pop(contig_desc))
		    break;

		agc_d.decompress_contig(contig_desc, zstd_ctx, ctg);
		agc_d.convert_to_alpha(ctg);

		out << ">" << contig_desc.name_range.str() << '\n';
		size_t nt_count = 0;
		for (auto nt : ctg) {
		    out << nt;
		    ++nt_count;
		    if (nt_count%80 == 0) {
			out << '\n';
		    }
		}
		out << std::endl;
	    }

	agc_d.pq_contigs_to_save->MarkCompleted();

	ZSTD_freeDCtx(zstd_ctx);
    });

    for (auto& task : v_tasks)
	agc_d.q_contig_tasks->Push(task, 0);
    agc_d.q_contig_tasks->MarkCompleted();

    agc_d.join_threads(v_threads);

    v_threads.clear();

    gio_thread.join();

    agc_d.q_contig_tasks.release();
    agc_d.pq_contigs_to_save.release();

    return std::istringstream(out.str());
}

int main() {
    std::string sample_name = "ref";

    std::istringstream in = open_agc_istream(sample_name);
    
    std::string line;
    size_t line_nr = 0;
    while (std::getline(in, line)) {
	std::cerr << line_nr << '\t' << line << std::endl;
	++line_nr;
    }

    return 1;
}
