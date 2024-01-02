#define protected public // I don't want to make yet another derived class from CAGCDecompressorLibrary

#include "zstd.h"

#include "agc_decompressor_lib.h" // CAGCDecompressorLibrary and sample_contig_data_t and contig_task_t
#include "collection_v3.h" // sample_desc_t
#include "queue.h" // CBoundedQueue and CPriorityQueue
#include "genome_io.h" // CGenomeIO
#include "defs.h" // contig_t

void start_decompressing_threads2(CAGCDecompressorLibrary &agc, std::vector<std::thread>& v_threads, const uint32_t n_t, bool converted_to_alpha, std::ostream &out)
{
    for (uint32_t i = 0; i < n_t; ++i)
	v_threads.emplace_back([&agc, i, converted_to_alpha, &out] {

	    ZSTD_DCtx *zstd_ctx = ZSTD_createDCtx();

	    std::vector<uint8_t> ctg;
	    CAGCDecompressorLibrary::contig_task_t contig_desc;

	    while (!agc.q_contig_tasks->IsCompleted())
		{
		    if (!agc.q_contig_tasks->Pop(contig_desc))
			break;

		    if (!agc.decompress_contig(contig_desc, zstd_ctx, ctg))
			continue;

		    agc.convert_to_alpha(ctg);

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

	    agc.pq_contigs_to_save->MarkCompleted();

	    ZSTD_freeDCtx(zstd_ctx);
	});
}

int main() {
    std::vector<std::string> sample_names = { "ref" };
    size_t no_threads = 1;
    size_t _line_length = 80;
    std::string _file_name;

    CAGCDecompressorLibrary agc_d(true);

    agc_d.Open("ref.agc", false);

    std::vector<sample_desc_t> v_sample_desc;

    for (const auto &s : sample_names)
	{
	    sample_desc_t sample_desc;

	    if (!agc_d.collection_desc->get_sample_desc(s, sample_desc))
		{
		    std::cerr << "There is no sample " << s << std::endl;

		    return false;
		}

	    v_sample_desc.emplace_back(sample_desc);
	}

    agc_d.q_contig_tasks = std::make_unique<CBoundedQueue<CAGCDecompressorLibrary::contig_task_t>>(1, 1);
    agc_d.pq_contigs_to_save = std::make_unique<CPriorityQueue<CAGCDecompressorLibrary::sample_contig_data_t>>(no_threads * v_sample_desc.size());

    std::vector<CAGCDecompressorLibrary::contig_task_t> v_tasks;
    std::vector<std::thread> v_threads;

    // Saving thread
    std::thread gio_thread([&] {
	CGenomeIO gio;
	CAGCDecompressorLibrary::sample_contig_data_t ctg;

	gio.Open(_file_name, true);

	while (!agc_d.pq_contigs_to_save->IsCompleted())
	    {
		if (!agc_d.pq_contigs_to_save->Pop(ctg))
		    break;

		gio.SaveContig(ctg.contig_name, ctg.contig_data, _line_length);
	    }

	gio.Close();
    });

    v_threads.clear();
    v_threads.reserve(no_threads);

    uint32_t global_i = 0;

    std::ostringstream text;

    for (auto& sample_desc : v_sample_desc)
	{
	    v_tasks.clear();
	    v_tasks.shrink_to_fit();
	    v_tasks.reserve(sample_desc.size());

	    for (uint32_t i = 0; i < sample_desc.size(); ++i, ++global_i)
		v_tasks.emplace_back(global_i, "", sample_desc[i].first, sample_desc[i].second);

	    std::sort(v_tasks.begin(), v_tasks.end(), [](auto& x, auto& y) {return x.segments.size() > y.segments.size(); });

	    agc_d.q_contig_tasks->Restart(1);

	    start_decompressing_threads2(agc_d, v_threads, 1, false, text);

	    for (auto& task : v_tasks)
		agc_d.q_contig_tasks->Push(task, 0);
	    agc_d.q_contig_tasks->MarkCompleted();

	    agc_d.join_threads(v_threads);

	    v_threads.clear();
	}

    gio_thread.join();

    agc_d.q_contig_tasks.release();
    agc_d.pq_contigs_to_save.release();

    std::cerr << text.str();

    return 1;
}
