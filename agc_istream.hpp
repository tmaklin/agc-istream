// agc_istream: istream [agc](https://github.com/refresh-bio/agc) archives.
//
// MIT License
//
// Copyright (c) 2023 Tommi Mäklin <tommi 'at' maklin.fi>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
#include "zstd.h"
#include "agc_decompressor_lib.h"

namespace agc {
class AgcStreamer : public CAGCDecompressorLibrary {
private:
    ZSTD_DCtx* zstd_ctx;
    std::string archive_path;

    sample_desc_t sample_desc;
    std::vector<std::vector<uint8_t>> contigs;
    std::vector<std::string> contig_names;

    size_t line_length = 80;

public:
    AgcStreamer(std::string _archive_path) : CAGCDecompressorLibrary(true) {
	this->archive_path = _archive_path;
	this->Open(this->archive_path);
	this->zstd_ctx = ZSTD_createDCtx();
    }

    ~AgcStreamer() {
	ZSTD_freeDCtx(this->zstd_ctx);
    }

    int extract(const std::string &sample_name) {
	if (!this->collection_desc->get_sample_desc(sample_name, this->sample_desc)) {
	    return 0;
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

	return 1;
    }

    size_t n_contigs() const { return this->contigs.size(); }
    const std::vector<std::vector<uint8_t>>& get_contigs() const { return this->contigs; }
    const std::vector<std::string>& get_contig_names() const { return this->contig_names; }

    std::vector<std::string> get_contig_strings() {
        std::vector<std::string> contig_strings(this->n_contigs());
	for (size_t i = 0; i < this->n_contigs(); ++i) {
	    contig_strings[i] += this->contig_names[i];
	    contig_strings[i] += '\n';
	    contig_strings[i] += std::string(this->contigs[i].begin(), this->contigs[i].end());
	    contig_strings[i] += '\n';
	}
	return contig_strings;
    }

};

class istream : public std::istringstream {
private:
    AgcStreamer stream;

public:
    istream(std::string _archive_path) : stream(_archive_path) {};

    void find(const std::string &sample_name) {
	int ret = this->stream.extract(sample_name);
	if (ret == 1) {
	    const std::vector<std::string> &contig_strings = this->stream.get_contig_strings();
	    const std::string &concat = std::accumulate(contig_strings.begin(), contig_strings.end(), std::string(""));
	    std::istringstream contigs(concat);
	    this->swap(contigs);
	} else {
	    this->setstate(std::ios_base::failbit);
	}
    }
};
}
