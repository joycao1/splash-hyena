#include "processreads.h"
#include "../../src/bkc/fq_reader.h"

#include <iostream>
#include <fstream>
#include <ostream>
#include <atomic>
#include <cmath>
#include <cstring>
#include <unordered_set>
#include <random>
#include <filesystem>

#include "../../libs/refresh/hash_tables/lib/murmur_hash.h"
#include "../../libs/refresh/sort/lib/pdqsort_par.h"
#include "../../shared/filters/illumina_adapters_static.h"
#include "../../shared/types/kmer.h"
// #include "/Users/joycao/Documents/salzmanwork/SPLASH/src/common/version.h"

#include "../../libs/mimalloc/include/mimalloc.h"

//#define AGGRESIVE_MEMORY_SAVING
#define USE_READ_COMPRESSION

// *********************************************************************************************
void CBarcodedCounter::enumerate_kmer_leaders_from_read(uint8_t* bases, vector<leader_t>& kmer_leaders)
{
	// cout<< leader_len << ", " << gap_len << ", " << follower_len << endl;
	CKmer leader(leader_len, kmer_mode_t::direct);
	CKmer follower(follower_len, kmer_mode_t::direct);

	int read_len = strlen((char*)bases);
	int follower_start_pos = leader_len + gap_len;

	if (leader_len + gap_len + follower_len > (uint32_t) read_len)
		return;

	for(uint32_t i = 0; i < leader_len-1; ++i)
	{
		uint64_t symbol = dna_code(bases[i]);
		if (symbol < 4)
			leader.insert(symbol);
		else
			leader.Reset();
	}

	for(uint32_t i = follower_start_pos; i < follower_start_pos + follower_len-1; ++i)
	{
		uint64_t symbol = dna_code(bases[i]);
		if (symbol < 4)
			follower.insert(symbol);
		else
			follower.Reset();
	}

	// leader and follower contain almost complete k-mers (without last symbols)

	for (int i = follower_start_pos + follower_len - 1; i < read_len; ++i)
	{
		uint64_t t_symbol = dna_code(bases[i]);
		uint64_t a_symbol = dna_code(bases[i - follower_len - gap_len]);

		if (t_symbol < 4)
			follower.insert(t_symbol);
		else
			follower.Reset();

		if (a_symbol < 4)
			leader.insert(a_symbol);
		else
			leader.Reset();

		if (leader.is_full() && follower.is_full())
		// using data_aligned_dir() instead of data_dir() to ensure the leader is in same format as dictionary
			kmer_leaders.emplace_back(leader.data_aligned_dir());
	}
}

// *********************************************************************************************
void CBarcodedCounter::enumerate_kmer_leaders_for_cbc(cbc_t cbc, vector<leader_t>& kmer_leaders)
{
	kmer_leaders.clear();

	uint64_t file_id;
	uint64_t read_id;

#ifdef USE_READ_COMPRESSION
	vector<uint8_t> decompressed_read;
#endif

	for (auto x : global_cbc_dict[cbc])
	{
		tie(file_id, read_id) = decode_read_id(x);

#ifdef USE_READ_COMPRESSION
		base_coding3.decode_bases(sample_reads[file_id][read_id], decompressed_read);
		enumerate_kmer_leaders_from_read(decompressed_read.data(), kmer_leaders);
#else
		enumerate_kmer_leaders_from_read(sample_reads[file_id][read_id], kmer_leaders);
#endif
	}
}

// *********************************************************************************************
void CBarcodedCounter::enumerate_kmer_pairs_from_read(uint8_t* bases, vector<leader_follower_t>& kmer_pairs)
{

	CKmer leader(leader_len, kmer_mode_t::direct);
	CKmer follower(follower_len, kmer_mode_t::direct);

	int read_len = (int) strlen((char*)bases);
	int follower_start_pos = leader_len + gap_len;

	if (leader_len + gap_len + follower_len > (uint32_t) read_len)
		return;

	for(uint32_t i = 0; i < leader_len-1; ++i)
	{
		uint64_t symbol = dna_code(bases[i]);
		if (symbol < 4)
			leader.insert(symbol);
		else
			leader.Reset();
	}

	for(uint32_t i = follower_start_pos; i < follower_start_pos + follower_len-1; ++i)
	{
		uint64_t symbol = dna_code(bases[i]);
		if (symbol < 4)
			follower.insert(symbol);
		else
			follower.Reset();
	}

	// leader and follower contain almost complete k-mers (without last symbols)

	for (int i = follower_start_pos + follower_len - 1; i < read_len; ++i)
	{
		uint64_t t_symbol = dna_code(bases[i]);
		uint64_t a_symbol = dna_code(bases[i - follower_len - gap_len]);

		if (t_symbol < 4)
			follower.insert(t_symbol);
		else
			follower.Reset();

		if (a_symbol < 4)
			leader.insert(a_symbol);
		else
			leader.Reset();

		if (leader.is_full() && follower.is_full() && (min_leader_count <= 1 || valid_leaders.count(leader.data())))
			kmer_pairs.emplace_back(leader.data_aligned_dir(), follower.data_aligned_dir());
	}
}

// *********************************************************************************************
void CBarcodedCounter::enumerate_kmer_pairs_for_cbc(cbc_t cbc, vector<leader_follower_t>& kmer_pairs)
{
	kmer_pairs.clear();

	uint64_t file_id;
	uint64_t read_id;

#ifdef USE_READ_COMPRESSION
	vector<uint8_t> decompressed_read;
#endif

	for (auto x : global_cbc_dict[cbc])
	{
		tie(file_id, read_id) = decode_read_id(x);

#ifdef USE_READ_COMPRESSION
		base_coding3.decode_bases(sample_reads[file_id][read_id], decompressed_read);

		enumerate_kmer_pairs_from_read(decompressed_read.data(), kmer_pairs);
#else
		enumerate_kmer_pairs_from_read(sample_reads[file_id][read_id], kmer_pairs);
#endif
	}

#ifdef AGGRESIVE_MEMORY_SAVING
	kmer_pairs.shrink_to_fit();
#endif
}

// *********************************************************************************************
void CBarcodedCounter::sort_and_gather_kmer_pairs_for_cbc(vector<leader_follower_t>& kmer_pairs, vector<leader_follower_count_t>& kmer_pair_counts)
{
//	std::sort(kmer_pairs.begin(), kmer_pairs.end());
	refresh::sort::pdqsort(kmer_pairs.begin(), kmer_pairs.end());

	kmer_pair_counts.clear();

	if (kmer_pairs.empty())
		return;

	kmer_pair_counts.emplace_back(kmer_pairs.front());

	for (int i = 1; i < (int) kmer_pairs.size(); ++i)
		if (kmer_pair_counts.back().equal_lf(kmer_pairs[i]))
			kmer_pair_counts.back().count++;
		else
		{
			// if (poly_ACGT_filter.IsPolyACGT(kmer_pair_counts.back().leader, leader_len) ||
			// 		artifacts_filter.ContainsArtifact(kmer_pair_counts.back().leader, leader_len))
			// 	kmer_pair_counts.pop_back();
			kmer_pair_counts.emplace_back(kmer_pairs[i]);
		}
	// if (poly_ACGT_filter.IsPolyACGT(kmer_pair_counts.back().leader, leader_len) ||
	// 		artifacts_filter.ContainsArtifact(kmer_pair_counts.back().leader, leader_len))
	// 	kmer_pair_counts.pop_back();

#ifdef AGGRESIVE_MEMORY_SAVING
	clear_vec(kmer_pairs);
	kmer_pair_counts.shrink_to_fit();
#else
	kmer_pairs.clear();
#endif
}

// *********************************************************************************************
void CBarcodedCounter::extract_fafq_style_anchor_target_pairs(
    const cbc_t& cbc,
    vector<leader_follower_t>& kmer_pairs)
{
    kmer_pairs.clear();

    vector<leader_t> kmer_leaders;
    enumerate_kmer_leaders_for_cbc(cbc, kmer_leaders);  // 1. get anchors (leaders)
	// std::cout<< "collecting accepted anchors" << endl;
	// for (int i = 0; i < 30; i++) {
	// 	cout << kmer_leaders[i] << endl;
	// }
	
    for (const auto& leader : kmer_leaders) {
        if (!accepted_anchors || accepted_anchors->IsAccepted(leader)) {
			// cout << "leader raw = " << leader << endl;
			// cout << "leader str = " << kmer_to_string(leader, leader_len) << endl;
			if (!accepted_anchors) {
				std::cout << "anchor list not used" << endl;
			}
			// if (accepted_anchors->IsAccepted(leader)) {
			// 	std::cout<< "Anchor " << leader << " is accepted" << endl;
			// }
            // Search for target(s) in reads belonging to this CBC
            for (const auto& read_id : global_cbc_dict[cbc]) {
                uint64_t file_id, local_read_id;
                tie(file_id, local_read_id) = decode_read_id(read_id);

#ifdef USE_READ_COMPRESSION
                vector<uint8_t> decompressed_read;
                base_coding3.decode_bases(sample_reads[file_id][local_read_id], decompressed_read);
                const uint8_t* bases = decompressed_read.data();
#else
                const uint8_t* bases = sample_reads[file_id][local_read_id];
#endif
                int read_len = (int)strlen((char*)bases);
				
                for (int i = 0; i <= read_len - (leader_len + gap_len + follower_len); ++i) {
                    // extract candidate leader
                    CKmer anchor_kmer(leader_len, kmer_mode_t::direct);
                    bool valid_leader = true;
                    for (int j = 0; j < leader_len; ++j) {
                        uint64_t sym = dna_code(bases[i + j]);
                        if (sym >= 4) { valid_leader = false; break; }
                        anchor_kmer.insert(sym);
                    }
                    if (!valid_leader || anchor_kmer.data_aligned_dir() != leader) continue;

                    // extract candidate follower
                    CKmer target_kmer(follower_len, kmer_mode_t::direct);
                    bool valid_follower = true;
                    for (int j = 0; j < follower_len; ++j) {
                        uint64_t sym = dna_code(bases[i + leader_len + gap_len + j]);
                        if (sym >= 4) { valid_follower = false; break; }
                        target_kmer.insert(sym);
                    }
                    if (!valid_follower) continue;

                    kmer_pairs.emplace_back(leader, target_kmer.data_aligned_dir());
                }
            }
        }
    }
	// std::cout<< "anchor list updated" <<endl;
}

// *********************************************************************************************
void CBarcodedCounter::count_kmer_pairs()
{
	atomic_int id{ 0 };
	atomic_uint64_t total_no_after_removal{ 0 };

	vector<thread> threads;
	vector<cbc_t> cbcs;
	
	total_no_kmer_pair_counts = 0;
	sum_kmer_pair_counts = 0;

	cbcs.reserve(global_cbc_dict.size());
	for (auto& x : global_cbc_dict)
		cbcs.emplace_back(x.first);

	threads.reserve(no_threads);

	for (int i = 0; i < no_threads; ++i)
		threads.emplace_back([&] {
		int curr_id = -1;

		vector<leader_follower_t> kmer_pairs;
		vector<leader_follower_count_t> kmer_pair_counts;

		vector<vector<bkc_record_t>> record_buffers;

		vector<uint8_t> packed_buffer;

		record_buffers.resize(no_splits);
		// cout<< leader_len<< ", " << gap_len << ", " << follower_len << ", " << zstd_level << endl;

//		zstd_in_memory zim{ (int) zstd_level };
//		vector<uint8_t> zstd_working_space;

		while (true)
		{
			curr_id = id.fetch_add(1);
			if (curr_id >= (int)cbcs.size())
				break;

			// enumerate_kmer_pairs_for_cbc(cbcs[curr_id], kmer_pairs);
			// std::cout << "extracting pairs fafq style" << endl;
			extract_fafq_style_anchor_target_pairs(cbcs[curr_id], kmer_pairs);
			// std::cout << "sort and gathering" << endl;
			sort_and_gather_kmer_pairs_for_cbc(kmer_pairs, kmer_pair_counts);
			// filter_rare_leader_sample_cbc(kmer_pair_counts);
			store_kmer_pairs(cbcs[curr_id], kmer_pair_counts, record_buffers);

			for (uint32_t i = 0; i < no_splits; ++i)
				if ((int) record_buffers[i].size() >= max_records_in_buffer)
				{
					pack_records(record_buffers[i], packed_buffer);
					bkc_files[i]->AddPacked(packed_buffer);

/*					zstd_working_space.resize(packed_buffer.size() + zim.get_overhead(packed_buffer.size()));
					auto packed_size = zim.compress(packed_buffer.data(), packed_buffer.size(), zstd_working_space.data(), zstd_working_space.size(), (int) zstd_level);
					zstd_working_space.resize(packed_size);
					bkc_files[i]->AddPacked(zstd_working_space);*/

					record_buffers[i].clear();
				}
		}

		for (uint32_t i = 0; i < no_splits; ++i)
		{
			pack_records(record_buffers[i], packed_buffer);
			bkc_files[i]->AddPacked(packed_buffer);
		}
		});

	join_threads(threads);

	if (verbosity_level >= 2)
	{
		std::cerr << "Total no. k-mer pair counts: " << total_no_kmer_pair_counts << endl;
		std::cerr << "Sum of k-mer pair counts: " << sum_kmer_pair_counts << endl;
	}
}

// *********************************************************************************************
void CBarcodedCounter::pack_records(vector<bkc_record_t>& records, vector<uint8_t>& packed_buffer)
{
	vector<uint8_t> rec_prev, rec_curr;

#ifdef COMPACT_ENCODING
	packed_buffer.clear();

	for (auto& x : records)
	{
		rec_curr.clear();

		append_int_msb(rec_curr, x.sample_id, sample_id_size_in_bytes);
		append_int_msb(rec_curr, x.barcode, barcode_size_in_bytes);
		append_int_msb(rec_curr, x.leader, leader_size_in_bytes);
		append_int_msb(rec_curr, x.follower, follower_size_in_bytes);
		append_int_msb(rec_curr, x.count, counter_size_in_bytes);

		encode_shared_prefix(packed_buffer, rec_prev, rec_curr);

		swap(rec_prev, rec_curr);
	}
#else
	save_int_lsb(sample_id, sample_id_size_in_bytes);
	save_int_lsb(barcode, barcode_size_in_bytes);
	save_int_lsb(leader, leader_size_in_bytes);
	save_int_lsb(follower, follower_size_in_bytes);
	save_int_lsb(count, counter_size_in_bytes);
#endif
}

// *********************************************************************************************
void CBarcodedCounter::store_kmer_pairs(cbc_t cbc, vector<leader_follower_count_t>& kmer_pair_counts, vector<vector<bkc_record_t>>& record_buffers)
{
	total_no_kmer_pair_counts += kmer_pair_counts.size();

	string cbc_str = base_coding4.decode_bases_2b(cbc, cbc_len);

	refresh::MurMur64Hash mh;

	uint64_t sum = 0;
	for (const auto& x : kmer_pair_counts)
	{
		uint64_t h = mh(x.leader) % no_splits;

		record_buffers[h].emplace_back(sample_id, cbc, x.leader, x.follower, x.count);
		sum += x.count;
	}
	
	sum_kmer_pair_counts += sum;

#ifdef AGGRESIVE_MEMORY_SAVING
	clear_vec(kmer_pair_counts);
#else
	kmer_pair_counts.clear();
#endif
}

// *********************************************************************************************
void CBarcodedCounter::store_kmers(cbc_t cbc, vector<kmer_count_t>& kmer_counts, vector<vector<bkc_record_t>>& record_buffers)
{
	total_no_kmer_counts += kmer_counts.size();

	string cbc_str = base_coding4.decode_bases_2b(cbc, cbc_len);

	refresh::MurMur64Hash mh;

	uint64_t sum = 0;
	for (const auto& x : kmer_counts)
	{
		uint64_t h = mh(x.kmer) % no_splits;

		record_buffers[h].emplace_back(sample_id, cbc, x.kmer, 0, x.count);
		sum += x.count;
	}
	
	sum_kmer_pair_counts += sum;

#ifdef AGGRESIVE_MEMORY_SAVING
	clear_vec(kmer_counts);
#else
	kmer_counts.clear();
#endif
}

// *********************************************************************************************
string CBarcodedCounter::kmer_to_string(uint64_t kmer, int len)
{
	string str;

	for (int i = 0; i < len; ++i)
	{
		auto c = kmer & 3;
		str.push_back("ACGT"[c]);
		kmer >>= 2;
	}

	reverse(str.begin(), str.end());

	return str;
}

// *********************************************************************************************
bool CBarcodedCounter::ProcessReads()
{
	set_read_file_names();

	std::cout<<"filenames set"<<endl;

	if (!no_threads || file_names.empty())
		return false;

	no_reading_threads = max(min(no_threads / 2, (int)file_names.size()), 1);

	if (verbosity_level >= 1)
		std::cerr << "Reads loading\n";

	reinit_queues();

	init_bkc_files();

	start_reading_threads();
	start_reads_loading_threads();

	join_threads(reading_threads);
	join_threads(reads_loading_threads);
	mi_collect(true);

	std::cout<<"threads read loaded & joined"<<endl;

	if (verbosity_level >= 2)
	{
		std::cerr << "Total no. of loaded reads: " << a_total_no_reads << endl;
		std::cerr << "Total len of loaded reads: " << a_total_read_len << endl;
	}

	times.emplace_back("Reads loading", high_resolution_clock::now());

#if 0			// Currently not used
	if (min_leader_count > 1)
	{
		if (verbosity_level >= 1)
			std::cerr << "Enumerating and counting leader k-mers\n";
		count_leaders();
		times.emplace_back("Enumerating and counting leader k-mers", high_resolution_clock::now());

		if (verbosity_level >= 1)
			std::cerr << "Determine valid leaders\n";
		determine_valid_leaders();
		times.emplace_back("Determine valid leaders", high_resolution_clock::now());
	}
#endif

	// if (counting_mode == counting_mode_t::single)
	// {
	// 	if (verbosity_level >= 1)
	// 		std::cerr << "Enumerating and counting k-mers\n";
	// 	count_kmers();
	// 	mi_collect(true);
	// 	times.emplace_back("Enumerating and counting k-mers", high_resolution_clock::now());
	// }
	// else
	// {
		if (verbosity_level >= 1)
			std::cerr << "Enumerating and counting leader-follower pairs\n";
		std::cout<<"counting kmer pairs"<<endl;
		count_kmer_pairs();
		mi_collect(true);
		times.emplace_back("Enumerating and counting leader-follower pairs", high_resolution_clock::now());
	// }

	mma.clear();

	bkc_files.clear();

	return true;
}

// EOF
