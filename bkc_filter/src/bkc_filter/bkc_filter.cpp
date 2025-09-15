/*
 * File: bkc.cpp
 * Description: This file contains the implementation of the BKC (Barcoded K-mer Counter) tool.
 *              It processes barcoded reads to count k-mers or k-mer pairs, applying various filters
 *              and corrections based on user-specified parameters.
 * 
 * Main Functions:
 * - parse_args: Parses command-line arguments and sets the parameters for processing.
 * - usage: Displays usage information and available options for the tool.
 * - load_predefined_cbc_visium: Loads predefined CBCs for Visium technology.
 * - load_predefined_cbc_plain: Loads predefined CBCs for plain technology.
 * - main: Entry point of the program, orchestrates the processing of barcoded reads.
 * 
 * Dependencies:
 * - Includes standard libraries for I/O, string manipulation, and regex.
 * - Includes custom headers like "kmer_counter.h" and "params.h" for specific functionality.
 * 
 * Author: [Joy Cao]
 * Date: [7/16/2025]
 * Version: [Insert Version]
 */

 #include <string>
 #include <iostream>
 #include <fstream>
 #include <cstdint>
 #include <regex>
 #include "processreads.h"
//  #include "/Users/joycao/Documents/salzmanwork/SPLASH/src/common/version.h"
 
 CParams params;
 
 shared_ptr<AcceptedAnchors> CBarcodedCounter::accepted_anchors;

array<uint8_t, 256> CBarcodedCounter::char2bits;

void CBarcodedCounter::init() {
    std::fill(char2bits.begin(), char2bits.end(), 4);
    char2bits['A'] = char2bits['a'] = 0;
    char2bits['C'] = char2bits['c'] = 1;
    char2bits['G'] = char2bits['g'] = 2;
    char2bits['T'] = char2bits['t'] = 3;
}

 void usage();
 bool parse_args(int argc, char** argv);
 bool load_predefined_cbc_visium();
 bool load_predefined_cbc_plain();
 
 // remember to include in class definition
 bool load_strings(vector<string>& vec, const string& fn)
 {
     ifstream ifs(fn);
     if (!ifs) {
         cerr << "Error: Cannot open: " << fn << endl;
         return false;
     }
 
     istream_iterator<string> ie;
     vec.assign(istream_iterator<string>(ifs), ie);
 
     return true;
 }
 
 bool CBarcodedCounter::prepare_anchor_dict(const vector<string>& vec)
 {
     vector<kmer_t> anchors;
     anchors.reserve(vec.size());
 
     for (const auto& s : vec)
     {
         if (s.size() != params.leader_len.get()) {
             cerr << "Error: Wrong anchor length: " << s << endl;
             return false;
         }
 
         kmer_t anchor = 0;
         for (auto c : s) {
			kmer_t x = char2bits[(uint8_t) c];
			if (x > 3)
			{
				cerr << "Error: Anchor contains strange symbols: " << s << endl;
				return false;
			}

			anchor <<= 2;
			anchor += x;
		}
         anchors.emplace_back(anchor);
     }

     for (size_t i = 0; i < anchors.size() && i < 30; i++) {
        cout << vec[i] << " => " << anchors[i] << endl;
    }
 
     accepted_anchors = make_shared<AcceptedAnchors>(anchors, params.leader_len.get());
     return true;
 }
 
 // *********************************************************************************************
 bool parse_args(int argc, char** argv)
 {
     // Parses command-line arguments and sets the parameters for processing.
     // Returns false if any argument is invalid or missing required values.
     string input_name;
     string dict_name;
 
     for (int i = 1; i < argc; ++i)
     {
         if (argv[i] == "--leader_len"s && i + 1 < argc)
         {
             if (!params.leader_len.set(atoi(argv[++i])))
             {
                 cerr << "Incorrect value for leader_len: " << argv[i] << endl;
                 return false;
             }
         }
         else if (argv[i] == "--follower_len"s && i + 1 < argc)
         {
             if (!params.follower_len.set(atoi(argv[++i])))
             {
                 cerr << "Incorrect value for follower_len: " << argv[i] << endl;
                 return false;
             }
         }
         else if (argv[i] == "--gap_len"s && i + 1 < argc)
         {
             if (!params.gap_len.set(atoi(argv[++i])))
             {
                 cerr << "Incorrect value for gap_len: " << argv[i] << endl;
                 return false;
             }
         }
         else if (string(argv[i])== "-d" && i + 1 < argc) {
             dict_name = argv[++i];
         }

 
         else if (argv[i] == "--n_splits"s && i + 1 < argc)
         {
             if (!params.no_splits.set(atoi(argv[++i])))
             {
                 cerr << "Incorrect value for n_splits: " << argv[i] << endl;
                 return false;
             }
         }
         else if (argv[i] == "--cbc_len"s && i + 1 < argc)
         {
             if (!params.cbc_len.set(atoi(argv[++i])))
             {
                 cerr << "Incorrect value for cbc_len: " << argv[i] << endl;
                 return false;
             }
         }
         else if (argv[i] == "--umi_len"s && i + 1 < argc)
         {
             if (!params.umi_len.set(atoi(argv[++i])))
             {
                 cerr << "Incorrect value for umi_len: " << argv[i] << endl;
                 return false;
             }
         }
         else if (argv[i] == "--soft_cbc_umi_len_limit"s && i + 1 < argc)
         {
             if (!params.soft_cbc_umi_len_limit.set(atoi(argv[++i])))
             {
                 cerr << "Incorrect value for soft_cbc_umi_len_limit: " << argv[i] << endl;
                 return false;
             }
         }
         else if (argv[i] == "--n_threads"s && i + 1 < argc)
         {
             if (!params.no_threads.set(atoi(argv[++i])))
             {
                 cerr << "Incorrect value for n_threads: " << argv[i] << endl;
                 return false;
             }
         }
         else if (argv[i] == "--zstd_level"s && i + 1 < argc)
         {
             if (!params.zstd_level.set(atoi(argv[++i])))
             {
                 cerr << "Incorrect value for zstd_level: " << argv[i] << endl;
                 return false;
             }
         }
         else if (argv[i] == "--verbose"s && i + 1 < argc)
         {
             if (!params.verbosity_level.set(atoi(argv[++i])))
             {
                 cerr << "Incorrect value for verbosity_level: " << argv[i] << endl;
                 return false;
             }
         }
         else if (argv[i] == "--leader_sample_counts_threshold"s && i + 1 < argc)
         {
             if (!params.rare_leader_thr.set(atoi(argv[++i])))
             {
                 cerr << "Incorrect value for leader_sample_counts_threshold: " << argv[i] << endl;
                 return false;
             }
         }
         else if (argv[i] == "--cbc_filtering_thr"s && i + 1 < argc)
         {
             if (!params.cbc_filtering_thr.set(atoi(argv[++i])))
             {
                 cerr << "Incorrect value for cbc_filtering_thr: " << argv[i] << endl;
                 return false;
             }
         }
         else if (argv[i] == "--max_count"s && i + 1 < argc)
         {
             if (!params.max_count.set(atoi(argv[++i])))
             {
                 cerr << "Incorrect value for max_count: " << argv[i] << endl;
                 return false;
             }
         }
         else if (argv[i] == "--sample_id"s && i + 1 < argc)
             params.sample_id = atoi(argv[++i]);
         else if (argv[i] == "--output_name"s && i + 1 < argc)
             params.out_file_name = argv[++i];
         else if (argv[i] == "--poly_ACGT_len"s && i + 1 < argc)
         {
             if (!params.poly_ACGT_len.set(atoi(argv[++i])))
             {
                 cerr << "Incorrect value for poly_ACGT_len: " << argv[i] << endl;
                 return false;
             }
         }
         else if (argv[i] == "--canonical"s)
             params.canonical_mode = true;
         else if (argv[i] == "--artifacts"s && i + 1 < argc)
             params.artifacts = argv[++i];
         else if (argv[i] == "--apply_filter_illumina_adapters"s)
             params.apply_filter_illumina_adapters = true;
         else if (argv[i] == "--apply_cbc_correction"s)
             params.apply_cbc_correction = true;
         else if (argv[i] == "--log_name"s && i + 1 < argc)
         {
             params.export_cbc_logs = true;
             params.cbc_log_file_name = argv[++i];
         }
         else if (argv[i] == "--filtered_input_path"s && i + 1 < argc)
             params.filtered_input_path = argv[++i];
         else if (argv[i] == "--export_filtered_input_mode"s && i + 1 < argc)
         {
             ++i;
             if (argv[i] == "none"s)
                 params.export_filtered_input = export_filtered_input_t::none;
             else if (argv[i] == "first"s)
                 params.export_filtered_input = export_filtered_input_t::first;
             else if (argv[i] == "second"s)
                 params.export_filtered_input = export_filtered_input_t::second;
             else if (argv[i] == "both"s)
                 params.export_filtered_input = export_filtered_input_t::both;
             else
             {
                 cerr << "Wrong value for filtered_input_mode: " << argv[i] << endl;
                 return false;
             }
         }
         else if (argv[i] == "--allow_strange_cbc_umi_reads"s)
             params.allow_strange_cbc_umi_reads = true;
         else if (argv[i] == "--input_name"s && i + 1 < argc)
             input_name = argv[++i];
         else if (argv[i] == "--technology"s && i + 1 < argc)
         {
             ++i;
             if (argv[i] == "10X"s || argv[i] == "10x"s)
                 params.technology = technology_t::ten_x;
             else if (argv[i] == "visium"s)
                 params.technology = technology_t::visium;
             else
             {
                 cerr << "Wrong value for technology: " << argv[i] << endl;
                 return false;
             }
         }
         else if (argv[i] == "--input_format"s && i + 1 < argc)
         {
             ++i;
             params.input_format = input_format_from_string(argv[i]);
             if (params.input_format == input_format_t::unknown)
             {
                 cerr << "Wrong value for input_format: " << argv[i] << endl;
                 return false;
             }
         }
         else if (argv[i] == "--output_format"s && i + 1 < argc)
         {
             ++i;
             params.output_format = output_format_from_string(argv[i]);
             if (params.output_format == output_format_t::unknown)
             {
                 cerr << "Wrong value for output_format: " << argv[i] << endl;
                 return false;
             }
         }
         else if (argv[i] == "--mode"s && i + 1 < argc)
         {
             ++i;
             params.counting_mode = counting_mode_from_string(argv[i]);
             if (params.counting_mode == counting_mode_t::unknown)
             {
                 cerr << "Wrong value for mode: " << argv[i] << endl;
                 return false;
             }
         }
         else if (argv[i] == "--predefined_cbc"s && i + 1 < argc)
             params.predefined_cbc_fn = argv[++i];
         else
         {
             cerr << "Unknown parameter: " << argv[i] << endl;
             return false;
         }
     }

     vector<string> ad;
         if (!dict_name.empty()) {
             if (!load_strings(ad, dict_name)) {
                 cerr << "Error: Could not load dictionary file\n";
                 return false;
             }
             if (!CBarcodedCounter::prepare_anchor_dict(ad)) {
                 cerr << "Error: Failed to process dictionary\n";
                 return false;
             }
         }
 
     // Validate input name and technology.
     if (input_name.empty())
     {
         cerr << "No input name provided\n";
         return false;
     }
 
     if (params.technology == technology_t::unknown)
     {
         cerr << "Unknown technology\n";
         return false;
     }
 
     // Process input file and load predefined CBCs if specified.
     ifstream ifs(input_name);
     string s1, s2, s;
 
     while (ifs >> s)
     {
         auto p = find(s.begin(), s.end(), ',');
         if (p == s.end())
         {
             cerr << "Wrong line in input name file: " << s << endl;
             return false;
         }
 
         params.cbc_file_names.emplace_back(s.begin(), p);
         params.read_file_names.emplace_back(p+1, s.end());
     }
 
     if (!params.predefined_cbc_fn.empty())
     {
         if (params.technology == technology_t::visium)
             load_predefined_cbc_visium();
         else if (params.technology == technology_t::ten_x)
             load_predefined_cbc_plain();
     }
 
     return true;
 }
 
 // *********************************************************************************************
 void usage()
 {
     // Displays usage information and available options for the tool.
     cerr
         << "BKC: Counter of k-mers or k-mer pairs in barcoded reads (v." << BKC_VERSION << " [" << BKC_DATE << "])" << endl;
     cerr
         << "Usage:\n"
         << "    bxc [options]\n"
         << "Options - main:\n"
         << "    --mode <single|pair> - single k-mers or pairs of k-mers (default: single)\n"
         << "    --cbc_len <int> - CBC len " << params.cbc_len.str() << endl
         << "    --umi_len <int> - UMI len " << params.umi_len.str() << endl
         << "    --leader_len <int> - leader_len " << params.leader_len.str() << endl
         << "    --follower_len <int> - follower len " << params.follower_len.str() << endl
         << "    --gap_len <int> - gap len " << params.gap_len.str() << endl
         << "    --n_threads <int> - no. threads " << params.no_threads.str() << endl
         << "    --canonical - turn on canonical k-mers (default: false); works only in single mode" << endl
         << "    --verbose <int> - verbosity level " << params.verbosity_level.str() << endl
         << "Options - input:\n"
         << "    -d <file_name> - file with accepted anchors (one k-mer per line or TSV with 'anchor' column)" << endl
         << "    --input_format <fasta|fastq> - input format (default: fastq)\n"
         << "    --input_name <file_name> - file name with list of pairs (comma separated) of barcoded files; 1st contains CBC+UMI\n"
         << "    --technology <10x|visium> - sequencing technology (default: " << technology_str(params.technology) << ")\n"
         << "    --soft_cbc_umi_len_limit <int> - tolerance of CBC+UMI len " << params.soft_cbc_umi_len_limit.str() << endl
         << "    --cbc_filtering_thr <int> - CBC filtering threshold (0 is for auto) " << params.cbc_filtering_thr.str() << endl
         << "    --allow_strange_cbc_umi_reads - use to prevent the application from crashing when the CBC+UMI read length is outside the acceptable range (either shorter than CBC+UMI or longer than CBC+UMI+soft_cbc_umi_len_limit) (default: " << params.allow_strange_cbc_umi_reads << ")\n"
         << "    --apply_cbc_correction - apply CBC correction (default: " << params.apply_cbc_correction << ")\n"
         << "Options - output:\n"
         << "    --output_format <bkc|splash> (default: " << to_string(params.output_format) << ")\n"
         << "    --output_name <file_name> - output file name (default: " << params.out_file_name << ")\n"
         << "    --sample_id <int> - sample id (default: " << params.sample_id << ")\n"
         << "    --n_splits <int> - no. splits " << params.no_splits.str() << endl
         << "    --log_name <file_name> - path to cbc log files (default: " << params.cbc_log_file_name << "); if not provided, log will not be produced\n"
         << "    --filtered_input_path <string> - path to filtered input files (default: " << params.filtered_input_path << ")\n"
         << "    --export_filtered_input_mode <none|first|second|both> - specifies which reads will be outputted (default: " << to_string(params.export_filtered_input) << ")\n"
         << "    --max_count <int> - max. counter value " << params.max_count.str() << endl
         << "    --zstd_level <int> - internal compression level " << params.zstd_level.str() << endl
         << "Options - filtering:\n"
         << "    --predefined_cbc <file_name> - path to file with predefined CBCs (default: " << params.predefined_cbc_fn << ")\n"
         << "    --poly_ACGT_len <int> - all leaders containing polyACGT of this length will be filtered out (0 means no filtering) " << params.poly_ACGT_len.str() << endl
         << "    --artifacts <file_name> - path to artifacts, each leader containing artifact will be filtered out\n"
         << "    --apply_filter_illumina_adapters - if used leaders containing Illumina adapters will be filtered out\n"					
         << "    --leader_sample_counts_threshold <int> - keep only leaders with counts > leader_sample_counts_threshold " << params.rare_leader_thr.str() << endl	// kmer_counts_in_cbc_thr ?
 //		<< "Options - other:\n"
         ;
 }
 
 // *********************************************************************************************
 bool load_predefined_cbc_visium()
 {
     ifstream ifs(params.predefined_cbc_fn);
 
     params.predefined_cbc.clear();
 
     regex re("([ACGT]+)-(.+),([0-9]+),[0-9]+,[0-9]+,[0-9]+,[0-9]+");
     smatch sm;
 
     string s;
 
     while (ifs >> s)
     {
         if (!regex_match(s, sm, re))
         {
             cerr << "Unknown trusted CBC description: " << s << endl;
             return false;
         }
 
         if (stoi(sm[3]) == 1)
             params.predefined_cbc.emplace_back(sm[1]);
     }
 
     return true;
 }
 
 // *********************************************************************************************
 bool load_predefined_cbc_plain()
 {
     ifstream ifs(params.predefined_cbc_fn);
 
     params.predefined_cbc.clear();
 
     string s;
 
     while (ifs >> s)
         params.predefined_cbc.emplace_back(s);
 
     return true;
 }
 
 // *********************************************************************************************
 int main(int argc, char **argv)
 {
    CBarcodedCounter::init();
    
    if (!parse_args(argc, argv))
     {
         usage();
         return 1;
     }
 
     CBarcodedCounter barcoded_counter;
     std::cout << "setting params" << std::endl;
     barcoded_counter.SetParams(params);
     
     std::cout << "processing CBC" << std::endl;
     barcoded_counter.ProcessCBC();
    
     std::cout << "checking for filter" << std::endl;
     if (((uint32_t)params.export_filtered_input) & (uint32_t)export_filtered_input_t::first)
         barcoded_counter.ProcessExportFilteredCBCReads();
 
     if (params.counting_mode == counting_mode_t::filter)
     {
         if (((uint32_t)params.export_filtered_input) & (uint32_t)export_filtered_input_t::second)
             barcoded_counter.ProcessExportFilteredReads();
         
         return 0;
     }
     std::cout << "processing reads" << std::endl;
     barcoded_counter.ProcessReads();
 
     barcoded_counter.ShowTimings();
 
     return 0;
 }
 