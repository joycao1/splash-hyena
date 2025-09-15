./bin/bkc_filter \
  --mode pair \
  --input_name fl.txt \
  -d 8mer_list_for_single_cell_testing.txt \
  --cbc_len 16 \
  --umi_len 12 \
  --leader_len 8 \
  --follower_len 31 \
  --gap_len 0 \
  --verbose 2 \
  --output_name ERR13720422_pair_1.bkc

  ./bin/bkc_dump --input_name ERR13720422_pair_1.bkc --output_name ERR13720422_pair_1.txt
