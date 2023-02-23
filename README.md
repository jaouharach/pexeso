# PEXESO: Efficient Joinable Table Discovery in Data Lakes: A High-Dimensional Similarity-Based Approach (Dong Yuyang, Takeoka Kunihiro et al.)
Paper: Efficient Joinable Table Discovery in Data Lakes: A High-Dimensional Similarity-Based Approach  

URL link for the PEXESO's paper (The original implementation was in Python):\
> https://arxiv.org/abs/2010.13273


This C version was developed by: Jaouhara Chanchaf on 02/22/2023
Copyright 2023 Mohammed 6 Polytechnic University. All rights reserved.\

> (!) Make sure that the c compiler support function isnan() in package math.h.
> To compile, go to code-directory/src and run:\
`gcc globals.c select_pivots.c gsl_matrix.c inv_index.c match_map.c  query_engine.c level.c file_buffer_manager.c file_buffer.c file_loader.c cell.c pexeso.c hgrid.c  main.c -o ../bin/pexeso -lm -g -lgsl -lgslcblas -lmcheck`


#### DATA PREPARATION
1. Download and extract tables from :
> https://data.dws.informatik.uni-mannheim.de/webtables/2015-07/englishCorpus/compressed/

2. Clean data: \
To remove numerical data go to data-preparation/clean and run:\
`python3 clean.py [directory_of_raw_tables] [directory_for_clean_tables]`

3. Generate embeddings: \
    i. download fasttext binary file from:\
        > https://fasttext.cc/docs/en/crawl-vectors.html \
    ii. then go to data-preparation/embed/fasttext.py and update the following line with the correct path to fasttext file, if you want to normalize vectors set normalize_vectors to True:\
        `FASTTEXT_PATH = '/home/jchanchaf/local/src/fasttext-en/cc.en.300.bin`\
        `normalize_vectors = False` \
    iii. Finally, run:\
        `fasttext_embed.py [directory_of_clean_tables] [directory_to_store_embbedings] [path_to_fasttext_binary] [embeddings_length]`\

> Note: you can use the same workflow in (3.) to generate embbeddings using glove (https://nlp.stanford.edu/data/glove.6B.zip).

4. Generate Binary files:\
Go to data-preparation/bin and run: \
`tobin.py [directory_of_table_embeddings] [directory_to_store_binary_files] [embeddings_length]`

#### EXPERIMENTS
1. Build the index:\
`bin/pexeso --work-dir  $HOME/Projects/pexeso/grids/100/ --dataset /data/real/jchanchaf/wdc-2015-en-full/subsets/100/ --dataset-size 100 --metric-space-dim 50 --pivot-space-dim 3 --num-levels 4 --leaf-size 1500 --metric-buffer-size 60 --qgrid-metric-buffer-size 5 --index-path $HOME/Projects/pexeso/grids/100/ --mode 0 --fft-scale 7 --best-fft 0`

2. Query the index:\
`bin/pexeso --work-dir  $HOME/Projects/pexeso/grids/100/ --dataset /data/real/jchanchaf/wdc-2015-en-full/subsets/50k/ --dataset-size 100 --queries /data/real/jchanchaf/wdc-2015-en-full/10m/query/10q-size100/  --num-query-sets 1 --min-query-set-size 100 --max-query-set-size 100 --metric-space-dim 50 --pivot-space-dim 3 --num-levels 4 --index-path $HOME/Projects/pexeso/grids/100/Hgrid/ --mode 1 --join-threshold 0.1 --dist-threshold 0.06`

3. Parameters:\
* `--work-dir [string]` : work directory.
* `--dataset [string]` : directory where data files are stored.--dataset-size
* `--dataset-size [x]` : nb of files (tables) in the dataset.
* `--queries [string]` : directory where query files are stored.
* `--num-query-sets [x]` : number of query column to execute.
* `--max-query-set-size [x]` : Maximum query column size.
* `--min-query-set-size [x]` : Minumum query column size.
* `--metric-space-dim [x]` : embeddings length.
* `--pivot-space-dim [x]` : number of pivots (|P|).
* `--num-levels [x]` : number of levels (m).
* `--leaf-size [x]` : max number of vectors that can be stored in a leaf cell (a cell at the bottom level).
* `--metric-buffer-size [x]` : memory size to store buffered (metric space) vectors in MBs.
* `--index-path [string]` : directory where pexeso grid is stored.
* `--mode [0/1/2]` : 0 for index building , 1 for querying, 2 for building and quering the index.
* `--fft-scale [x]` : a constant used to perform pivot selection the higher the value the longer it takes to build the index (the recommended value is in  [7 - 30]).
* `--best-fft [0/1]` : set to 1 if you want to search for the best pivots using a range of fft_scale values (test values in [fft_scale - max-fft]), if  `--best-fft` is set to one `--max-fft [x]` must be defined.
* `--max-fft [x]` : max fft_scale value to test (default value is 30).
* `--join-threshold [x]` : T.
* `--dist-threshold [x]` : tau.

