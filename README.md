# ncbi_genomes_utils

Various utilities to ease the specification, downloading and verification of NCBI Genome data files, for multiple genomes.

## mk_get_genomes_script.pl

Usage summary:

>mk_get_genomes_script.pl <local_dir> <list_file> <file_extension> [ <file_index> ] > get_my_genomes.bash

### Background

The officially-recommended approach to obtaining data files of multiple genomes is to automate the process with scripts which
run rsync (or a web or ftp client, such as wget, curl, ncftpget, etc) (https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#protocols).

This script is a 'meta-script' - its purpose is not to do the above, but to generate a script which will in turn do the above.
The latter script should then be run multiple times (unless the number of files to download is small). When it is run the first
time, the practical reality is that it is unlikely that all files will download completely. Therefore there should be at least one
repeat run (after the first one has finished and exited).

The scripts created by mk_get_genomes_script.pl use ncftpget to download the files. rsync would probably be preferable, but like
rsync, ncftpget is efficient in terms of checking first whether there is a complete local copy of the remote file, before
attempting to download it; if only an incomplete copy is present, then ncftpget is efficient in terms of partially downloading the
remote copy of that file instead of starting from scratch.

### Input to mk_get_genomes_script.pl

Input is the contents of the file named in the second argument (see above). This file should ideally be in "assembly summary" format as
used by NCBI GenBank and NCBI RefSeq to specify genome data. Some example files:

* ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt

* ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt

The default value of <file_index> (optional fourth argument, see above) is 19, meaning that field #19 (the 20th field) contains the URLs of each remote directory
(containing the data files for one genome).

Note that the <file_extension> argument is intentionally mandatory, because there are various extensions which are likely to be of
primary interest (e.g. .fna, .faa, .gff) and it would not be adviseable to download all files by default (as that would be a large
amount of data). If multiple types of file are needed, then the mk_get_genomes_script.pl script should be run once for each; it is
wise to keep the downloads for each file extension separate, for clarity and ease of management.

One approach to focussing on particular organisms is to download the whole assembly summary file (e.g. for all bacteria) and then
produce a smaller version, filtered by organism name using e.g. `grep`. Other useful screens can be done using other columns, e.g.
the degree of completeness of the assembly, the ASM name, the strain name, etc.

### Output of mk_get_genomes_script.pl

Principal output is sent to standard output, and this should be redirected into a new (Bash script) file with a suitable name.
When that script is run in turn, then each genome-specific directory in the input list file will be checked locally to determine if it
exists, and it will be created if not. Then, those remote data files which match the filename extension (third argument) will be
downloaded, into a local directory of the same name as the remote one (i.e., specific to a single genome); the local parent of
this directory will be that specified by the first argument, <local_dir>.

Additional output from mk_get_genomes_script.pl is sent to standard error, consisting of brief commentary about every directory
in the input list file.


