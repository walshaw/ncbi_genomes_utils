#!/usr/bin/perl

use strict;
use warnings;
use autodie;

use Carp;

my $usage = qq{$0 COMPR_DIR UNCOMPR_DIR LIST_FILE FILE_EXT [ FIELD_INDEX [ COMPR_EXT ] ]

	COMPR_DIR	parent dir containing per-genome dirs containing compressed files
	UNCOMPR_DIR	parent dir containing per-genome dirs to contain uncompressed files
	LIST_FILE	e.g. ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
	FILE_EXT	e.g. '.faa' or '.fna' or '.gff'
	FIELD_INDEX	refers to tab-delimited columns (first is 0) in LIST_FILE;
			must specify the ftp path column; default is 19
	COMPR_EXT	compressed file extension; default is '.gz'

};

my $base_local_compr_dir = shift or croak $usage;
my $base_local_unc_dir   = shift or croak $usage;

# listing file e.g.
# ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
# or
# ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt

my $listing_file  = shift or croak $usage;

my $file_ext      = shift || croak $usage;

my $field_index   = shift || 19;

my $cmpr_file_ext = shift || '.gz';

my $file_ext_qm      = quotemeta $file_ext;
my $cmpr_file_ext_qm = quotemeta $cmpr_file_ext;

print STDERR qq{filename extension is : '$file_ext'\n}; #; quotemeta'd: '$file_ext_qm'\n};

open my $fh, '<', $listing_file;

my $n_genome_dirs         = 0; # input, i.e. compressed
my $n_missing_genome_dirs = 0; # ditto
my $n_processed_genome_dirs = 0; # no. of input dirs processed = no. of output
                                 # dirs 'to be created'
###my $n_present_genome_dirs   = 0; # output, i.e. already present; had already been created

my $n_genome_files        = 0; # total 'to be created'; same as no. of input files
###my $n_present_genome_files = 0; # number that were already present; had already been created
###my $n_created_genome_files = 0;

# Note that the concept is to run this script once, to produce bash
# script that could well be run more than once; for example if one
# run of the bash script was interrupted, it could be re-run and it
# checks to see what's already been done, and so skips on to the
# next yet to be done tasks (N.B. if interrupted while uncompressing
# one file, then the partial uncompressed file should be deleted
# before rerunning - otherwise it will be left as-is).

print qq{N_PRESENT_OUTPUT_DIRS=0
N_CREATED_OUTPUT_DIRS=0
N_PRESENT_OUTPUT_FILES=0
N_CREATED_OUTPUT_FILES=0

};

LINE:
while (defined (my $line = <$fh>)) {

    next if $line =~ m{ \A [#] }xms;

    chomp $line;
    my @field = split m{ \t }xms, $line;
    my $url = $field[$field_index];

    croak qq{unexpected URL: '$url'}
        if $url !~ m{ \A ftp[:][/]{2} ([^/]+) [/] (.*) \z }ixms;

    my ($host, $path) = ($1, $2);
    (my $dir_name = $path)     =~ s{ \A .* [/] ([^/]+) \z }{$1}xms;
    my $local_compressed_dir   = qq{$base_local_compr_dir/$dir_name};
    my $local_uncompressed_dir = qq{$base_local_unc_dir/$dir_name};

    $n_genome_dirs++;

    if (! -d $local_compressed_dir) {
        print STDERR qq{WARNING: no local dir $local_compressed_dir\n};
        $n_missing_genome_dirs++;
        next LINE;
    }

    $n_processed_genome_dirs++;

    opendir(my $dh, $local_compressed_dir);
    my @compressed_files = grep { m{ $file_ext_qm $cmpr_file_ext_qm \z }xms } readdir $dh ;
    closedir $dh;

    print qq{if [[ ! -d $local_uncompressed_dir ]]; then
    MKDIR="mkdir $local_uncompressed_dir"
    echo \$MKDIR
    eval \$MKDIR
    N_CREATED_OUTPUT_DIRS=\$(( \$N_CREATED_OUTPUT_DIRS + 1 ));
else
    echo "$local_uncompressed_dir already exists"
    N_PRESENT_OUTPUT_DIRS=\$(( \$N_PRESENT_OUTPUT_DIRS + 1 ));

fi

};

    for my $compr_file (@compressed_files) {

        $n_genome_files++;
        my $local_cmpr_path   = qq{$local_compressed_dir/$compr_file};
        (my $local_uncmpr_path = qq{$local_uncompressed_dir/$compr_file})
            =~ s{ $cmpr_file_ext_qm \z }{}xms;
        print qq{if [[ -f $local_uncmpr_path ]]; then
    echo "$local_uncmpr_path already exists"
    N_PRESENT_OUTPUT_FILES=\$(( \$N_PRESENT_OUTPUT_FILES + 1 ))
else
    GUNZIP='gunzip -c $local_cmpr_path > $local_uncmpr_path'
    echo \$GUNZIP
    eval \$GUNZIP
    N_CREATED_OUTPUT_FILES=\$(( \$N_CREATED_OUTPUT_FILES + 1 ))
fi
\n};

    }

    print qq{\n\n};

}

close $fh;

print qq{
echo "\$N_PRESENT_OUTPUT_DIRS output (uncompressed) dirs were already present"
echo "\$N_CREATED_OUTPUT_DIRS output (uncompressed) dirs were created"
echo "\$N_PRESENT_OUTPUT_FILES output (uncompressed) files were already present"
echo "\$N_CREATED_OUTPUT_FILES output (uncompressed) files were created"

};

print STDERR qq{Read $n_genome_dirs genome lines from $listing_file
locally, $n_missing_genome_dirs expected genome dirs were missing and
$n_processed_genome_dirs were present and processed;
those $n_processed_genome_dirs dirs contained $n_genome_files ${file_ext}$cmpr_file_ext files
which were uncompressed.\n};

