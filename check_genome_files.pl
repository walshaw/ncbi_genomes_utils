#!/usr/bin/perl

use strict;
use warnings;
use Carp;

use List::Util qw( min sum );
use Bio::SeqIO;

my $root_dir = shift || '.';

my $reqd_ext = shift || '.faa';

croak qq{$root_dir does not exist or is not a directory} if (! -d $root_dir);

my $seq_len_interval = 100;
my $file_seqs_interval = 10;

opendir my $dh, $root_dir or croak qq{couldn't get listing of directory $root_dir};

my @dirs = grep { (! m{ \A [.] {1,2} \z }xms) && (-d qq{$root_dir/$_} ) } readdir $dh;

closedir $dh;

##my @total_seq_len_frq; # indices represent sequence lengths; values are frequencies of files of those lengths

my @file_count_freq;

for my $genome_dir (sort @dirs) {

    my $genome_dir_path = qq{$root_dir/$genome_dir};

    opendir my $dh, $genome_dir_path or croak qq{couldn't get listing of directory $genome_dir_path};
###print qq{reading '$genome_dir_path'\n};
    my @files = grep { (-f qq{$genome_dir_path/$_}) } (readdir $dh);
#    my @files = readdir $dh;

    closedir $dh;
###print qq{no. of files is }, scalar @files, qq{\n};
    my %file_tally; # keys are filename extensions, e.g. faa, gff
    my $n_zero_byte_files = 0;
    my $n_files_in_dir    = 0;
    my @reqdext_files;
    my @dir_seq_len_freq;   # indices represent seq lengths; values are frequencies of sequences of that length (irrespective
                            # of how those sequences are distributed amongst different .faa files)
    my @dir_seq_count_freq; # indices represent numbers of seqs in a .faa file); values are frequencies of files with that number
    my $n_seqs_in_dir = 0;
    my $total_aa_in_dir = 0;

    for my $file_nopath (@files) {
        my $file = qq{$genome_dir_path/$file_nopath};
###print qq{TESTING FILE '$file'\n\n};
        my @fstrings = split m{ [.] }xms, $file_nopath;

        my $ext = $fstrings[-1];

        $file_tally{$ext}++;
        $n_files_in_dir++;

        my @stat = stat $file;

        my $size = $stat[7];
###print qq{\t\t$file\t$size\n};
        if (!$size) { $n_zero_byte_files++; }

        if ($ext eq $reqd_ext) {
            push @reqdext_files, $file;
            my $seqio = Bio::SeqIO->new(-file => $file, -format => 'Fasta');
            my $n_seqs_in_file;
            while (my $seq = $seqio->next_seq()) {
                $n_seqs_in_file++;
                my $seq_length = $seq->length();
                $dir_seq_len_freq[$seq_length]++;
                $total_aa_in_dir += $seq_length;
            }
            $dir_seq_count_freq[$n_seqs_in_file]++;
            $n_seqs_in_dir += $n_seqs_in_file;
        }
    }

    my $n_reqdext_files = $file_tally{$reqd_ext} || 0;

    my $mean_seqs_per_file;
    if ($n_reqdext_files) {
        $mean_seqs_per_file = $n_seqs_in_dir   / $n_reqdext_files;
    }
    my $mean_aa_per_seq;
    if ($n_seqs_in_dir) {
       $mean_aa_per_seq = $total_aa_in_dir / $n_seqs_in_dir;
    }

    print qq{\n\n$genome_dir\t$n_reqdext_files $reqd_ext file},
              ($n_reqdext_files       == 1) ? qq{} : qq{s},
              qq{;\t$n_files_in_dir file},
              ($n_files_in_dir    == 1) ? qq{} : qq{s},
              qq{ in all;\t$n_zero_byte_files zero-byte file},
              ($n_zero_byte_files == 1) ? qq{} : qq{s},
              ($n_zero_byte_files)      ? qq{ WARNING} : qq{}, qq{\n};

    if ($n_seqs_in_dir) { printf qq{\ttotal %4d sequences; total %7d aa; mean %6.1f sequences per file; mean %7.1f aa per sequence\n},
                                 $n_seqs_in_dir, $total_aa_in_dir, $mean_seqs_per_file, $mean_aa_per_seq; }

    if ($n_reqdext_files) { print qq{\n\tsequences per file:\n}; }
    if ($dir_seq_count_freq[0]) {
        printf qq{\t%3d file(s) have 0 sequences\n}, $dir_seq_count_freq[0];
    }
    for (my $interval_start = 1; $interval_start <= $#dir_seq_count_freq; $interval_start += $file_seqs_interval) {
        my $max_c = $interval_start + $file_seqs_interval - 1;
        $max_c = min $max_c, $#dir_seq_count_freq;
        my $interval_total_files = 0;
        for my $c ( $interval_start .. $max_c) {
            if ($dir_seq_count_freq[$c]) { $interval_total_files += $dir_seq_count_freq[$c]; }
        }
        next if (!$interval_total_files);
        printf qq{\t%3d file(s) have %3d - %3d sequences\n}, $interval_total_files, $interval_start, $max_c;
    }

    if ($n_seqs_in_dir) { print qq{\n\taa per sequence:\n}; }
    if ($dir_seq_len_freq[0]) {
        printf qq{\t%3d sequence(s) have 0 aa\n}, $dir_seq_len_freq[0];
    }
    for (my $interval_start = 1; $interval_start <= $#dir_seq_len_freq; $interval_start += $seq_len_interval) {
        my $max_l = $interval_start + $seq_len_interval - 1;
        $max_l = min $max_l, $#dir_seq_len_freq;
        my $interval_total_seqs = 0;
        for my $l ( $interval_start .. $max_l ) {
            if ($dir_seq_len_freq[$l]) { $interval_total_seqs += $dir_seq_len_freq[$l]; }
        }
        next if (!$interval_total_seqs);
        printf qq{\t%4d sequence(s) have %4d - %4d aa\n}, $interval_total_seqs, $interval_start, $max_l;

    }

    $file_count_freq[$n_reqdext_files]++;

}

print qq{\n\nAll }, scalar @dirs, qq{ genome dirs checked.\n\n};

for my $f (0 .. $#file_count_freq) {

    if ($file_count_freq[$f]) {
        printf qq{\t%4d dirs have %4d $reqd_ext file%s\n}, $file_count_freq[$f], $f, ($f == 1) ? qq{} : qq{s};
    }
}

