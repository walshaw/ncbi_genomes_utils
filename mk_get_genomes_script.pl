#!/usr/bin/perl

use strict;
use warnings;
use autodie;

use Carp;

my $base_local_dir = shift or croak;

# listing file e.g.
# ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
# or
# ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt

my $listing_file = shift or croak;

my $file_ext = shift || croak;

my $field_index = shift || 19;

my %host_count; # should really only be one key

my $file_ext_qm = quotemeta $file_ext;

print STDERR qq{filename extension is : '$file_ext'\n}; #; quotemeta'd: '$file_ext_qm'\n};

open my $fh, '<', $listing_file;

my $n_lines = 0;

while (defined (my $line = <$fh>)) {

    next if $line =~ m{ \A [#] }xms;

    chomp $line;
    my @field = split m{ \t }xms, $line;
    my $url = $field[$field_index];

    croak qq{unexpected URL: '$url'}
        if $url !~ m{ \A ftp[:][/]{2} ([^/]+) [/] (.*) \z }ixms;

    my ($host, $path) = ($1, $2);
    (my $dir_name = $path) =~ s{ \A .* [/] ([^/]+) \z }{$1}xms;
    my $local_dir = qq{$base_local_dir/$dir_name};

    print qq{if [[ ! -d $local_dir ]]; then
    MKDIR="mkdir $local_dir"
    echo \$MKDIR
    eval \$MKDIR
else
    echo "$local_dir already exists"
fi

NCFTPGET="ncftpget $host $local_dir $path/*${file_ext}*"
echo \$NCFTPGET
eval \$NCFTPGET

};
    print STDERR qq{$dir_name\t(path: $path)\n};

    
    $host_count{$host}++;
    $n_lines++;
}

close $fh;

print STDERR qq{read $n_lines lines from $listing_file\n};


