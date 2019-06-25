#!/usr/bin/perl

use strict;
use warnings;
use autodie;
use Carp;

use Getopt::Long;
use Data::Dumper;
use List::Util qw( first );
use List::MoreUtils qw( none uniq );
use File::List;
use File::Basename;

my $default_parent_dir = q{.};
# only one of $parent_dir or $file_path_list can be specified by the user
my $parent_dir;
my $file_path_list;               # file containing list of paths
my $file_ext           = q{.gff};
my $output_dest; # a dir name or 'STDOUT'; if undefined, then same as 'STDOUT'
my $verbosity          = 0;
my $max_n_buffer_lines = 100000;
my $outfile_key_length =      4;
my $tmp_basedir        = q{/tmp};

my $command = qq{$0 }.join(q{ },@ARGV);

my $usage = qq{Usage:

    $0 [OPTIONS]

Options:

    -parentdir	STRING	name of dir representing top of tree in which to
                        search for files; default '$default_parent_dir' ;
                        see -filelist which is an alternative, and mutually
                        exclusive

    -filelist   STRING  name of file containing paths of all files to input;
                        it's an alternative to -parentdir, and avoids having
                        to do the 'find' operation which -parentdir
                        necessitates, and which can take a long time;
                        -filelist is also more suitable for parallelisation,
                        i.e. multiple invocations of the script with different
                        sub-lists as input; because the order of the input
                        files may not be guaranteed to be identical each time;
                        there is no default value

    -ext	STRING	filename extension specifying files to search for;
                        default '$file_ext'

    -output	STRING	name of directory in which to create output files;
    			special value is 'STDOUT', which is also the default;
			in which case, no files are written, and all output
			is to STDOUT

    -buffer	INTEGER	maximum number of lines in buffer, before the next
    			write to output file(s) (or STDOUT); default $max_n_buffer_lines

    -keylength	INTEGER	number of characters to use in output filenames;
    			the first N characters of the protein IDs are
			used as the filenames, and all lines with matching
			protein IDs (first field) are written to the
			corresponding file; default $outfile_key_length ;
			in the few cases where the protein ID may be
			missing (will be replaced by '[undefined]') then
			the output file is named '_undefined_'

    -tmpdir	STRING	name of temporary dir used to do some file
    			manipulations; for example, if two or more GFF files
			are encountered with the same name, in different parts
			of the traversed tree, then these will be compared and
			if non-identical, a more detailed comparison is done,
			which involves making temporary new versions of them
			(e.g. when known 'null' strings are stripped out, the
			resulting files are often the same)

    -verbosity	INTEGER	verbosity level of commentary sent to STDERR

};

my %optctl = (

    "parentdir=s" => \$parent_dir,
    "filelist=s"  => \$file_path_list,
    "ext=s"       => \$file_ext,
    "output=s"    => \$output_dest,
    "buffer=i"    => \$max_n_buffer_lines,
    "keylength=i" => \$outfile_key_length,
    "tmpdir=s"    => \$tmp_basedir,
    "verbosity=i" => \$verbosity,

);

croak qq{$usage} if !(&GetOptions(%optctl));

$parent_dir         ||= shift;
$file_ext           ||= shift;
$output_dest        ||= shift;
$max_n_buffer_lines ||= shift;
$outfile_key_length ||= shift;
$tmp_basedir        ||= shift;
if (!defined $verbosity) { # as 0 is a legitimate value
    $verbosity = shift || 0;
}

my $parent_dir_or_list_file;

if ($parent_dir) {
    croak qq{specify either -parentdir or -filelist but not both}
        if $file_path_list;
    $parent_dir =~ s{ [\/]+ \z }{}xms;
    croak qq{$parent_dir does not exist or is not a directory}
      if ! -d $parent_dir;
    $parent_dir_or_list_file = $parent_dir;
}
elsif ($file_path_list) {
    croak qq{$file_path_list does not exist or is not a file}
        if ! -f $file_path_list;
    $parent_dir_or_list_file = $file_path_list;
}
else {
    $parent_dir = $default_parent_dir;
    $parent_dir_or_list_file = $parent_dir;
}

print STDERR qq{$command},
    ($parent_dir) ? qq{

    Searching for GFF format files with filename extension $file_ext,
    in directory tree '$parent_dir/'} : qq{

    Reading paths of GFF format files with filename extension $file_ext
    from file '$file_path_list'},
    qq{ . Buffer of protein_id-keyed lines
    has a maximum size of $max_n_buffer_lines; lines will be written
    and the buffer emptied whenever the buffer size is exceeded.

};


if ((!defined $output_dest) || (uc $output_dest eq 'STDOUT')) {
    print STDERR qq{    Writing protein_id - keyed lines to STDOUT\n};
    $output_dest = \*STDOUT;
}
else {
    croak qq{'$output_dest' does not exist or is not a directory}
      if (!-d $output_dest);
    print STDERR qq{    Writing protein_id-keyed lines to files in directory
                        $output_dest ; names of files are the first $outfile_key_length
			characters of the protein_ids of the lines they contain.\n};
}

print STDERR qq{    Verbosity level of commentary is $verbosity\n\n};

process_gff_dirs($parent_dir_or_list_file, $file_ext, $output_dest,
                 $max_n_buffer_lines, $outfile_key_length,
		 $verbosity);


exit 0;

sub process_gff_dirs {

    my $parent_dir_or_list_file    = shift ||    q{.};
    my $file_ext           = shift || q{.gff};
    my $output_dest        = shift ||    q{.};
    my $max_n_buffer_lines = shift ||  100000;
    my $outfile_key_length = shift ||       4; # 4 characters, representing first 4 chars
                                     # of the protein IDs, and also the names
				     # of the output files
    my $verbosity          = shift ||       0;

    my $n_processed_lines  =                0;
    my $n_processed_files  =                0;
    my $n_ignored_files    =                0;

    my @output_buffer;

    my $file_ext_re  = qr{$file_ext};
    my $file_ext_end = $file_ext_re . q{$};

###    print qq{$file_ext_re\n};

    my @gff_files;

    if ( -f $parent_dir_or_list_file ) {
        my $list_file = $parent_dir_or_list_file;
        open my $fh, '<', $list_file;
        @gff_files = grep { m{ $file_ext_end }xms } (<$fh>);
        @gff_files = map { chomp; $_ } @gff_files;
        close $fh;
    }
    else {

        (my $parent_dir = $parent_dir_or_list_file) =~ s{ [/] \z }{}xms;

        croak qq{Directory $parent_dir does not exist or is not a directory}
            if (! -d $parent_dir);

        my $file_list = File::List->new($parent_dir);
        @gff_files = @{ $file_list->find($file_ext_end) };
    }

    # To ensure that the same file isn't parsed twice or more, and
    # to examine the file/parent dir if it is; e.g. the parent
    # directory could have been encountered twice, once as the
    # original dir and once as a symlink.

    my %location_of_gff_file;
    my %inode_of_gff_file;


    print STDERR scalar @gff_files, qq{ GFF files to process\n\n};

    ##print STDERR join(qq{\n},@gff_files), qq{\n\n};

    GFF_FILE:
    for my $gff_file_path (@gff_files) {

        my $gff_file  = basename($gff_file_path);
	my $gff_dir   = dirname($gff_file_path);
        my $gff_inode = (stat $gff_file_path)[1];
	if ($verbosity > 3) {
	    print STDERR qq{$gff_file is in $gff_dir\n};
        }

        if ($location_of_gff_file{$gff_file}) {
	    my $existing_path = $location_of_gff_file{$gff_file}.qq{\/$gff_file};
	    print STDERR qq{$gff_file already encountered:\n\t};
	    if ($gff_inode == $inode_of_gff_file{$gff_file}) {
	        print STDERR
		  qq{$gff_file_path is the same file as\n\t$existing_path (same inode)\n};
		$n_ignored_files++;
		next GFF_FILE;
	    }
	    else {
	        print STDERR qq{$gff_file_path is different to\n\t$existing_path\n};
		my $diff_command = qq{diff --brief $existing_path $gff_file_path | wc -l};
#		open my $fh, '-|', $diff_command or print STDERR qq{couldn't do: $diff_command\n};
#		my $diff_result = <$fh>;
#		close $fh or croak STDERR qq{couldn't complete: $diff_command\n};
                my $diff_result = `$diff_command`;
		chomp $diff_result;
		print STDERR qq{\tcomparing content of the two files : },
		    ($diff_result) ? qq{DIFFERENT CONTENT - }
		                   : qq{identical - },
		    qq{IGNORING $gff_file_path\n\n};
		if ($diff_result) {
		    my $strip_string = qq{RefSeq\tSTS};
		    my $diff2_result
		      = check_stripped_content($gff_file_path, $existing_path, $strip_string);
		    print STDERR qq{\tContent of both $gff_file files has been compared: },
		      ($diff2_result)
		        ? qq{could not detect a simple reason for differences; PLEASE CHECK}
			: qq{content is identical when all lines containing '$strip_string' are ignored}
			, qq{\n\n};
                    if ($diff2_result && ($verbosity > 2)) {
		        my $diff3_command = qq{diff $existing_path $gff_file_path};
		        print STDERR qq{$diff3_command\n}, `$diff3_command`, qq{\n\n};
		    }
		}
		$n_ignored_files++;
		next GFF_FILE;
	    }
	}
	else {

	    $location_of_gff_file{$gff_file} = $gff_dir;
	    $inode_of_gff_file{$gff_file} = $gff_inode;

	    my @output_lines = process_gff_file($gff_file_path,$verbosity);
	    push @output_buffer, @output_lines;
	    $n_processed_lines += scalar @output_lines;

	    if (@output_buffer >= $max_n_buffer_lines) {
	        print STDERR qq{$n_processed_lines lines in $n_processed_files files processed }
		             .qq{($n_ignored_files files ignored)\n};
	        write_buffer(\@output_buffer, $outfile_key_length, $output_dest);
		@output_buffer = ();
	    }
	}

        $n_processed_files++;
    }

   # Write the remaining lines from the buffer
   if (@output_buffer) {
       print STDERR qq{$n_processed_lines lines in $n_processed_files files processed }
		    .qq{($n_ignored_files files ignored)\n};
       write_buffer(\@output_buffer, $outfile_key_length, $output_dest);
       @output_buffer = ();

   }

}

sub check_stripped_content {

    # Compares the stripped-down content of two files whose original content
    # is known to be different. The stripped-down versions are the result
    # of removing all lines which contain any of the strings (not regexes)
    # specified in the remove_strings list.

    my $file1          = shift;
    my $file2          = shift;
    my @remove_strings = @_;

    my $tmp_dir = qq{$tmp_basedir/tmp_gp$$};
    my @tmp_files;

    for my $file ($file1, $file2) {
        croak qq{file '$file1' does not exist or is not a file}
	  if (! -f $file);
    }

    mkdir $tmp_dir or croak qq{couldn't create temporary directory $tmp_dir};

    my $i = 0;
    for my $file ($file1, $file2) {
        $i++;
        my $tmp_file = qq{$tmp_dir/}. basename($file1) .qq{.$i.tmp};
	push @tmp_files, $tmp_file;
	# probably easier doing it the long-winded way than trying
	# to do an egrep on multiple strings (which potentially will
	# contain characters which need to be escaped)
        open my $fh,  '<', $file     or croak qq{couldn't read file $file};
	my @lines = <$fh>;
	close $fh  or croak qq{couldn't close file $file};

        my @screened_lines = grep { my $line = $_;
	                            none { my $re = qr($_); $line =~ m{ $re }xms } @remove_strings
				  } @lines;

	open my $ofh, '>', $tmp_file or croak qq{couldn't create temporary file $tmp_file};
	print $ofh @screened_lines;
	close $ofh or croak qq{couldn't close temporary file $tmp_file};
    }
    my $diff_command = q{diff --brief }. join(q{ }, @tmp_files) . q{ | wc -l};
    ###print STDERR qq{$diff_command\n};
    my $diff_result = `$diff_command`;
    chomp $diff_result; # otherwise the returned value, even if zero, is the STRING "0\n"
                        # - which evaluates as true
    unlink @tmp_files or croak qq{couldn't delete temporary files }.join(q{,},@tmp_files);
    rmdir $tmp_dir        or croak qq{couldn't delete temporary dir $tmp_dir};

    return $diff_result;
}

sub process_gff_file {

    my $gff_file = shift;
    my $verbosity = shift;

    my $null_value = '[undefined]';

    my ($feature_list_href, $n_lines) = read_gff_file($gff_file, $verbosity);

    my %feature_list  = %{$feature_list_href};

    my @cds_features;
    my @gene_features;
    my @mrna_features;

    my @output_lines;

    for my $feat_type (qw( cds gene )) {
        if (!defined $feature_list{$feat_type}) {
            print STDERR qq{!!! WARNING: no $feat_type features were obtained from $gff_file\n\n};
        }
    }

    if (defined $feature_list{cds}) {
        @cds_features  = @{$feature_list{cds}}; # list of hrefs
    }
    if (defined $feature_list{gene}) {
        @gene_features = @{$feature_list{gene}}; # list of hrefs
    }
    if (defined $feature_list{mrna}) {
        @mrna_features = @{$feature_list{mrna}}; # list of hrefs
    }

    if ($verbosity) {
        print STDERR qq{$gff_file\n\t},
	  scalar @cds_features , qq{ CDS, },
	  scalar @gene_features, qq{ gene, },
	  scalar @mrna_features, qq{ mRNA features; $n_lines lines\n\n},
    
    }

    if ($verbosity > 4) { print STDERR qq{sorting\n}; }

#    my @sorted_cds_features = sort { $a->{protein_id} cmp $b->{protein_id}  } @cds_features;
    # The PRIMARY value will be the same as protein_id if defined, but
    # otherwise '[undefined]' (in which case, non-unique)
    my @sorted_cds_features = sort { $a->{PRIMARY} cmp $b->{PRIMARY}  } @cds_features;

    CDS:
    for my $cds_feature_href (@sorted_cds_features) {
###print STDERR Data::Dumper->Dump([$cds_feature_href]);
        my $cds_parent;
	my $parent_gene_href;

        my $cds_id = $cds_feature_href->{ID};
	if (!defined $cds_id) {
	    print STDERR qq{CDS feature without an ID!\n};
	    next CDS;
	}

        for my $attrib ('Name', 'product', 'protein_id',
	                'sequenceID',  'begin', 'end') {
	    if (!defined $cds_feature_href->{$attrib}) {
	        print STDERR qq{No $attrib defined for CDS feature $cds_id\n};
            }
	}

        my $cds_name        = $cds_feature_href->{Name}        || $null_value;
        my $cds_sequence_id = $cds_feature_href->{sequenceID}  || $null_value;
        my $cds_begin       = $cds_feature_href->{begin}       || $null_value;
        my $cds_end         = $cds_feature_href->{end}         || $null_value;
        my $cds_product     = $cds_feature_href->{product}     || $null_value;
        my $cds_protein_id  = $cds_feature_href->{protein_id}  || $null_value;

        if (defined $cds_feature_href->{Parent}) {
	    $cds_parent = $cds_feature_href->{Parent};
	    $parent_gene_href = first { $_->{ID} eq $cds_parent } @gene_features;

	    if (!defined $parent_gene_href) {

	        print STDERR
		    qq{CDS feature ($cds_id ("$cds_name"): could not find Parent gene ($cds_parent)\n};

		# In some cases the gene is the 'grandparent'; as the Parent of
		# the CDS is an mRNA, and the parent of the mRNA is the gene

		my $parent_mrna_href = first { $_->{ID} eq $cds_parent } @mrna_features;

		if (defined $parent_mrna_href) {
		    print STDERR qq{Found parent mRNA ($cds_parent)\n};

                    if (defined $parent_mrna_href->{Parent}) {

		        my $mrna_parent = $parent_mrna_href->{Parent};

		        print STDERR
			    qq{Parent ($cds_parent) is mRNA; which has Parent $mrna_parent which };

		        $parent_gene_href = first { $_->{ID} eq $mrna_parent } @gene_features;

		        print STDERR qq{}, (defined $parent_gene_href)
		              ? qq{is a gene; will use grandparent $mrna_parent as the Parent of $cds_id\n\n}
			      : qq{could not be found}
			      , qq{\n\n};
		    }
		    else {
		        print STDERR qq{mRNA $cds_parent has no defined Parent\n};
		    }
		    if (!defined $parent_gene_href) {
		        print STDERR qq{Failed to find gene grandparent of cds $cds_id\n\n};
		    }
		}
	    }
	}
	else {
	    print STDERR qq{CDS feature ($cds_id ("$cds_name") has no defined Parent\n};
	}
	my $gene_id        = $parent_gene_href->{ID}        || $null_value;
	my $gene_name      = $parent_gene_href->{Name}      || $null_value;
	my $gene_begin     = $parent_gene_href->{begin}     || $null_value;
	my $gene_end       = $parent_gene_href->{end}       || $null_value;
	my $gene_strand    = $parent_gene_href->{strand}    || $null_value;
	my $gene_locus_tag = $parent_gene_href->{locus_tag} || $null_value;

        # Note the name (includes relative path) of the GFF file, written as
	# the last field; this enables parsers of the output to (for example)
	# easily read neighbouring gene data from the file quickly, without
	# having to determine what the original GFF file was.

	my $output_line = join(qq{\t}, $cds_protein_id, $cds_sequence_id, $cds_begin,
	                   $cds_end       , $gene_id        , $gene_name     ,
			   $gene_begin    , $gene_end       , $gene_strand   ,
			   $cds_product   , $cds_id         , $gene_locus_tag,
			   $gff_file      ,)
          .qq{\n};

        push @output_lines, $output_line;
    }

###print Data::Dumper->Dump(\@sorted_cd_features);


#proteinID	sequenceID	begin	end	parentID	parentName	parentBegin	parentEnd	product

    return @output_lines;
}



sub read_gff_file {

    my $file      = shift or croak qq{supply filename};
    my $verbosity = shift || 0;

    croak qq{file '$file' doesn't exist or is not a file}
        if (! -e $file);

    open my $fh, '<', $file or croak qq{couldn't read file '$file'};

    my ($feature_list_href, $n_lines) = read_gff_fh($fh, $verbosity);

    close $fh or croak qq{couldn't close file '$file'};

    return ($feature_list_href, $n_lines);
}

sub read_gff_fh {

    my $fh        = shift or croak qq{supply filehandle};
    my $verbosity = shift || 0;
    my $header    = shift || q{##gff-version 3};

    my %get_attributes = (
        gene => {
                  attributes => [ 'ID', 'Name', 'locus_tag' ],
                  fields     => { begin => 3, end => 4, strand => 6}, # field numbering starts at 0
                },
        cds  => {
                  attributes => [ 'ID', 'Name', 'Parent', 'product', 'protein_id' ],
                  fields     => { sequenceID => 0, begin => 3, end => 4},
		  key        => 'protein_id', 
                },
        mrna => {
                  attributes => [ 'ID', 'Parent' ],
                  fields     => {  },
                },
    );

    my @get_types = (keys %get_attributes);

    my %feature_list;

    my $found_header = 0;

    my $n_lines      = 0;

    LINE:
    while (defined (my $line = <$fh>)) {
        $line =~ s{ \A \s* ( \S .* \S ) \s* \n? \z }{$1}xms;

        if ($line =~ m{ \A [#] }xms) {
            next LINE if ($found_header);
            if ($line =~ m{ \A $header \z }xms) {
                $found_header = 1;
                next LINE;
            }
        }

        if (!$found_header) {
            print STDERR qq{WARNING: missing header (expecting: "$header")\n\n};
            $found_header = -1; # prevents error being repeated every line
        }

        my @fields = split qq{\t}, $line;
        if ($verbosity > 3) {
            for my $i (0 .. $#fields) { printf qq{\t(%d) %s\n},$i+1,$fields[$i]; }
            print qq{\n};
        }

        my ($sequenceID, $source, $type, $begin, $end, $score, $strand, $phase,
            $attribute_list) = @fields;

        my @attributes = split m{ [;] }xms, $attribute_list;
        my %attribute;

        for my $attrib_pair (@attributes) {
            my ($key, $value) = split m{ [=] }xms, $attrib_pair;
            while ($value =~ m{ \A ( [^%]* ) [%] ([0-9A-F]{2} ) ( .* ) \z}xms) {
                my ($prefix, $hexstr, $suffix) = ($1, $2, $3);
                my $hexno = eval "0x$hexstr";
                $value = $prefix . chr($hexno). $suffix;
            }
            $attribute{$key} = $value;
        }

        if ($verbosity > 4) {
            for my $attrib_name (sort keys %attribute) {
                print STDERR qq{\t\t$attrib_name\t=\t$attribute{$attrib_name}\n};
            }
            print STDERR qq{\n};
        }

        if (first { lc $_ eq lc $type } @get_types) {

            #my $type = $attribute{$type};
            my %this_feature;

            for my $attrib_name (@{$get_attributes{lc $type}->{attributes}}) {
                my $value = undef;
                if (defined $attribute{$attrib_name}) {
                    $value = $attribute{$attrib_name};
		    $this_feature{$attrib_name} = $value;
                }
                

            }
            for my $field_name (keys %{$get_attributes{lc $type}->{fields}}) {
                my $field_no = $get_attributes{lc $type}->{fields}->{$field_name};
                my $value = undef;
                if (defined $fields[$field_no]) {
                    $value = $fields[$field_no];
		    $this_feature{$field_name} = $value;
                }
                
            }

            my $use_key;
            if (defined $get_attributes{lc $type}->{key}) {
		$use_key = $get_attributes{lc $type}->{key};
		if (!defined $this_feature{$use_key}) {
		    $this_feature{PRIMARY} = q{[undefined]};
		}
		else {
		    $this_feature{PRIMARY} = $this_feature{$use_key};
		}
	    }

            if ($verbosity > 4) {
                print STDERR qq{\trecorded feature attributes:\n\n};
		my @check_attribs = (sort keys %this_feature);
		# this fills in the (non-unique) primary key as '[undefined]' in the
		# event that the key field is undefined
		if (defined $use_key) {
		    @check_attribs = uniq (@check_attribs, $use_key);
		}

                for my $f_name (@check_attribs) {
                    print qq{\t\t\t$f_name\t=\t},
		        defined ($this_feature{$f_name})
			    ? qq{$this_feature{$f_name}}
			    : qq{[undefined]}
			    , qq{\n};
			
                }
                print STDERR qq{\n};
            }

            push @{$feature_list{lc $type}}, \%this_feature;

        }
	$n_lines++;
    }

    return (\%feature_list, $n_lines);

}

sub write_buffer {

    my $buffer_lref = shift;
    my $key_length  = shift || 4;
    my $output_dest = shift || \*STDOUT;

    my $n_written_lines;

    DESTINATION: {
        (ref $output_dest eq 'GLOB') && do {
	    # could be to STDOUT for example
	    $n_written_lines = write_buffer_to_fh($buffer_lref, $key_length, $output_dest);
            last DESTINATION;
        };
        (-d $output_dest) && do {
	    $n_written_lines = write_buffer_to_files($buffer_lref, $key_length, $output_dest);
            last DESTINATION;
        };
	croak qq{$output_dest does not exist, or is neither a directory nor a filehandle};
    }

    return $n_written_lines;
}

sub write_buffer_to_fh {

    my $buffer_lref = shift;
    my $key_length  = shift || 4;
    my $output_fh   = shift || \*STDOUT;

    my @buffer_lines = sort @{$buffer_lref};

    print $output_fh @buffer_lines or croak qq{couldn't write to filehandle $output_fh};

    my $n_written_lines = scalar @buffer_lines;
    return $n_written_lines

}

sub write_buffer_to_files {

    my $buffer_lref = shift;
    my $key_length  = shift || 4;
    my $output_dir  = shift || q{.};

    my @buffer_lines = @{$buffer_lref};

    my $n_written_lines = 0;

    my $null_protein_id      = q{[undefined]};
    my $null_protein_id_key  = substr($null_protein_id, 0, $key_length);
    my $null_protein_id_file = q{_undefined_};

    # ensure that the abbreviation (length determined by $key_length) does
    # not extend beyond the first field, i.e. include whitespace and possibly
    # some of the second field as well
    my @unique_keys = uniq (map {
        my $abbrev = substr($_,0,$key_length);
        $abbrev =~ s{ \s .* $ }{}xms;
        $abbrev;
    } @buffer_lines);

    print STDERR qq{writing }, scalar @buffer_lines, qq{ lines to },
                 scalar @unique_keys, qq{ different files\n\n};

    for my $key (sort @unique_keys) {
###print STDERR qq{key is '$key'\n};
        # qr() will fall over with some argument strings. E.g. if
	# it is given '[und' then it falls foul of the opening bracket.
	# So it is necessary to apply quotemeta first. However, you
	# can't do this:
        # my $key_regexp = qr(quotemeta($key));
	# because then the 'quotemeta' keyword is treated as a literal string,
	# which qr regexp-izes.
	my $key_quotemeta = quotemeta($key);
###print STDERR qq{key quotemeta'd is '$key_quotemeta'\n};
	my $key_regexp = qr{$key_quotemeta};
###print STDERR qq{key regexp is '$key_regexp'\n};
        my @outlines = grep { m{ \A $key_regexp }xms } @buffer_lines;

        my $key_file = $key;
	# take care of the few cases of CDS with no defined protein ID;
	# these will have ended up with '[undefined]' as the protein ID
	# string
	if ($key eq $null_protein_id_key) {
	    $key_file = $null_protein_id_file;
	}
        my $outfile = qq{$output_dir/$key_file};
	my $fh;

	if (-f $outfile) {
	    open $fh, '>>', $outfile or croak qq{couldn't append to file $outfile\n};
	    print STDERR qq{\tappending };
	}
	else {
	    open $fh, '>', $outfile or croak qq{couldn't create file $outfile\n};
	    print STDERR qq{\tcreating and writing };
	}
	print STDERR scalar @outlines, qq{ lines to $outfile\n};
        print $fh @outlines;
	$n_written_lines += scalar @outlines; # sanity check
	close $fh or croak qq{\tcouldn't close file $outfile\n};
    }
    print STDERR qq{\twrote $n_written_lines lines\n\n};
    return $n_written_lines;
}
