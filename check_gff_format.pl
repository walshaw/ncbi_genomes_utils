#!/usr/bin/perl

use strict;
use warnings;
use autodie;

use Carp;

use List::Util qw{ first sum };
###use List::MoreUtils qw{ all any each_array first_index mesh pairwise uniq };

my $gff_file = shift;

my @field_names = qw(
altitude
anticodon
bio-material
biotype
biovar
chromosome
codons
collected-by
collection-date
country
culture-collection
description
Dbxref
end_range
environmental-sample
exception
gbkey
gb-synonym
gene
gene_biotype
gene_synonym
genome
genotype
group
ID
identified-by
Is_circular
is_ordered
isolate
isolation-source
lat-lon
locus_tag
metagenome-source
metagenomic
mol_type
nat-host
Name
ncrna_class
Note
old-lineage
old_locus_tag
old-name
Parent
part
partial
pathovar
plasmid-name
product
protein_id
pseudo
pseudogene
serogroup
serotype
serovar
start_range
strain
sub-species
subgroup
substrain
subtype
transl_except
transl_table
type
);


my $fname_regexp = join( ' | ', @field_names );

my %field_tally;

open my $fh, '<', $gff_file;

while (defined (my $line=<$fh>) ) {

    next if $line =~ m{ \A \s* [#] }xms;
    next if $line !~ m{ \S }xms;

    chomp $line;

    my @fields = split m{ \t }xms, $line;

    my ($sequenceID, $source, $type, $begin, $end, $score, $strand, $phase,
        $attribute_list) = @fields;

    pos $attribute_list = 0;

    my $last_match_pos = 0;

    my %data;

    while ($attribute_list =~ m{ \G ( [^=]+ ) ( [=] [^;]+ ) (?: [;] | \z ) }gxms) {
#####        while ($attribute_list =~ m{ \G ( $fname_regexp ) [=] [^;]+ (?: [;] | \z ) }gxms) {
        my ($field_name, $field_value) = ($1, $2);
        $last_match_pos = pos $attribute_list;
        $field_tally{$field_name}++;

        croak qq{unknown field name '$field_name'}
            if (!first { lc $_ eq lc $field_name } @field_names);

        croak qq{(file $gff_file):problem with value '$field_value' (field $field_name)\n}.
              qq{attribute list string:\n(file $gff_file):$attribute_list\n}.
              qq{in line:\n(file $gff_file):$line} if check_parentheses($field_value);
    }

#print qq{length }, length $attribute_list, qq{; pos }, (defined (pos $attribute_list)) ?
 #pos $attribute_list : q{undefined}, qq{; last match pos $last_match_pos\n};

    croak qq{(file $gff_file):unexpected attribute list format:\n$attribute_list\n}.
          qq{\n(file $gff_file):problem at/after position $last_match_pos (length of string is }.
          length($attribute_list).qq{, last match at $last_match_pos}
        if ($last_match_pos != length $attribute_list);

}

close $fh;

#for my $fname (keys %field_tally) {

#    print qq{$fname\t$field_tally{$fname}\n};

#}

print qq{OK	file $gff_file\n};

exit 0;

sub check_parentheses {

    my $string = shift;

    my $n_opening_parentheses = ($string =~ s{ [(] }{}gxms);
    my $n_closing_parentheses = ($string =~ s{ [)] }{}gxms);
    my $n_opening_brackets    = ($string =~ s{ \[  }{}gxms);
    my $n_closing_brackets    = ($string =~ s{ \]  }{}gxms);
    my $n_opening_braces      = ($string =~ s{ \{  }{}gxms);
    my $n_closing_braces      = ($string =~ s{ \}  }{}gxms);
    
    my $parenthesis_problems = $n_opening_parentheses - $n_closing_parentheses;
    my $bracket_problems     = $n_opening_brackets    - $n_closing_brackets;
    my $brace_problems       = $n_opening_braces      - $n_closing_braces;

    my @is_problem = map { ($_ != 0 ) || 0 }
                  ($parenthesis_problems, $bracket_problems, $brace_problems);

    my $n_problems = sum @is_problem;

    return ($n_problems);
}

=pod
    pos $string = 0;

    my $nested = 0;

    CHAR:
    while ( $string =~ m{ \G ( . | \z ) }gxms ) {

        my $char = $1;

        ($char eq '(') && do {
            $nested++;
            $gene_label .= $char;
                            next CHAR;
        };

                        ($char eq ')') && do {
                            croak qq{unmatched closing parenthesis:\n$string}
                                if --$nested < 0;
                            $gene_label .= $char;
                            next CHAR;
                        };

                        ($char eq ';') && do {
                            if ($nested) {
                                $gene_label .= $char;
                                next CHAR;
                            }
                            push @gene_labels, $gene_label;
                            last CHAR;
                        };

                        if (!length $char) { push @gene_labels, $gene_label; last LABEL;}
                        $gene_label .= $char;

                    } # end CHAR block
                }



}

=cut

