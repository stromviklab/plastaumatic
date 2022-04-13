#!/usr/bin/env perl
#David A. Parry, 5th November 2015

use warnings;
use strict;
use Getopt::Long;
my %codons =
(
    AAA=>"K", AAC=>"N", AAG=>"K", AAT=>"N",
    ACA=>"T", ACC=>"T", ACG=>"T", ACT=>"T",
    AGA=>"R", AGC=>"S", AGG=>"R", AGT=>"S",
    ATA=>"I", ATC=>"I", ATG=>"M", ATT=>"I",
    CAA=>"Q", CAC=>"H", CAG=>"Q", CAT=>"H",
    CCA=>"P", CCC=>"P", CCG=>"P", CCT=>"P",
    CGA=>"R", CGC=>"R", CGG=>"R", CGT=>"R",
    CTA=>"L", CTC=>"L", CTG=>"L", CTT=>"L",
    GAA=>"E", GAC=>"D", GAG=>"E", GAT=>"D",
    GCA=>"A", GCC=>"A", GCG=>"A", GCT=>"A",
    GGA=>"G", GGC=>"G", GGG=>"G", GGT=>"G",
    GTA=>"V", GTC=>"V", GTG=>"V", GTT=>"V",
    TAA=>"*", TAC=>"Y", TAG=>"*", TAT=>"Y",
    TCA=>"S", TCC=>"S", TCG=>"S", TCT=>"S",
    TGA=>"*", TGC=>"C", TGG=>"W", TGT=>"C",
    TTA=>"L", TTC=>"F", TTG=>"L", TTT=>"F"
);

my %opts = ();
GetOptions(
    \%opts,
    'i|input=s',        #dna input
    'o|output=s',       #optional output file
    'l|line_length=i',  #no. characters per line in output
    'f|frame=i',        #frame of translation (1 to 3 or -1 to -3)
    'n|numbers',        #flag if true than add numbers to beginning of each line
    'd|dna',            #output DNA as well
    's|stop',           #stop at first termination codon
    'h|?|help',
) or usage("Syntax error.\n");

usage() if $opts{h};

usage("-i/--input is required" ) if (not $opts{i});

my $line_length = $opts{l} ? $opts{l} : 60;
usage("-l/--line_length must be greater than 0") if $line_length < 0;

my $frame = defined $opts{f} ? $opts{f} : 1;
if (abs($frame) > 3 or not $frame){
    usage("-f/--frame must be between 1 and 3 or between -1 and -3\n");
}

my $IN;
if (-f $opts{i} or -p $opts{i} or $opts{i} eq '-'){
    open ($IN, $opts{i}) or die "Could not open $opts{i} for reading: $!\n";
}else{
    $opts{i} =~ s/['"]//g;
    if ($opts{i} =~ /^[\sacgtuACGTU]+$/){
        $opts{i} =~ s/\s+/\n/g;
        open ($IN, "<", \$opts{i})
            or die "Error getting filehandle from input string: $!\n";
    }else{
        die "ERROR: '$opts{i}' does not exist as a file and does not look ".
            "like a DNA sequence.\n";
    }
}

my $OUT = \*STDOUT; # main output defaults to STDOUT
if ($opts{o}){
    open ($OUT, ">", $opts{o}) or die "Can't open $opts{o} for writing: $!\n";
}
my $prev_header = '';
my $dna = '';
while (my $line = <$IN>){
    chomp $line;
    if ($line =~ /^>(.*)/){
        my $header = $1;
        translateDna($dna, $prev_header) if $dna;
        $prev_header = $header;
        $dna = '';
    }else{
        $dna .= $line;
    }
}
translateDna($dna, $prev_header) if $dna;

#########################################################
sub translateDna{
    my ($dna, $name) = @_;
    my $cdna = '';
    my $protein = '';
    ($dna, $cdna) = getDnaInFrame($dna);
    while ($cdna =~ /(...)/g){
        my $codon = uc($1);
        my $aa = '?';
        if (exists $codons{$codon}){
            $aa = $codons{$codon};
        }
        $protein .= $aa;
        last if ($opts{s} and $aa eq '*');
    }
    outputTranslation($dna, $protein, $name);
}

#########################################################
sub outputTranslation{
    my ($dna, $protein, $name) = @_;
    if ($opts{d}){
        outputDnaAndProtein($dna, $protein, $name);
    }else{
        outputTranslationOnly($protein, $name);
    }
    print $OUT "\n";
}

#########################################################
sub outputTranslationOnly{
    my ($protein, $name) = @_;
    print $OUT ">$name\n" if ($name);
    my $num_length = length(length($protein));
    for (my $i = 0; $i < length($protein); $i += $line_length){
        my $l = $line_length;
        if ($i + $line_length > length($protein)){
            $l = length($protein) - $i ;
        }
        my $p = substr($protein, $i, $l);
        if ($opts{n}){
            my $n = $i + 1; #protein position
            $p = sprintf("%${num_length}d: %s", $n, $p);
        }
        print $OUT "$p\n";
    }
}

#########################################################
sub outputDnaAndProtein{
    my ($dna, $protein, $name) = @_;
    my $max_protein = length($protein);
    my $num_length = length(length($dna));
    $protein = join('', map { "-$_-" } split('', $protein) );
    $protein = '-' x (abs($frame) -1)  . $protein;
    print $OUT ">$name\n\n" if ($name);
    for (my $i = 0; $i < length($dna); $i += $line_length){
        my $l = $line_length;
        if ($i + $line_length > length($dna)){
            $l = length($dna) - $i ;
        }
        my $pl = $line_length;
        if ($i + $line_length > length($protein)){
            $pl = length($protein) - $i ;
        }
        my $d = substr($dna, $i, $l);
        my $p = '';
        if ($i < length($protein)){
            $p = substr($protein, $i, $pl);
        }
        if ($opts{n}){
            my $n = $i + 1; # dna position
            $d = sprintf("%${num_length}d: %s", $n, $d);
            my $cn = $i - (abs($frame) - 1); # 0-based coding DNA position of first letter of line
            my $pn = 1; # protein position...
            if ($i > 0){
                $pn = int($cn/3) + 1 ;
                $pn++ if ($cn % 3 == 2); # if first nucleotide of line is the last letter of
                                         # codon we represent the next protein 'letter'
            }
            $pn = $pn <= $max_protein ? $pn : $max_protein;
            $p = sprintf("%${num_length}d: %s", $pn, $p);
        }
        print "$d\n$p\n\n";
    }
}

#########################################################
sub getDnaInFrame{
    my $dna = shift;
    if ($frame < 0){
        $dna = revComp($dna);
    }else{
        $dna =~ tr/uU/tT/;
    }
    my $f = abs($frame);
    return ($dna, substr($dna, $f - 1) );
}

#########################################################
sub revComp{
    my $dna = shift;
    $dna =~ tr/uU/tT/;
    $dna =~ tr/acgtACGT/tgcaTGCA/;
    return reverse($dna);
}

#########################################################
sub usage{
    my $msg = shift;
    print STDERR "ERROR: $msg\n" if $msg;
    print <<EOT

    usage: $0 -i <dna_input.fa>  [options]
    usage: $0 -i ATGCCGCTACGC... [options]

    Options:

    -i, --input FILE/STRING
        Either an input file containing a single DNA sequence to be
        translated (FASTA header optional) or several FASTA formatted sequences
        OR
        a single DNA sequence (containing only the standard DNA/RNA letters
        and whitespace) to be translated. If your sequence contains any spaces
        you must enclose it with quotes (e.g. "ATG CCC GGG").

    -o, --ouptut FILE
        Output file. Optional. Default is STDOUT.

    -l, --line_length INT
        Number of letters per line. Default = 60.

    -f, --frame INT
        Frame of translation. Valid values are 1 to 3 for translation
        starting at the first to third letter of the DNA input or -1 to -3
        for translation of the reverse complement. Default = 1.

    -n, --numbers
        Use this flag to number the DNA and peptide output.

    -d, --dna
        Output DNA lined up with protein.

    -s, --stop
        Stop translation at first termination codon. Default behaviour is
        to translate all codons in input.

    -?, -h, --help
        Show this help message.


    Copyright (C) 2015  David A. Parry

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
    MA  02110-1301, USA.

EOT
;
    exit 1 if $msg;
    exit;
}
