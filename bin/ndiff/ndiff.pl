#! /usr/bin/perl
# 
# Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
#
# Please make a following refer to Flow123d on your project site if you use the program for any purpose,
# especially for academic research:
# Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
#
# This program is free software; you can redistribute it and/or modify it under the terms
# of the GNU General Public License version 3 as published by the Free Software Foundation.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program; if not,
# write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
#
# $Id$
# $Revision$
# $LastChangedBy$
# $LastChangedDate$
#

use File::Basename;


# BEGIN block is executed when perl compiles the file
BEGIN {
	$scriptdir = &File::Basename::dirname($0);
	# May want to use shift to put it at the beginning
	push @INC, $scriptdir;
}


#
# Usage:
# ndiff [-w] [-r rel_tol] [-a abs_tol] [-o log_file] file1 file2
#
# Compare files 'file1' 'file2' like the diff system command, but
# try to detect numbers on matching places of both files and 
# treat them as equal if their absolute difference or relative difference 
# is less then given values.
#
# -w 			Ignore diference in white space material.
# -r rel_tol	Set tolerance for relative difference. (default: 0.01)
# -a abs_tol	Set tolerance for absolute difference. (default: 0.0001)
# -o log_file   Description of differences is redirected to the given file
#				(given filename beginning with 'greater then' like ">logfile"
#				allow appending to the file)
#		default: ndiff.log
#
# Program reports total number of numerical diferences and max. norm of all
# absolute and all relative differences.
#
#
# TODO : 
# - gruped discards, line number info on the separate line
# - faster initial parsing

use lib "$scriptdir";
use Diff;
use Scalar::Util;

# global variables
# constants, setting
my $machine_prec = 1.0E-30;
my $rel_tol=0.01;
my $abs_tol=0.0001;
my $log_name="-"; # dafult is write to the output
my $inputfilter;

# work space
my @left;
my @right;
my @left_discard;
my @right_discard;
my $abs_norm=0;			# Linf norm of absolute differences
my $rel_norm=0;   		# Linf norm of relative differences
my $num_of_diffs=0;		# total number of value differences


# ---------------------------------- MAIN PROGRAM

# read arguments
while ($_=shift @ARGV) {
	if (/-r/) { 
		$rel_tol=shift @ARGV;
	} elsif (/-a/) {
		$abs_tol=shift @ARGV;
	} elsif (/-w/) { # ignore white space == all white spece are equal to one space
		$inputfilter='s/[ \t\n]*/ /g';
	} elsif (/-o/) {
		$log_name=shift @ARGV;
	} else {
		unshift(@ARGV,$_);
		last;
	}
}
#print "arg: @ARGV\n";

# read the files
@left=readfile(shift @ARGV);	# reads the left file
@right=readfile(shift @ARGV);	# reads the right file

# open log file for output of differences
open( OUT_LOG, ">$log_name");

# ------- find LCS of both lines lists
Diff::traverse_sequences(\@left,\@right,
        {   MATCH=> \&match,
        	DISCARD_A => \&discard_left,
            DISCARD_B => \&discard_right
        },
        \&keygen_line       
    ); 

print "\n";
print "Total num. of diff.      : $num_of_diffs\n";
print "L-inf norm of abs. diff. : $abs_norm\n";
print "L-inf norm of rel. diff. : $rel_norm\n";
close(OUT_LOG);
if ($num_of_diffs > 0) {exit 1;}    
    
# ------------------------------------------------------------ SUBROUTINES    
    
#####################    
# given element of the lines list, we return the filtered line string    
sub keygen_line {
	return $_[0]->{text_only};
}    

#####################
# push discarded elements to the bag    
sub discard_left {
	my ($ia, $ib)=@_;
	my $aln=$left[$ia];
	my $bln=$right[$ib]; 
        my $b_line_num=$bln->{line_no} || "eof";
	print OUT_LOG  "$aln->{line_no},$b_line_num < ".$aln->{line};
	# push(@left_discard,$ia);	
        $num_of_diffs ++;
}    
sub discard_right {
	my ($ia, $ib)=@_;
	my $aln=$left[$ia];
	my $bln=$right[$ib]; 
        my $a_line_num=$aln->{line_no} || "eof";
 
	print OUT_LOG  "$a_line_num,$bln->{line_no} > ".$bln->{line};
	# push(@right_discard,$ib);	
        $num_of_diffs ++;
}    

#################################
# min, max
sub max {
	return ( $_[0]<$_[1] ? $_[1] : $_[0] );
}
sub min {
	return ( $_[0]>$_[1] ? $_[1] : $_[0] );
}

#####################
# if sturcture of lines match, we:
# 1) empty discard bags,
# 2) check if the values match too
sub match {
	my ($ia, $ib)=@_;
	my $aln=$left[$ia];
	my $bln=$right[$ib];
	my @messages=();
	
	# empty the bags (future - gruped discards)
	
	# check lines
#	print "$aln->{reals}, $bln->{reals}\n";
	my $l_ref=$aln->{reals};
	my $r_ref=$bln->{reals};
	
	if (scalar(@{$l_ref}) != scalar(@{$r_ref}) ) {
                $num_of_diffs++;
                push(@messages,"#floats differs, ");
		#die( "PrgErr: lines match but have different number of floats:\n".$aln->{line}.$bln->{line} );
	}		 
		
	$field_num=0;	
	while ( scalar(@{$l_ref}) &&  scalar(@{$r_ref}) ) {
                $lfloat=shift @{$l_ref};
                $rfloat=shift @{$r_ref}; 
		$field_num++;
		$max = max( abs($lfloat),abs($rfloat) );
		if ($max < $machine_prec ) {next;}	# skip too small numbers
		$adiff=abs($lfloat-$rfloat);
		$abs_norm=max($abs_norm,$adiff);
		$rdiff=$adiff/$max;
		$rel_norm=max($rel_norm,$rdiff);
		if ( $adiff < $abs_tol ) {next;}		
		if ( $rdiff < $rel_tol ) {next;}
		
		$num_of_diffs++;
		# report bigger of differences
		if ($adiff > $rdiff) {
			push(@messages,"A$field_num $adiff ");	
		} else {
			push(@messages,"R$field_num $rdiff ");
		}
	}	
	
	if (scalar(@messages)) {
		print OUT_LOG ">ndiff< @messages\n";
		print OUT_LOG "$aln->{line_no},$bln->{line_no} < ".$aln->{line};
		print OUT_LOG "$aln->{line_no},$bln->{line_no} > ".$bln->{line};
                print OUT_LOG "--\n"
	}		
}

# read the file and retrieve the floats
sub readfile {
	$fname=shift @_;
	my $line_no;	
	my $float;
        my $token;
        my $tmp;
	my @file_list;
        my @line_list;
	my $float_regexp = qr/([+-.0-9]+(?:[Ee][+-]?\d+)?)\b/;
	
	open FILE, $fname  or die "$!  \"$fname\" "; # or return the system error message
	# read all lines
	$line_no=1;
	while($_=<FILE>) {
		# cereate a new line item
		$line_ref={
			line=>$_,
			line_no=>$line_no
		};
#		print "line $line_no: $_";
		
		# substitute all floats by '0' and store them
		$reals=[];
		
# TODO : faster detection of real numbers
# this match only real number format		 
#		while ( s/([+-]?(?=\d|\.\d)\d*(?:\.\d*)?(?:[Ee](?:[+-]?\d+))?)\b// ) {
# this match more strings, but it is faster	
#		/([+-.0-9]+(?:[Ee][+-]?\d+)?)\b/
#
		#while ( s/$float_regexp// ) {
		#	#print "match: <$1> line: $_";
		#	$float=$1+0;	# cast to float -  does perl remember that we store number not string			
		#	push(@{$reals},$float);
		#}
		
		# this is little bit faster, but there is some bug  
		@line_list=split;
                foreach my $token (@line_list) {  
                    #$save=$token."0";
                    #$tmp=$save;  
                    #$diff = abs( (++$tmp) - 1 - $save );
                    #print "dgb: /$token/$save/$tmp/$diff/\n";
                    #if ($diff < 1.0e-15) {
                    if (Scalar::Util::looks_like_number($token) && $token=~/[.eE]/) {
                        $float = $token + 0;
                        $token='0';
			push(@{$reals},$float);
		    }
                }
		$_=join " ", @line_list;
                #print "$fname : $_\n";
#                foreach my $i (@{$reals}) {print "$i ";}
#                print "\n";

		s/\r//;	# remove CR (for windows compatibility)
		eval $input_filter; # input filter command s/// ...
		
		# finish the line item 
		$line_ref->{text_only}=$_;
		$line_ref->{reals}=$reals;
		push(@file_list,$line_ref);
		
		$line_no++;
	}
	close FILE;

        # add \n to the last line if it is not present
        my $end_ref = $file_list[-1];
        $end_ref->{line} =~ s/[^\n]$/\n/;
	return @file_list;
}

