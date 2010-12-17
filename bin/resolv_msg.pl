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

$msg_fname="messages.h";
$err_fname="error_messages.h";
$verb_treshold=3; # greater numbers leads to MsgVerb

my %msg_hash; # hash of format strings for xprintf
my %err_hash; # hash of type of error and format string 'type, "string"'
my %verb_hash; # count of verbosity types  

sub read_msg {
	open(IN,$msg_fname);
	$line=0;
	while (<IN>) {
		$line++;
		if (/^ *#define ([^ ]*) (".*")/) {
			$msg_hash{$1}=$2;
		} else {
			print "warning: line $line does not match msg file format.\n";
		}	
	}
	close(IN);	
#	print "--------------- messages: \n";
#	while ( my ($key, $value) = each(%msg_hash) ) {
#        print "$key => $value\n";
#    }
}

sub read_err {
	open(IN,$err_fname);
	print "Error file $err_fname. Ingoring lines:\n";
	$line=0;
	while (<IN>) {
		# print $_;
		$line++;
		if (/^\s*{ *([-0-9]*), "(.*)", *(USER|PGMR)/) {
		 	#print "match '$1' '$2' '$3'\n";
			if ($3 eq "USER") {
				$err_hash{$1}="UsrErr, \"(OLD): $2\"";
			} elsif ($3 eq PGMR) {
				$err_hash{$1}="PrgErr, \"(OLD): $2\"";
			} else {
				print "!! wrong error type on line $line\n";
			}
		} else {
			print "$_";
		}	
	}		
	close(IN);
#	print "---------------------- errors:\n";
#	while ( my ($key, $value) = each(%err_hash) ) {
#        print "$key => $value\n";
#    }
}

# ---------------------------------------- MAIN

read_msg();
read_err();

while ($fname=shift @ARGV) {
	print "file $fname --------------------------\n";
	open(IN,"<${fname}");
	open(OUT,">${fname}.subst");
	$line=0;
	while (<IN>) {
		$line++;
		# substitute all mprintf
		while (/myprintf/) {
			if (! /myprintf\( *([0-9]*) *,([^,]*)(,.*)?\)/) {
				print "line $line: wrong myprintf call\n";
				s|myprintf|xprintf/* TODO specify parameters*/|;
				next;
			}
			$verb_hash{$1}++; # update statistics
			$subst=$msg_hash{$2} or	$subst=$2; # substitute by message or let it be
			# select verbosity
			if ($1>$verb_treshold) {
				$subst=" MsgVerb, $subst".$3;
			} else {
				$subst=" Msg, $subst".$3;
			}
			print "l: $line myprintf '$2'\n";
			s/myprintf\( *([0-9]*) *,([^,]*)(,.*)?\)/xprintf($subst)/;									
		}
		# substitute all terminate
		while (/terminate/) {
			if (! /terminate\([^,]*,[^,]*, *([0-9]+) *(,.*)?\)/) {
				print "line $line: wrong terminate call\n";
				s|terminate|xprintf/* TODO specify parameters*/|;
				next;
			}
			if (! ( $subst = $err_hash{$1}) ) {
				print "line $line: wrong error number\n$_";
				$subst=$err_has{"-1"};
			}
			$subst.=$2; # add optional parameters
			print "l: $line terminate '$1'\n";
			s/terminate\([^,]*,[^,]*,([ 0-9]*)(,.*)?\)/xprintf($subst)/;						
		}
	}
	print "verbostity statistics:\n";
	while ( my ($key, $value) = each(%verb_hash) ) {
        print "$key => $value\n";
    }
	close(IN);
	close(OUT);
}