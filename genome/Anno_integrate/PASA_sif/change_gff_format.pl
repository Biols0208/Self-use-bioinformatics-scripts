#!/usr/bin/perl
use warnings;
use strict;
@ARGV|| die "Usage: perl  $0 file1 out \n";
my ($file1,$out)=@ARGV;
#my %hash;
my $i;
open FILE1,"<$file1";
open OUT,">$out";
while(<FILE1>){
	chomp;
	my@data=split/[\t||=||;]/,$_;
	if($data[2]=~/mRNA/){
 		$i=0;
		print OUT "$data[0]\t$data[1]\tgene\t$data[3]\t$data[4]\t$data[5]\t$data[6]\t$data[7]\t$data[8]=gene_$data[9];\n";
		print OUT "$data[0]\t$data[1]\tmRNA\t$data[3]\t$data[4]\t$data[5]\t$data[6]\t$data[7]\t$data[8]=$data[9];Parent=gene_$data[9];\n";
#		print $_."_"."mRNA;Parent=$data[9]_gene;\n";
		
	}
	if($data[2]=~/CDS/){
		print OUT "$data[0]\t$data[1]\texon\t$data[3]\t$data[4]\t$data[5]\t$data[6]\t$data[7]\tID=$data[9]_exon$i;Parent=$data[9];\n";
		print OUT "$data[0]\t$data[1]\tCDS\t$data[3]\t$data[4]\t$data[5]\t$data[6]\t$data[7]\tID=cds_$data[9];Parent=$data[9];\n";
		#print $_."_"."cds;Parent=$data[9];\n";
	}
		$i++;
}
close FILE1;
