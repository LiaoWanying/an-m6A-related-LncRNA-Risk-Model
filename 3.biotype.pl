use strict;
use warnings;


my $gtfFile="human.gtf";
my $expFile="symbol.txt";
my $outmRNAFile="mRNA.txt";
my $outlncRNAFile="lncRNA.txt";
my %hash=();


open(RF,"$gtfFile") or die $!;
while(my $line=<RF>)
{
	chomp($line);
	if($line=~/gene_id \"(.+?)\"\;.+gene_name "(.+?)"\;.+gene_biotype \"(.+?)\"\;/)
	{
	    my $symbol=$2;
	    my $biotype=$3;
        if($biotype eq "protein_coding"){
			$hash{$symbol}=$biotype;
		}
		unless(exists $hash{$symbol}){
		    $hash{$symbol}=$biotype;
		}
	}
}
close(RF);			  

			  
open(RF,"$expFile") or die $!;
open(WF,">$outmRNAFile") or die $!;
open(ZF,">$outlncRNAFile") or die $!;
while(my $line=<RF>)
{
	if($.==1){
		print WF $line;
		print ZF $line;
		next;
	}
	chomp($line);
	my @arr=split(/\t/,$line);
	if(exists $hash{$arr[0]}){
	  	if($hash{$arr[0]} eq "protein_coding"){
		    print WF $line,"\n";	 
		}
		elsif($hash{$arr[0]} eq "lncRNA"){
		    print ZF $line,"\n";		
		}
	}
}
close(WF); 
close(RF);	
close(ZF);			  


			  
	  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
