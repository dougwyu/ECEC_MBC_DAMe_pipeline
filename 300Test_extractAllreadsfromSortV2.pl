#!/usr/local/bin/perl -w
#
# Writer:         wangxy <wxyang1988@126.com>
# Program Date:
# Modifier:       wangxy <wxyang1988@126.com>
# Last Modified:
##########################################################

my $ver="1.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

#########################################################

my %opts;
GetOptions(\%opts,"id=s","o=s","h" );

#&help()if(defined $opts{h});
if(!defined($opts{id}) || !defined($opts{o}) ||defined($opts{h}))
{
	print <<"	Usage End.";
	Description:

		Version: $ver

	Usage:

		Please put all libraries (A-F) into 1 folder

		-id     Input folder which contains A-F                must be given

		-o      Output otu-id with taxon-info                 must be given

		-h    Help document

		example: 

		# cd 300Test/data  # this folder holds the seqs/ folder, which holds Folder{A,B,C,D,E,F}/
		# perl 300Test/scripts/300Test_extractAllreadsfromSort.pl -id seqs/ -o sep_pools

	Usage End.

	exit;
}

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
################
my $programe_dir=basename($0);
my $ind=$opts{id};
my $out=$opts{o};

my @allfiles = glob "$ind/*/pool*/Tag*.txt";
my $count=0;

open OUT,">$out" || die "Can't create $out,$!\n";

for (my $i=0;$i<@allfiles;$i++)
{
	my @a = split/\//,$allfiles[$i];
	my $rep = $a[-3];
	my $num = (split/pool/,$a[-2])[1];
	my $tag_name = (split/\./,$a[-1])[0];
	my ($t1,$t2) = split/_/,$tag_name;
	if ($t1 eq $t2)
	{
		open IN,"$allfiles[$i]";
		while (<IN>)
		{
			chomp;
			$count++;
			my @b = split/\t/,$_;
			print OUT ">$rep","$num"," ","$t1-$t2:$count ","count=$b[3]\n$b[4]\n";
			#if ($b[3] == 1) {$count++;print OUT ">$rep","$num","_","$t1-$t2","_","$count\n$b[4]\n";}
			#else {for (my $k=1;$k<=$b[3];$k++) {$count++;print OUT ">$rep","$num","_","$t1-$t2","_","$count\n$b[4]\n";}}
		}
	}



}

close (OUT);
###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";

###############Subs
sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
