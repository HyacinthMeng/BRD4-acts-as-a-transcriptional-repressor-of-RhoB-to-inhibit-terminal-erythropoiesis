use strict;
use warnings;

my %hash=();
open(RF,"id.txt") or die $!;
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/,$line);
	$hash{$arr[$#arr]}="$arr[0]";
}
close(RF);

open(KEGG,"KEGG.txt") or die $!;
open(WF,">KEGG.xls") or die $!;
while(my $line=<KEGG>){
	if($.==1){
		print WF $line;
		next;
	}
	chomp($line);
	my @arr=split(/\t/,$line);
	my @idArr=split(/\//,$arr[$#arr-1]);
	my @symbols=();
	foreach my $id(@idArr){
		if(exists $hash{$id}){
		  push(@symbols,$hash{$id});
		}
	}
	$arr[$#arr-1]=join("/",@symbols);
	print WF join("\t",@arr) . "\n";
}
close(WF);
close(KEGG);



###Video source: http://study.163.com/u/biowolf
######Video source: https://shop119322454.taobao.com
######速科生物: http://www.biowolf.cn/
######作者邮箱：2740881706@qq.com
######作者微信: seqBio
######QQ群:  259208034
