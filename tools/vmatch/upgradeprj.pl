#!/usr/bin/env perl

# read old versions of Vmatch prj file and output new format.

use strict;
use warnings;

my $argcount = scalar @ARGV;
my $prjfile;

if ($argcount eq 1)
{
  $prjfile = $ARGV[0];
} else
{
  print "Usage: $0 <prjfile>\n";
  exit 1;
}

unless ( -e $prjfile)
{
  print STDERR "file \"$prjfile\" does not exist\n";
  exit 1;
}

unless(open(FP, $prjfile))
{
  print STDERR "Cannot open file $prjfile\n";
  exit 1;
}

while(my $line = <FP>)
{
  if($line =~ m/^specialcharacters/)
  {
    my @linearray = split(' ',$line);
    printf("%s\n",$linearray[0]);
    printf("specialranges=%s\n",$linearray[1]);
    printf("lengthofspecialprefix=%d\n",$linearray[2]);
    printf("lengthofspecialsuffix=%d\n",$linearray[3]);
  } else
  {
    print $line;
  }
}
