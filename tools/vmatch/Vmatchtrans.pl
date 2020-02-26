#!/usr/bin/env perl

# transform Vmatch output format
# Stefan Kurtz, Feb 16, 2005.

use strict;
use warnings;

my $usage = 'chainer|open|posdist [inputfile]';

my $numofargs = scalar @ARGV;
my $filename = '';
my $format;

if($numofargs eq 1)
{
  $format = $ARGV[0];
} elsif($numofargs eq 2)
{
  $format = $ARGV[0];
  $filename = $ARGV[1];
} else
{
  print STDERR "Usage: $0 $usage\n";
  exit 1;
}

my $chainernum = 0;
my $opennum = 1;
my $posdist = 2;
my $formatnum;

if($format eq 'chainer')
{
  $formatnum = $chainernum;
  print ">CHA 2\n";
} elsif($format eq 'open')
{
  $formatnum = $opennum;
} elsif($format eq 'posdist')
{
  $formatnum = $posdist;
} else
{
  print STDERR "$0: format \"$format\" not yet implemented\n";
  exit 1;
}

if($filename ne '')
{
  unless(open(FILEPTR,$filename))
  {
    print STDERR "$0: cannot open $filename\n";
    exit 1
  }
  while(my $line = <FILEPTR>)
  {
    transforminputline($formatnum,$line);
  }
} else
{
  while(my $line = <STDIN>)
  {
    transforminputline($formatnum,$line);
  }
}

sub transforminputline
{
  my($formatnum,$line) = @_; 

  if($line =~ m/^\#/)
  {
    if($formatnum eq $opennum)
    {
      print $line;
    }
  } else
  {
    my @linearray = split(' ',$line);
    my $length1 = $linearray[0];
    my $start1 = $linearray[2];
    my $length2 = $linearray[4];
    my $start2 = $linearray[6];
    my $weight = $linearray[9];
    if($formatnum eq $chainernum)
    {
      printf("# %d\n",$weight);
      printf("[%d,%d] [%d,%d]\n",
             $start1,$start1+$length1-1,
             $start2,$start2+$length2-1);
    } elsif($formatnum eq $opennum)
    {
      printf("%7d\t%7d\t%7d\t%7d\t%7d\n",
             $start1,$start1+$length1-1,
             $start2,$start2+$length2-1,
             $weight);
    } elsif($formatnum eq $posdist)
    {
      my $end1 = $start1+$length1-1;
      my $gap = $start2-($end1+1);
      printf("%d %d %d %d %d\n",
             $start1,$end1,
             $start2,$start2+$length2-1,
             $gap);
    }
  }
}
