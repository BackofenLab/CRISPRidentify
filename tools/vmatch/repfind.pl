#!/usr/bin/env perl

use strict;
use warnings;

# run repfind via Vmatch
# Stefan Kurtz, March, 6 2005

emulaterepfind(\@ARGV);

sub emulaterepfind
{
  my($argvptr) = @_;
  my @argv = @$argvptr;
  my $argcount = scalar @argv;

  if($argcount == 0)
  {
    printhelpbasic($0);
    exit 1;
  }

  if ($argcount == 1 && $argv[0] eq '-help')
  {
    printhelpall();
    exit 0;
  }

  if ($argcount == 1 && $argv[0] eq '-v')
  {
    printversion($0);
    exit 0;
  }

  my $vmatchoptionsref = analyzerepfindargs($0,\@argv,$argcount);
  my @vmatchoptions = @$vmatchoptionsref;
  my $inputfile = $argv[$argcount-1];
  my $indexname = `basename $inputfile`;
  chomp($indexname);
  callmkvtree($inputfile,$indexname);
  push(@vmatchoptions,$indexname);
  makesystemcall("vmatch " . join(' ',@vmatchoptions));
}

sub printhelpbasic
{
  my($programname) = @_;
  print STDERR "$programname: Missing Arguments\n";
  print STDERR "Usage: $programname [options] filename\n";
  print STDERR "try repfind.pl -help\n";
}

sub printhelpall
{
print <<HEREDOC;
-f           compute maximal forward repeats
-p           compute maximal palindromes
-l           specify that repeats must have the given length
-h           search for repeats up to the given hamming distance
-e           search for repeats up to the given edit distance
-seedsize    set the seed size
-allmax      show all maximal repeats in the order of their computation
-best        show the repeats with smallest E-value (default best 50)
-s           show the string content of the maximal repeats
-lw          format string output to given linewidth
-iub         print pair of different residues in IUB format
-nodistance  do not show distance values
-noevalue    do not compute evalues
-i           give info about number of different repeats
-v           show program version
-help        this option
HEREDOC
}

sub printversion
{
  my($programname) = @_;
  print "this is $programname,\n";
  print "a perl script emulating the options of the C-program repfind\n";
  print "by calling mkvtree and vmatch\n";
  print "copyright Stefan Kurtz 2000-2005, kurtz\@zbh.uni-hamburg.de\n";
}

sub checkdbfile
{
  my($inputfile,$prjfile) = @_;
  my $dbfile = '';

  if(-e $prjfile)
  {
    unless(open(PRJFILEPTR,$prjfile))
    {
      return 0;
    }
    while(my $line = <PRJFILEPTR>)
    {
      if($line =~ /^dbfile=(\S+) (\d+)/)
      {
        if($dbfile eq '')
        {
          $dbfile = $1;
          my $dbfilesize = $2;
          if($dbfile eq $inputfile)
          {
            if(my ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,
                   $atime,$mtime,$ctime,$blksize,$blocks) = stat($dbfile)) 
            {
              if($size == $dbfilesize)
              {
                close PRJFILEPTR;
                return 1;
              }
            }
          } 
        }
      }
      close PRJFILEPTR;
      return 0;
    }
  }
  return 0;
}

sub makesystemcall
{
  my ($argstring) = @_;
  my @args = split(' ',$argstring);
  my($retcode) = system($argstring);
  $retcode = $? >> 8;
  if($retcode ne 0)
  {
    print STDERR "failure: \"$argstring\", errorcode $?\n";
    exit 1;
  }
  print STDERR "# $argstring\n";
}

sub callmkvtree
{
  my ($inputfile,$indexname) = @_;

  if(not (checkdbfile($inputfile,"${indexname}.prj")))
  {
    makesystemcall("mkvtree -db $inputfile " .
                   "-dna -pl -lcp -suf -tis -ois -bwt -bck -sti1");
  }
}

sub analyzerepfindargs
{
  my ($program,$argvref,$argcount) = @_;
  my @vmatchoptions = ();
  my %notsupported = 
  (
    "-r" => 1,
    "-c" => 1,
    "-hrate" => 1,
    "-erate" => 1,
    "-o" => 1,
    "-b" => 1,
    "-warn" => 1,
    "-iw" => 1,
    "-mem" => 1
  );

  my @argv = @$argvref;
  my $argnum;
  my $stringoption = 0;
  my $linewidth = 0;
  my $doiub = 0;
  my $bestoption = 0;
  my $allmaxoption = 0;
  for($argnum=0; $argnum<$argcount-1; $argnum++)
  {
    my $currentarg = $argv[$argnum];
    if($currentarg eq '-f')
    {
      push(@vmatchoptions,"-d");
    } elsif($currentarg eq '-p')
    {
      push(@vmatchoptions,"-p");
    } elsif($currentarg eq '-l' || 
            $currentarg eq '-seedsize' ||
            $currentarg eq '-best')
    {
      if($currentarg eq '-seedsize')
      {
        push(@vmatchoptions,"-seedlength");
      } else
      {
        push(@vmatchoptions,$currentarg);
      }
      $argnum++;
      if($argnum >= $argcount-1 || $argv[$argnum] =~ m/^\-/)
      {
        print STDERR "$program: missing argument for option \"$currentarg\"\n";
        exit 1;
      }
      if($currentarg eq '-best')
      {
        $bestoption = 1;
      }
      push(@vmatchoptions,$argv[$argnum]);
    } elsif($currentarg eq '-lw')
    {
      $argnum++;
      if($argnum >= $argcount-1 || $argv[$argnum] =~ m/^\-/)
      {
        print STDERR "$program: missing argument for option \"$currentarg\"\n";
        exit 1;
      }
      $linewidth = $argv[$argnum];
      if($linewidth <= 0)
      {
        print STDERR "$program: illegal argument \"$linewidth\" ";
        print STDERR "to option \"-lw\"\n";
        exit 1;
      }
    } elsif($currentarg eq '-h' || $currentarg eq '-e')
    {
      push(@vmatchoptions,$currentarg);
      $argnum++;
      if($argnum >= $argcount-1 || $argv[$argnum] =~ m/^-/)
      {
        push(@vmatchoptions,"4");
      } else
      {
        push(@vmatchoptions,$argv[$argnum]);
      }
    } elsif($currentarg eq '-allmax')
    {
      $allmaxoption = 1;
      push(@vmatchoptions,"-allmax");
    } elsif($currentarg eq '-s')
    {
      $stringoption = 1;
    } elsif($currentarg eq '-iub')
    {
      $doiub = 1;
    } elsif($currentarg eq '-nodistance')
    {
      push(@vmatchoptions,"-nodist");
    } elsif($currentarg eq '-noevalue' || $currentarg eq '-i')
    {
      push(@vmatchoptions,$currentarg);
    } else
    {
      if(exists $notsupported{$currentarg})
      {
        print STDERR "$program: repfind option \"$currentarg\" is ";
        print STDERR "not supported\n";
      } else
      {
        print STDERR "$program: illegal option \"$currentarg\"\n";
      }
      exit 1;
    }
  }
  if($argnum == $argcount-1)
  {
    if($argv[$argnum] =~ m/^-/)
    {
      print STDERR "$program: last argument must be filename, ";
      print STDERR "not beginning with \"-\"\n";
      exit 1;
    } 
  } elsif($argnum > $argcount-1)
  {
    print STDERR "$program: missing last argument\n";
    exit 1;
  } elsif($argnum < $argcount-1)
  {
    print STDERR "$program: superfluous argument $argv[$argnum]\n";
    exit 1;
  }
  if(not @vmatchoptions)
  {
    print STDERR "$program: at least one option is required\n";
    exit 1;
  }
  if((not $bestoption) && (not $allmaxoption))
  {
    push(@vmatchoptions,"-best");
    push(@vmatchoptions,"50");
  }
  if($stringoption)
  {
    push(@vmatchoptions,"-s");
    if($linewidth > 0)
    {
      push(@vmatchoptions,$linewidth);
    }
    if($doiub)
    {
      push(@vmatchoptions,"abbreviub");
    }
  }
  push(@vmatchoptions,"-noscore");
  push(@vmatchoptions,"-noidentity");
  push(@vmatchoptions,"-absolute");
  return \@vmatchoptions;
}

#-r           compute maximal reverse repeats
#-c           compute maximal complemented repeats
#-hrate       search for repeats up to the given hamming error rate
#-erate       search for repeats up to the given edit error rate
#-o           specify the maximal output size
#-b           binary output
#-warn        warn about non-base symbols in input file
#-iw          ignore wildcard characters
#-mem         show memory usage
