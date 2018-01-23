#!/usr/bin/perl
##############################################################################
# SCRIPT NAME:	g2gl_ct.pl
# DESCRIPTION:	NMDP HapLogic file harness genotype input to genotype-list
#
# DATE WRITTEN: 2017-08-09
# WRITTEN BY:   Loren Gragert
#
##############################################################################
use strict;    # always
use warnings;  # or else
# use MAC;
use Getopt::Std;

my %opts;
getopts('i:o:t:', \%opts);
die "$0: -i inputfile -o outputfile" unless defined $opts{i} && defined $opts{o} && defined $opts{t};

# set to absolute path of your installation
my $top = "../..";


# set to input genotypes file
my $gfile = $opts{i};
my $ofile = $opts{o};
my $subject_type = $opts{t};

# load AC

  my $ac_file = "./data/alpha.v3.txt";
  my %AC;
  open (ACFILE,"$ac_file") or die "$!: Missing $ac_file";
  while(<ACFILE>) {
    chomp;
    my ($asterisk, $code, $ac_list) = split /\t/,$_;
    $AC{$code} = $ac_list;
  }
  close ACFILE;


##############################################################################
# parse input file
##############################################################################
open OFILE, ">$ofile" or die "$!: $ofile";
open GFILE, $gfile or die "$!: $gfile";
while(<GFILE>) {
  chomp;
  #D000001%A*11:XX+A*02:XX^C*UUUU+C*UUUU^B*13:XX+B*44:XX^DRB1*07:XX+DRB1*01:XX^DQB1*UUUU+DQB1*UUUU
  my ($recip_id, $recip_phen_seq, $recip_detail_race, $recip_ethnicity, 
    $recip_broad_race, $recip_A_dna, $recip_A_1, $recip_A_2, $recip_B_dna, 
    $recip_B_1, $recip_B_2, $recip_DRB1_dna, $recip_DRB1_1, $recip_DRB1_2, 
    $recip_C_dna, $recip_C_1, $recip_C_2, $recip_DQB1_dna, $recip_DQB1_1, 
    $recip_DQB1_2, $recip_DRB3_dna, $recip_DRB3_1, $recip_DRB3_2, 
    $recip_DRB4_dna, $recip_DRB4_1, $recip_DRB4_2, $recip_DRB5_dna, 
    $recip_DRB5_1, $recip_DRB5_2, $donor_id,  
    $donor_detail_race, $donor_ethnicity, $donor_broad_race, 
    $donor_A_dna, $donor_A_1, $donor_A_2, $donor_A_pdtl, 
    $donor_B_dna, $donor_B_1, $donor_B_2, $donor_B_pdtl, 
    $donor_DRB1_dna, $donor_DRB1_1, $donor_DRB1_2, $donor_DRB1_pdtl,
    $donor_C_dna, $donor_C_1, $donor_C_2, $donor_C_pdtl, 
    $donor_DQB1_dna, $donor_DQB1_1, $donor_DQB1_2, $donor_DQB1_pdtl,
    $donor_DRB3_dna, $donor_DRB3_1, $donor_DRB3_2, $donor_DRB3_pdtl, 
    $donor_DRB4_dna, $donor_DRB4_1, $donor_DRB4_2, $donor_DRB4_pdtl, 
    $donor_DRB5_dna, $donor_DRB5_1, $donor_DRB5_2, $donor_DRB5_pdtl) 
    = split /\t/;
  
    # print "A1: $donor_A_1 A2: $donor_A_2 DNA: $donor_A_dna\n";

    my @gl_array = ();
    my ($gl_A,$gl_C,$gl_B,$gl_DRB1,$gl_DQB1);
    if ($subject_type eq "R") {
      $gl_A = loc_gl($recip_id,"A",$recip_A_1,$recip_A_2,$recip_A_dna);
      $gl_C = loc_gl($recip_id,"C",$recip_C_1,$recip_C_2,$recip_C_dna);
      $gl_B = loc_gl($recip_id,"B",$recip_B_1,$recip_B_2,$recip_B_dna);
      $gl_DRB1 = loc_gl($recip_id,"DRB1",$recip_DRB1_1,$recip_DRB1_2,$recip_DRB1_dna);
      $gl_DQB1 = loc_gl($recip_id,"DQB1",$recip_DQB1_1,$recip_DQB1_2,$recip_DQB1_dna);
    } else {
      $gl_A = loc_gl($donor_id,"A",$donor_A_1,$donor_A_2,$donor_A_dna);
      $gl_C = loc_gl($donor_id,"C",$donor_C_1,$donor_C_2,$donor_C_dna);
      $gl_B = loc_gl($donor_id,"B",$donor_B_1,$donor_B_2,$donor_B_dna);
      $gl_DRB1 = loc_gl($donor_id,"DRB1",$donor_DRB1_1,$donor_DRB1_2,$donor_DRB1_dna);
      $gl_DQB1 = loc_gl($donor_id,"DQB1",$donor_DQB1_1,$donor_DQB1_2,$donor_DQB1_dna);     
    }

    if ($gl_A ne "SKIP") { push @gl_array, $gl_A; }
    if ($gl_C ne "SKIP") { push @gl_array, $gl_C; }
    if ($gl_B ne "SKIP") { push @gl_array, $gl_B; }
    if ($gl_DRB1 ne "SKIP") { push @gl_array, $gl_DRB1; }
    if ($gl_DQB1 ne "SKIP") { push @gl_array, $gl_DQB1; }



    my $gl_string = join ('^', @gl_array);
    $gl_string=~s/HLA\-//g;
    # print "$id$gl_string\n";
    if ($subject_type eq "R") {
      print OFILE join ('%', $recip_id, $gl_string, "\n");
    }
    else {
      print OFILE join ('%', $donor_id, $gl_string, "\n");
    }
}
close OFILE;
close GFILE;

exit 0;


# expandTyp is too slow because it hits MAC service over the web
sub expandTyp {
  my ($id, $t) = @_;
  my $who_typ = "HLA-".$t;
  my @al;

  # expand
  my $rc = MAC::decode($who_typ, \@al);
  if ($rc ne 200) {
    warn "decode failed with code: $rc for $who_typ for id: $id";
    next;
  }
  return @al;
}

sub loc_gl {
  my ($id, $loc, $typ1, $typ2, $dna) = @_;
  # print "loc_gl $id $loc $typ1 $typ2 $dna\n";
  if ($typ2 eq "") { $typ2 = $typ1; }
  if ($typ1 eq "" || $typ2 eq "" || 
    $typ1 eq "NEW" || $typ2 eq "NEW" || $dna ne "1") {
    return "SKIP"; 
  }
  my $loctyp1 = join ('*',$loc,$typ1);
  my $loctyp2 = join ('*',$loc,$typ2);
  # my @al1 = expandTyp($id, $loctyp1);
  # my @al2 = expandTyp($id, $loctyp2);
  my @al1 = expand_AC(\%AC, $loctyp1);
  my @al2 = expand_AC(\%AC, $loctyp2);
  my $al1 = join ('/',@al1);
  my $al2 = join ('/',@al2);
  my $gl = join ('+', $al1, $al2);
  return $gl;
}

sub expand_AC {
  my ($AC,$loctyp) = @_;
  my ($loc,$typ) = split /\*/,$loctyp;
  my ($family,$allele) = split /:/,$typ;
  my $MAC;
  my @MAC_loctyp;
  if (exists $AC{$allele}) { 
    $MAC = $AC{$allele};
  } else {
    push @MAC_loctyp, $loctyp;
    return @MAC_loctyp; # high res allele
  }
  my @MAC = split /\//,$MAC;

  foreach my $allele (@MAC) {
    # check to see if MAC is 4-digit
    $loctyp = "";
    if ($allele =~ m/\:/) {
      $loctyp = $loc . "*". $allele;
    }
    else {
      $loctyp = $loc . "*". $family . ":" . $allele;
    }
    push @MAC_loctyp, $loctyp;
  }
  return @MAC_loctyp;
}

