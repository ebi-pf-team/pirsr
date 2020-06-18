#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

use Smart::Comments;

use FindBin qw($Bin);
use lib "$Bin";


use PIRSR;





use Smart::Comments;

my $rules_file;
my $fasta_file;
my $target_file;
# my $hmm_folder;


GetOptions(
  'help'     => sub { pod2usage( -verbose => 1 ) },
  'man'      => sub { pod2usage( -verbose => 2 ) },

  'rules=s'  => \$rules_file,
  'fasta=s'  => \$fasta_file,
  # 'hmms=s'   => \$hmm_folder,

  'target=s' => \$target_file,

) or pod2usage(2);

### $rules_file



my $fasta = PIRSR::read_template_fasta($fasta_file);
## $fasta
my $rules = PIRSR::get_rules($rules_file);
## $rules


my $align = PIRSR::align_template($rules, $fasta);
## $align
