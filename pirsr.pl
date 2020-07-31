#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

use Smart::Comments;


use File::Copy;
use FindBin qw($Bin);
use lib "$Bin";


use PIRSR;






use Smart::Comments;

my $rules_file;
my $fasta_file;
my $target_file;


my $data_folder;

my $hmm_folder;
my $template_folder;
my $rule_folder;



my $query_file;


GetOptions(
  'help'     => sub { pod2usage( -verbose => 1 ) },
  'man'      => sub { pod2usage( -verbose => 2 ) },



  'data=s'   => \$data_folder,

  'hmms=s'   => \$hmm_folder,
  'templates=s' => \$template_folder,
  'rules=s'  => \$rule_folder,



  'rules=s'  => \$rules_file,
  'fasta=s'  => \$fasta_file,


  'query=s'  => \$query_file,



  'target=s' => \$target_file,

) or pod2usage(2);

## $rules_file

$data_folder =~ s/\/$//;
### $data_folder

# set the default hmm_folder
if (!$hmm_folder) {
    $hmm_folder = "${data_folder}/sr_hmm";
}

# set the default template_folder
if (!$template_folder) {
    $template_folder = "${data_folder}/sr_tp";
}

# set the default rule_folder and move the uru file there
if (!$rule_folder) {
    $rule_folder = "${data_folder}/sr_uru";
    if (!-d $rule_folder) {
        mkdir($rule_folder);
    }

    if (!-e "${rule_folder}/PIRSR.uru") {
        move("${data_folder}/PIRSR.uru", "${rule_folder}/PIRSR.uru");
    }
}


my $pirsr = PIRSR->new(
    template_folder => $template_folder,
    hmm_folder      => $hmm_folder,
    rule_folder     => $rule_folder,
);

### $pirsr

my $bla = $pirsr->process_data();
### $bla

# my $template_result = $pirsr->process_template_folder();
# # $template_result

# my $hmm_result = $pirsr->process_hmm_folder();
# # $hmm_result


# my $rule_result = $pirsr->process_rule_folder();
# ## $rule_result





my $search = $pirsr->run_query($query_file);
### $search





# my $fasta = PIRSR::read_template_fasta($fasta_file);
# ## $fasta


# my $rules = PIRSR::get_rules($rules_file);
# ## $rules


# my $align = PIRSR::align_template($rules, $template_folder);
# ### $align


# ## $result




# my $search = PIRSR::run_search($query_file, $hmm_folder);
# ### $search