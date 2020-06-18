use strict;
use warnings;
use Test::More;
use Test::Exception;
use Test::Warnings;
use Test::Warn;

use FindBin '$Bin';
use lib "$Bin/../";

use Smart::Comments;

use_ok('PIRSR');




# fasta reading and processing step

my $sample_fasta = "$Bin/data/fasta/sr_tp.seq.min";

my $fasta = PIRSR::read_template_fasta($sample_fasta);
## $fasta

is(scalar keys %{$fasta} , 16,  "16 sequences found on fasta file");

my $example_seq_obj = $fasta->{'P19235'};

my $expected_P19235 = {
  file => "$Bin/data/fasta/P19235.fa",
  hmm => 'PIRSR001959-2',
  seq => 'MDHLGASLWPQVGSLCLLLAGAAWAPPPNLPDPKFESKAALLAARGPEELLCFTERLEDLVCFWEEAASAGVGPGNYSFSYQLEDEPWKLCRLHQAPTARGAVRFWCSLPTADTSSFVP#LELRVTAASGAPRYHRVIHINEVVLLDAPVGLVARLADESGHVVLRWLPPPETPMTSHIRYEVDVSAGNGAGSVQRVEILEGRTECVLSNLRGRTRYTFAVRARMAEPSFGGFWSAWSEPVSLLTPSDLDPLILTLSLILVVILVLLTVLALLSHRRALKQKIWPGIPSPESEFEGLFTTHKGNFQLWLYQNDGCLWWSPCTPFTEDPPASLEVLSERCWGTMQAVEPGTDDEGPLLEPVGSEHAQDTYLVLDKWLLPRNPPSEDLPGPGGSVDIVAMDEGSEASSCSSALASKPSPEGASAASFEYTILDPSSQLLRPWTLCPELPPTPPHLKYLYLVVSDSGISTDYSSGDSQGAQGGLSDGPYSNPYENSLIPAAEPLPPSYVACS'
};

is_deeply($example_seq_obj, $expected_P19235, "P19235 fasta sequence data is correct");



# rules file reading and processing step

my $sample_rule = "$Bin/data/test1.uru";

my $rules = PIRSR::get_rules($sample_rule);
## $rules

is(scalar keys %{$rules} , 13,  "13 rules found on rules file");


my $expected_PIRSR001341_5 = {
  AC => 'PIRSR001341-5',
  DC => 'Domain',
  Feature => {
               from => 'P00862',
               model => 'SRHMM001341-5'
             },
  Groups => {
              1 => [
                     {
                       condition => 'S',
                       desc => 'Histidine decarboxylase beta chain.',
                       end => '82',
                       group => '1',
                       label => 'CHAIN',
                       start => 'Nter'
                     }
                   ],
              2 => [
                     {
                       condition => 'S',
                       desc => 'Histidine decarboxylase alpha chain.',
                       end => 'Cter',
                       group => '2',
                       label => 'CHAIN',
                       start => '83'
                     }
                   ]
            },
  Related => 'None',
  Scope => [
             'Bacteria'
           ],
  Size => 'unlimited',
  TR => 'PIRSF; PIRSF001341; -; 1'
};

is_deeply($rules->{'PIRSR001341-5'}, $expected_PIRSR001341_5, "Rule PIRSR001341-5 imported correctly");


my $aligned_rules;
warning_like{ $aligned_rules = PIRSR::align_template($rules, $fasta) } qr/rule PIRSR001341-5/, 'warning for Nter/Cter start/end position';
## $aligned_rules

# my $aligned_rules = PIRSR::align_template($rules, $fasta);
# ## $aligned_rules

is(scalar keys %{$aligned_rules} , 12,  "12 rules now found on rules file");

ok(!defined $aligned_rules->{'PIRSR001341-5'} ,  "rule PIRSR001341-5 has been dropped");


my $expected_PIRSR000008_2 = {
  AC => 'PIRSR000008-2',
  DC => 'Domain',
  Feature => {
               from => 'P29899',
               model => 'SRHMM000008-2'
             },
  Groups => {
              1 => [
                     {
                       condition => 'H',
                       desc => 'Iron (heme axial ligand).',
                       end => 83,
                       group => '1',
                       hmmEnd => 87,
                       hmmStart => 87,
                       label => 'METAL',
                       start => 83
                     }
                   ]
            },
  Related => 'None',
  Scope => [
             'Bacteria'
           ],
  Size => 'unlimited',
  TR => 'PIRSF; PIRSF000008; -; 1'
};

is_deeply($aligned_rules->{'PIRSR000008-2'}, $expected_PIRSR000008_2, "Rule PIRSR000008-2 aligned correctly");




done_testing();








# my $dat_file = "$Bin/../../../../data/pirsf/3.10/pirsf.dat";

# my ($pirsf_data, $children) = PIRSF::read_pirsf_dat($dat_file);
# my $pirsf_count = scalar keys %{$pirsf_data};

# is(scalar keys %{$pirsf_data} , 4,  "4 PIRSF ids found on pirsf_dat");

# my $example_pirsf_dat->{'PIRSF500175'} = $pirsf_data->{'PIRSF500175'};

# my $expected_ex_pirsf_dat = {
#     blast => 0,
#     meanL => '427.5',
#     meanS => '570.346666666667',
#     minS => '442.9',
#     name => 'Glutamyl-tRNA(Gln) amidotransferase, subunit D',
#     stdL => '16.8926631697177',
#     stdS => '46.7185242248408'
# };

# is_deeply($example_pirsf_dat->{'PIRSF500175'}, $expected_ex_pirsf_dat, "Random PIRSF id data is correct");
# is($pirsf_data->{'PIRSF500175'}->{'children'}, undef, "No children for PIRSF500175");

# my $expected_children = {
#     PIRSF500175 => 1,
#     PIRSF500176 => 1
# };

# is_deeply($pirsf_data->{'PIRSF001220'}->{'children'}, $expected_children, "Correct children for PIRSF001220");


# my $input = "$Bin/data/test.fasta";
# my $single_input = "$Bin/data/UPI000000078D.fasta";
# my $dominput = "$Bin/data/UPI000000078D.hmmscan.domtblout";
# my $sf_hmm = "$Bin/../../../../data/pirsf/3.10/sf_hmm_all";
# my $hmmer_path = "$Bin/../../../hmmer/hmmer3/3.1b1";
# # my $cpus = 1;
# my $mode = "hmmscan";
# my $i5_tmpdir = '/tmp';
# my $pirsf_bin = "$Bin/../";


# my @targets = sort keys %{PIRSF::read_fasta($dat_file)};
# my @expected_targets = qw/PIRSF001220 PIRSF001789 PIRSF500175 PIRSF500176/;

# is( @targets, @expected_targets, "Correct targets found on input fasta file");

# my $dom_results = PIRSF::read_dom_input($dominput);
# my $exp_results = ['NGFB                 PIRSF001789   253 UPI000000078D        -            226   1.8e-94  304.0   2.9   1   1   6.5e-95     2e-94  303.9   2.9    36   252     5   225     1   226 0.94 -'];

# is($dom_results->[0], $exp_results->[0], "Correct results read from dom file");


# my $UPI000000078D_i5_expected = <<UPI000000078D_i5;
# 226	1	PIRSF001789	UPI000000078D	1.8e-94	..	225	5	303.9	304.0	2.9	6.5e-95	2e-94	226	1	0.94	2.9
# UPI000000078D_i5

# my $UPI000000078D_pi_expected = <<UPI000000078D_pirsf;
# Query sequence: UPI000000078D matches PIRSF001789: Nerve growth factor, subunit beta
#    #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
#    1 !  303.9   2.9   6.5e-95     2e-94      36     252 ..       5     225 ..       1     226 [] 0.94
# UPI000000078D_pirsf


# my $single_prot_i5_cmd = "perl ${pirsf_bin}/pirsf.pl --fasta $single_input --domtbl $dominput -path $hmmer_path --hmmlib $sf_hmm -dat $dat_file --mode hmmscan --outfmt i5";
# my $single_prot_pi_cmd = "perl ${pirsf_bin}/pirsf.pl --fasta $single_input --domtbl $dominput -path $hmmer_path --hmmlib $sf_hmm -dat $dat_file --mode hmmscan --outfmt pirsf";

# my $single_prot_i5 = qx/$single_prot_i5_cmd/;
# my $single_prot_pi = qx/$single_prot_pi_cmd/;

# is($single_prot_i5, $UPI000000078D_i5_expected, "expected i5 output for UPI000000078D");
# is($single_prot_pi, $UPI000000078D_pi_expected, "expected pirsf output for UPI000000078D");

# my $UPI000000078D_hmmsearch_expected = <<UPI000000078D_hmmsearch;
# Query sequence: UPI000000078D matches PIRSF001789: Nerve growth factor, subunit beta
#    #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
#    1 !  303.9   2.9   6.5e-95   6.5e-95      36     252 ..       5     225 ..       1     226 [] 0.94
# UPI000000078D_hmmsearch

# my $single_prot_hmmsearch_cmd = "perl ${pirsf_bin}/pirsf.pl --fasta $single_input -path $hmmer_path --hmmlib $sf_hmm -dat $dat_file --mode hmmsearch --outfmt pirsf";
# my $single_prot_hmmsearch = qx/$single_prot_hmmsearch_cmd/;

# is($single_prot_hmmsearch, $UPI000000078D_hmmsearch_expected, "expected pirsf output for UPI000000078D");


# my $test_prot_cmd = "perl ${pirsf_bin}/pirsf.pl --fasta $input -path $hmmer_path --hmmlib $sf_hmm -dat $dat_file --mode hmmscan --outfmt pirsf";
# my $test_prot_pirsf = qx/$test_prot_cmd/;

# my $georgetown_testing = `for file in data/test/*
#                           do
#                             perl pirsf.2017.pl "\$file"
#                           done`;

# is($test_prot_pirsf, $georgetown_testing, "results of test.fasta match those of georgetown pirsf.2017.pl");


# done_testing();
