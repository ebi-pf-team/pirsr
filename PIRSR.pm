package PIRSR;

use strict;
use warnings;

use JSON;

use Moose;
use namespace::autoclean;

use Bio::Pfam::HMM::HMMResultsIO;
use File::Temp qw/tempfile/;


use Smart::Comments;



has 'template_folder' => (
    is       => 'rw',
    isa      => 'Str',
    required => 1
);

has 'hmm_folder' => (
    is       => 'rw',
    isa      => 'Str',
    required => 1
);

has 'rule_folder' => (
    is       => 'rw',
    isa      => 'Str',
    required => 1
);



sub process_data {
    my ($self) = @_;

    $self->process_template_folder();
    $self->process_hmm_folder();
    $self->process_rule_folder();

    return;
}




sub process_template_folder {
    my ($self) = @_;

    my $template_folder = $self->template_folder;

## $template_folder

    # parse the sr_tp.seq a seq block at a time
    do {
        local $/ = ">";
        open (my $in, '<', "$template_folder/sr_tp.seq") or die "Failed top open $template_folder/sr_tp.seq file: $!\n";
        while (my $block = <$in>) {

            $block =~ s/\n?>?\z//;
            next if !$block;


            # the first line (seq identifier) has the format: 'Q59771   PIRSR000188-2'
            if ($block =~ /\A((\w+)\s+\w+-\d+)\n/) {
                # the whole line is the $header, the first bit is the $prot_id
                my $header = $1;
                my $prot_it = $2;

                # now we have $header and $prot_it, get rid of the first line and grab the seq
                $block =~ s/\A(\w+)\s+(\w+-\d+)\n//;
                $block =~ s/\n//g;

                my $seq = $block;
## $seq
                # $data->{$prot_it}->{'seq'} = $block;

                open (my $out, '>', "${template_folder}/${prot_it}.fa" ) or die "Failed top open ${prot_it}.fa file: $!\n";

                print $out ">${header}\n${seq}\n";

                close($out) or die "Failed to close ${prot_it}.fa\n";

                # $data->{$prot_it}->{'file'} = "${template_folder}/${prot_it}.fa";

            } else {
                die "ERROR: failed to parse fasta block: \"$block\"\n";
            }
            ## $block
            # my $stuff = _process_rules_block($block);
            ## $stuffs
        }
        close($in) or die "Failed to close $template_folder/sr_tp.seq file\n";
    };



    return;
}


sub process_hmm_folder {
    my ($self) = @_;

    my $hmm_folder = $self->hmm_folder;
## $hmm_folder

    # opendir(my $dir, $hmm_folder) or die "Could not open folder '$hmm_folder': $!\n";

    open(my $out_fh, '>', "$hmm_folder/sr_hmm.hmm") or die "Can't open '$hmm_folder/sr_hmm.hmm': $!\n";

    foreach my $hmm_file ( glob("$hmm_folder/PIRSR*.hmm") ) {
## $hmm_file
        open(my $in_fh, '<', "$hmm_file") or die "Can't open '$hmm_file': $!\n";
        # slurp individual hmm file
        my $hmm = do { local $/; <$in_fh> };
        close($in_fh) or die "Can't close '$hmm_file' after reading: $!\n";

        # fix hmm name
        my $hmm_name = $hmm_file;
        # remove path
        $hmm_name =~ s/.*\///;
        # remove file type
        $hmm_name =~ s/\.hmm//;
        # replace with hmm/rule name
        $hmm =~ s/NAME\s+.*?\n/NAME  ${hmm_name}\n/;

        # print it out to library file
        print $out_fh $hmm;

    }

    close($out_fh) or die "Can't close '$hmm_folder/sr_hmm.hmm' after writing: $!\n";

    # create auxfiles
    my $cmd = "hmmpress ${hmm_folder}/sr_hmm.hmm";

    foreach my $ext (qw(h3p h3m h3f h3i)) {
        if (!-e "${hmm_folder}/sr_hmm.hmm.${ext}") {
            # Looks like the hmm database is not pressed
            warn "Running hmmpress on ${hmm_folder}/sr_hmm.hmm\n";
            system("hmmpress ${hmm_folder}/sr_hmm.hmm") and die "Could not run hmmpress: $!\n";
            last;
        }
    }

    return;
}


sub process_rule_folder {
    my ($self) = @_;
## $self
    my $rule_folder = $self->rule_folder;
    ## $rule_folder

    do {
        local $/ = "//\n";
        open (my $in, '<', "${rule_folder}/PIRSR.uru") or die "Failed top open ${rule_folder}/PIRSR.uru file: $!\n";
        while (my $block = <$in>) {
            my $rule_hash = _parse_rules($block);
            ## $rule_hash

            $rule_hash = $self->align_template($rule_hash);

            my $rule_acc = $rule_hash->{'AC'};
## $rule_acc

            open (my $out, '>', "${rule_folder}/${rule_acc}.json" ) or die "Failed top open ${rule_acc}.json file: $!\n";

            my $json_uru = to_json( $rule_hash, { pretty => 1 } );
            print $out $json_uru;

            close($out) or die "Failed to close ${rule_acc}.json\n";
        }
        close($in) or die "Failed to close ${rule_folder}/PIRSR.uru: $!\n";
    };



    return;
}






sub _parse_rules {
    my ($block) = @_;
# print ">>> ENTERING _parse_rules\n";
    ## $block
    my @rule = split(/\n/,$block);
    ## @rule

    my $rule_hash = {};

    while (my $line = shift @rule) {
# print "\tON WHILE: \"$line\"\n";

        if ($line =~ /\A(\*\*|Comments|XX|end case|\/\/)/ ) {
# print " skipping \"$line\"\n";
        } elsif ($line =~ /\AAC\s+(.+);/ ) {
            $rule_hash->{'AC'} = $1;
        } elsif ($line =~ /\ADC\s+(.+);/ ) {
            $rule_hash->{'DC'} = $1;
        } elsif ($line =~ /\ATR\s+(.+);/ ) {
            $rule_hash->{'TR'} = $1;
        } elsif ($line =~ /\ASize:\s+(.+);/ ) {
            $rule_hash->{'Size'} = $1;
        } elsif ($line =~ /\ARelated:\s+(.+);/ ) {
            $rule_hash->{'Related'} = $1;
        } elsif ($line =~ /\AScope:/ ) {
            while ($rule[0] =~ /\A\s+(.*)/ ) {
                push( @{$rule_hash->{'Scope'}}, $1);
                shift @rule;
            }
        } elsif ($line =~ /\Acase\s\<Feature:(\S+)\>\z/) {
            my $feature->{'model'} = $1;
# print "\tENTERING case \"$line\"\n";
            my $group;
            while ($rule[0] !~ /\Aend case/ ) {
                $line = shift @rule;
# print "\t\t \"$line\"\n";


                if ($line =~ /\AFT\s+From\: (\S+)\z/) {
                    $feature->{'from'} = $1;
# print "\t\t\t feature from\n";#
                } elsif ($line =~ /\AFT\s{3}(\S+)\s+(\w+)\s+(\w+)(\s+(.*))?\z/) {
# print "\t\t\t start of group\n";#
                        undef $group;
                        $group->{'label'} = $1;
                        $group->{'start'} = $2;
                        $group->{'end'} = $3;
                        $group->{'desc'} = $5;

                } elsif ($line =~ /\AFT\s+Group\:\s+(\d+)\;\s+Condition\:\s+(.*)\z/) {
# print "\t\t\t group conditio \n";#
                    $group->{'group'} = $1;
                    $group->{'condition'} = $2;
## $group
                    push( @{$rule_hash->{'Groups'}->{$group->{'group'}}}, $group);

                } elsif ($line =~ /\AFT\s+(.*)\z/) {
# print "\t\t\t desc continuation\n";#
                    $group->{'desc'} .= " $1";


                } else {
                    die "ERROR: failed to parse Group line: \"$line\"\n";
                }


            }
            $rule_hash->{'Feature'} = $feature;

        } elsif ($line =~ /\Acase\s\<Feature:(\S+)\>.+/){
# print " ignoring \"$line\"\n";
            while ($rule[0] !~ /\Aend case/ ) {
                $line = shift @rule;
# print " ignoring \"$line\"\n";
            }
        } else {
            die "ERROR: failed to parse rules line: \"$line\"\n";
        }
    }

## $rule_hash

    return $rule_hash;
}






sub align_template {
    my ($self, $rule) = @_;


    my $prot_id = $rule->{'Feature'}->{'from'};
    ## $prot_id

    my $fasta_file = "$self->{template_folder}/${prot_id}.fa";
## $fasta_file
    open (my $in, '<', "$fasta_file") or die "Failed top open $fasta_file file: $!\n";



    my ($prot_model, $prot_seq);

    my $fasta = <$in>;
    if ($fasta =~ m/.*\t(.*)\n(.*)/) {
        $prot_model = $1;
        $prot_seq = $2;
    }
    ## $prot_model
    ## $prot_seq

# die;

    my $stockholm = `hmmalign data/sr_hmm/${prot_model}.hmm ${fasta_file}`;
## $stockholm
    my $alignment_str;
    while ($stockholm =~ /\n${prot_id}\s+([^\n]*)/g) {
        $alignment_str .= $1;
    }
    ## $alignment_str


    my $alignment = align_map($alignment_str);
    ## @map;
## $alignment



## $rules


    foreach my $grp (keys %{$rule->{'Groups'}}) {
            ## $grp

        for my $pos (0 .. $#{$rule->{'Groups'}->{$grp}} ) {
            ## $pos
            # my $both = 0;


            # process start
            if ($rule->{'Groups'}->{$grp}->[$pos]->{'start'} eq 'Nter') {
                warn "rule $rule->{AC}, group $grp, pos $pos: Start is Nter.\n";
                $rule->{'Groups'}->{$grp}->[$pos]->{hmmStart} = $rule->{'Groups'}->{$grp}->[$pos]->{'start'};
            } else {
                # if ($alignment->[$rule->{'Groups'}->{$grp}->[$pos]->{'start'}] eq "-") {
                #     warn $prot_model." ".$prot_id."/".$rule->{'Groups'}->{$grp}->[$pos]->{'start'}." is non-match state.\n";
                #     $both++;
                #     ## second_if: $rule->{'Groups'}->{$grp}->[$pos]
                # }
                $rule->{'Groups'}->{$grp}->[$pos]->{hmmStart} = $alignment->[$rule->{'Groups'}->{$grp}->[$pos]->{'start'}];
            }


            # process end
            if ($rule->{'Groups'}->{$grp}->[$pos]->{'end'} eq 'Cter') {
                warn "rule $rule->{AC}, group $grp, pos $pos: End is Cter.\n";
                $rule->{'Groups'}->{$grp}->[$pos]->{hmmEnd} = $rule->{'Groups'}->{$grp}->[$pos]->{'end'};
                ## $rule
            } else {
                # if ($alignment->[$rule->{'Groups'}->{$grp}->[$pos]->{'end'}] eq "-") {
                #     warn $prot_model." ".$prot_id."/".$rule->{'Groups'}->{$grp}->[$pos]->{'end'}." is non-match state.\n";
                #     $both++;

                #     ## third_if: $rule->{'Groups'}->{$grp}->[$pos]
                # }
                $rule->{'Groups'}->{$grp}->[$pos]->{hmmEnd} = $alignment->[$rule->{'Groups'}->{$grp}->[$pos]->{'end'}];
            }



            #Some of the disulphide bridges lack on of the positions. However, the two current failing
            #rules means that we are completely dependent on a length.    This is not ideal, but we
            #can try.
            # if ($both == 2) {
            #     warn "No anchor point for $prot_model ".$prot_id."/".$rule->{'Groups'}->{$grp}->[$pos]->{'start'}."-".$rule->{'Groups'}->{$grp}->[$pos]->{'end'}."\n";
            # }


            ## $pos
            ## $rule

        }

## $rules

    }


# $rules




    return $rule;
}





# returns array of seq->hmm positions for a given alignment string
sub align_map {
    my ($alignment) = @_;

    my $seq_pos = 1;
    my $hmm_pos = 1;
    my @map;

    foreach my $pos (split(//, $alignment)) {
        #print STDERR "$pos";
        if ($pos eq "-") {
            #delete position
            $hmm_pos++;
        } elsif ($pos eq ".") {
            #Alignment insert, skip
        } elsif ($pos =~ /[A-Z]/) {
            #Match position
            $map[$seq_pos] = $hmm_pos;
            $seq_pos++;
            $hmm_pos++;
        } elsif ($pos =~ /[a-z]/) {
            #hmm insert pos

            # non-match state
            # $map[$seq_pos] = '-';
            # Use valid hmm pos
            $map[$seq_pos] = $hmm_pos;
            $seq_pos++;
        }
    }

    return \@map;
}







sub run_query {
    my ($self, $query_file) = @_;

## $query_file

    my %query_rules;


    my (undef, $out) = tempfile(
        DIR => '/tmp',
        UNLINK => 1
    );

    ## $out

    # my $out = 'result.out';

    my $hmm_library = $self->{hmm_folder} . '/sr_hmm.hmm';
    ## $hmm_library

    my $cmd = "hmmscan --notextw -o $out $hmm_library $query_file";
    ## $cmd
    system($cmd) && die qq(Failed to run "$cmd");
    ## $out
# my $hmmRes = Bio::Pfam::HMM::HMMResultsIO->new;
# ## $hmmRes
    my $res_obj = Bio::Pfam::HMM::HMMResultsIO->new->parseMultiHMMER3($out);
## $res_obj



    # print "Sequence acc\tSeq start\tSeq end\tRule acc\tFeature group\tTrigger\tModel\tTemplate\tTemplate start\tTemplate end\tTaxonomic scope\tFeature type\tFeature description\n";
    foreach my $query_match (@{$res_obj}) {
## $query_match
        my $query_id = $query_match->{'seqName'};
        next if !$query_id;
## $query_id
        $query_rules{$query_id} = {} unless $query_rules{$query_id};

        # loop over HMMs - i.e. Unirule profiles
        # my $rule_id = $query_match->{hmmName};
        # ## $rule_id
        # $rule_id =~ s{\.}{\-}; # FIXME - check what format we are using for HMMs - may not need this


        foreach my $target_match (@{$query_match->{'units'}}) {
            ## $target_match
            # loop over alignments
            my $rule_id = $target_match->{'name'};
            if ($query_rules{$query_id}{$rule_id}) {
                next;
            };
            ## $rule_id
            open(my $in, '<', "$self->{rule_folder}/${rule_id}.json") or die "Failed to open $self->{rule_folder}/${rule_id}.json file: $!\n";
            my $json_string = do { local $/; <$in> };
            ## $json_string
            close($in) or die "Failed to close /${rule_id}.json: $!\n";

            my $rule = from_json($json_string);
            ## $rule

            my $tname = $target_match->{name}; # the sequence ID/accession
            my $hmm_seq = $target_match->{hmmalign}{hmm};
            my $query_seq = $target_match->{hmmalign}{seq};
            my $hmm_from = $target_match->{hmmFrom};
            my $seq_from = $target_match->{seqFrom};
            my $map = map_hmm_to_seq($hmm_from, $hmm_seq, $query_seq);
            ## $map
            # while (my ($group_id, $group_rules) = each %{$rules->{$rule_id}}) {
                # Loop over the groups in the rule - check and print
                # print_group_matches($rule_id, $group_id, $tname, $query_seq, $seq_from, $map, $group_rules, $fh, 0);
            # }


            my $pass = 0;

            foreach my $grp (sort keys %{$rule->{'Groups'}}) {
## $grp

                my $pass_count = 0;  # pass

                foreach my $pos (0 .. $#{$rule->{'Groups'}->{$grp}} ) {

                    ## $pos



                    my $condition = $rule->{'Groups'}->{$grp}->[$pos]->{'condition'};
                    ## $condition

                    $condition =~ s/-//g;
                    $condition =~ tr/\(\)x/\{\}\./;
                    ## $condition

                    my $condition_regex = qr/\A${condition}\z/;
                    ## $condition_regex





                    my $rule_hmm_start = $rule->{'Groups'}->{$grp}->[$pos]->{'hmmStart'};
                    my $rule_hmm_end = $rule->{'Groups'}->{$grp}->[$pos]->{'hmmEnd'};


                    # Fix for single residue Nter/Cter location (there can be only one at a time)
                    if ($rule_hmm_start eq 'Nter') {
                        my $cond_match = $condition;
                        $rule_hmm_start = $rule_hmm_end - get_ter_offset($cond_match);
                    }
                    if ($rule_hmm_end eq 'Cter') {
                        my $cond_match = $condition;
                        $rule_hmm_end = $rule_hmm_start + get_ter_offset($cond_match);
                    }
                    ## $rule_hmm_start
                    ## $rule_hmm_end



                    # my $target_seq = substr($query_seq, $map->[$rule_hmm_start], $map->[$rule_hmm_start] - $map->[$rule_hmm_start] + 1);
                    # ## $target_seq

# my $target_seq = substr($query_seq, $map->[$rule_hmm_start], $rule_hmm_end - $rule_hmm_start + 6);

                    my $seq_start = $map->[$rule_hmm_start];
                    my $seq_end = $map->[$rule_hmm_end];
                    ## $seq_start
                    ## $seq_end


                    $query_seq =~ s/-//g;
                    my $target_seq = '';
                    $target_seq = substr($query_seq, $seq_start, $seq_end - $seq_start + 1) unless (!defined $seq_start || !defined $seq_end);


                    # debug failures
                    # if ($target_seq =~ /${condition_regex}/) {
                    #     ### TARGET MATCH
                    #     my $condition = $rule->{'Groups'}->{$grp}->[$pos]->{'condition'};
                    #     ### $condition
                    #     ### $condition_regex
                    #     ### $target_seq
                    # }
                    # else {
                    #     ### TARGET NO MATCH
                    #     ### $condition
                    #     ### $condition_regex
                    #     ### $target_seq

                    # }



                    $pass_count += ($target_seq =~ /${condition_regex}/);
                    ## $pass_count

                }


                # my $pos_count = @{$rule->{'Groups'}->{$grp}};
                ## $pos_count
                ## $pass_count


                if (@{$rule->{'Groups'}->{$grp}} == $pass_count) {
                    ## WE HAVE A PASS!
                    $pass = 1;
                }
                # else {
                #     ## NO PASS!
                #     ## $rule
                #     # die;
                # }


            }
## END OF GROUP

        $query_rules{$query_id}{$rule_id} = $pass;


    

        }

    
    }

### %query_rules
}



# map base positions from alignment, from query HMM coords to (ungapped) target sequence coords
sub map_hmm_to_seq {
    my ($hmm_pos, $hmm, $seq) = @_;

    # so we can extract residues direct from the alignment; need to add the offset later
    my $seq_pos = 0;
    my @map;

    for (my $i = 0; $i < length $hmm; $i++) {
        $map[$hmm_pos] = $seq_pos;
        $hmm_pos++ if substr($hmm, $i, 1) ne '.';
        $seq_pos++ if substr($seq, $i, 1) ne '-';
    }

    return \@map;
}


# When hmmStart/hmmEnd is Nter/Cter, calculate offset for termination based on how many bases long the condition is
sub get_ter_offset {
    my ($cond_match) = @_;

    my $offset = -1;

    # transform the condition in a string of char counts
    $cond_match =~ s/\[\w+\]/1/;
    $cond_match =~ s/[A-Z]/1/g;
    $cond_match =~ s/\.\{//g;
    $cond_match =~ s/\}//g;

    foreach my $char (split //, $cond_match) {
      $offset += $char;
    }

    return $offset;

}











## DEPRECATED CODE

# sub get_rules {
#     my ($rules_file) = @_;



#     my %data;

#     do {
#         local $/ = "//\n";
#         open (my $in, '<', $rules_file) or die "Failed top open $rules_file file: $!\n";
#         while (my $block = <$in>) {
#             my $rule_hash = _parse_rules($block);
#             ## $rule_hash
#             $data{$rule_hash->{'AC'}} = $rule_hash;
#         }
#         close($in) or die "Failed to close $rules_file\n";
#     };

# ## %data
#     #return data structure
#     return \%data;
# }





# sub align_template {
#     my ($rules, $template_folder) = @_;
# ## $rules
#     RULE:
#     foreach my $acc (keys %{$rules}) {
#         ## $acc

#         my $cur_rule = $rules->{$acc};
# ## $cur_rule

#         my $prot_id = $cur_rule->{'Feature'}->{'from'};
#         ## $prot_id

#         my $fasta_file = "${template_folder}/${prot_id}.fa";
# ## $fasta_file
#         open (my $in, '<', "$fasta_file") or die "Failed top open $fasta_file file: $!\n";


#         my $prot_model = <$in>;
#         $prot_model =~ s/.*\t//;
#         $prot_model =~ s/\n//;
#         ## $prot_model

#         my $prot_seq = <$in>;
#         $prot_seq =~ s/\n//;
#         ## $prot_seq

#         my $stockholm = `hmmalign data/sr_hmm/${prot_model}.hmm ${fasta_file}`;
# ## $stockholm
#         my $alignment_str;
#         while ($stockholm =~ /\n${prot_id}\s+([^\n]*)/g) {
#             $alignment_str .= $1;
#         }
#         ## $alignment_str


#         my $alignment = align_map($alignment_str);
#         ## @map;

#         foreach my $grp (keys %{$cur_rule->{'Groups'}}) {
#                 ## $grp

#             for my $pos (0 .. $#{$cur_rule->{'Groups'}->{$grp}} ) {
#                 ## $pos
#                 my $both = 0;

#                 if ($cur_rule->{'Groups'}->{$grp}->[$pos]->{'start'} eq 'Nter' || $cur_rule->{'Groups'}->{$grp}->[$pos]->{'end'} eq 'Cter') {
#                     warn "rule $acc, group $grp, pos $pos: Start/End is Nter/Cter";
#                     delete $rules->{$acc};
#                     # undef $cur_rule;
#                     next RULE;
#                     # if (length ($cur_rule->{'Groups'}->{$grp}->[$pos]->{'condition'}) == 1) {
#                     #     $cur_rule->{'Groups'}->{$grp}->[$pos]->{'start'} = $cur_rule->{'Groups'}->{$grp}->[$pos]->{'end'};
#                     # } else {
#                     #     die "Nter with range not supported. Fixme";
#                     # }
#                 }

#                 if ($alignment->[$cur_rule->{'Groups'}->{$grp}->[$pos]->{'start'}] eq "-") {
#                     warn $prot_model." ".$prot_id."/".$cur_rule->{'Groups'}->{$grp}->[$pos]->{'start'}." is non-match state.\n";
#                     $both++;
#                     ## second_if: $cur_rule->{'Groups'}->{$grp}->[$pos]
#                 }
#                 if ($alignment->[$cur_rule->{'Groups'}->{$grp}->[$pos]->{'end'}] eq "-") {
#                     warn $prot_model." ".$prot_id."/".$cur_rule->{'Groups'}->{$grp}->[$pos]->{'end'}." is non-match state.\n";
#                     $both++;
#                     ## third_if: $cur_rule->{'Groups'}->{$grp}->[$pos]
#                 }
#                 #Some of the disulphide bridges lack on of the positions. However, the two current failing
#                 #rules means that we are completely dependent on a length.    This is not ideal, but we
#                 #can try.
#                 if ($both == 2) {
#                     warn "No anchor point for $prot_model ".$prot_id."/".$cur_rule->{'Groups'}->{$grp}->[$pos]->{'start'}."-".$cur_rule->{'Groups'}->{$grp}->[$pos]->{'end'}."\n";
#                 }

#                 $cur_rule->{'Groups'}->{$grp}->[$pos]->{hmmStart} = $alignment->[$cur_rule->{'Groups'}->{$grp}->[$pos]->{'start'}];
#                 $cur_rule->{'Groups'}->{$grp}->[$pos]->{hmmEnd} = $alignment->[$cur_rule->{'Groups'}->{$grp}->[$pos]->{'end'}];


#                 ## $pos
#                 ## $cur_rule

#             }

#         }

#     }

#     return $rules;
# }





# sub read_template_fasta {
#     my ($fasta_file) = @_;

#     my ($path) = $fasta_file =~ /(.+)\/[^\/+]/;

#     my $data;

#     do {
#         local $/ = ">";
#         open (my $in, '<', $fasta_file) or die "Failed top open $fasta_file file: $!\n";
#         while (my $block = <$in>) {

#             $block =~ s/\n?>?\z//;
#             next if !$block;



#             if ($block =~ /\A(\w+)\s+(\w+-\d+)\n/) {
#                 my $prot_it = $1;
#                 my $hmm_id = $2;

#                 $data->{$prot_it}->{'hmm'} = $hmm_id;

#                 $block =~ s/\A(\w+)\s+(\w+-\d+)\n//;
#                 $block =~ s/\n//g;

#                 $data->{$prot_it}->{'seq'} = $block;

#                 open (my $out, '>', "${path}/${prot_it}.fa" ) or die "Failed top open ${prot_it}.fa file: $!\n";

#                 print $out ">${prot_it}\n$data->{$prot_it}->{'seq'}\n";

#                 close($out) or die "Failed to close ${prot_it}.fa\n";

#                 $data->{$prot_it}->{'file'} = "${path}/${prot_it}.fa";

#             } else {
#                 die "ERROR: failed to parse fasta block: \"$block\"\n";
#             }
#             ## $block
#             # my $stuff = _process_rules_block($block);
#             ## $stuffs
#         }
#         close($in) or die "Failed to close $fasta_file\n";
#     };

#     #return data structure
#     return $data;
# }







# sub run_search {
#     my ($query_file, $hmm_folder) = @_;

# ### $query_file
# ### $hmm_folder


#     my $out = 'result.out';

#     my $cmd = "hmmscan --notextw -o $out ${hmm_folder}/sr_hmm.hmm $query_file";
#     ### $cmd
#     system($cmd) && die qq(Failed to run "$cmd");
#     ### $out
# # my $hmmRes = Bio::Pfam::HMM::HMMResultsIO->new;
# # ### $hmmRes
#     my $res_obj = Bio::Pfam::HMM::HMMResultsIO->new->parseMultiHMMER3($out);
# ## $res_obj



#     # print "Sequence acc\tSeq start\tSeq end\tRule acc\tFeature group\tTrigger\tModel\tTemplate\tTemplate start\tTemplate end\tTaxonomic scope\tFeature type\tFeature description\n";
#     foreach my $query_match (@{$res_obj}) {
# ## $query_match
#         my $query_id = $query_match->{'seqName'};
#         ## $query_id


#         # loop over HMMs - i.e. Unirule profiles
#         # my $rule_id = $query_match->{hmmName};
#         # ### $rule_id
#         # $rule_id =~ s{\.}{\-}; # FIXME - check what format we are using for HMMs - may not need this


#         foreach my $target_match (@{$query_match->{'units'}}) {
#             ## $target_match
#             # loop over alignments
#             my $rule_id = $target_match->{'name'};

#             ## $rule_id

#             my $tname = $target_match->{name}; # the sequence ID/accession
#             my $hmm_seq = $target_match->{hmmalign}{hmm};
#             my $target_seq = $target_match->{hmmalign}{seq};
#             my $hmm_from = $target_match->{hmmFrom};
#             my $seq_from = $target_match->{seqFrom};
#             # my $map = map_hmm_to_seq($hmm_from, $hmm_seq, $target_seq);
#             # while (my ($group_id, $group_rules) = each %{$rules->{$rule_id}}) {
#                 # Loop over the groups in the rule - check and print
#                 # print_group_matches($rule_id, $group_id, $tname, $target_seq, $seq_from, $map, $group_rules, $fh, 0);
#             # }
#         }
#     }





# }


__PACKAGE__->meta->make_immutable;

1;

__END__