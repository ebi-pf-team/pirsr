package PIRSR;

use strict;
use warnings;

use Smart::Comments;






sub read_template_fasta {
    my ($fasta_file) = @_;

    my ($path) = $fasta_file =~ /(.+)\/[^\/+]/;

    my $data;

    do {
        local $/ = ">";
        open (my $in, '<', $fasta_file) or die "Failed top open $fasta_file file: $!\n";
        while (my $block = <$in>) {

            $block =~ s/\n?>?\z//;
            next if !$block;



            if ($block =~ /\A(\w+)\s+(\w+-\d+)\n/) {
                my $prot_it = $1;
                my $hmm_id = $2;

                $data->{$prot_it}->{'hmm'} = $hmm_id;

                $block =~ s/\A(\w+)\s+(\w+-\d+)\n//;
                $block =~ s/\n//g;

                $data->{$prot_it}->{'seq'} = $block;

                open (my $out, '>', "${path}/${prot_it}.fa" ) or die "Failed top open ${prot_it}.fa file: $!\n";

                print $out ">${prot_it}\n$data->{$prot_it}->{'seq'}\n";

                close($out) or die "Failed to close ${prot_it}.fa\n";

                $data->{$prot_it}->{'file'} = "${path}/${prot_it}.fa";

            } else {
                die "ERROR: failed to parse fasta block: \"$block\"\n";
            }
            ## $block
            # my $stuff = _process_rules_block($block);
            ## $stuffs
        }
        close($in) or die "Failed to close $fasta_file\n";
    };

    #return data structure
    return $data;
}






sub get_rules {
    my ($rules_file) = @_;



    my %data;

    do {
        local $/ = "//\n";
        open (my $in, '<', $rules_file) or die "Failed top open $rules_file file: $!\n";
        while (my $block = <$in>) {
            my $rule_hash = _parse_rules($block);
            ## $rule_hash
            $data{$rule_hash->{'AC'}} = $rule_hash;
        }
        close($in) or die "Failed to close $rules_file\n";
    };

    #return data structure
    return \%data;
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
                } elsif ($line =~ /\AFT\s+(\S+)\s+(\w+)\s+(\w+)(\s+(.*))?\z/) {
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
    my ($rules, $fasta) = @_;
## $rules
    RULE:
    foreach my $acc (keys %{$rules}) {
        ## $acc

        my $cur_rule = $rules->{$acc};
## $cur_rule

        my $prot_id = $cur_rule->{'Feature'}->{'from'};
        ## $prot_id

        if (!$fasta->{$prot_id}) {
            warn "FASTA for template sequence $prot_id not available.\n";
            next;
        }

        my $prot_model = $fasta->{$prot_id}->{'hmm'};
        ## $prot_model

        my $prot_seq = $fasta->{$prot_id}->{'seq'};
        ## $prot_seq

        my $fasta_file = $fasta->{$prot_id}->{'file'};
        ## $fasta_file



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


# $rules->{$acc}->{'Groups'}->{$grp}->{$pos}

        foreach my $grp (keys %{$cur_rule->{'Groups'}}) {
                ## $grp

            for my $pos (0 .. $#{$cur_rule->{'Groups'}->{$grp}} ) {
                ## $pos
                my $both = 0;

                if ($cur_rule->{'Groups'}->{$grp}->[$pos]->{'start'} eq 'Nter' || $cur_rule->{'Groups'}->{$grp}->[$pos]->{'end'} eq 'Cter') {
                    warn "rule $acc, group $grp, pos $pos: Start/End is Nter/Cter";
                    delete $rules->{$acc};
                    # undef $cur_rule;
                    next RULE;
                    # if (length ($cur_rule->{'Groups'}->{$grp}->[$pos]->{'condition'}) == 1) {
                    #     $cur_rule->{'Groups'}->{$grp}->[$pos]->{'start'} = $cur_rule->{'Groups'}->{$grp}->[$pos]->{'end'};
                    # } else {
                    #     die "Nter with range not supported. Fixme";
                    # }
                }

                if ($alignment->[$cur_rule->{'Groups'}->{$grp}->[$pos]->{'start'}] eq "-") {
                    warn $prot_model." ".$prot_id."/".$cur_rule->{'Groups'}->{$grp}->[$pos]->{'start'}." is non-match state.\n";
                    $both++;
                    ## second_if: $cur_rule->{'Groups'}->{$grp}->[$pos]
                }
                if ($alignment->[$cur_rule->{'Groups'}->{$grp}->[$pos]->{'end'}] eq "-") {
                    warn $prot_model." ".$prot_id."/".$cur_rule->{'Groups'}->{$grp}->[$pos]->{'end'}." is non-match state.\n";
                    $both++;
                    ## third_if: $cur_rule->{'Groups'}->{$grp}->[$pos]
                }
                #Some of the disulphide bridges lack on of the positions. However, the two current failing
                #rules means that we are completely dependent on a length.    This is not ideal, but we
                #can try.
                if ($both == 2) {
                    warn "No anchor point for $prot_model ".$prot_id."/".$cur_rule->{'Groups'}->{$grp}->[$pos]->{'start'}."-".$cur_rule->{'Groups'}->{$grp}->[$pos]->{'end'}."\n";
                }

                $cur_rule->{'Groups'}->{$grp}->[$pos]->{hmmStart} = $alignment->[$cur_rule->{'Groups'}->{$grp}->[$pos]->{'start'}];
                $cur_rule->{'Groups'}->{$grp}->[$pos]->{hmmEnd} = $alignment->[$cur_rule->{'Groups'}->{$grp}->[$pos]->{'end'}];


                ## $pos
                ## $cur_rule

            }

## $rules

        }


## $rules

# die;

    }



    return $rules;
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
            $map[$seq_pos] = '-';
            $seq_pos++;
        }
    }

    return \@map;
}






1;