#!/usr/bin/perl

#Give Perl the modules to use. If one of these is not installed, the script will fail. However, Perl should print to the screen which one is missing.
use Bio::SeqIO;
use strict;
use warnings;

#Usage statement
my $usage = "Takes single GenBank file and produces two input files for vConTACT2 (proteins and mapping).\nUsage: perl single_GenBank_to_vConTACT_inputs.pl input.gbk";

#Tell Perl where the genbank file is
my $gb = $ARGV[0] or die "$usage\n";

#Define an output for the proteins
my $out = Bio::SeqIO->new(-file => ">vConTACT2_proteins.faa", 
                          -format => 'fasta');

#Define an output for the csv
open(OUT,">vConTACT2_gene_to_genome.csv") or die;

#Read it in as a Bio::SeqIO object...
my $seq_in = Bio::SeqIO->new(-file => "$gb",
                             -format => "genbank");

#Go through the Genbank file one entry at a time
while (my $inseq_obj = $seq_in->next_seq) {

    #Get the name of the sequence
    my $acc =  $inseq_obj->accession_number;
    chomp $acc;
    
    #Get features
    for my $feat_object ($inseq_obj->get_SeqFeatures) {
    
        #Is it a CDS?
        if ($feat_object->primary_tag eq 'CDS') {
        
            #Just double check it has a locus tag and translation, that it's nothing weird
            if ($feat_object->has_tag('locus_tag') && $feat_object->has_tag('translation')) {
                    
                #Get the ID
                my ($protein_id) = $feat_object->get_tag_values('locus_tag');
                chomp $protein_id;
           
                #Print the things we want
                print OUT "$protein_id",",","$acc",",","none","\n";
                                       
                #Get the translation of the CDS
                my ($trans) = $feat_object->get_tag_values('translation');
                chomp $trans;
                
                #Define a new sequence to be written out
                my $protein = Bio::Seq->new(-display_id => $protein_id, 
                                            -seq => $trans,
                                            -alphabet => 'protein');
                
                #Write the sequence to output
                $out->write_seq($protein);   
            }
        }
    }  
}