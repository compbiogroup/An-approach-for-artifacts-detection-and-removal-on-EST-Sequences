#!/usr/bin/perl -w

use strict;

#############################################################################
# BD2006 - Trimming procedure                                               #
# Version 1.1                                                               #
# 2007-09-30                                                                #
# Perl script                                                               #
#                                                                           #
# Usage:                                                                    #
#   perl bd2006trimmer.pl <fastafile> <qualfile>                            #
#                                                                           #
# WARNING: <fastafile> and <qualfile> must be in the same order.            #
#                                                                           #
#          The result will be written to:                                   #
#          - <fastafile>.trimmed                                            #
#          - <fastafile>.trimmed.qual                                       #
#                                                                           #
# ------------------------------------------------------------------------- #
#  Trimming Procedures Set Reference:                                       #
#   C. Baudet and Z. Dias. New EST Trimming Procedure Applied to SUCEST     #
#   Sequences. In Marie-France Sagot and Maria Emilia M. T. Walter,	    #
#   editors, Advances in Bioinformatics and Computational Biology,	    #
#   volume 4643 of LNBI, pages 57-68. Springler 2007.			    #
# 									    #
#  \@InProceedings{BD2007,						    #
#    editor    = "Marie-France Sagot and Maria Emilia M.T. Walter",	    #
#    booktitle = "Advances in Bioinformatics and Computational Biology",    #
#    publisher = "Springer",						    #
#    location  = "Heidelberg",						    #
#    series    = "LNBI",						    #
#    volume    = "4643",						    #
#    year      = "2007",						    #
#    isbn      = "978-3-540-73730-8",					    #
#    author    = "C. Baudet and Z. Dias",				    #
#    title     = "{N}ew {EST} {T}rimming {P}rocedure {A}pplied to {SUCEST}  #
#                 {S}equences",                                             #
#    pages     = "57--68"                                                   #
#  }                                                                        #
# ------------------------------------------------------------------------- #
#                                                                           #
# ------------------------------------------------------------------------- #
#  Slippage Methods Reference:                                              #
#   C. Baudet and Z. Dias. Analysis of slipped sequences in EST projects.   #
#   Genetics and Molecular Research, 5(1):169-181,2006.			    #
# 									    #
#  \@Article{BD2006a,							    #
#    author    = "C. Baudet and Z. Dias",				    #
#    title     = "{A}nalysis of slipped sequences in {EST} projects",	    #
#    journal   = "Genetics and Molecular Research",			    #
#    volume    = "5",							    #
#    number    = "1",							    #
#    pages     = "169--181",						    #
#    year      = "2006"							    #
#  }                                                                        #
# ------------------------------------------------------------------------- #
#                                                                           #
# ------------------------------------------------------------------------- #
#  Master Thesis Reference:                                                 #
#   C. Baudet. Uma abordagem para deteccao e remocao de artefatos em	    #
#   sequencias ESTs. Master's thesis, University of Campinas, Brazil,	    #
#   2006. In Portuguese.						    #
# 									    #
#  \@MastersThesis{Baudet2006,						    #
#    author    = "C. Baudet",						    #
#    title     = "Uma abordagem para detec{\c{c}}{\~{a}}o e		    #
#                  remo{\c{c}}{\~{a}}o de artefatos em			    #
#                  seq{\"{u}}{\^{e}}ncias {EST}s",			    #
#    school    = "University of Campinas, Brazil",			    #
#    year      = "2006",						    #
#    note      = "In Portuguese"					    #
#  }                                                                        #
# ------------------------------------------------------------------------- #
#############################################################################

#############################################################################
#                                                                           #
# CHANGELOG:                                                                #
#                                                                           #
# Version 1.1: 2007-09-30                                                   #
#   - Vector removal based on Univec and VecScreen added.                   #
#   - New references added                                                  #
#                                                                           #
# Version 1.0: 2006-12-21                                                   #
#   - First version: Methods developed during my master thesis              #
#                                                                           #
#############################################################################

#############################################################################
# GLOBAL VARIABLES AND CONSTANTS SECTION BEGIN                              #
#############################################################################

# Script configuration file
my $CONFIGURATION_FILE = "./bd2006trimmer.conf";

# Configuration Hash : Holds all configuration parameters
my %configuration;

# FASTA file: Fasta with all sequences that will be processed
my $fastaFile;

# QUAL file: Quality file (must be in the same order of the fasta file)
my $qualFile;

# Number of sequences found in the FASTA file
my $numberOfSequences;

# Number of FASTA files
my $numberOfFASTAs;

# Hash that maps sequence names to sequence ids
my %nameToId;

# Hash that maps sequence ids to sequence names
my %idToName;

# Hash that maps sequence ids to sequence size
my %sequenceSizes;

# Array that holds error probability values defined by the index
# Example: $errorProbabilityArray[10] = error probability value for quality 10
my @errorProbabilityArray;

# Temporary Files : Holds all temporary files
my %temporaryFiles;

# Ribosomal artifacts hash: Holds all ribosomal artifacts
my %ribosomalArtifacts;

# Vector artifacts hash: Holds all vector artifacts
my %vectorArtifacts;

# Univec vector artifacts hash: Holds all univec vectors artifacts
my %univecArtifacts;

# Adaptor artifacts hash: Holds all adaptor artifacts
my %adaptorArtifacts;

# Poly-A artifacts hash: Holds all poly-A artifacts
my %polyAArtifacts;

# Poly-T artifacts hash: Holds all poly-T artifacts
my %polyTArtifacts;

# Slippage artifacts hash: Holds all slippage artifacts
my %slippageArtifacts;

# Low Quality artifacts hash: Holds all low quality artifacts
my %lowQualityArtifacts;


#############################################################################
# GLOBAL VARIABLES AND CONSTANTS SECTION END                                #
#############################################################################



#############################################################################
# MAIN PROGRAM SECTION BEGIN                                                #
#############################################################################

printReference();

open(SCRIPTLOG, ">./BD2006TRIMMER.LOG");
$temporaryFiles{"bd2006trimmer.err"} = "bd2006trimmer.err";

configureScript();
checkParameters();

open(TRIMMEDFASTA, ">$fastaFile.trimmed");

if ($configuration{"PERFORM_LOW_QUALITY_DETECTION"} == 1) {
    open(TRIMMEDQUAL, ">$fastaFile.trimmed.qual");
}

# First step - Clear input files
print SCRIPTLOG "\nChecking input files:\n";
clearFastaFile();
clearQualFile();
print SCRIPTLOG "Checking input files done.\n";


# Second step - Run all auxiliar softwares
print SCRIPTLOG "\nRunning auxiliar programs:\n";
if ($configuration{"PERFORM_RIBOSOMAL_DETECTION"} == 1) {
    runRibosomalBLAST();
}
if ($configuration{"PERFORM_VECTOR_DETECTION"} == 1) {
    if ($configuration{"VECTOR_DETECTION_METHOD"} eq "cm") {
	runVectorCrossMatch();
    } elsif ($configuration{"VECTOR_DETECTION_METHOD"} eq "uv") {
	runUnivecBLAST();
    } else {
	runVectorCrossMatch();
	runUnivecBLAST();
    }
}
if ($configuration{"PERFORM_ADAPTOR_DETECTION"} == 1) {
    runAdaptorSwat();
}
if ($configuration{"PERFORM_POLY-A_DETECTION"} == 1) {
    createPolyTemplate($configuration{"POLYA_TEMPLATE_SIZE"}, "a");
    runPolyASwat();
}
if ($configuration{"PERFORM_POLY-T_DETECTION"} == 1) {
    createPolyTemplate($configuration{"POLYT_TEMPLATE_SIZE"}, "t");
    runPolyTSwat();
}
print SCRIPTLOG "Running auxiliar programs done.\n\n";


# Third step - Process auxiliar programs output to extract artifacts
print SCRIPTLOG "\nReading output files to extract artifacts:\n";
if ($configuration{"PERFORM_RIBOSOMAL_DETECTION"} == 1) {
    processRibosomalBLASTOutput();
}
if ($configuration{"PERFORM_VECTOR_DETECTION"} == 1) {
    if ($configuration{"VECTOR_DETECTION_METHOD"} eq "cm") {
	processVectorCrossMatchOutput();
    } elsif ($configuration{"VECTOR_DETECTION_METHOD"} eq "uv") {
	processUnivecBLASTOutput();
    } else {
	processVectorCrossMatchOutput();
	processUnivecBLASTOutput();
    }
}
if ($configuration{"PERFORM_ADAPTOR_DETECTION"} == 1) {
    processAdaptorSwatOutput();
}
if ($configuration{"PERFORM_POLY-A_DETECTION"} == 1) {
    processPolyASwatOutput();
}
if ($configuration{"PERFORM_POLY-T_DETECTION"} == 1) {
    processPolyTSwatOutput();
}
print SCRIPTLOG "Reading output files to extract artifacts done.\n\n";


# Fourth step - Process all sequences to extract low quality and
# slippage artifacts
print SCRIPTLOG "\nDetecting slippage and low quality and identifying insert:\n\n";

for (my $i = 0; $i < $numberOfFASTAs; $i++) {

    my $key1 = "cleared_fasta_$i";
    my $key2 = "cleared_qual_$i";

    open(FASTA, $temporaryFiles{$key1});

    if ($configuration{"PERFORM_LOW_QUALITY_DETECTION"} == 1) {
	initErrorProbabilityArray();
	open(QUAL, $temporaryFiles{$key2});
    }

    while (my $line = <FASTA>) {
	
	if ($line =~ /^>(\d+)$/) {

	    chop($line);

	    my $sequenceId = $1;

	    print SCRIPTLOG "[$sequenceId]\t$idToName{$sequenceId}\n";

	    my $sequence = <FASTA>;
	    chop($sequence);

	    my @qualities;
	    if ($configuration{"PERFORM_LOW_QUALITY_DETECTION"} == 1) {
		my $qualityString = <QUAL>; #Ignore header
		$qualityString = <QUAL>;
		chop($qualityString);
		@qualities = split(" ", $qualityString);
	    }

	    if (! defined $ribosomalArtifacts{$sequenceId}) {

		# Detects slippage and low quality artifacts
		if ($configuration{"PERFORM_SLIPPAGE_DETECTION"} == 1) {
		    detectSlippageArtifacts($sequenceId, $sequence);
		}
		if ($configuration{"PERFORM_LOW_QUALITY_DETECTION"} == 1) {
		    qualityAnalysis($sequenceId, @qualities);
		}
		
		# Mask all artifacts found
		my $sequenceSize = $sequenceSizes{$sequenceId};
		my $maskedSequence = $sequence;
		
		if ($configuration{"PERFORM_VECTOR_DETECTION"} == 1 &&
		    defined $vectorArtifacts{$sequenceId}) {
		    $maskedSequence = 
			maskArtifacts($maskedSequence, 
				      $sequenceSize, 
				      $vectorArtifacts{$sequenceId});
		    print SCRIPTLOG "[VE]\t$vectorArtifacts{$sequenceId}\n";
		
		}
		if ($configuration{"PERFORM_VECTOR_DETECTION"} == 1 &&
		    defined $univecArtifacts{$sequenceId}) {
		    mergeUnivecArtifacts($sequenceId);
		    $maskedSequence = 
			maskArtifacts($maskedSequence, 
				      $sequenceSize, 
				      $univecArtifacts{$sequenceId});
		    print SCRIPTLOG "[UV]\t$univecArtifacts{$sequenceId}\n";
		}
		if ($configuration{"PERFORM_ADAPTOR_DETECTION"} == 1 &&
		    defined $adaptorArtifacts{$sequenceId}) {
		    $maskedSequence = 
			maskArtifacts($maskedSequence, 
				      $sequenceSize,
				      $adaptorArtifacts{$sequenceId});
		    print SCRIPTLOG "[AD]\t$adaptorArtifacts{$sequenceId}\n";
		
		}
		if ($configuration{"PERFORM_POLY-A_DETECTION"} == 1 &&
		    defined $polyAArtifacts{$sequenceId}) {
		    $maskedSequence = 
			maskArtifacts($maskedSequence, 
				      $sequenceSize, 
				      $polyAArtifacts{$sequenceId});
		    print SCRIPTLOG "[PA]\t$polyAArtifacts{$sequenceId}\n";
		}
		if ($configuration{"PERFORM_POLY-T_DETECTION"} == 1 &&
		    defined $polyTArtifacts{$sequenceId}) {
		    $maskedSequence = 
			maskArtifacts($maskedSequence, 
				      $sequenceSize, 
				      $polyTArtifacts{$sequenceId});
		    print SCRIPTLOG "[PT]\t$polyTArtifacts{$sequenceId}\n";
		}
		if ($configuration{"PERFORM_SLIPPAGE_DETECTION"} == 1 &&
		    defined $slippageArtifacts{$sequenceId}) {
		    $maskedSequence = 
			maskArtifacts($maskedSequence, 
				      $sequenceSize, 
				      $slippageArtifacts{$sequenceId});
		    print SCRIPTLOG "[SL]\t$slippageArtifacts{$sequenceId}\n";
		}
		if ($configuration{"PERFORM_LOW_QUALITY_DETECTION"} == 1 &&
		    defined $lowQualityArtifacts{$sequenceId}) {
		    $maskedSequence = 
			maskArtifacts($maskedSequence, 
				      $sequenceSize, 
				      $lowQualityArtifacts{$sequenceId});
		    print SCRIPTLOG "[LQ]\t$lowQualityArtifacts{$sequenceId}\n";
		}

		# Finally, identifies the insert in the masked sequence
		my ($insertBegin, $insertEnd) = 
		    extractInsert($maskedSequence, $sequenceSize, @qualities);
	
		if ($insertBegin == -1) {
		    print SCRIPTLOG "[DISCARDED]\t$idToName{$sequenceId}\n";
		}
		else {
		    my $goodSequence = substr($sequence, $insertBegin, 
					      $insertEnd - $insertBegin);

		    print TRIMMEDFASTA "$idToName{$sequenceId}\n";
		    print TRIMMEDFASTA "$goodSequence\n";

		    if ($configuration{"PERFORM_LOW_QUALITY_DETECTION"} == 1) {
			my @goodQuality = 
			    @qualities[$insertBegin..($insertEnd-1)];
			print TRIMMEDQUAL "$idToName{$sequenceId}\n";
			print TRIMMEDQUAL "@goodQuality\n";
		    }
		}
	    }
	    else {
		# Ribosomal sequence is fully discarded
		print SCRIPTLOG "[RI]\t$ribosomalArtifacts{$sequenceId}\n";
		print SCRIPTLOG "[DISCARDED]\t$idToName{$sequenceId}\n";
	    }
	}# if ($line =~ /^>(\d+)$/) {...}
	
    }# while (my $line = <FASTA>) {...}

    close(FASTA);
    if ($configuration{"PERFORM_LOW_QUALITY_DETECTION"} == 1) {
	close(QUAL);
    }

} # for (my $i = 0; $i < $numberOfFastas; $i++) {...}

close(TRIMMEDFASTA);
if ($configuration{"PERFORM_LOW_QUALITY_DETECTION"} == 1) {
    close(TRIMMEDQUAL);
}

print SCRIPTLOG "\nDetecting slippage and low quality and identifying insert done.\n\n";

close(SCRIPTLOG);

# Finally... clearing up the mess
removeTemporaryFiles();


#############################################################################
# MAIN PROGRAM SECTION END                                                  #
#############################################################################



#############################################################################
# SUBROUTINES SECTION BEGIN                                                 #
#############################################################################

#############################################################################
# Verifies script's parameters
sub checkParameters {

    # Read the parameters
    ($fastaFile, $qualFile) = @ARGV;

    if (! defined $fastaFile) {
	print STDERR "ERROR: You must specify the FASTA file.\n";
	print SCRIPTLOG "ERROR: You must specify the FASTA file.\n";
	printUsage();
    }

    if (! -e $fastaFile) {
	print STDERR "ERROR: Could not find the FASTA file: $fastaFile.\n";
	print SCRIPTLOG "ERROR: Could not find the FASTA file: $fastaFile.\n";
	printUsage();
    }

    if (! defined $qualFile) {
	$configuration{"PERFORM_LOW_QUALITY_DETECTION"} = "0";
    }
    else {
	if (! -e $qualFile) {
	    print STDERR "ERROR: Could not find the QUAL file: $qualFile.\n";
	    print SCRIPTLOG "ERROR: Could not find the QUAL file: $qualFile.\n";
	    printUsage();
	}
    }
}


#############################################################################
# Loads script's configuration
sub configureScript {

    # Filling the hash with default values
    %configuration = (
		      "BLASTALL_PATH" => "/usr/local/blast-2.2.11/bin/blastall",
		      "CROSS_MATCH_PATH" => "/usr/local/bin/cross_match",
		      "SWAT_PATH" => "/usr/local/bin/swat",

		      "SWAT_MATRIX" => "./auxiliarFiles/MATRIX",

		      "NUMBER_OF_SEQUENCES_IN_AUXILIAR_FASTA" => "50000",
		      
		      "PERFORM_RIBOSOMAL_DETECTION" => "1",
		      "PERFORM_LOW_QUALITY_DETECTION" => "1",
		      "PERFORM_VECTOR_DETECTION" => "1",
		      "PERFORM_ADAPTOR_DETECTION" => "1",
		      "PERFORM_POLY-A_DETECTION" => "1",
		      "PERFORM_POLY-T_DETECTION" => "1",
		      "PERFORM_SLIPPAGE_DETECTION" => "1",

		      "RIBOSOMAL_DATABASE_PATH" => "./ribosomal/ribosomal_sucest",
		      "RIBOSOMAL_BLAST_PARAMETERS" => "-p blastn -b 1 -v 1 -m 8",
		      "RIBOSOMAL_MAXIMUM_EVALUE" => "1e-10",
		      
		      "VECTOR_DETECTION_METHOD" => "cm",
		      "VECTOR_FILE_PATH" => "./vector/pSport1.fasta",
		      "VECTOR_CROSSMATCH_PARAMETERS" => "-minmatch 12 -minscore 20",
		      "UNIVEC_DATABASE_PATH" => "./univec/UniVec_Core",
		      "UNIVEC_BLAST_MISMATCH" => "-5",
		      "UNIVEC_BLAST_OPENGAP" => "3",
		      "UNIVEC_BLAST_EXTENDEDGAP" => "3",
		      "UNIVEC_BLAST_SEARCHSPACE" => "1.75e12",
		      "UNIVEC_BLAST_EVALUE" => "700",
		      "UNIVEC_BLAST_FILTER" => "m D",
		      "UNIVEC_REMOVAL_CRITERIA" => "suspect",
		      "UNIVEC_TERMINAL_DISTANCE_CRITERIA" => "25",
		      "UNIVEC_TERMINAL_MINIMUM_SCORE_STRONG" => "24",
		      "UNIVEC_TERMINAL_MINIMUM_SCORE_MODERATE" => "19",
		      "UNIVEC_TERMINAL_MINIMUM_SCORE_WEAK" => "16",
		      "UNIVEC_NONTERMINAL_MINIMUM_SCORE_STRONG" => "30",
		      "UNIVEC_NONTERMINAL_MINIMUM_SCORE_MODERATE" => "25",
		      "UNIVEC_NONTERMINAL_MINIMUM_SCORE_WEAK" => "23",
		      "UNIVEC_SUSPECT_TO_END_DISTANCE" => "50",
		      "UNIVEC_SUSPECT_INTER_FRAGMENTS_DISTANCE" => "50",
		      
		      "ADAPTOR_FILE_PATH" => "./adaptor/pSport1_adaptors.fasta",
		      "ADAPTOR_SWAT_PARAMETERS" => "-gap_init -5 -gap_ext -5 -ins_gap_ext -5 -del_gap_ext -5 -end_gap -5",
		      "ADAPTOR_SUBTRACT_FROM_ADAPTOR_LENGTH" => "4",

		      "POLYA_SWAT_PARAMETERS" => "-gap_init -5 -gap_ext -5 -ins_gap_ext -5 -del_gap_ext -5 -end_gap -5 -minscore 10",
		      "POLYA_MINIMUM_SCORE" => "10",
		      "POLYA_TEMPLATE_SIZE" => "600",
		      
		      "POLYT_SWAT_PARAMETERS" => "-gap_init -5 -gap_ext -5 -ins_gap_ext -5 -del_gap_ext -5 -end_gap -5 -minscore 10",
		      "POLYT_MINIMUM_SCORE" => "10",
		      "POLYT_TEMPLATE_SIZE" => "600",
		      
		      "SLIPPAGE_TRIMMING_ALGORITHM" => "ec",
		      "SLIPPAGE_TRIMMING_STRATEGY" => "subsequence",
		      "SLIPPAGE_MINIMUM_NUMBER_OF_ECHOES" => "8",
		      "SLIPPAGE_MINIMUM_ECHO_SIZE" => "5",
		      "SLIPPAGE_THRESHOLD_VALUE" => "0.25",
		      
		      "LOW_QUALITY_TRIMMING_ALGORITHM" => "ms",
		      "LOW_QUALITY_SEARCH_FOR_ISLANDS" => "1",
		      "MAXIMUM_SUBSEQUENCE_MINIMUM_QUALITY" => "11",
		      "SLIDDING_WINDOW_WINDOW_SIZE" => "10",
		      "SLIDDING_WINDOW_QUALITY_THRESHOLD" => "16",
		      "SLIDDING_WINDOW_BAD_BASES_THRESHOLD" => "3",
		      "LOW_QUALITY_ISLAND_WINDOW_SIZE" => "10",
		      "LOW_QUALITY_ISLAND_MINIMUM_ERROR_PROBABILITY" => "0.2",
		      
		      "MINIMUM_SEQUENCE_SIZE" => "100"
		      );

    # Reading configuration file (if it exists)
    if (-e $CONFIGURATION_FILE) {
	
	print SCRIPTLOG "\nReading configuration file...";
	
	open(CONFIG, "$CONFIGURATION_FILE");
	while (my $line = <CONFIG>) {
	
	    chop($line);
	
	    if ($line =~ /^[^\#]/ && !($line =~ /^\s+$/)) {
		
		my ($key, $value) = split("=", $line);
		
		if (! defined $key || ! defined $value) {
		    print STDERR "\nERROR: Invalid pair key/value -> $line\n";
		    exit;
		}
		else {
		
		    if (defined $configuration{$key}) {
			
			$value =~ s/^\s+//g;
			$value =~ s/\s+$//g;
			$configuration{$key} = $value;
		    }
		    else {
			print STDERR "\nERROR Invalid configuration key -> $key\n";
			exit;
		    }
		}
	    }
	}
	close(CONFIG);
	
	print SCRIPTLOG "done.\n";
    }
}



#############################################################################
# Prints script's usage information
sub printUsage {
    print "\n-------------------------------------------------------------------------\n";
    print "  BD2006 - Trimming procedure\n";
    print "  Version 1.1\n";
    print "  2007-09-30\n";
    print "  Perl script\n";
    print "\n";
    print "  Usage:\n";
    print "    perl bd2006trimmer.pl <fastafile> <qualfile>\n";
    print "\n";
    print "  WARNING:\n";
    print "    <fastafile> and <qualfile> must be in the same order.\n";
    print "    <qualfile> is optional.\n";
    print "\n";
    print "  The result will be written to:\n";
    print "    + <fastafile>.trimmed\n";
    print "    + <fastafile>.trimmed.qual\n";
    print "-------------------------------------------------------------------------\n";
    print "\n";
    exit;
}


#############################################################################
# Prints script's reference information
sub printReference {
    print "-------------------------------------------------------------------------\n";
    print " Trimming Procedures Set Reference:\n\n";
    print "  C. Baudet and Z. Dias. New EST Trimming Procedure Applied to SUCEST\n";
    print "  Sequences. In Marie-France Sagot and Maria Emilia M. T. Walter,\n";
    print "  editors, Advances in Bioinformatics and Computational Biology,\n";
    print "  volume 4643 of LNBI, pages 57-68. Springler 2007.\n";
    print "\n";
    print " \@InProceedings{BD2007,\n";
    print "   editor    = \"Marie-France Sagot and Maria Emilia M.T. Walter\",\n";
    print "   booktitle = \"Advances in Bioinformatics and Computational Biology\",\n";
    print "   publisher = \"Springer\",\n";
    print "   location  = \"Heidelberg\",\n";
    print "   series    = \"LNBI\",\n";
    print "   volume    = \"4643\",\n";
    print "   year      = \"2007\",\n";
    print "   isbn      = \"978-3-540-73730-8\",\n";
    print "   author    = \"C. Baudet and Z. Dias\",\n";
    print "   title     = \"{N}ew {EST} {T}rimming {P}rocedure {A}pplied to {SUCEST}\n";
    print "                {S}equences\",\n";
    print "   pages     = \"57--68\"\n";
    print " }\n";
    print "-------------------------------------------------------------------------\n";
    print "\n";
    print "-------------------------------------------------------------------------\n";
    print " Slippage Methods Reference:\n\n";
    print "  C. Baudet and Z. Dias. Analysis of slipped sequences in EST projects.\n";
    print "  Genetics and Molecular Research, 5(1):169-181,2006.\n";
    print "\n";
    print " \@Article{BD2006a,\n";
    print "   author    = \"C. Baudet and Z. Dias\",\n";
    print "   title     = \"{A}nalysis of slipped sequences in {EST} projects\",\n";
    print "   journal   = \"Genetics and Molecular Research\",\n";
    print "   volume    = \"5\",\n";
    print "   number    = \"1\",\n";
    print "   pages     = \"169--181\",\n";
    print "   year      = \"2006\"\n";
    print " }\n";
    print "-------------------------------------------------------------------------\n";
    print "\n";
    print "-------------------------------------------------------------------------\n";
    print " Master Thesis Reference:\n\n";
    print "  C. Baudet. Uma abordagem para deteccao e remocao de artefatos em\n";
    print "  sequencias ESTs. Master's thesis, University of Campinas, Brazil,\n";
    print "  2006. In Portuguese.\n";
    print "\n";
    print " \@MastersThesis{Baudet2006,\n";
    print "   author    = \"C. Baudet\",\n";
    print "   title     = \"Uma abordagem para detec{\c{c}}{\~{a}}o e\n";
    print "                 remo{\c{c}}{\~{a}}o de artefatos em\n";
    print "                 seq{\"{u}}{\^{e}}ncias {EST}s\",\n";
    print "   school    = \"University of Campinas, Brazil\",\n";
    print "   year      = \"2006\",\n";
    print "   note      = \"In Portuguese\"\n";
    print " }\n";
    print "\n";
    print "-------------------------------------------------------------------------\n\n";
}


#############################################################################
# Creates a random name for auxitiar files
sub createTemporaryFileName {
    my $aux = time().rand();
    return "./$aux.tmp";
}


#############################################################################
# Runs BLAST against the ribosomal database
sub runRibosomalBLAST {
    print SCRIPTLOG "Running BLAST for ribosomal detection...";

    my $blast_output = createTemporaryFileName();
    $temporaryFiles{"ribosomal_blast"} = $blast_output;
    $temporaryFiles{"blast_error_log"} = "error.log";
    system("rm -f $blast_output");

    for (my $i = 0; $i < $numberOfFASTAs; $i++ ) {
	my $key = "cleared_fasta_$i";
	my $commandLine = "$configuration{'BLASTALL_PATH'} $configuration{'RIBOSOMAL_BLAST_PARAMETERS'} -d $configuration{'RIBOSOMAL_DATABASE_PATH'} -i $temporaryFiles{$key} >> $blast_output 2>> bd2006trimmer.err";
	system("$commandLine");
    }

    print SCRIPTLOG "done.\n";
}


#############################################################################
# Runs cross_match to detect vector regions
sub runVectorCrossMatch {
    print SCRIPTLOG "Running cross_match for vector detection...";

    my $cross_match_output = createTemporaryFileName();
    $temporaryFiles{"vector_crossmatch"} = $cross_match_output;
    system("rm -f $cross_match_output");

    for (my $i = 0; $i < $numberOfFASTAs; $i++ ) {
	my $key1 = "cleared_fasta_$i";
	my $key2 = "vector_error_log_$i";
	$temporaryFiles{$key2} = "$temporaryFiles{$key1}.log";
	my $commandLine = "$configuration{'CROSS_MATCH_PATH'} $temporaryFiles{$key1} $configuration{'VECTOR_FILE_PATH'} $configuration{'VECTOR_CROSSMATCH_PARAMETERS'} >> $cross_match_output 2>> bd2006trimmer.err";
	system("$commandLine");
    }

    print SCRIPTLOG "done.\n";
}


#############################################################################
# Runs BLAST to detect vector regions using Univec
sub runUnivecBLAST {
    print SCRIPTLOG "Running BLAST against Univec for vector detection...";

    my $blast_univec_output = createTemporaryFileName();
    $temporaryFiles{"blast_univec_output"} = $blast_univec_output;
    system("rm -f $blast_univec_output");
    
    for (my $i = 0; $i < $numberOfFASTAs; $i++ ) {
	my $key1 = "cleared_fasta_$i";
	my $commandLine = "$configuration{'BLASTALL_PATH'} -q $configuration{'UNIVEC_BLAST_MISMATCH'} -G $configuration{'UNIVEC_BLAST_OPENGAP'} -E $configuration{'UNIVEC_BLAST_EXTENDEDGAP'} -F \"$configuration{'UNIVEC_BLAST_FILTER'}\" -e $configuration{'UNIVEC_BLAST_EVALUE'} -Y $configuration{'UNIVEC_BLAST_SEARCHSPACE'} -p blastn -m 8 -i $temporaryFiles{$key1} -d $configuration{'UNIVEC_DATABASE_PATH'} >> $blast_univec_output 2>> bd2006trimmer.err";
	system("$commandLine");
    }

    print SCRIPTLOG "done.\n";
}


#############################################################################
# Runs swat to detect adaptor regions
sub runAdaptorSwat {
    print SCRIPTLOG "Running swat for adaptor detection...";

    my $swat_output = createTemporaryFileName();
    $temporaryFiles{"adaptor_swat"} = $swat_output;
    system("rm -f $swat_output");

    for (my $i = 0; $i < $numberOfFASTAs; $i++ ) {
	my $key = "cleared_fasta_$i";
	my $commandLine = "$configuration{'SWAT_PATH'} $configuration{'ADAPTOR_FILE_PATH'} $temporaryFiles{$key} $configuration{'ADAPTOR_SWAT_PARAMETERS'} -raw -N 60000000 -M $configuration{'SWAT_MATRIX'} >> $swat_output 2>> bd2006trimmer.err";
	system("$commandLine");
    }

    print SCRIPTLOG "done.\n";
}


#############################################################################
# Runs swat to detect poly-A regions
sub runPolyASwat {
    print SCRIPTLOG "Running swat for poly-a detection...";

    my $polyA_output = createTemporaryFileName();
    $temporaryFiles{"polya_swat"} = $polyA_output;
    system("rm -f $polyA_output");

    for (my $i = 0; $i < $numberOfFASTAs; $i++ ) {
	my $key = "cleared_fasta_$i";
	my $commandLine = "$configuration{'SWAT_PATH'} $temporaryFiles{'poly_a_template'} $temporaryFiles{$key} $configuration{'POLYA_SWAT_PARAMETERS'} -raw -N 60000000 -M $configuration{'SWAT_MATRIX'} >> $polyA_output 2>> bd2006trimmer.err";
	system("$commandLine");
    }

    print SCRIPTLOG "done.\n";
}


#############################################################################
# Runs swat to detect poly-T regions
sub runPolyTSwat {
    print SCRIPTLOG "Running swat for poly-t detection...";

    my $polyT_output = createTemporaryFileName();
    $temporaryFiles{"polyt_swat"} = $polyT_output;
    system("rm -f $polyT_output");

    for (my $i = 0; $i < $numberOfFASTAs; $i++ ) {
	my $key = "cleared_fasta_$i";
	my $commandLine = "$configuration{'SWAT_PATH'} $temporaryFiles{'poly_t_template'} $temporaryFiles{$key} $configuration{'POLYT_SWAT_PARAMETERS'} -raw -N 60000000 -M $configuration{'SWAT_MATRIX'} >> $polyT_output 2>> bd2006trimmer.err";
        system("$commandLine");
    }

    print SCRIPTLOG "done.\n";
}


#############################################################################
# Deletes all temporary files
sub removeTemporaryFiles {
    foreach my $key (keys %temporaryFiles) {
	my $file = $temporaryFiles{$key};
	if (defined $file && -e $file) {
	    system("rm -f $file");
	}
    }
    exit;
}


#############################################################################
# Builds the ribosomal artifact list
sub processRibosomalBLASTOutput {

    print SCRIPTLOG "Processing ribosomal BLAST output...";
    open(IN, $temporaryFiles{"ribosomal_blast"});

    while (my $line = <IN>) {

	chop($line);
	
	my ($queryId,
	    $subjectId,
	    $identity,
	    $alignmentLength,
	    $mismatches,
	    $gapopenings,
	    $queryStart,
	    $queryEnd,
	    $subjectStart,
	    $subjectEnd,
	    $eValue,
	    $bitScore) = split("\t",$line);

	if ($eValue =~ /^e/) {
	    $eValue = "1$eValue";
	}

	if ($eValue <= $configuration{"RIBOSOMAL_MAXIMUM_EVALUE"}) {
	
	    my $begin = $queryStart - 1;
	    my $end = $queryEnd;

	    my $aux = $ribosomalArtifacts{$queryId};
	    if (defined $aux) {
		$aux .= "\t$begin\t$end";
	    }
	    else {
		$aux = "$begin\t$end";
	    }
	    $ribosomalArtifacts{$queryId} = $aux;
	}
    }
    close(IN);
    print SCRIPTLOG "done.\n";
}


#############################################################################
# Builds the vector artifact list
sub processVectorCrossMatchOutput {

    print SCRIPTLOG "Processing vector cross_match output...";
    open(IN, $temporaryFiles{"vector_crossmatch"});
    while (my $line = <IN>) {
	
	if ($line =~ /^\s+\d+\s+\d\.\d\d\s+\d\.\d\d\s+\d\.\d\d\s+(\d+)\s+(\d+)\s+(\d+)\s+\(\d+\)/) {
	    my $seqName = $1;
	    my $begin = $2 - 1;
	    my $end = $3;
	    my $aux = $vectorArtifacts{$seqName};
	    if (defined $aux) {
		$aux .= "\t$begin\t$end";
	    }
	    else {
		$aux = "$begin\t$end";
	    }
	    $vectorArtifacts{$seqName} = $aux;
	}

    }
    close(IN);
    print SCRIPTLOG "done.\n";
}


#############################################################################
# Builds the vector artifact list 
sub processUnivecBLASTOutput {

    print SCRIPTLOG "Processing BLAST against Univec output...";

    open(IN, $temporaryFiles{"blast_univec_output"});
    while (my $line = <IN>) {

	chomp($line);

	my ($query, $subject, $identity, $length, $mismatches,
	    $gapOpenings, $queryStart, $queryEnd, $subjectStart,
	    $subjectEnd, $eValue, $bitScore) = split("\t", $line); 

	my $sequenceSize = $sequenceSizes{$query};

	my $numberOfQueryBases = $queryEnd - $queryStart + 1;
	my $numberOfGaps = $length - $numberOfQueryBases;
	my $matches = $numberOfQueryBases - $mismatches;

	my $score = $matches + 
	    ($mismatches * $configuration{"UNIVEC_BLAST_MISMATCH"}) - 
	    ($gapOpenings * $configuration{'UNIVEC_BLAST_OPENGAP'}) - 
	    ($numberOfGaps * $configuration{'UNIVEC_BLAST_EXTENDEDGAP'});

	# None     = -1
	# Strong   = 0
	# Moderate = 1
	# Weak     = 2

	my $category = -1; # None
	
	my $terminalDistanceCriteria = 
	    $configuration{"UNIVEC_TERMINAL_DISTANCE_CRITERIA"};

	my $terminalStrong = 
	    $configuration{"UNIVEC_TERMINAL_MINIMUM_SCORE_STRONG"};
	my $terminalModerate = 
	    $configuration{"UNIVEC_TERMINAL_MINIMUM_SCORE_MODERATE"};
	my $terminalWeak = 
	    $configuration{"UNIVEC_TERMINAL_MINIMUM_SCORE_WEAK"};
	
	my $nonTerminalStrong = 
	    $configuration{"UNIVEC_NONTERMINAL_MINIMUM_SCORE_STRONG"};
	my $nonTerminalModerate = 
	    $configuration{"UNIVEC_NONTERMINAL_MINIMUM_SCORE_MODERATE"};
	my $nonTerminalWeak = 
	    $configuration{"UNIVEC_NONTERMINAL_MINIMUM_SCORE_WEAK"};
	
	my $criteria = $configuration{"UNIVEC_REMOVAL_CRITERIA"};
	
	if ($queryStart <= $terminalDistanceCriteria || 
	    ($sequenceSize - $queryEnd + 1) <= $terminalDistanceCriteria) {
	    
	    # Found a terminal match
	    if ($score >= $terminalStrong) {
		$category = 0; # Strong
	    } elsif ($score >= $terminalModerate) {
		$category = 1; # Moderate
	    } elsif ($score >= $terminalWeak) {
		$category = 2; # Weak
	    }

	} else {

	    # Found an internal match
	    if ($score >= $nonTerminalStrong) {
		$category = 0; # Strong
	    } elsif ($score >= $nonTerminalModerate) {
		$category = 1; # Moderate
	    } elsif ($score >= $nonTerminalWeak) {
		$category = 2; # Weak
	    }
	    
	}
	
	my $addArtifact = 0;
	
	if ($category == 0) {
	    # Strong is always removed
	    $addArtifact = 1;
	} elsif ($category == 1 && $criteria ne "strong") {
	    # Moderate is removed only when criteria is not strong
	    $addArtifact = 1;
	} elsif ($category == 2 && $criteria ne "strong" && 
		 $criteria ne "moderate") {
	    # Weak is removed only when criteria is not strong or moderate
	    $addArtifact = 1;
	}
	
	if ($addArtifact == 1) {

	    my $aux = $univecArtifacts{$query};
	    
	    my $begin = $queryStart - 1;
	    my $end = $queryEnd;
	    my $artifact = "$begin\t$end";
	    if (defined $aux) {
		if (!($aux =~ /$artifact/)) {
		    $aux .= "\t$artifact";
		}
	    } else {
		$aux = $artifact;
	    }
	    $univecArtifacts{$query} = $aux;
	}
    } 
    close(IN); 
    print SCRIPTLOG "done.\n"; 
}

#############################################################################
# Merge artifact that show overlaps, find suspect regions
sub mergeUnivecArtifacts {
    my ($sequenceId) = @_;

    my $artifacts = $univecArtifacts{$sequenceId};
    my $sequenceSize = $sequenceSizes{$sequenceId};

    my @array = split("\t", $artifacts);
    my $arraySize = @array;
    my $nArtifacts = $arraySize / 2;

    # Someday I have to do a better job here....
    # A small macarronic code to sort the artifacs without implement
    # a sort algorithm to order pairs of values....
    my @sortedArray;
    for (my $i = 0; $i < $arraySize; $i+=2) {
	$sortedArray[@sortedArray] = sprintf("%15d", $array[$i]).
	    "\t".sprintf("%15d", $array[$i+1]);
    }
    @sortedArray = sort @sortedArray;
    
    my $suspect = 0;
    if ($configuration{"UNIVEC_REMOVAL_CRITERIA"} eq "suspect") {
	$suspect = 1;
    }

    my $toEndDistance = 
	$configuration{"UNIVEC_SUSPECT_TO_END_DISTANCE"};
    my $interFragmentsDistance = 
	$configuration{"UNIVEC_SUSPECT_INTER_FRAGMENTS_DISTANCE"};

    
    my $fragment = $sortedArray[0];
    $fragment =~ s/ //g;

    my ($previousBegin, $previousEnd) = split("\t", $fragment);

    if ($suspect == 1 && $previousBegin < $toEndDistance) {
	$previousBegin = 0;
    }

    my $newArtifacs;

    for (my $i = 1; $i < $nArtifacts; $i++) {

	$fragment = $sortedArray[$i];
	$fragment =~ s/ //g;
	my ($currentBegin, $currentEnd) = split("\t", $fragment);

	if ($currentBegin >= $previousBegin && $currentBegin < $previousEnd && 
	    $currentEnd > $previousEnd) {

	    # Overlap found: Merge both artifacts
	    $previousEnd = $currentEnd;

	} elsif ($currentBegin >= $previousEnd) {

	    # There is not a overlap, we must check the distance if
	    # suspect is turned on, otherwise we must register the
	    # artifact and start a new one....

	    # Decision [0|1] => 0 = Merge, 1 = Start a new one
	    my $decision = 1;

	    if ($suspect == 1) {

		if ($currentEnd - $previousEnd <= $interFragmentsDistance) {
		    $decision = 0;
		} 
		
	    }
	    
	    if ($decision == 1) {
		
		# Register and start a new one
		if (! defined $newArtifacs) {
		    $newArtifacs = "$previousBegin\t$previousEnd";
		} else {
		    $newArtifacs .= "\t$previousBegin\t$previousEnd";
		}
		$previousBegin = $currentBegin;
		$previousEnd = $currentEnd;

	    } else {

		# Merge
		$previousEnd = $currentEnd;
		
	    }

	} 

    } # for (my $i = 1; $i < $nArtifacts; $i++) {...}

    if ($suspect == 1 && 
	(($sequenceSize - $previousEnd + 1) <= $toEndDistance)) {
	$previousEnd = $sequenceSize;
    }

    # Register the last artifact...
    if (!defined $newArtifacs) {
	$newArtifacs = "$previousBegin\t$previousEnd";
    } else {
	$newArtifacs .= "\t$previousBegin\t$previousEnd";
    }
   
    $univecArtifacts{$sequenceId} = $newArtifacs;

}

#############################################################################
# Builds the adaptor artifact list
sub processAdaptorSwatOutput {

    print SCRIPTLOG "Processing adaptor swat output...";

    open(IN, $temporaryFiles{"adaptor_swat"});

    my ($sequenceId, $adaptorLength, $subjectBegin, $subjectEnd);

    while (my $line = <IN>) {

	if ($line =~ /^(\d+)\s+Length:\s+\d+/) {
	
	    if (defined $sequenceId) {
		
		my $threshold = $adaptorLength - $configuration{"ADAPTOR_SUBTRACT_FROM_ADAPTOR_LENGTH"};
		my $alignmentSize = $subjectEnd - $subjectBegin;
		
		if ($alignmentSize >= $threshold) {
		    my $aux = $adaptorArtifacts{$sequenceId};
		    if (! defined $aux) {
			$aux = "$subjectBegin\t$subjectEnd";
		    }
		    else {
			$aux .= "\t$subjectBegin\t$subjectEnd";
		    }
		    $adaptorArtifacts{$sequenceId} = $aux;
		}
	    }
	    $sequenceId = $1;
	    undef($subjectBegin);
	}
	elsif ($line =~ /^Subject\s+(\d+)\s+[ATCGNX\s\-]+(\d+)$/) {
	    if (! defined $subjectBegin) {
		$subjectBegin = $1 - 1;
	    }
	    $subjectEnd = $2;
	}
	elsif ($line =~ /^Query:.+Length:\s+(\d+)\s+residues$/) {
	    $adaptorLength = $1;
	}
    }
    close(IN);

    if (defined $sequenceId) {

	my $threshold = $adaptorLength - $configuration{"ADAPTOR_SUBTRACT_FROM_ADAPTOR_LENGTH"};
	my $alignmentSize = $subjectEnd - $subjectBegin;
	
	if ($alignmentSize >= $threshold) {
	    my $aux = $adaptorArtifacts{$sequenceId};
	    if (! defined $aux) {
		$aux = "$subjectBegin\t$subjectEnd";
	    }
	    else {
		$aux .= "\t$subjectBegin\t$subjectEnd";
	    }
	    $adaptorArtifacts{$sequenceId} = $aux;
	}
    }

    print SCRIPTLOG "done.\n";
}


#############################################################################
# Builds the poly-A artifact list
sub processPolyASwatOutput {
    print SCRIPTLOG "Processing poly-A swat output...";

    open(IN, $temporaryFiles{"polya_swat"});

    my ($sequenceId, $score, $subjectBegin, $subjectEnd);

    while (my $line = <IN>) {
	
	if ($line =~ /^(\d+)\s+Length:\s+\d+/) {
	    if (defined $sequenceId) {
		if ($score >= $configuration{"POLYA_MINIMUM_SCORE"}) {
		    my $aux = $polyAArtifacts{$sequenceId};
		    if (! defined $aux) {
			$aux = "$subjectBegin\t$subjectEnd";
		    }
		    else {
			$aux .= "\t$subjectBegin\t$subjectEnd";
		    }
		    $polyAArtifacts{$sequenceId} = $aux;
		}
	    }
	    $sequenceId = $1;
	    undef($subjectBegin);
	}
	elsif ($line =~ /^\s+Score:\s+(\d+)/) {
	    $score = $1;
	}
	elsif ($line =~ /^Subject\s+(\d+)\s+[ATCGNX\s\-]+(\d+)$/) {
	    if (! defined $subjectBegin) {
		$subjectBegin = $1 - 1;
	    }
	    $subjectEnd = $2;
	}
    }
    close(IN);

    if (defined $sequenceId) {
	if ($score >= $configuration{"POLYA_MINIMUM_SCORE"}) {
	    my $aux = $polyAArtifacts{$sequenceId};
	    if (! defined $aux) {
		$aux = "$subjectBegin\t$subjectEnd";
	    }
	    else {
		$aux .= "\t$subjectBegin\t$subjectEnd";
	    }
	    $polyAArtifacts{$sequenceId} = $aux;
	}
    }

    print SCRIPTLOG "done.\n";
}


#############################################################################
# Builds the adaptor artifact list
sub processPolyTSwatOutput {

    print SCRIPTLOG "Processing poly-T swat output...";

    open(IN, $temporaryFiles{"polyt_swat"});

    my ($sequenceId, $score, $subjectBegin, $subjectEnd);

    while (my $line = <IN>) {
	
	if ($line =~ /^(\d+)\s+Length:\s+\d+/) {
	    if (defined $sequenceId) {
		if ($score >= $configuration{"POLYT_MINIMUM_SCORE"}) {
		    my $aux = $polyTArtifacts{$sequenceId};
		    if (! defined $aux) {
			$aux = "$subjectBegin\t$subjectEnd";
		    }
		    else {
			$aux .= "\t$subjectBegin\t$subjectEnd";
		    }
		    $polyTArtifacts{$sequenceId} = $aux;
		}
	    }
	    $sequenceId = $1;
	    undef($subjectBegin);
	}
	elsif ($line =~ /^\s+Score:\s+(\d+)/) {
	    $score = $1;
	}
	elsif ($line =~ /^Subject\s+(\d+)\s+[ATCGNX\s\-]+(\d+)$/) {
	    if (! defined $subjectBegin) {
		$subjectBegin = $1 - 1;
	    }
	    $subjectEnd = $2;
	}
    }
    close(IN);

    if (defined $sequenceId) {
	if ($score >= $configuration{"POLYT_MINIMUM_SCORE"}) {
	    my $aux = $polyTArtifacts{$sequenceId};
	    if (! defined $aux) {
		$aux = "$subjectBegin\t$subjectEnd";
	    }
	    else {
		$aux .= "\t$subjectBegin\t$subjectEnd";
	    }
	    $polyTArtifacts{$sequenceId} = $aux;
	}
    }

    print SCRIPTLOG "done.\n";
}



#############################################################################
# Clears the fasta file and build auxiliar hashes
sub clearFastaFile {

    my $clearedFasta = createTemporaryFileName();
    my $key = "cleared_fasta_0";
    $temporaryFiles{$key} = $clearedFasta;

    $numberOfFASTAs = 1;

    my $n = 0;
    my $sequenceId = 0;
    my $sequenceName;
    my $sequence = "";

    open(IN, $fastaFile);
    open(OUT, ">$clearedFasta");

    while (my $line = <IN>) {
	
	chop($line);
	
	if ($line =~ /^>/) {

	    if (defined $sequenceName) {

		if ($sequence eq "") {
		    print STDERR "ERROR: File has empty sequence $sequenceName.\n\n";
		    print SCRIPTLOG "ERROR: File has empty sequence $sequenceName.\n\n";
		    close(OUT); close(IN);
		    close(SCRIPTLOG);
		    removeTemporaryFiles();
		}
		
		$sequence =~ s/\s//g;
		$sequence =~ tr/A-Z/a-z/;

		$sequenceSizes{$sequenceId} = length($sequence);
		
		print OUT ">$sequenceId\n$sequence\n";
	    }
	
	    $sequence = "";

	    $line =~ s/^(> ?[^ ]+).*$/$1/;
	    $sequenceName = $line;

	    $sequenceId++;
	    $n++;

	    if ($n > $configuration{"NUMBER_OF_SEQUENCES_IN_AUXILIAR_FASTA"}) {
		close(OUT);
		$clearedFasta = createTemporaryFileName();
		$key = "cleared_fasta_$numberOfFASTAs";
		$temporaryFiles{$key} = $clearedFasta;
		open(OUT, ">$clearedFasta");
		$numberOfFASTAs++;
		$n = 1;
	    }
	
	    $nameToId{$sequenceName} = $sequenceId;
	    $idToName{$sequenceId} = $sequenceName;
	}
	else {
	    $sequence .= $line;
	
	}
    }

    close(IN);

    if (defined $sequenceName) {
	
	if ($sequence eq "") {
	    print STDERR "ERROR: File has empty sequence $sequenceName.\n\n";
	    print SCRIPTLOG "ERROR: File has empty sequence $sequenceName.\n\n";
	    close(OUT);
	    close(SCRIPTLOG);
	    removeTemporaryFiles();
	}

	$sequence =~ s/\s//g;
	$sequence =~ tr/A-Z/a-z/;
	
	$sequenceSizes{$sequenceId} = length($sequence);
	
	print OUT ">$sequenceId\n$sequence\n";
    }

    $numberOfSequences = $sequenceId;

    close(OUT);
}



#############################################################################
# Clears the qual file and verifies it compatibility with the fasta file
sub clearQualFile {

    if (defined $qualFile) {

	my $clearedQual = createTemporaryFileName();
	my $key = "cleared_qual_0";
	$temporaryFiles{$key} = $clearedQual;

	my $nQualFile = 1;
	my $n = 0;
	my $index = 0;
	my $sequenceName;
	my $sequence = "";

	open(IN, $qualFile);
	open(OUT, ">$clearedQual");

	while (my $line = <IN>) {

	    chop($line);

	    if ($line =~ /^>/) {
		
		if (defined $sequenceName) {
		
		    if ($sequence eq "") {
			print STDERR "ERROR: File has empty quality sequence: $sequenceName.\n\n";
			print SCRIPTLOG "ERROR: File has empty quality sequence: $sequenceName.\n\n";
			close(OUT); close(IN);
			close(SCRIPTLOG);
			removeTemporaryFiles();
		    }
		
		    $sequence =~ s/\s+/ /g;
		    $sequence =~ s/^ //;
		    $sequence =~ s/ $//;
		
		    if ($sequence =~ /[^ 0-9]/) {
			print STDERR "ERROR: Quality sequence has invalid characters: $sequenceName.\n\n";
			print SCRIPTLOG "ERROR: Quality sequence has invalid characters: $sequenceName.\n\n";
			close(OUT); close(IN);
			close(SCRIPTLOG);
			removeTemporaryFiles();
		    }
		
		    my @qual = split(" ", $sequence);
		    my $length = @qual;

		    if ($length != $sequenceSizes{$index}) {
			print STDERR "ERROR: Fasta and Qual sequences have different lengths: $sequenceName.\n\n";
			print SCRIPTLOG "ERROR: Fasta and Qual sequences have different lengths: $sequenceName.\n\n";
			close(OUT); close(IN);
			close(SCRIPTLOG);
			removeTemporaryFiles();
		    }
		
		    print OUT ">$index\n$sequence\n";

		}
		
		$sequence = "";

		$line =~ s/^(> ?[^ ]+).*$/$1/;
		$sequenceName = $line;

		my $id = $nameToId{$sequenceName};

		$index++;
		
		if (! defined $id) {
		    print STDERR "ERROR: Sequence found in Qual file was not found in Fasta file: $sequenceName.\n\n";
		    print SCRIPTLOG "ERROR: Sequence found in Qual file was not found in Fasta file: $sequenceName.\n\n";
		    close(OUT); close(IN);
		    close(SCRIPTLOG);
		    removeTemporaryFiles();
		}
		
		if ($id != $index) {
		    print STDERR "ERROR: Fasta and Qual files have different orders or different numbers of sequences.\n\n";
		    print SCRIPTLOG "ERROR: Fasta and Qual files have different orders or different numbers of sequences.\n\n";
		    close(OUT); close(IN);
		    close(SCRIPTLOG);
		    removeTemporaryFiles();
		}
		
		$n++;

		if ($n > $configuration{"NUMBER_OF_SEQUENCES_IN_AUXILIAR_FASTA"}) {
		    close(OUT);
		    $clearedQual = createTemporaryFileName();
		    $key = "cleared_qual_$nQualFile";
		    $temporaryFiles{$key} = $clearedQual;
		    open(OUT, ">$clearedQual");
		    $nQualFile++;
		    $n = 1;
		}
	    }
	    else {
		$sequence .= " $line";
	    } # if ($line =~ /^>/) {...} else {...}
	
	} # while (my $line = <IN>) {...}
	
	close(IN);
	
	if (defined $sequenceName) {
	
	    if ($sequence eq "") {
		print STDERR "ERROR: File has empty quality sequence $sequenceName.\n\n";
		print SCRIPTLOG "ERROR: File has empty quality sequence $sequenceName.\n\n";
		close(OUT);
		close(SCRIPTLOG);
		removeTemporaryFiles();
	    }
	
	    $sequence =~ s/\s+/ /g;
	    $sequence =~ s/^ //;
	    $sequence =~ s/ $//;
	
	    if ($sequence =~ /[^ 0-9]/) {
		print STDERR "ERROR: Quality sequence has invalid characters $sequenceName.\n\n";
		print SCRIPTLOG "ERROR: Quality sequence has invalid characters $sequenceName.\n\n";
		close(OUT);
		close(SCRIPTLOG);
		removeTemporaryFiles();
	    }
	
	    my @qual = split(" ", $sequence);
	    my $length = @qual;

	    if ($length != $sequenceSizes{$index}) {
		print STDERR "ERROR: Fasta and Qual sequences have different lengths: $sequenceName.\n\n";
		print SCRIPTLOG "ERROR: Fasta and Qual sequences have different lengths: $sequenceName.\n\n";
		close(OUT);
		close(SCRIPTLOG);
		removeTemporaryFiles();
	    }
	
	    print OUT ">$index\n$sequence\n";
	}
	
	close(OUT);
	
	if ($index != $numberOfSequences) {
	    print STDERR "ERROR: Fasta and Qual files have different numbers of sequences.\n\n";
	    print SCRIPTLOG "ERROR: Fasta and Qual files have different numbers of sequences.\n\n";
	    close(SCRIPTLOG);
	    removeTemporaryFiles();
	}
    }
}



#############################################################################
# Detects slippage artifacts in the sequence passed by parameter
sub detectSlippageArtifacts {

    my ($sequenceId, $sequence) = @_;

    my $sequenceSize = length($sequence);
	
    my $minimum_echo_size = $configuration{"SLIPPAGE_MINIMUM_ECHO_SIZE"};
    my $minimum_number_of_echoes = $configuration{"SLIPPAGE_MINIMUM_NUMBER_OF_ECHOES"};
    my $threshold_value = $configuration{"SLIPPAGE_THRESHOLD_VALUE"};
    my $algorithm = $configuration{"SLIPPAGE_TRIMMING_ALGORITHM"};
    my $strategy = $configuration{"SLIPPAGE_TRIMMING_STRATEGY"};
    
    if ($strategy eq "suffix") {
	# Find the longest slipped suffix
 	my ($score, $initialPosition) = 
	    findLongestSlippedSuffix($algorithm, $minimum_echo_size,
				     $minimum_number_of_echoes, $threshold_value,
				     $sequence);
	if ($score > 0) {
 	    $slippageArtifacts{$sequenceId} = "$initialPosition\t$sequenceSize";
 	}
    }
    else {
	
	my $originalSequence = $sequence;
	my $initialPositionCorrection = 0;
	my $currentInitialPosition = -1;
	my $currentRegionSize      = 0;
	my $currentEndPosition     = 0;
	
	while (1) {
	    
	    # Find the first slipped region in the sequence
	    my ($score, $initialPosition, $regionSize) = 
		findFirstSlippedRegion($algorithm, $minimum_echo_size,
				       $minimum_number_of_echoes, $threshold_value, 
				       $sequence);
	    
	    # There are no echoed groups
	    if($initialPosition == -1){
		last;
	    }

	    # Correcting the initial position
	    my $newSubsequenceInitialPosition = $initialPositionCorrection + $initialPosition;
	    
	    # Removing non-echoed region at sequence begin
	    $sequence = substr($originalSequence, $newSubsequenceInitialPosition,
			       $sequenceSize - $newSubsequenceInitialPosition + 1);
	
	    # Removing the echoed group to accelerate the process
	    $sequence =~ s/(a{1,}|t{1,}|c{1,}|g{1,}|n{1,}|x{1,})//;
	    
	    # Update the position correction factor
	    $initialPositionCorrection = $newSubsequenceInitialPosition + length($1);
	    
	    # A slipped region was found
	    if ($score > 0) {
		
		if ($currentInitialPosition == -1) {

		    # We have a new slipped region
		    $currentInitialPosition = $newSubsequenceInitialPosition;
		    $currentRegionSize = $regionSize;
		    $currentEndPosition = $currentInitialPosition + $currentRegionSize;
		    
		} else {
		    
		    # We had already found a region. We must check for overlapping occurence
		    if ($currentEndPosition < $newSubsequenceInitialPosition) {
			
			# No overlapping
			# Adds the previous artifact in the list and starts a new one
			
			my $finalPosition = $currentInitialPosition + $currentRegionSize;
			
			my $aux = $slippageArtifacts{$sequenceId};
			if (! defined $aux) {
			    $aux = "$currentInitialPosition\t$finalPosition";
			}
			else {
			    $aux .= "\t$currentInitialPosition\t$finalPosition";
			}
			$slippageArtifacts{$sequenceId} = $aux;
			
			# New slipped region
			$currentInitialPosition = $newSubsequenceInitialPosition;
			$currentRegionSize = $regionSize;
			$currentEndPosition = $currentInitialPosition + $currentRegionSize;
		    }
		    else {
			
			# Overlapping
			my $newEndPosition = $newSubsequenceInitialPosition + $regionSize;
			
			# Checking if the new slipped region end is
			# farter than the old slipped region end
			
			if ($currentEndPosition != $newEndPosition) {
			    # The slipped region grows
			    $currentEndPosition = $newEndPosition;
			    $currentRegionSize = $currentEndPosition - $currentInitialPosition;
			}
		    }

		} # if ($currentInitialPosition == -1) {...} else {...}
		
	    } # if ($score > 0) {...}
	    
	} # while (1) {...}
	
	if ($currentInitialPosition != -1) {
	    
	    my $finalPosition = $currentInitialPosition + $currentRegionSize;
	    my $aux = $slippageArtifacts{$sequenceId};
	    if (! defined $aux) {
		$aux = "$currentInitialPosition\t$finalPosition";
	    }
	    else {
		$aux .= "\t$currentInitialPosition\t$finalPosition";
	    }
	    $slippageArtifacts{$sequenceId} = $aux;
	}
	
    } # if ($strategy eq "suffix") {...} else {...}
}


#############################################################################
# Detects the first slipped region and returns an array with score, 
# region first base index and region size. Returns (0, -1, 0)
# when nothing is found.
sub findFirstSlippedRegion {

    my ($algorithm, $minimum_echo_size, 
	$minimum_number_of_echoes, $threshold_value,
	$sequence) = @_;
    
    my $initialPosition = -1;

    # Searching for the first echoed group of the sequence
    # Note that the first non-echoed groups are ignored
    # (We are searching subsequences that aren't anchored to 
    # the sequence extremity)
    my $auxInitialPosition = 0;

    while ($sequence =~ /(a{1,}|t{1,}|c{1,}|g{1,}|n{1,})/g) {
	
	my $group = $1;
	my $groupSize = length($group);

	if ( ($groupSize >= $minimum_echo_size) && !($group =~ m/n+/) ) {
	    $initialPosition = $auxInitialPosition;
	    last;
	}
	$auxInitialPosition += $groupSize;
    }

    # If an echoed group has been found ...
    if ($initialPosition != -1) {

	# Removing non-echoed region at sequence begin
	$sequence = substr($sequence, $initialPosition,
			   length($sequence) - $initialPosition + 1);

	my ($score, $regionSize);
	
	if ($algorithm eq "am") {
	    ($score, $regionSize) = 
		findLongestPrefix_ArithmeticMean($minimum_echo_size, 
						 $minimum_number_of_echoes, 
						 $threshold_value,
						 $sequence);
	}
	elsif ($algorithm eq "gm") {	
	    ($score, $regionSize) = 
		findLongestPrefix_GeometricMean($minimum_echo_size, 
						$minimum_number_of_echoes, 
						$threshold_value,
						$sequence);
	}
	else {
	    ($score, $regionSize) = 
		findLongestPrefix_EchoCoverage($minimum_echo_size, 
					       $minimum_number_of_echoes, 
					       $threshold_value,
					       $sequence);
	}
	
	return ($score, $initialPosition, $regionSize);
	
    } # if ($initialPosition != -1) {...}

    #Nao encontrou nenhum grupo ecoado
    return (0, -1, 0);
}


#############################################################################
# Detects the longest slipped sequence suffix and returns an array with
# score and the region first base index. Returns zero score 
# when nothing is found.
sub findLongestSlippedSuffix {
    
    my ($algorithm, $minimum_echo_size, 
	$minimum_number_of_echoes, $threshold_value,
	$sequence) = @_;

    my $sequenceSize = length($sequence);
    
    $sequence = reverse($sequence);

    # Note that, in this case, the non-echoed bases at the sequence begin 
    # are considered (different of findFirstSlippedRegion())
    # (We are searching subsequences that are anchored to 
    # the sequence extremity)

    my ($score, $regionSize);

    if ($algorithm eq "am") {
	($score, $regionSize) = 
	    findLongestPrefix_ArithmeticMean($minimum_echo_size, 
					     $minimum_number_of_echoes, 
					     $threshold_value,
					     $sequence);
    }
    elsif ($algorithm eq "gm") {	
	($score, $regionSize) = 
	    findLongestPrefix_GeometricMean($minimum_echo_size, 
					    $minimum_number_of_echoes, 
					    $threshold_value,
					    $sequence);
    }
    else {
	($score, $regionSize) = 
	    findLongestPrefix_EchoCoverage($minimum_echo_size, 
					   $minimum_number_of_echoes, 
					   $threshold_value,
					   $sequence);
    }
    
    my $initialPosition = $sequenceSize - $regionSize;
    return ($score, $initialPosition);
}

#############################################################################
# Detects the longest slipped sequence prefix and returns an array with
# score and slipped region size. Returns zero score when nothing is found.
# Implements Echo Coverage Method
sub findLongestPrefix_EchoCoverage {

    my ($minimum_echo_size, $minimum_number_of_echoes, 
	$threshold_value, $sequence) = @_;

    my $regionSize     = 0;
    my $numberOf1s     = 0;
    my $numberOfGroups = 0;
    my $score           = 0;
    my $scoreRegionSize = 0;

    while ($sequence =~ /(a{1,}|t{1,}|c{1,}|g{1,}|n{1,})/g) {
	
	my $group = $1;
	my $groupSize = length($group);
	
	$regionSize += $groupSize;
	
	$numberOfGroups++;
	
	if (($groupSize >= $minimum_echo_size) && !($group =~ m/n+/)) {
	    
	    $numberOf1s++;
	    
	    if ($numberOf1s >= $minimum_number_of_echoes) {
		
		my $auxScore = $numberOf1s / $numberOfGroups;
		if ($auxScore >= $threshold_value) {
		    $score = $auxScore;
		    $scoreRegionSize = $regionSize;
		}
	    }
	}
    }
    return ($score, $scoreRegionSize);
}

#############################################################################
# Detects the longest slipped sequence prefix and returns an array with
# score and slipped region size. Returns zero score when nothing is found.
# Implements Arithmetic Mean Method
sub findLongestPrefix_ArithmeticMean {

    my ($minimum_echo_size, $minimum_number_of_echoes, 
	$threshold_value, $sequence) = @_;

    my $regionSize     = 0;
    my $numberOfGroups = 0;
    my $echoedGroupSum = 0;
    my $numberOfEchoedGroups = 0;
    
    my $score = 0;
    my $scoreRegionSize = 0;
    
    while ($sequence =~ /(a{1,}|t{1,}|c{1,}|g{1,}|n{1,})/g) {
	
	my $group = $1;
	my $groupSize = length($group);
	
	$regionSize += $groupSize;
	
	$numberOfGroups++;
	
	if (($groupSize >= $minimum_echo_size) && !($group =~ m/n+/)) {

	    $echoedGroupSum += $groupSize;
	    $numberOfEchoedGroups++;
	    
	    if ($numberOfEchoedGroups >= $minimum_number_of_echoes) {
		
		my $auxScore = $echoedGroupSum / $numberOfGroups;
		if ($auxScore >= $threshold_value) {
		    $score = $auxScore;
		    $scoreRegionSize = $regionSize;
		}
	    }
	}
    }
    return ($score, $scoreRegionSize);
}


#############################################################################
# Detects the longest slipped sequence prefix and returns an array with
# score and slipped region size. Returns zero score when nothing is found.
# Implements Geometric Mean Method
sub findLongestPrefix_GeometricMean {

    my ($minimum_echo_size, $minimum_number_of_echoes, 
	$threshold_value, $sequence) = @_;
    
    my $regionSize     = 0;
    my $numberOfGroups = 0;
    my $echoedGroupProduct = 1;
    my $numberOfEchoedGroups = 0;
    
    my $score = 0;
    my $scoreRegionSize = 0;
    
    while ($sequence =~ /(a{1,}|t{1,}|c{1,}|g{1,}|n{1,})/g) {
	
	my $group = $1;
	my $groupSize = length($group);
	
	$regionSize += $groupSize;
	
	$numberOfGroups++;
	
	if (($groupSize >= $minimum_echo_size) && !($group =~ m/n+/)) {

	    $echoedGroupProduct *= $groupSize;
	    $numberOfEchoedGroups++;
	    
	    if ($numberOfEchoedGroups >= $minimum_number_of_echoes) {
		
		my $auxScore = $echoedGroupProduct**(1/$numberOfGroups);
		if ($auxScore >= $threshold_value) {
		    $score = $auxScore;
		    $scoreRegionSize = $regionSize;
		}
	    }
	}
    }
    return ($score, $scoreRegionSize);
}


######################################################################
sub qualityAnalysis {

    my ($sequenceId, @qualities) =  @_;
    
    my $minimum_sequence_size = $configuration{"MINIMUM_SEQUENCE_SIZE"};
    my $low_quality_trimming_algorithm = $configuration{"LOW_QUALITY_TRIMMING_ALGORITHM"};

    # 1st step - Search maximum subsequence
    

    my $sequenceSize = @qualities;

    my ($subsequenceSize_1, $begin5_1, $end5_1, $begin3_1, $end3_1);

    if ($low_quality_trimming_algorithm eq "ms" || $low_quality_trimming_algorithm eq "MS") {
	($subsequenceSize_1, $begin5_1, $end5_1, $begin3_1, $end3_1) =
	    maximumSubsequenceMethod(0, $sequenceSize, @qualities);
    }
    else {
	($subsequenceSize_1, $begin5_1, $end5_1, $begin3_1, $end3_1) =
	    sliddingWindowMethod(0, $sequenceSize, @qualities);
    }

    if ($configuration{"LOW_QUALITY_SEARCH_FOR_ISLANDS"} == 1) {
	
	# 2nd step - Low quality island detection
	if ($subsequenceSize_1 < $minimum_sequence_size) {
	    
	    # Too small to check for low quality island
	    my $aux;
	    if ($end5_1 ne "--") {
		$aux = "$begin5_1\t$end5_1";
	    }
	    if ($begin3_1 ne "--") {
		if (defined $aux) {
		    $aux .= "\t$begin3_1\t$end3_1";
		}
		else {
		    $aux = "$begin3_1\t$end3_1";
		}
	    }
	    if (defined $aux) {
		$lowQualityArtifacts{$sequenceId} = $aux;
	    }
	    
	}
	else{

	    my $subSeqBegin = 0;
	    my $subSeqEnd = $sequenceSize;
	    
	    if ($end5_1 ne "--") {
		$subSeqBegin = $end5_1;
	    }
	    
	    if ($begin3_1 ne "--") {
		$subSeqEnd = $begin3_1;
	    }
	
	    # Getting the low quality island cut positions
	    my @indexes = badWindowDetection($subSeqBegin, $subSeqEnd, @qualities);
	    
	    my $numberOfIndexes = @indexes;
	
	    if ($numberOfIndexes == 0) {
		
		# No low quality island has been found
		my $aux;
		if ($end5_1 ne "--") {
		    $aux = "$begin5_1\t$end5_1";
		}
		if ($begin3_1 ne "--") {
		    if (defined $aux) {
			$aux .= "\t$begin3_1\t$end3_1";
		    }
		    else {
			$aux = "$begin3_1\t$end3_1";
		    }
		}
		if (defined $aux) {
		    $lowQualityArtifacts{$sequenceId} = $aux;
		}
	    }
	    else {
		
		# At least one low quality island has been found
		my $newGoodBegin = $subSeqBegin;
		my $newGoodEnd = $subSeqEnd;
		my $newGoodSize = 0;
		
		@indexes = ($subSeqBegin, @indexes, $subSeqEnd);
		$numberOfIndexes += 2;

		# Processes the all fragments and preserves only the longest
		for (my $i = 1; $i < $numberOfIndexes; $i++) {
		    
		    my $begin = $indexes[$i-1];
		    my $end = $indexes[$i];
		    
		    my ($subsequenceSize_2, $begin5_2, $end5_2, $begin3_2, $end3_2);
		    if ($low_quality_trimming_algorithm eq "ms" || 
			$low_quality_trimming_algorithm eq "MS") {
			my $minimum_quality = $configuration{"MAXIMUM_SUBSEQUENCE_MINIMUM_QUALITY"};
			my $maximumErrorProbability = $errorProbabilityArray[$minimum_quality];
			($subsequenceSize_2, $begin5_2, $end5_2, $begin3_2, $end3_2) =
			    maximumSubsequenceMethod($begin, $end, @qualities);
		    }
		    else {
			($subsequenceSize_2, $begin5_2, $end5_2, $begin3_2, $end3_2) =
			    sliddingWindowMethod($begin, $end, @qualities);
		    }
		    
		    if ($end5_2 ne "--") {
			# 5' end was trimmed. Update sequence extremity
			$begin = $end5_2;
		    }
		    if ($begin3_2 ne "--") {
			# 3' end was trimmed. Update sequence extremity
			$end = $begin3_2;
		    }
		    
		    # If this fragment is greater than the previous ones,
		    # it will be the new good fragment
		    if ($subsequenceSize_2 >= $newGoodSize) {
			$newGoodSize = $subsequenceSize_2;
			$newGoodBegin = $begin;
			$newGoodEnd = $end;
		    }
		}
		
		my $aux;
		if ($newGoodBegin != 0) {
		    $aux = "0\t$newGoodBegin";
		}
		
		if ($newGoodEnd != $sequenceSize) {
		    if (defined $aux) {
			$aux .= "\t$newGoodEnd\t$sequenceSize";
		    }
		    else {
			$aux = "$newGoodEnd\t$sequenceSize";
		    }
		}
		if (defined $aux) {
		    $lowQualityArtifacts{$sequenceId} = $aux;
		}
		
	    } # if ($numberOfIndexes == 0) {...} else {...}
	    
	}# if ($subsequenceSize_1 < $minimum_sequence_size) {...} else {...}
	
    }
    else {
	
	my $aux;
	if ($end5_1 ne "--") {
	    $aux = "$begin5_1\t$end5_1";
	}
	if ($begin3_1 ne "--") {
	    if (defined $aux) {
		$aux .= "\t$begin3_1\t$end3_1";
	    }
	    else {
		$aux = "$begin3_1\t$end3_1";
	    }
	}
	if (defined $aux) {
	    $lowQualityArtifacts{$sequenceId} = $aux;
	}
	
    } # if ($configuration{"LOW_QUALITY_SEARCH_FOR_ISLANDS"} == 1) {...} else {...}
    
}



######################################################################
# Maximum subsequence low quality detection method
sub maximumSubsequenceMethod {

    my ($indexStart, $indexFinish, @qualities) = @_;

    my $minimum_quality = $configuration{"MAXIMUM_SUBSEQUENCE_MINIMUM_QUALITY"};
    my $maximumErrorProbability = $errorProbabilityArray[$minimum_quality];


    my $maxGlobal = 0;
    my $maxGlobalBegin = $indexStart;
    my $maxGlobalEnd = $indexStart;

    my $maxSuffix = 0;
    my $maxSuffixBegin = $indexStart;
    my $auxSum = 0;

    for (my $i = $indexStart; $i < $indexFinish; $i ++) {
	
	$auxSum = $maxSuffix +  $maximumErrorProbability - $errorProbabilityArray[$qualities[$i]];
	
	if ( $auxSum > $maxGlobal ) {
	
	    $maxSuffix = $auxSum;
	    $maxGlobal = $maxSuffix;
	    $maxGlobalBegin = $maxSuffixBegin;
	    $maxGlobalEnd = $i + 1;
	
	}
	else {
	
	    if ($auxSum <= 0) {
		$maxSuffix = 0;
		$maxSuffixBegin = $i + 1;
	    }
	    else {
		$maxSuffix = $auxSum;
	    }
	}
    }

    my $left1  = "--";
    my $right1 = "--";
    my $left2  = "--";
    my $right2 = "--";
    my $subsequenceSize = 0;

    if ($maxGlobalBegin == $indexStart && $maxGlobalEnd == $indexStart) {

	$left1 = $indexStart;
	$right1 = $indexFinish;
	
    }
    else {
	
	$subsequenceSize = $maxGlobalEnd - $maxGlobalBegin;
	if ($maxGlobalBegin > $indexStart) {
	    $left1 = $indexStart;
	    $right1 = $maxGlobalBegin;
	}
	if ($maxGlobalEnd < $indexFinish ) {
	    $left2 = $maxGlobalEnd;
	    $right2 = $indexFinish;
	}
    }
    return ($subsequenceSize, $left1, $right1, $left2, $right2);
}



######################################################################
# Slidding window low quality detection method
sub sliddingWindowMethod {

    my ($indexStart, $indexFinish, @qualities) =  @_;

    my $window_size = $configuration{"SLIDDING_WINDOW_WINDOW_SIZE"};
    my $quality_threshold = $configuration{"SLIDDING_WINDOW_QUALITY_THRESHOLD"};
    my $bad_bases_threshold = $configuration{"SLIDDING_WINDOW_BAD_BASES_THRESHOLD"};

    my $qualitiesSize = @qualities;
    
    my $numberOfBadBases = 0;
    
    for (my $i = $indexStart; $i < $window_size && $i < $indexFinish; $i++) {
	if($qualities[$i] < $quality_threshold) {
	    $numberOfBadBases++;
	}
    }

    my $leftCut = $indexStart;
    my $rightCut = $indexFinish - 1;

    if ($window_size <= $indexFinish) {

	my $limit = $indexFinish - $window_size;
	for (my $i = $indexStart; $i < $limit && $numberOfBadBases > $bad_bases_threshold; $i ++) {
	    $leftCut++;
	    if ($qualities[$i] < $quality_threshold) {
		$numberOfBadBases--;
	    }
	    if ($qualities[$i + $window_size] < $quality_threshold) {
		$numberOfBadBases++;
	    }
	}

	if ($leftCut < $limit) {
    
	    $numberOfBadBases = 0;
	    
	    for (my $i = $indexFinish - 1; $i >= $limit; $i--) {
		if ($qualities[$i] < $quality_threshold) {
		    $numberOfBadBases++;
		}
	    }

	    $limit = $leftCut + $window_size - 1;

	    for (my $i = $indexFinish - 1; 
		 $i > $limit &&  $numberOfBadBases > $bad_bases_threshold && $rightCut >= $leftCut; 
		 $i--) {

		$rightCut--;
		if ($qualities[$i] < $quality_threshold) {
		    $numberOfBadBases--;
		}
		if ($qualities[$i - $window_size] < $quality_threshold) {
		    $numberOfBadBases++;
		}
	    }
	}
	else {
	    if ($numberOfBadBases <= $bad_bases_threshold) {
		$rightCut = $indexFinish - 1;
	    }
	    else {
		$rightCut = - 1;
	    }
	} # if ($leftCut < $limit) {...} else {...}
    }
    else {
	
	if ($indexFinish - $numberOfBadBases >= $window_size - $bad_bases_threshold) {
	    $leftCut = $indexStart;
	    $rightCut = $indexFinish - 1;
	}
	else {
	    $rightCut = - 1;
	}
    } # if ($window_size <= $indexFinish) {...} else {...}
    
    my $subsequenceSize = 0;
    my $left1 = "--";
    my $right1 = "--";

    my $left2 = "--";
    my $right2 = "--";

    if ($leftCut > $rightCut) {
	$left1 = $indexStart;
	$right1 = $indexFinish;

    }
    else {
	if($leftCut > $indexStart) {
	    $left1 = $indexStart;
	    $right1 = $leftCut;
	}
	if($rightCut < $indexFinish - 1) {
	    $left2 = $rightCut + 1;
	    $right2 = $indexFinish;
	}
    }
    
    my $leftSize = 0;
    my $rightSize = 0;

    if ($right1 ne "--") {

	$leftSize = $right1 - $left1;

	if ($left2 ne "--") {
	    $rightSize = $right2 - $left2;
	    $subsequenceSize = $left2 - $right1;
	}
	else {
	    $subsequenceSize = $indexFinish  - $right1;
	}
    } 
    else {
	if($left2 ne "--") {
	    $rightSize = $right2 - $left2;
	    $subsequenceSize = $left2;
	}
	else {
	    $subsequenceSize = $indexFinish;
	}
    }
    
    return ($subsequenceSize, $left1, $right1, $left2, $right2);
}



######################################################################
# Analyses the quality sequence and returns the low quality cut
# positions
sub badWindowDetection {

    my ($indexStart, $indexFinish, @qualities) = @_;

    my $sumWindow = 0;
    my $averageWindow = 0;
    my $windowBegin = 0;
    my $windowEnd = 0;

    my @begin;
    my @end;
    my @indexes;

    my $window_size = $configuration{"LOW_QUALITY_ISLAND_WINDOW_SIZE"};
    my $minimum_error_probability = $configuration{"LOW_QUALITY_ISLAND_MINIMUM_ERROR_PROBABILITY"};


    # Calculates the average error probability for the first window
    my $limit = $indexStart + $window_size;

    for (my $i = $indexStart; $i < $limit; $i++) {
	$sumWindow += $errorProbabilityArray[$qualities[$i]];
    }
    $averageWindow = $sumWindow / $window_size;


    # If the average is greater than ou equal to the threshold value,
    # a new island is started
    if ($averageWindow >= $minimum_error_probability) {
	$windowBegin = $indexStart;
	$windowEnd = $limit;
    }

    # Searchs for a new low quality island or extends the current one
    $limit = $indexFinish - $window_size;

    for (my $i = $indexStart + 1; $i <= $limit; $i++) {

	my $i1 = $i - 1;
	my $i2 = $i + $window_size;

	$sumWindow -= $errorProbabilityArray[$qualities[$i1]];
	$sumWindow += $errorProbabilityArray[$qualities[$i2 - 1]];
	$averageWindow = $sumWindow / $window_size;

	if ($averageWindow >= $minimum_error_probability) {
	    # Extends or starts a new island

	    if ($windowBegin == 0) {
		# Starts a new one
		$windowBegin = $i;
	    }

	    $windowEnd = $i2;
	}
	else {
	    # The window has average lower than the threshold value

	    if ($windowBegin != 0) {
		# There is a island
		# Registering in the array its limits
		
		my $index = @begin;
		$begin[$index] = $windowBegin;
		$end[$index] = $windowEnd;
		
		$windowBegin = 0;
		$windowEnd = 0;
		
		if ($i2 >=  $limit) {
		    # There is no space left to build a new island
		    last;
		}
		else {

		    # Calculate the average error probability for the
		    # new window
		    $i = $i2;
		    $sumWindow = 0;
		
		    for (my $j = $i2; $j < $i2 + $window_size; $j ++) {
			$sumWindow += $errorProbabilityArray[$qualities[$j]];
		    }
		    $averageWindow = $sumWindow / $window_size;
		    if ($averageWindow >= $minimum_error_probability) {
			$windowBegin = $i2;
			$windowEnd = $i2 + $window_size;
		    }
		
		} # if ($i2 >=  $limit) {...} else {...}
		
	    } # if ($windowBegin != 0) {...}
	
	} # if ($averageWindow >= $minimum_error_probability) {...} else {...}
	
    } # for (my $i = $indexStart + 1; $i <= $limit; $i++) {...}


    if ($windowBegin != 0) {
	# There is a island
	# Registering in the array its limits
	my $index = @begin;
	$begin[$index] = $windowBegin;
	$end[$index] = $windowEnd;
    }

    # Searches the minimum quality point of each island
    my $n = @begin;
    for (my $i = 0; $i < $n; $i++) {
	my $min = $qualities[$begin[$i]];
	my $index = $begin[$i];
 	for (my $j = $begin[$i] + 1; $j < $end[$i]; $j++) {
 	    if ($min > $qualities[$j]) {
 		$min = $qualities[$j];
 		$index = $j;
 	    }
 	}
 	$indexes[@indexes] = $index;
    }
    return @indexes;
}



######################################################################
# Fills the error probability array
sub initErrorProbabilityArray {
    for(my $i = 0; $i < 99; $i++) {
	$errorProbabilityArray[$i] = errorProbabilityValue($i);
    }
}



######################################################################
# Calculates the error probability for the quality value specified by
# parameter
sub errorProbabilityValue {
    my ($quality) = @_;
    return 10**(-($quality/10));
}


######################################################################
# Receives a sequence and a list of artifacts and returns the sequence
# with all artifact regions masked
sub maskArtifacts {

    my ($sequence, $sequenceSize, $artifacts) = @_;

    my $maskedSequence = $sequence;

    my @fragments = split("\t", $artifacts);

    for (my $i = 0; $i < @fragments; $i+=2) {
	
	my $begin = $fragments[$i];
	my $end = $fragments[$i + 1];
		
	my $fragmentSize = $end - $begin;
	my $lastFragmentSize = $sequenceSize - $end;
		
	$maskedSequence =~ /^([atcgnx]{0,$begin})([atcgnx]{0,$fragmentSize})([atcgnx]{0,$lastFragmentSize})$/;
	my $fragment1 = $1;
	my $fragment2 = $2;
	my $fragment3 = $3;

	$fragment2 =~ s/./x/g;
	$maskedSequence = "$fragment1$fragment2$fragment3";
    }

    return $maskedSequence;
}



######################################################################
# Analyses the masked sequence and returns the insert
sub extractInsert {

    my ($maskedSequence, $sequenceSize, @qualities) = @_;
    my $minimum_sequence_size = $configuration{"MINIMUM_SEQUENCE_SIZE"};

    my $cutPosition = 0;
    my @fragments;
    do {
	if ($maskedSequence =~ /^([atcgn]+)/) {
	    my $begin = $cutPosition;
	    my $end = $cutPosition + length($1);
	    $cutPosition = $end;
	    $maskedSequence =~ s/^$1//;
	    $fragments[@fragments] = $begin;
	    $fragments[@fragments] = $end;
	} elsif ($maskedSequence =~ /^(x+)/) {
	    my $begin = $cutPosition;
	    my $end = $cutPosition + length($1);
	    $cutPosition = $end;
	    $maskedSequence =~ s/^$1//;
	}
    } while($cutPosition < $sequenceSize);

    my $numberOfFragments = @fragments;

    if ($numberOfFragments > 0) {

	my $goodBegin = -1;
	my $goodEnd = -1;
	my $goodQualitySum = 0;
	my $goodSize = 0;
			
	for (my $i = 0; $i < $numberOfFragments; $i += 2) {
	
	    my $begin = $fragments[$i];
	    my $end = $fragments[$i + 1];
	
	    my $size = $end - $begin;
	    if ($size >= $minimum_sequence_size) {

		my $sum = 0;
		if ($configuration{"PERFORM_LOW_QUALITY_DETECTION"} == 1) {
		    for (my $j = $begin; $j < $end; $j++) {
			$sum += $qualities[$j];
		    }
		}

		if (($sum > $goodQualitySum) ||
		    ($sum == $goodQualitySum && $goodSize < $size)) {
		    $goodQualitySum = $sum;
		    $goodBegin = $begin;
		    $goodEnd = $end;
		    $goodSize = $size;
		}
	    }
	}
	return ($goodBegin, $goodEnd);
    }
    return (-1, -1);
}

######################################################################
sub createPolyTemplate {

    my ($size, $letter) = @_;
    my $key = "poly_".$letter."_template";
    my $poly_file = createTemporaryFileName();
    $temporaryFiles{$key} = $poly_file;
    
    my $string = "";
    for(my $i = 0; $i < $size; $i++) {
	$string .= $letter;
    }
    
    open(FILE, ">$poly_file");
    print FILE ">poly-$letter\n$string\n";
    close(FILE);
}


#############################################################################
# SUBROUTINES SECTION END                                                   #
#############################################################################
