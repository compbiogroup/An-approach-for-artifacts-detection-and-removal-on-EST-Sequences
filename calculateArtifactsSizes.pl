#!/usr/bin/perl -w



my $DEFAULT_FILENAME = "./BD2006TRIMMER.LOG";

my $file = $ARGV[0];

printUsage();

if (! defined $file) {
    $file = $DEFAULT_FILENAME;
}

if (! -e $file) {
    print STDERR "\nERROR: $file not found.\n\n";
    exit;
}

my @vector;
my @univec;
my @lowQuality5;
my @lowQuality3;
my @adaptor;
my @ribosomal;
my @polyA;
my @polyT;
my @slippage;

my $sum_vector = 0;
my $sum_univec = 0;
my $sum_lowQuality5 = 0;
my $sum_lowQuality3 = 0;
my $sum_adaptor = 0;
my $sum_ribosomal = 0;
my $sum_polyA = 0;
my $sum_polyT = 0;
my $sum_slippage = 0;

my $num_seq_vector = 0;
my $num_seq_univec = 0;
my $num_seq_lowQuality5 = 0;
my $num_seq_lowQuality3 = 0;
my $num_seq_adaptor = 0;
my $num_seq_ribosomal = 0;
my $num_seq_polyA = 0;
my $num_seq_polyT = 0;
my $num_seq_slippage = 0;

my $number_of_sequences;
my $number_of_discarded_sequences;

open(IN, "$file");

while (my $line =<IN>) {

    chop($line);
    
    if ($line =~ /^\[\d+\]/) {
	$number_of_sequences++;

    } elsif ($line =~ /^\[AD\]/) {
	$num_seq_adaptor++;
	my @aux = split("\t", $line);
	for (my $i = 1; $i < @aux; $i +=2) {
	    my $size = $aux[$i+1] - $aux[$i];
	    $adaptor[@adaptor] = $size;
	    $sum_adaptor += $size;
	}

    } elsif ($line =~ /^\[LQ\]/) {
	my @aux = split("\t", $line);
	if (! defined $aux[3]) {
	    my $size = $aux[2] - $aux[1];
	    if ($aux[1] == 0) {
		$num_seq_lowQuality5++;
		$lowQuality5[@lowQuality5] = $size;
		$sum_lowQuality5 += $size;
	    } else {
		$num_seq_lowQuality3++;
		$lowQuality3[@lowQuality3] = $size;
		$sum_lowQuality3 += $size;
	    }
	} else {
	    $num_seq_lowQuality5++;
	    $num_seq_lowQuality3++;
	    
	    my @aux = split("\t", $line);
	    
	    my $size1 = $aux[2] - $aux[1];
	    my $size2 = $aux[4] - $aux[3];
	    
	    $lowQuality5[@lowQuality5] = $size1;
	    $lowQuality3[@lowQuality3] = $size2;
	    
	    $sum_lowQuality5 += $size1;
	    $sum_lowQuality3 += $size2;
	}

    } elsif ($line =~ /^\[PA\]/) {
	$num_seq_polyA++;
	my @aux = split("\t", $line);
	for (my $i = 1; $i < @aux; $i +=2) {
	    my $size = $aux[$i+1] - $aux[$i];
	    $polyA[@polyA] = $size;
	    $sum_polyA += $size;
	}

    } elsif ($line =~ /^\[PT\]/) {
	$num_seq_polyT++;
	my @aux = split("\t", $line);
	for (my $i = 1; $i < @aux; $i +=2) {
	    my $size = $aux[$i+1] - $aux[$i];
	    $polyT[@polyT] = $size;
	    $sum_polyT += $size;
	}

    } elsif ($line =~ /^\[RI\]/) {
	$num_seq_ribosomal++;
	my @aux = split("\t", $line);
	
	for (my $i = 1; $i < @aux; $i +=2) {
	    my $size = $aux[$i+1] - $aux[$i];
	    $ribosomal[@ribosomal] = $size;
	    $sum_ribosomal += $size;
	}

    } elsif ($line =~ /^\[SL\]/) {
	$num_seq_slippage++;
	my @aux = split("\t", $line);
	for (my $i = 1; $i < @aux; $i +=2) {
	    my $size = $aux[$i+1] - $aux[$i];
	    $slippage[@slippage] = $size;
	    $sum_slippage += $size;
	}

    } elsif ($line =~ /^\[UV\]/) {
	$num_seq_univec++;
	my @aux = split("\t", $line);
	for (my $i = 1; $i < @aux; $i +=2) {
	    my $size = $aux[$i+1] - $aux[$i];
	    $univec[@univec] = $size;
	    $sum_univec += $size;
	}

    } elsif ($line =~ /^\[VE\]/) {
	$num_seq_vector++;
	my @aux = split("\t", $line);
	for (my $i = 1; $i < @aux; $i +=2) {
	    my $size = $aux[$i+1] - $aux[$i];
	    $vector[@vector] = $size;
	    $sum_vector += $size;
	}
    } elsif ($line =~ /^\[DISCARDED\]/) {
	$number_of_discarded_sequences++;
    }
}
close(IN);

print "NUMBER OF PROCESSED SEQUENCES : $number_of_sequences\n";
print "NUMBER OF DISCARDED SEQUENCES : $number_of_discarded_sequences\n\n";


my $n_adaptor = @adaptor;
my $average_adaptor = 0;
my $sd_adaptor = 0;
if ($n_adaptor > 0) {
    $average_adaptor = $sum_adaptor / $n_adaptor;
    $sd_adaptor = calcSTDEV($average_adaptor, $n_adaptor, @adaptor);
    $average_adaptor = sprintf("%0.3f", $average_adaptor);
    $sd_adaptor = sprintf("%0.3f", $sd_adaptor);
}

my $n_lowQuality5 = @lowQuality5;
my $average_lowQuality5 = 0;
my $sd_lowQuality5 = 0;
if ($n_lowQuality5 > 0) {
    $average_lowQuality5 = $sum_lowQuality5 / $n_lowQuality5;
    $sd_lowQuality5 = calcSTDEV($average_lowQuality5, $n_lowQuality5, @lowQuality5);
    $average_lowQuality5 = sprintf("%0.3f", $average_lowQuality5);
    $sd_lowQuality5 = sprintf("%0.3f", $sd_lowQuality5);
}

my $n_lowQuality3 = @lowQuality3;
my $average_lowQuality3 = 0;
my $sd_lowQuality3 = 0;
if ($n_lowQuality3 > 0) {
    $average_lowQuality3 = $sum_lowQuality3 / $n_lowQuality3;
    $sd_lowQuality3 = calcSTDEV($average_lowQuality3, $n_lowQuality3, @lowQuality3);
    $average_lowQuality3 = sprintf("%0.3f", $average_lowQuality3);
    $sd_lowQuality3 = sprintf("%0.3f", $sd_lowQuality3);
}

my $n_polyA = @polyA;
my $average_polyA = 0;
my $sd_polyA = 0;
if ($n_polyA > 0) {
    $average_polyA = $sum_polyA / $n_polyA;
    $sd_polyA = calcSTDEV($average_polyA, $n_polyA, @polyA);
    $average_polyA = sprintf("%0.3f", $average_polyA);
    $sd_polyA = sprintf("%0.3f", $sd_polyA);
}

my $n_polyT = @polyT;
my $average_polyT = 0;
my $sd_polyT = 0;
if ($n_polyT > 0) {
    $average_polyT = $sum_polyT / $n_polyT;
    $sd_polyT = calcSTDEV($average_polyT, $n_polyT, @polyT);
    $average_polyT = sprintf("%0.3f", $average_polyT);
    $sd_polyT = sprintf("%0.3f", $sd_polyT);
}

my $n_ribosomal = @ribosomal;
my $average_ribosomal = 0;
my $sd_ribosomal = 0;
if ($n_ribosomal > 0) {
    $average_ribosomal = $sum_ribosomal / $n_ribosomal;
    $sd_ribosomal = calcSTDEV($average_ribosomal, $n_ribosomal, @ribosomal);
    $average_ribosomal = sprintf("%0.3f", $average_ribosomal);
    $sd_ribosomal = sprintf("%0.3f", $sd_ribosomal);
}

my $n_slippage = @slippage;
my $average_slippage = 0;
my $sd_slippage = 0;
if ($n_slippage > 0) {
    $average_slippage = $sum_slippage / $n_slippage;
    $sd_slippage = calcSTDEV($average_slippage, $n_slippage, @slippage);
    $average_slippage = sprintf("%0.3f", $average_slippage);
    $sd_slippage = sprintf("%0.3f", $sd_slippage);
}

my $n_vector = @vector;
my $average_vector = 0;
my $sd_vector = 0;
if ($n_vector > 0) {
    $average_vector = $sum_vector / $n_vector;
    $sd_vector = calcSTDEV($average_vector, $n_vector, @vector);
    $average_vector = sprintf("%0.3f", $average_vector);
    $sd_vector = sprintf("%0.3f", $sd_vector);
}

my $n_univec = @univec;
my $average_univec = 0;
my $sd_univec = 0;
if ($n_univec > 0) {
    $average_univec = $sum_univec / $n_univec;
    $sd_univec = calcSTDEV($average_univec, $n_univec, @univec);
    $average_univec = sprintf("%0.3f", $average_univec);
    $sd_univec = sprintf("%0.3f", $sd_univec);
}


print "\nLEGEND:\n";
print "[1]\tArtifact type\n";
print "[2]\tNumber of sequences that showed at least one artifact\n";
print "[3]\tNumber of observed artifacts\n";
print "[4]\tMean artifact size\n";
print "[5]\tStandard deviation\n\n";

print "[1]               \t[2]\t[3]\t[4]\t[5]\n";
print "ADAPTOR           \t$num_seq_adaptor\t$n_adaptor\t$average_adaptor\t$sd_adaptor\n";

print "LOW QUALITY 5' end\t$num_seq_lowQuality5\t$n_lowQuality5\t$average_lowQuality5\t$sd_lowQuality5\n";

print "LOW QUALITY 3' end\t$num_seq_lowQuality3\t$n_lowQuality3\t$average_lowQuality3\t$sd_lowQuality3\n";

print "POLY-A             \t$num_seq_polyA\t$n_polyA\t$average_polyA\t$sd_polyA\n";

print "POLY-T             \t$num_seq_polyT\t$n_polyT\t$average_polyT\t$sd_polyT\n";

print "RIBOSOMAL          \t$num_seq_ribosomal\t$n_ribosomal\t$average_ribosomal\t$sd_ribosomal\n";

print "SLIPPAGE           \t$num_seq_slippage\t$n_slippage\t$average_slippage\t$sd_slippage\n";

print "VECTOR             \t$num_seq_vector\t$n_vector\t$average_vector\t$sd_vector\n";

print "UNIVEC             \t$num_seq_univec\t$n_univec\t$average_univec\t$sd_univec\n";

print "\n";



############################################################################
# Calculates the standard deviation
sub calcSTDEV {
    
    my ($mean, $size, @values) = @_;
    my $sum = 0;
    for (my $i = 0; $i < $size; $i++) {
	$sum += ($mean - $values[$i])**2;
    }
    return sqrt($sum/$size);
}

############################################################################
# Prints the usage
sub printUsage {
    print "\n----------------------------------------------------------------\n";
    print "Usage: \n";
    print "  perl calculateArtifactMeanSizes.pl <bd2006trimmer_log_file>\n\n";
    print "Parameters:\n";
    print "  bd2006trimmer_log_file : (optional) bd2006trimmer.pl log file\n";
    print "                           if no file is specified, the script\n";
    print "                           will search for BD2006TRIMMER.LOG.\n";
    print "----------------------------------------------------------------\n\n\n";
}
