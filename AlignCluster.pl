#!/usr/bin/perl

use strict;
use warnings;
use Benchmark;
# install: sudo cpan Parallel::ForkManager
use Parallel::ForkManager;
# The firstidx function from List::MoreUtils
use List::MoreUtils qw(firstidx);

# Numero de procesos
my $nprocess          = 100;
#
my $numAtoms;
#
my %Atomic_number = ( '89'  => 'Ac', '13'  => 'Al', '95'  => 'Am', '51'  => 'Sb',
                      '18'  => 'Ar', '33'  => 'As', '85'  => 'At', '16'  => 'S',
                      '56'  => 'Ba', '4'   => 'Be', '97'  => 'Bk', '83'  => 'Bi',
                      '107' => 'Bh', '5'   => 'B', 	'35'  => 'Br', '48'  => 'Cd',
                      '20'  => 'Ca', '98'  => 'Cf',	'6'   => 'C',  '58'  => 'Ce',
                      '55'  => 'Cs', '17'  => 'Cl',	'27'  => 'Co', '29'  => 'Cu',
                      '24'  => 'Cr', '96'  => 'Cm', '110' => 'Ds', '66'  => 'Dy',
                      '105' => 'Db', '99'  => 'Es', '68'  => 'Er', '21'  => 'Sc',
                      '50'  => 'Sn', '38'  => 'Sr', '63'  => 'Eu', '100' => 'Fm',
                      '9'   => 'F',  '15'  => 'P',  '87'  => 'Fr', '64'  => 'Gd',
                      '31'  => 'Ga', '32'  => 'Ge', '72'  => 'Hf', '108' => 'Hs',
                      '2'   => 'He', '1'   => 'H',  '26'  => 'Fe', '67'  => 'Ho',
                      '49'  => 'In', '53'  => 'I',  '77'  => 'Ir', '70'  => 'Yb',
                      '39'  => 'Y',  '36'  => 'Kr', '57'  => 'La', '103' => 'Lr',
                      '3'   => 'Li', '71'  => 'Lu', '12'  => 'Mg', '25'  => 'Mn',
                      '109' => 'Mt', '101' => 'Md', '80'  => 'Hg', '42'  => 'Mo',
                      '60'  => 'Nd', '10'  => 'Ne', '93'  => 'Np', '41'  => 'Nb',
                      '28'  => 'Ni', '7'   => 'N',  '102' => 'No', '79'  => 'Au',
                      '76'  => 'Os', '8'   => 'O', 	'46'  => 'Pd', '47'  => 'Ag',
                      '78'  => 'Pt', '82'  => 'Pb',	'94'  => 'Pu', '84'  => 'Po',
                      '19'  => 'K',  '59'  => 'Pr', '61'  => 'Pm', '91'  => 'Pa',
                      '88'  => 'Ra', '86'  => 'Rn', '75'  => 'Re', '45'  => 'Rh',
                      '37'  => 'Rb', '44'  => 'Ru', '104' => 'Rf', '62'  => 'Sm',
                      '106' => 'Sg', '34'  => 'Se', '14'  => 'Si', '11'  => 'Na',
                      '81'  => 'Tl', '73'  => 'Ta', '43'  => 'Tc', '52'  => 'Te',
                      '65'  => 'Tb', '22'  => 'Ti', '90'  => 'Th', '69'  => 'Tm',
                      '112' => 'Uub','116' => 'Uuh','111' => 'Uuu','118' => 'Uuo',
                      '115' => 'Uup','114' => 'Uuq','117' => 'Uus','113' => 'Uut',
                      '92'  => 'U',  '23'  => 'V',  '74'  => 'W',  '54'  => 'Xe',
                      '30'  => 'Zn', '40'  => 'Zr' );

###################################
# Read files
sub read_file {
	# filename
	my ($input_file) = @_;
	my @array        = ();
	# open file
	open(FILE, "<", $input_file ) || die "Can't open $input_file: $!";
	while (my $row = <FILE>) {
		chomp($row);
		push (@array,$row);
	}
	close (FILE);
	# return array
	return @array;
}
###################################
# Reference file XYZ
sub format_xyz {
	my ($input_file) = @_;
	#
	my @array_coord = ();
	#
	$numAtoms = @$input_file[0];
	#
	my $tam = scalar (@{$input_file});
	for ( my $i = 2 ; $i < $tam ; $i = $i + 1 ){
		if ( length(@$input_file[$i]) > 2) {
			my @array_tabs  = split (/\s+/,@$input_file[$i]);
			my $radii_val;
			if ( exists $Atomic_number{$array_tabs[0]} ) {
				# exists
				$radii_val = $Atomic_number{$array_tabs[0]};
			} else {
				# not exists
				$radii_val = $array_tabs[0] ;
			}
			my $strong = "$radii_val\t$array_tabs[1]\t$array_tabs[2]\t$array_tabs[3]";
			push (@array_coord,$strong);
		}
	}
	return @array_coord;
}
###################################
# Fisher-Yates Algorithm
sub fisher_yates_shuffle {
	my $array = shift;
	my $i = @$array;
	while ( --$i ) {
		my $j = int rand( $i+1 );
		@$array[$i,$j] = @$array[$j,$i];
	}
	return @$array;
}

#sub tridimensional_structure_cluster {
#	my ($input_file) = @_;
	#
#	my @array_coord = ();
	#
#	$numAtoms = @$input_file[0];
	#
#	my $tam = scalar (@{$input_file});
#}

###################################
# Delete repeat data
sub uniq {
	my %seen;
	grep !$seen{$_}++, @_;
}

# # # # # # # # # # # # # # # # # #
# MAIN
#
#print_logo();
# Funcion para el tiempo de ejecucion del programa
#
my ($fileXYZ,$numbCombination) = @ARGV;
if (not defined $fileXYZ) {
	die "\nAlignCluster must be run with:\n\nUsage:\n\tperl AlignCluster [XYZ-file] [Numb-of-Combinations]\n\n\n";
	exit(1);
}
# Numero de combinaciones
if (not defined $numbCombination) {
	die "\nAlignCluster must be run with:\n\nUsage:\n\tperl AlignCluster [XYZ-file] [Numb-of-Combinations]\n\n\n";
	exit(1);
}
# Initial Time
my $tiempo_inicial  = new Benchmark;
my $datestringStart = localtime();
# Read and parse XYZ format
my @data_file    = read_file($fileXYZ);
my @coords_file  = format_xyz (\@data_file);
#
my @elements = ();
my @axis     = ();
#
my @ClutterArray = ();
my @originalPos  = ();
my $count        = 0;
#
# Valor Total
my $totalValue = 0;
#
my $seqOriginal;
foreach my $data_coords (@coords_file) {
	my @Cartesians              = split '\s+', $data_coords;
	my ($Atom_label, @orig_xyz) = @Cartesians;
	$seqOriginal.="$Atom_label\t";
	#
  push (@elements,$Atom_label);
	push (@axis,"$orig_xyz[0]\t$orig_xyz[1]\t$orig_xyz[2]");
	push (@ClutterArray,$count);
	push (@originalPos,$count);
	#
	$totalValue = $totalValue + $count;
	#
	$count++;
}
#
my @seqElements = ();
my @seqIndex    = ();
for ( my $i = 0 ; $i < $numbCombination ; $i = $i + 1 ){
	my @shuffled_atoms = fisher_yates_shuffle(\@ClutterArray);
	#
	my $stringElem;
	my $stringIndex;
	for ( my $j = 0 ; $j < scalar (@shuffled_atoms) ; $j = $j + 1 ){
		$stringElem.="$elements[$shuffled_atoms[$j]]\t";
		$stringIndex.="$shuffled_atoms[$j]\t";
	}
	push (@seqElements,$stringElem);
	push (@seqIndex,$stringIndex);
}
#
# Multiple alignment of chemical element sequences
# Score
my @totalScore = ();
# Porcentaje de Coincidencia
my @matchRate = ();
my @signsRate = ();
for ( my $i = 0 ; $i < scalar (@seqElements) ; $i = $i + 1 ){
		my @seq1 = split ("\t",$seqOriginal); # Mi secuencia original
		my @seq2 = split ("\t",$seqElements[$i]); # Secuencias generadas
		#
		my $score = 0;
		for ( my $j = 0 ; $j < scalar (@seq1) ; $j = $j + 1 ){
			if (($seq1[$j] eq $seq2[$j])){
				#print "$seq1[$k] eq $seq2[$k]\n";
				$score++;
			}
		}
		#
		my @pos1 = @originalPos ; # Mi posicion original
		my @pos2 = split ("\t",$seqIndex[$i]); # Posiciones generadas
		#
		my $matrat = 0;
		my $stringSign;
		for ( my $j = 0 ; $j < scalar (@pos1) ; $j = $j + 1 ){
			if (($pos1[$j] eq $pos2[$j])){
				$matrat = $matrat + $pos1[$j];
				$stringSign.="*\t";
			} else {
				$stringSign.="-\t";
			}
		}
		my $porcent = ( ($matrat * 100) / $totalValue );
		my $value = sprintf '%.2f', $porcent;
		#print " $value \n";
		push (@signsRate,$stringSign);
		push (@totalScore,$score);
		push (@matchRate,$value);
}
#
#
my @scoreSort     = ();
my @rateSort      = ();
my @sequencesSort = ();
my @indexSort     = ();
my @signsSort     = ();
my @idx = sort { $matchRate[$b] <=> $matchRate[$a] } 0 .. $#matchRate;
@scoreSort     = @totalScore[@idx];
@rateSort      = @matchRate[@idx];
@sequencesSort = @seqElements[@idx];
@indexSort     = @seqIndex[@idx];
@signsSort     = @signsRate[@idx];
#
# Remove repeated item sequences
my @filtered   = uniq(@sequencesSort);
my @arrayIndex = ();
foreach my $i (@filtered) {
	my $indexSeq = firstidx { $_ eq $i } @sequencesSort;
	push (@arrayIndex,$indexSeq);
}
#
my $poblTotal = scalar (@sequencesSort);
my $sinRepeat = scalar (@filtered);
#
print "\n\n";
print "    * Unique sequences = $sinRepeat\n";
print "    * Total sequences  = $poblTotal\n";
#
# Mi secuencia original
open(FILE, ">MultAlign.txt");
my @seqOrig = split ("\t",$seqOriginal);
#
for ( my $i = 0 ; $i < scalar (@seqOrig) ; $i = $i + 1 ){
	print FILE "$seqOrig[$i]$originalPos[$i] ";
}
print FILE "\n\n";
#
for ( my $i = 0 ; $i < scalar (@arrayIndex) ; $i = $i + 1 ){
	my $indexTmp = $arrayIndex[$i];
	# Signos de coincidencia
	my @signsGen = split ("\t",$signsSort[$indexTmp]);
	# Secuencias generadas
	my @seqGen   = split ("\t",$sequencesSort[$indexTmp]);
	# Indexacion atomos
	my @seqIndex = split ("\t",$indexSort[$indexTmp]);
	# Escribir archivo
	for ( my $j = 0 ; $j < scalar (@seqGen) ; $j = $j + 1 ){
			print FILE "$seqGen[$j]$seqIndex[$j] ";
	}
	print FILE ": Matching rate = $rateSort[$indexTmp]  Score = $scoreSort[$indexTmp]\n";
	for ( my $j = 0 ; $j < scalar (@seqGen) ; $j = $j + 1 ){
			print FILE "  $signsGen[$j] ";
	}
	print FILE "\n";
}
close (FILE);
#
# Estructuras 3-D
print "\n\n";
print "    $numAtoms\n";
print "    Original Structure\n";
for ( my $i = 0 ; $i < scalar (@ClutterArray) ; $i = $i + 1 ){
	print "    $elements[$ClutterArray[$i]]\t$axis[$ClutterArray[$i]]\n";
}
#
#
open(FILE, ">Coords.xyz");
for ( my $i = 0 ; $i < scalar (@arrayIndex) ; $i = $i + 1 ){
	my $indexTmp = $arrayIndex[$i];
	#
	print FILE "$numAtoms\n";
	print FILE "Matching rate = $rateSort[$indexTmp] Score = $scoreSort[$indexTmp]  \n";
	my @dataIndex = split ("\t",$indexSort[$indexTmp]); # Secuencias generadas
	my $count = 0;
	for ( my $j = 0 ; $j < scalar (@dataIndex) ; $j = $j + 1 ){
			print FILE "$elements[$dataIndex[$j]]\t$axis[$count]\n";
			$count++;
	}
}
close (FILE);
#
# Time in console is printed
my $tiempo_final  = new Benchmark;
my $datestringEnd = localtime();
my $tiempo_total  = timediff($tiempo_final, $tiempo_inicial);
my $timeT         = timestr($tiempo_total);
#
print "\n\n";
#
print "    * N° Process : $nprocess\n";
# Start of the program
print "    * Start : $datestringStart\n";
# Final execution of the programme
print "    * End : $datestringEnd\n";
# Execution time
print "    * Execution time :$timeT\n\n";
#
print "\n\n";
