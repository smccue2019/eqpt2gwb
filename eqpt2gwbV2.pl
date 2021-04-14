#!/usr/bin/perl
# -d:ptkdb

# EQPT2GWB, a script to convert an EQ3/6 thermodynamic database
# to a Geochemists's Workbench (GWB) thermodynamic database.
#
# Scott McCue
# Woods Hole Oceanographic Institution
# Dept of Marine Chemistry and Geochemistry
# January 2007
#
# Usage:
#        prompt> eqpt2gwbVx.pl eqptdbfilename gwbdbfilename
#                         or
#        prompt> eqpt2gwbVx.pl eqptdbfilename > gwbdbfilename
#                         or
#        prompt> cat eqptdbfilename | eqpt2gwbVx.pl > gwbdbfilename
#
# Requires:
#        * a working installation of the PERL scripting language.
#        * the POSIX PERL module, installed with above.
#        * this main script
#        * the accompanying subroutine module, XXXXXXXXXXXX.pl


# For a complete description of the EQ3/6 data0 file format, see
# _EQPT, A Data File Preprocessor for the EQ3/6 Software Package:
# User's Guide and Related Cocumentation_. The version used for
# development of this product is Version 7.0 by Stephanie A.
# Daveler and Thomas J. Woldery, Dec 17 1992.
#
# In particular, refer to chapter 4 of that document, _Data File
# Contents and Structures. pp 27-42. Many comments herein have
# been copied or paraphrased from that chapter.


###################################################################
#################### Initialization Section #######################
###################################################################


use POSIX;
use Tie::IxHash;

$TRUE = 1;
$FALSE = 0;

# The outfile is generally for debugging
#$outfile="blah.txt";
#open OUTF, ">$outfile";
################## Define data0 file #############################
#
#    my $infile="./DATA0RSF.250";
##    my $infile="./data0.250";
#
######### Define I/O filehandles from command line args ##########
#
$NUMARGS = $#ARGV + 1;
if ($NUMARGS == 2) {
    $infile = $ARGV[0];
    $outfile = $ARGV[1];
    open INF, "<$infile" or die "Couldnt open $infile\n";
    open OUTF, ">$outfile" or die "Couldnt open $outfile\n";
} elsif ($NUMARGS ==1 ) {
    $infile = $ARGV[0];
    open INF, "<$infile" or die "Couldnt open $infile\n";
    open OUTF, "|STDOUT" or die "Couldnt open STDOUT\n";
} else {
    open INF, "STDIN|" or die "Couldnt open STDIN\n";
    open OUTF, "|STDOUT" or die "Couldnt open STDOUT\n";
}  

# Set a default output filehandle and corresponding format.
# Important for the 'write' section in the subroutine
# 'write_formatted_el_mwgt';
# Filehandle of the select and format statements must match the
# output filehandle, defined above.
select OUTF; 
format OUTF =
@<<<<<<<<<<<<<< (@<)          mole wt.=   @<<<<<<
$elementname,    $e,                      $mwgt
.
# That dot isnt a speck on your monitor. It ends the format definition.

############ Read the input thermodatabase into a block ##########
    my @fulltext=();
      @fulltext = <INF>;
    close INF;
###################################################################
# Define a correspondence between element symbols and names       #
###################################################################

# 'Sred'- a special specie of Sulphur used by MKT to manipulate
#  GWB and EQPT.
# "Sred","ReducedSulphur",
# 

%sp_abbr2nameh = (
"Sred","ReducedSulphur",
"O","Oxygen",         
"Ag","Silver",
"Al","Aluminum",
"Am","Americium",
"Ar","Argon",
"Au","Gold",
"B","Boron",
"Ba","Barium",
"Be","Beryllium",
"Br","Bromine",
"Ca","Calcium",
"Cd","Cadmium",
"Ce","Cerium",
"Cl","Chlorine",
"Co","Cobalt",
"Cr","Chromium",
"Cs","Cesium",
"Cu","Copper",
"Dy","Dysprosium",
"Er","Erbium",
"Eu","Europium",
"F","Fluorine",
"Fe","Iron",
"Ga","Gallium",
"Gd","Gadolinium",
"H","Hydrogen",
"As","Arsenic",
"C","Carbon",
"P","Phosphorus",
"He","Helium",
"Hf","Hafnium",
"Hg","Mercury",
"Ho","Holmium",
"I","Iodine",
"In","Indium",
"K","Potassium",
"Kr","Krypton",
"La","Lanthanum",
"Li","Lithium",
"Lu","Lutetium",
"Mg","Magnesium",
"Mn","Manganese",
"Mo","Molybdenum",
"N","Nitrogen",
"Na","Sodium",
"Nd","Neodymium",
"Ne","Neon",
"Ni","Nickel",
"Np","Neptunium",
"Pb","Lead",
"Pd","Palladium",
"Pm","Prometium",
"Pr","Praseodymium",
"Pu","Plutonium",
"Ra","Radium",
"Rb","Rubidium",
"Re","Rhenium",
"Rn","Radon",
"Ru","Ruthenium",
"S","Sulfur",
"Sb","Antimony",
"Sc","Scandium",
"Se","Selenium",
"Si","Silicon",
"Sm","Samarium",
"Sn","Tin        ",
"Sr","Strontium",
"Tb","Terbium",
"Tc","Technetium",
"Th","Thorium",
"Ti","Titanium",
"Tl","Thallium",
"Tm","Thulium",
"U","Uranium",
"V","Vanadium",
"W","Tungsten",
"Xe","Xenon",
"Y","Yttrium",
"Yb","Ytterbium",
"Zn","Zinc",
"Zr","Zirconium",
);

###################################################################
############ Input Section- read from source database #############
###################################################################


################### Get data0 parameters ############################

($block0, $bt, $bb) = 
 capture_block("temperature limits", "temperatures",0, \@fulltext);
$d0_temp_limits = nums_from_block($block0);

($block1, $bt, $bb) = 
 capture_block("temperatures", "pressures",0, \@fulltext);
$d0_temps = nums_from_block($block1);

($block2, $bt, $bb) = 
 capture_block("pressures", "debye huckel a",0, \@fulltext);
$d0_pressures = nums_from_block($block2);

($block3, $bt, $bb) = 
 capture_block("debye huckel a", "debye huckel b",0, \@fulltext);
$d0_adh = nums_from_block($block3);

($block4, $bt, $bb) = 
 capture_block("debye huckel b", "bdot",0, \@fulltext);
$d0_bdh = nums_from_block($block4);

($block5, $bt, $bb) = 
 capture_block("bdot", "cco2",0, \@fulltext);
$d0_bdot = nums_from_block($block5);

($block6, $bt, $bb) = 
 capture_block("cco2", "log k for eh reaction",0, \@fulltext);
$d0_cco2 = nums_from_block($block6);

($block7, $bt, $bb) = 
 capture_block("log k for eh reaction", '\+-------------', 0, \@fulltext);
$d0_logk_eh = nums_from_block($block7);

################## Get bdot parameters ############################
# Section is bounded by lines containing begin: "bdot parameters"
# and end: "elements". 
# There are a couple lines of form "+--------------------" that
# should be ignored.

($bdotparblock, $bt, $bb) = 
 capture_block("bdot parameters", "elements", 0, \@fulltext);
$bdothash = parsebdothash($bdotparblock);

########################### Get elements ##############################

# From EQ3/6 User's Guide:
# The chemical elements block consists of a block header followed by the names
# of the chemical elements (represented by standard symbols) and their
# atomic weights (g/mole).

($elementsblock, $bt, $bb) = 
 capture_block("elements", "basis species", 0, \@fulltext);
($elementshash, $elementscount) = parse_elements($elementsblock);

######### Get aqueous species superblock and subblocks ##################

# The aqueous species superblock is comprised of a data block for each
# aqueous species. Structure is complicated somewhat in that these data
# blocks are organized into three "sub-superblocks", the first for strict
# basis species, the second for auxiliary basis species, and the third
# for non-basis aqueous species.

# Read in the superblock and then parse the sub-sbs? Or
# capture and parse the subs at this level? Decision: do the
# subs separately at this level.
#
# Each sub-block ends up saved into a hash table, from which the GWB
# formatted information will be written.
#
# delimiters:
# Subblock1: (strict) "basis species" -> "auxiliary basis species"
# Subblock2: (auxiliary) "auxiliary basis species" -> "aqueous species" 
# Subblock3: non-basis aqueous "aqueous species" -> "solids"
#
# Some blocks are of similar format and can be read using the same
# procedure. Reuse subroutines when possible.
#
#  Block                        Subroutine      GWB Type
#===================================================================
#  "basis species"              parse_basis     "Basis"
#  "auxiliary basis species"    parse_aa        "Redox"
#  "aqueous species"            parse_aa        "Aqueous"
#  "solids"                     parse_solids    "Minerals"
#  "liquids"                    parse_lg        "Minerals"
#  "gases"                      parse_lg        "Gases"

($strictblock, $bt, $bb) =
  capture_block("basis species","auxiliary basis species", 0, \@fulltext);
($stricthash, $strictcount)=parse_basis($strictblock);

($auxblock, $bt, $bb) = 
  capture_block("auxiliary basis species","aqueous species", 0, \@fulltext);
($auxhash, $auxcount) = parse_aa($auxblock);

($nbablock, $bt, $bb) = 
  capture_block("aqueous species","solids", 0, \@fulltext);
($nbahash, $nbacount) = parse_aa($nbablock);

######################## Get solids ##################################

($solidsblock, $bt, $bb) =
   capture_block("solids","liquids", 0, \@fulltext);
($solidshash, $solidscount) = parse_solids($solidsblock);

######################## Get liquids #################################

($liquidsblock, $bt, $bb) =
  capture_block("liquids", "gases", 0, \@fulltext);
($liquidshash, $liquidscount) = parse_lg($liquidsblock);

######################## Get liquids #################################

($gasesblock, $bt, $bb) =
  capture_block("gases", "solid solutions", 0, \@fulltext);
($gaseshash, $gasescount) = parse_lg($gasesblock);

# Completes the rending of the input database into hash tables, arrays,
# etc. It's now time to write out the info that's stored in these
# structures to get a thermodynamic database thats compatible with
# Geochemists WorkBench.
#

###################################################################
####################### New Devel Section #########################
###################################################################

#$testh=test_specie_dup_filter($auxhash);

$auxhash2=filter_specie_from_reaction_pair($auxhash);
$nbahash2=filter_specie_from_reaction_pair($nbahash);
$solidshash2=filter_specie_from_reaction_pair($solidshash);
$liquidshash2=filter_specie_from_reaction_pair($liquidshash);
$gaseshash2=filter_specie_from_reaction_pair($gaseshash);

#pause(1);
###################################################################
###################### Output Section #############################
###################################################################

# Put the writing of header info into a subroutine... Text was
# captured from the example GWB database.
top_of_output();

# Write out the basic params in 4 column format. p4 is a local function.
# OK, so the name isnt very imaginative. I didnt want to type a long name
# a bunch of times.
#
print OUTF "* temperatures\n";
p4(@{$d0_temps}[0 .. 3]);
p4(@{$d0_temps}[4 .. 7]);
print OUTF "* pressures\n";
p4(@{$d0_pressures}[0 .. 3]);
p4(@{$d0_pressures}[4 .. 7]);
print OUTF "* debye huckel a (adh)\n";
p4(@{$d0_adh}[0 .. 3]);
p4(@{$d0_adh}[4 .. 7]);
print OUTF "* debye huckel b (bdh)\n";
p4(@{$d0_bdh}[0 .. 3]);
p4(@{$d0_bdh}[4 .. 7]);
print OUTF "* bdot\n";
p4(@{$d0_bdot}[0 .. 3]);
p4(@{$d0_bdot}[4 .. 7]);
insert_c_co2_n2();

# The following is subsumed by insert_c_co2_n2
#print OUTF "* log k for eh reaction\n";
#p4(@{$d0_logk_eh}[0 .. 3]);
#p4(@{$d0_logk_eh}[4 .. 7]);

# The elements: names, symbols, and molecular weights
print OUTF "\n  $elementscount elements\n";
while (($element, $mwgt) = each %{$elementshash}) {
 undef $elementname;
 $elementname = $sp_abbr2nameh{$element};
 if (defined($elementname)) {
     write_formatted_el_mwgt($elementname, $element, $mwgt);
  } else {
    print STDERR "Missing element name for $element\n";
  }
}

print OUTF "\n-end-\n\n";

############################################################
# Generate the compound superblocks
############################################################

# The basis species type
print OUTF "  $strictcount basis species\n\n";

foreach $st_sp (keys %{$stricthash}) {
    $azero = $bdothash->{$st_sp}->{ionsize};
    print OUTF "\n\n$stricthash->{$st_sp}->{name}\n";
    print OUTF "charge= $stricthash->{$st_sp}->{chargeval}\tion size= $azero A\t mole wt.= $stricthash->{$st_sp}->{moleweight} g\n";
    print OUTF "$stricthash->{$st_sp}->{elcount} elements in species\n";
    while (($e, $ec) = each %{$stricthash->{$st_sp}->{elements}}) {
	printf OUTF "\t%s %s", $ec, $e;
    }
    print OUTF "\n";
    print OUTF "*===== Specie block from EQPT source file =====\n";
    foreach $line (@{$stricthash->{$st_sp}->{comments}}) {
        $chopped = chop($line);
	print OUTF "*$line\n";
    }
}

print OUTF "\n-end-\n\n";

# The redox block, which come from the 'auxiliary' type in EQPT.

print OUTF "  $auxcount redox couples\n\n";

foreach $aux_sp (keys %{$auxhash}) {
    $azero = $bdothash->{$aux_sp}->{ionsize};
    if (!defined($azero)) { $azero = $auxhash->{$aux_sp}->{this_azero}; }
#    $azero = $auxhash->{$aux_sp}->{this_azero};
    print OUTF "\n\n$auxhash->{$aux_sp}->{name}\n";
    print OUTF "charge= $auxhash->{$aux_sp}->{chargeval}\tion size= $azero A\t mole wt.= $auxhash->{$aux_sp}->{moleweight} g\n";
#    print OUTF "$auxhash->{$aux_sp}->{reactionSpecies} species in reaction\n";
    print OUTF "$auxhash->{$aux_sp}->{reactionSpeciesNewCount} species in reaction\n";
#    p3p($auxhash->{$aux_sp}->{reactionPairs});
    p3p($auxhash->{$aux_sp}->{reactionPairsNewList});
    p4(@{$auxhash->{$aux_sp}->{log_grid_vals}}[0 .. 3]);
    p4(@{$auxhash->{$aux_sp}->{log_grid_vals}}[4 .. 7]);
    print OUTF "\n";
    print OUTF "*===== Specie block from EQPT source file =====\n";
    foreach $line (@{$auxhash->{$aux_sp}->{comments}}) {
        $chopped = chop($line);
	print OUTF "*$line\n";
    }
}

print OUTF "\n-end-\n\n";

# Aqueous ('nba' for non-basis aqueous) goes to aqueous

print OUTF "  $nbacount aqueous species\n\n";

foreach $nba_sp (keys %{$nbahash}) {
    $azero = $bdothash->{$nba_sp}->{ionsize};
    if (!defined($azero)) { $azero = $nbahash->{$nba_sp}->{this_azero}; }
#    $azero = $auxhash->{$aux_sp}->{this_azero};
    print OUTF "\n\n$nbahash->{$nba_sp}->{name}\n";
    print OUTF "charge= $nbahash->{$nba_sp}->{chargeval}\tion size= $azero A\t mole wt.= $nbahash->{$nba_sp}->{moleweight} g\n";
#    print OUTF "$nbahash->{$nba_sp}->{reactionSpecies} species in reaction\n";
    print OUTF "$nbahash->{$nba_sp}->{reactionSpeciesNewCount} species in reaction\n";
#    p3p($nbahash->{$nba_sp}->{reactionPairs});
    p3p($nbahash->{$nba_sp}->{reactionPairsNewList});
    p4(@{$nbahash->{$nba_sp}->{log_grid_vals}}[0 .. 3]);
    p4(@{$nbahash->{$nba_sp}->{log_grid_vals}}[4 .. 7]);
    print OUTF "\n";
    print OUTF "*===== Specie block from EQPT source file =====\n";
    foreach $line (@{$nbahash->{$nba_sp}->{comments}}) {
        $chopped = chop($line);
	print OUTF "*$line\n";
    }
}

print OUTF "\n-end-\n\n";

# Solids and liquids from EQPT are minerals in GWB

$mineralscount = $solidscount + $liquidscount;
print OUTF "\n  $mineralscount minerals\n\n";

#Ahlfeldite                      type=
#     formula= NiSeO3:2H2O
#     mole vol.=   63.160 cc      mole wt.=  221.6788 g
#     3 species in reaction
#  1.0000 Ni++                1.0000 SeO3--              2.0000 H2O
#        -4.1213   -4.4894   -5.0144   -5.5657
#        -6.2112   -6.8762  500.0000  500.0000
#*    gflag = 1 [reported delG0f used]
#*    extrapolation algorithm: constant enthalpy approximatio
#*    reference-state data source = 74nau/ryz
#*         delG0f =   -218.800 kcal/mol
#*         delH0f =   -265.070 kcal/mol
#*         S0PrTr =     47.100 cal/(mol*K)




foreach $solids_cn (keys %{$solidshash}) {
    my $solids_cnud = upper_down($solids_cn);
    print OUTF "\n\n$solids_cnud\n";
    print OUTF "formula = $solidshash->{$solids_cn}->{name}\n";
    print OUTF "mole vol.= $solidshash->{$solids_cn}->{prtrval} cc\tmole wt.= $solidshash->{$solids_cn}->{moleweight} g\n";
#    print OUTF "$solidshash->{$solids_cn}->{reactionSpecies} species in reaction\n";
    print OUTF "$solidshash->{$solids_cn}->{reactionSpeciesNewCount} species in reaction\n";
#    p3p($solidshash->{$solids_cn}->{reactionPairs});
    p3p($solidshash->{$solids_cn}->{reactionPairsNewList});
    p4(@{$solidshash->{$solids_cn}->{log_grid_vals}}[0 .. 3]);
    p4(@{$solidshash->{$solids_cn}->{log_grid_vals}}[4 .. 7]);
    print OUTF "\n";
    print OUTF "*===== Specie block from EQPT source file =====\n";
    foreach $line (@{$solidshash->{$solids_cn}->{comments}}) {
        $chopped = chop($line);
	print OUTF "*$line\n";
    }
}

foreach $liquids_sp (keys %{$liquidshash}) {
    my $liquids_cnud = upper_down($solids_cn);
    print OUTF "\n\n$liquids_cnud\n";
    print OUTF "\n\nformula= $liquidshash->{$liquids_sp}->{name}\n";
    print OUTF "mole vol.= $liquidshash->{$liquids_sp}->{prtrval} cc\tmole wt.= $liquidshash->{$liquids_sp}->{moleweight} g\n";
#    print OUTF "$liquidshash->{$liquids_sp}->{reactionSpecies} species in reaction\n";
    print OUTF "$liquidshash->{$liquids_sp}->{reactionSpeciesNewCount} species in reaction\n";
#    p3p($liquidshash->{$liquids_sp}->{reactionPairs});
    p3p($liquidshash->{$liquids_sp}->{reactionPairsNewList});
    p4(@{$liquidshash->{$liquids_sp}->{log_grid_vals}}[0 .. 3]);
    p4(@{$liquidshash->{$liquids_sp}->{log_grid_vals}}[4 .. 7]);
    print OUTF "\n";
    print OUTF "*===== Specie block from EQPT source file =====\n";
    foreach $line (@{$liquidshash->{$liquids_sp}->{comments}}) {
        $chopped = chop($line);
	print OUTF "*$line\n";
    }
}

print OUTF "\n-end-\n\n";

# Gases to gases, dust to dust

foreach $gases_sp (keys %{$gaseshash}) {
    print OUTF "\n\n$gaseshash->{$gases_sp}->{name}\n";
    print OUTF "mole wt.= $gaseshash->{$gases_sp}->{moleweight} g\n";
#    print OUTF "$gaseshash->{$gases_sp}->{reactionSpecies} species in reaction\n";
    print OUTF "$gaseshash->{$gases_sp}->{reactionSpeciesNewCount} species in reaction\n";
#    p3p($gaseshash->{$gases_sp}->{reactionPairs});
    p3p($gaseshash->{$gases_sp}->{reactionPairsNewList});
    p4(@{$gaseshash->{$gases_sp}->{log_grid_vals}}[0 .. 3]);
    p4(@{$gaseshash->{$gases_sp}->{log_grid_vals}}[4 .. 7]);
    print OUTF "\n";
    print OUTF "*===== Specie block from EQPT source file =====\n";
    foreach $line (@{$gaseshash->{$gases_sp}->{comments}}) {
        $chopped = chop($line);
	print OUTF "*$line\n";
    }
}

print OUTF "\n-end-\n\n";

# The eqpt database doesnt provide information for an oxides section,
# which GWB needs. So, an oxides section from a GWB database was copied.
# Write out that oxides info into the database being created by this script.
insert_oxides_from_gwb_db();

####################### Tidy up and finish ############################

close(OUTF);

######################## END OF SCRIPT ################################

########################################################################
############################## Subroutines #############################
########################################################################

###
sub top_of_output {

print OUTF <<"EOB";
dataset of thermodynamic data for gwb programs
dataset format: oct94
activity model: debye-huckel
EOB
}

###
sub p4 {

    $item0=$_[0];
    $item1=$_[1];
    $item2=$_[2];
    $item3=$_[3];

    printf OUTF "%8.4f\t%8.4f\t%8.4f\t%8.4f\n",
    $item0,$item1,$item2,$item3;

}

###
sub p3p {

# Print hash pairs  in six-column (3 pairs) mode.
# Be forgiving if the pair count is less than 3,
# and print the lesser count in the same columns
# as if there were 3.

    my $inhash = $_[0];
    my $itemcount = 0;
   while  (($a1, $a2) = each %{$inhash}) {
       if (($itemcount%3) == 0) {print OUTF "\n";}
       printf OUTF "\t%s %s", $a2, $a1;
       $itemcount++;
   }
   print OUTF "\n";
}

###
sub write_formatted_el_mwgt {

# Use the forms printing functionality of perl to make the
# elements->mole weight table. The format statement is at the
# top of the script, where both out filehandle and format
# are defined with the OUTF label.

# MKT's use of "Reduced Sulphur", a pseudo-element that's Sulphur
# missing an electron, messes up the formatting a bit since it's
# a 4 character label and not one or two like everything else.
# Rather than screw up what works nicely for the usual, create
# an exception for element abbreviations that are longer than
# two characters.
   $elementname = $_[0];
   $e = $_[1];
   $mwgt = $_[2];

   if (length $e > 2) {
       printf OUTF "%s\t(%s)\t\tmole wt.=   %s\n", $elementname,$e,$mwgt;
   } else {
     write;
   }
}
###
sub insert_c_co2_n2 {

# The lines have been captured from a GWB database and will be put
# into the outfile verbatim. To quote MTK:
# "...the c co2 lines are used when doing a specific type of calculation
# that I never do.  So that is why they are there, but it is also why I
# don;t need to care whether the data within these lines are correct."

print OUTF << "EOB";
* c co2 1

         0.1224    0.1127    0.0934    0.0802

         0.0843    0.0989    0.1371    0.1967

* c co2 2

        -0.0047   -0.0105   -0.0036   -0.0015

        -0.0118   -0.0104   -0.0071   -0.0181

* c co2 3

        -0.0004    0.0015    0.0001    0.0005

         0.0031    0.0014   -0.0029   -0.0025

* c co2 4

         0.0000    0.0000    0.0000    0.0000

         0.0000    0.0000    0.0000    0.0000

* c h2o 1

      0.500000E+03   0.145397E+01   0.500000E+03   0.155510E+01

      0.162250E+01   0.500000E+03   0.500000E+03   0.500000E+03

* c h2o 2

      0.500000E+03   0.223570E-01   0.500000E+03   0.364780E-01

      0.458910E-01   0.500000E+03   0.500000E+03   0.500000E+03

* c h2o 3

      0.500000E+03   0.938040E-02   0.500000E+03   0.643660E-02

      0.452210E-02   0.500000E+03   0.500000E+03   0.500000E+03

* c h2o 4

      0.500000E+03  -0.536200E-03   0.500000E+03  -0.713200E-03

     -0.831200E-03   0.500000E+03   0.500000E+03   0.500000E+03

* log k for eh reaction

       -91.0448  -83.1049  -74.0534  -65.8641

       -57.8929  -51.6848  -46.7256  -42.6828

* log k for o2 gas solubility

        -2.6610   -2.8990   -3.0580   -3.1250

        -3.0630   -2.9140   -2.6600   -2.4100

* log k for h2 gas solubility

        -3.0240   -3.1120   -3.1440   -3.1120

        -3.0520   -2.9420   -2.7300   -2.4400

* log k for n2 gas solubility

        -2.9740   -3.1830   -3.3210   -3.3330

        -3.1740   -2.9960   -2.7480   -2.4530
EOB
}

###
sub insert_oxides_from_gwb_db {
# The block contained below was copied from the GWB database file
# "PaulThermo_gwb".


print OUTF <<"EOB";
  42 oxides

Ag2O
     mole wt.=   231.736 g
     3 species in reaction
 -2.0000 H+                  1.0000 H2O                 2.0000 Ag+


As2O5
     mole wt.=   229.840 g
     3 species in reaction
 -3.0000 H2O                 2.0000 H+                  2.0000 H2AsO4-

B2O3
     mole wt.=    69.620 g
     2 species in reaction
 -3.0000 H2O                 2.0000 B(OH)3(aq)

BaO
     mole wt.=   153.326 g
     3 species in reaction
 -2.0000 H+                  1.0000 Ba++                1.0000 H2O

Cs2O
     mole wt.=   281.810 g
     3 species in reaction
 -2.0000 H+                  1.0000 H2O                 2.0000 Cs+

Eu2O3(cubic)
     mole wt.=   351.928 g
     3 species in reaction
 -6.0000 H+                  2.0000 Eu+++               3.0000 H2O

FeO
     mole wt.=    71.846 g
     3 species in reaction
 -2.0000 H+                  1.0000 Fe++                1.0000 H2O

HBr
     mole wt.=    80.912 g
     2 species in reaction
  1.0000 Br-                 1.0000 H+

HF
     mole wt.=    20.006 g
     2 species in reaction
  1.0000 F-                  1.0000 H+

HI
     mole wt.=   127.912 g
     2 species in reaction
  1.0000 H+                  1.0000 I-

K2O
     mole wt.=    94.196 g
     3 species in reaction
 -2.0000 H+                  1.0000 H2O                 2.0000 K+

Li2O
     mole wt.=    29.881 g
     3 species in reaction
 -2.0000 H+                  1.0000 H2O                 2.0000 Li+

Na2O
     mole wt.=    61.979 g
     3 species in reaction
 -2.0000 H+                  1.0000 H2O                 2.0000 Na+

P2O5
     mole wt.=   141.945 g
     3 species in reaction
 -3.0000 H2O                 2.0000 HPO4--              4.0000 H+

SO3
     mole wt.=    80.064 g
     3 species in reaction
 -1.0000 H2O                 1.0000 SO4--               2.0000 H+

SrO
     mole wt.=   103.619 g
     3 species in reaction
 -2.0000 H+                  1.0000 H2O                 1.0000 Sr++

Al2O3
     mole wt.=   101.961 g
     3 species in reaction
 -6.0000 H+                  2.0000 Al+++               3.0000 H2O

CaO
     mole wt.=    56.077 g
     3 species in reaction
 -2.0000 H+                  1.0000 Ca++                1.0000 H2O

Cr2O3
     mole wt.=   151.990 g
     4 species in reaction
 -2.0000 H2O                -1.5000 O2(aq)              2.0000 CrO4--
  4.0000 H+

Cu2O
     mole wt.=   143.091 g
     3 species in reaction
 -2.0000 H+                  1.0000 H2O                 2.0000 Cu+

CuO
     mole wt.=    79.545 g
     3 species in reaction
 -2.0000 H+                  1.0000 Cu++                1.0000 H2O

Fe2O3
     mole wt.=   159.692 g
     3 species in reaction
 -6.0000 H+                  2.0000 Fe+++               3.0000 H2O

HgO
     mole wt.=   216.589 g
     3 species in reaction
 -2.0000 H+                  1.0000 H2O                 1.0000 Hg++

MgO
     mole wt.=    40.304 g
     3 species in reaction
 -2.0000 H+                  1.0000 H2O                 1.0000 Mg++

MnO
     mole wt.=    70.937 g
     3 species in reaction
 -2.0000 H+                  1.0000 H2O                 1.0000 Mn++

NiO
     mole wt.=    74.689 g
     3 species in reaction
 -2.0000 H+                  1.0000 H2O                 1.0000 Ni++

PbO
     mole wt.=   223.199 g
     3 species in reaction
 -2.0000 H+                  1.0000 H2O                 1.0000 Pb++

SeO2
     mole wt.=   110.959 g
     3 species in reaction
 -1.0000 H2O                 1.0000 SeO3--              2.0000 H+

SiO2
     mole wt.=    60.084 g
     1 species in reaction
  1.0000 SiO2(aq)

SnO2
     mole wt.=   150.709 g
     4 species in reaction
 -2.0000 H+                  0.5000 O2(aq)              1.0000 H2O
  1.0000 Sn++

V2O5
     mole wt.=   181.880 g
     3 species in reaction
 -2.0000 H+                  1.0000 H2O                 2.0000 VO2+

ZnO
     mole wt.=    81.389 g
     3 species in reaction
 -2.0000 H+                  1.0000 H2O                 1.0000 Zn++

-end-
EOB
################### Removed from block above ##################
# The block above, copied from a GWB database and used here to complete
# the eqpt -> gwb conversion, contains some elements that are not routinely
# found in eqpt databases. Moved to below so that they aren't printed out.

#NpO2
#     mole wt.=   269.047 g
#     3 species in reaction
# -4.0000 H+                  1.0000 Np++++              2.0000 H2O
#
#PuO2
#     mole wt.=   275.999 g
#     3 species in reaction
# -4.0000 H+                  1.0000 Pu++++              2.0000 H2O
#
#
#RaO
#     mole wt.=   242.024 g
#     3 species in reaction
# -2.0000 H+                  1.0000 H2O                 1.0000 Ra++
#
#ThO2
#     mole wt.=   264.037 g
#     3 species in reaction
# -4.0000 H+                  1.0000 Th++++              2.0000 H2O
#
#TiO2
#     mole wt.=    79.879 g
#     2 species in reaction
# -2.0000 H2O                 1.0000 Ti(OH)4(aq)
#
#Am2O3
#     mole wt.=   533.998 g
#     3 species in reaction
# -6.0000 H+                  2.0000 Am+++               3.0000 H2O
#
#HCl
#     mole wt.=    36.461 g
#     2 species in reaction
#  1.0000 Cl-                 1.0000 H+
#
#RuO2
#     mole wt.=   133.069 g
#     2 species in reaction
# -2.0000 H+                  1.0000 Ru(OH)2++
#
#TcO2
#     mole wt.=   129.999 g
#     4 species in reaction
# -0.7500 O2(aq)             -0.5000 H2O                 1.0000 H+
#  1.0000 TcO4-
#
#U3O8
#     mole wt.=   842.082 g
#     4 species in reaction
#-12.0000 H+                  1.0000 O2(aq)              3.0000 U++++
#  6.0000 H2O


}

####
sub nums_from_block {

    my $block = $_[0];
    my @stack = ();
    my @locblock = ();

    # Make a local copy of the input block
    foreach $inline (@{$block}) {
	push @locblock, $inline;
    }

    foreach $l (@locblock) {

     # There are at most 4 columns of numbers

	($num1, $num2, $num3, $num4) = split(" ", $l);

        if(defined($num1)) { push @stack, $num1;}
	if(defined($num2)) { push @stack, $num2;}
	if(defined($num3)) { push @stack, $num3;}
        if(defined($num4)) { push @stack, $num4;}
    }

    return \@stack;
}

####
sub parsebdothash {

    $bdblock = $_[0];
    %bdh = ();

    foreach $line (@{$bdblock}) {

        ($specie_eqpt, $ionsize, $iontype) = split(" ", $line);

        $specie_gwb = expand_charge_in_name($specie_eqpt);
        $bdh{$specie_gwb}{eqpt_name} = $specie_eqpt;
	$bdh{$specie_gwb}{ionsize} = $ionsize;
        $bdh{$specie_gwb}{iontype} = $iontype;
    }
    return \%bdh;
}

####
sub parse_elements {

    my $elblock = $_[0];
    my %elh = ();
    my @locblock = ();
    my $cnt = 0;

    # Make a local copy of the input block
    foreach $inline (@{$elblock}) {
	push @locblock, $inline;
    }

    foreach $line (@locblock) {
      $count = sprintf("%2.2d",$cnt);
      $line =~ m/^([a-zA-Z0-9_+\-,()]*)\s+(\d+\.\d+)/;
      if (defined($1)) {
        $elh{$1} = $2;
        $cnt++;
      }
    }

    return(\%elh, $cnt);
}   

####
sub capture_block {

    my $regexp1 =$_[0];
    my $regexp2 =$_[1];
    my $startline =$_[2];
    my $sourceblockin = $_[3];
    my @sourceblock = ();
    my $inline;
    

# Make a local copy of the input block
    foreach $inline (@{$sourceblockin}) {
	push @sourceblock, $inline;
    }

    my @blockbuffer = ();

    my $linecnt = $startline;
    my $line = $sourceblock[$linecnt];

    until ($line =~ m/$regexp1/s) {
     # Skip to next line, continue until
     # we get to the top of the target block.
     $linecnt++;
     $line = $sourceblock[$linecnt];
    }

    # If we're to here then we're positioned at the
    # top of the target block.
    $linecnt++;
    $blocktop = $linecnt;     # line num of block top with $sourceblock
    $line = $sourceblock[$linecnt];

    until ( $line =~ m/$regexp2/s) {
        push @blockbuffer, $sourceblock[$linecnt];
        $linecnt++;
	$blockbottom = $linecnt;
        $line = $sourceblock[$linecnt];
    }

    return (\@blockbuffer, $blocktop, $blockbottom);
}

###
sub expand_charge_in_name {

  # EQPT naming uses the format of, e.g., 'Al+3". GWB expects a
  # format of, e.g., 'Al+++'. This routine makes the transformation.

    my $eqpt_str = $_[0];
    my $down_cnt;
    my $alt_str = "";
    # in case the instring isn't changed.
    my $new_str = sprintf("%s", $eqpt_str);   

    $_ = $eqpt_str;

    /\w+([+-])(\d)/ && do {
        my($specie_root, $therest) = split(/[+-]/, $eqpt_str);
	my $charge_symbol = $1;
        my $expand_count = $2;

        for ($down_cnt = $expand_count; $down_cnt > 0; $down_cnt--) {
	    $alt_str = sprintf("%s%s",$alt_str, $charge_symbol);
        }
        $new_str = sprintf("%s%s",$specie_root, $alt_str);
    };

    return($new_str);
}

####
sub parse_lg {

#   'lg' to indicate that this is for liquids and gases.
    my $inblock = $_[0];    # the block to be parsed
    my @lblock = ();

    my $endofblock = $FALSE;
    my $blocktopline = 0;
    my $specieblock;
    my $specieblocksize;
    my $blockbottomline;
    my $blockcnt = 0;

    tie my %comph, 'Tie::IxHash';
    %comph = ();

# Make a local copy of the input block
    foreach $inline (@{$inblock}) {
	push @lblock, $inline;
    }

    until ($endofblock == $TRUE) {

      ($specieblock, $blocktopline, $blockbottomline) =
          capture_block('\+--------------','\+--------------',
          $blocktopline, \@lblock);

        my $specieblocksize= @{$specieblock};
       if ($specieblocksize > 1) {
        # Turn the species block into a hash table
         my ($comname, $sp_eqpt, $bt, $prtrval, $mw, $ec, $eh, $rc, $rh, $lg) = parse_lg_block($specieblock);
        
        # Now add the particulars for this species into the hash
        #for all (sub)species.
	$sp = expand_charge_in_name($sp_eqpt);
	$sp = comma2paren($sp);
        $comph{$sp}{name} = $sp;
        $comph{$sp}{eqpt_name} = $sp_eqpt;
        $comph{$sp}{keytype} = $bt;
        $comph{$sp}{prtrval}=$prtrval;
        $comph{$sp}{moleweight}=$mw;
        $comph{$sp}{elcount}=$ec;
        $comph{$sp}{elements}=$eh;
        $comph{$sp}{reactionSpecies} = $rc;
        $comph{$sp}{reactionPairs} = $rh;
        $comph{$sp}{log_grid_vals} = $lg;
        $comph{$sp}{comments}=$specieblock;


	$blockcnt++;

        if ($blockbottomline >= $#lblock) {$endofblock = $TRUE;}
    }
    }

    return (\%comph, $blockcnt);

}

####
sub parse_lg_block {

#Rn,g
#    date last revised =  23-may-1988
# keys   = gas              refstate         active
#       v0prtr = 24465.000 cm**3/mol  (source = supcrt92                )
#*      mwt    =   222.00000 g/mol
#     1 chemical elements =
#      1.0000 Rn
#     2 species in data0 reaction
#    -1.0000  Rn,g                         1.0000  Rn,aq
#*    LOG K (NEW) GRID (25-50-75-100/175-250-300-350)=Rn,g
#        -2.2780   -2.5020   -2.6100   -2.6390
#        -2.4360   -1.9890   -1.6070   -1.1320

    my $passedblock = $_[0];
    my $line;
    my @comments = ();
    my $elflag = $FALSE;
    my $els_list;
    my @el_els = ();
    my %el_pairs = ();
    my %rs_pairs = ();
    my @rs_els = ();
    my $rs_list;
    my @inblock=();
    my $i;
    my $rsflag=$FALSE;
    my $rslinecount=0;
    my $rslinestodo=0;
    my $log_grid_flag = $FALSE;
    my $lg_list;
    my @lg_nums;
    my @lg_stack = ();
    my $log_grid_lines_read = 0;
    my ($namefield1,$namefield2,$unkcnt);
    my ($specie,$keytype,$prtrval,$moleweight,$elcount,$reaction_species);
    
    undef $namefield1;
    undef $namefield2;
    undef $common_name;
    my($namefield1, $namefield2, $common_name);

# Make a local copy of the input block
    foreach $passedline (@{$passedblock}) {
	push @inblock, $passedline;
    }

    foreach $line (@inblock) {

	$_ = $line;

      LGCASE: {

          /^([a-zA-Z0-9_+\-,()]*)\s+([a-zA-Z0-9_+\-,():]*)\s+$/ && do {
              $namefield1 = $1;
              $namefield2 = $2;
              last LGCASE;
          };

          /keys\s+=\s+(\w+)/ && do {
              $keytype = $1;
              last LGCASE;
	  };
 
          /v0prtr\s+=\s+(-{0,1}\d+.\d+) cm\*\*3\/mol/ && do {
	      $prtrval = $1;
              last LGCASE;
	  };

	  /(\d) chemical elements\s+=/ && do {
              $elcount = $1;
              $elflag = $TRUE;
              last LGCASE;
	  };

          /mwt\s+=\s+(\d+.\d+)\s+g\/mol/ && do {
              $moleweight = $1;
              last LGCASE;
          };

          $elflag && do {
              $els_list = $_;
	      @el_els = split(' ', $els_list);
              for ($i=1; $i<=$#el_els; $i+=2) {$el_pairs{$el_els[$i]}=$el_els[$i-1];}
              $elflag = $FALSE;
              last LGCASE;
          };

          /(\d+)\s+species in data0 reaction/ && do {
              $reaction_species = $1;
              $rsflag = $TRUE;
              $rslinestodo = ceil($reaction_species/2);
              last LGCASE;
          };

          $rsflag && do {
	    $rs_list = $_;
	    # On the flipping of $rsflag $rslinestodo lines should be
            # read.
	    @rs_els = split(' ', $rs_list);
            for ($i=1; $i<=$#rs_els; $i+=2) {
              $rs_pairs{expand_charge_in_name($rs_els[$i])}=$rs_els[$i-1];
            }
            $rslinecount++;
            if ($rslinecount == $rslinestodo) {$rsflag = $FALSE;}
	    last LGCASE;
	  };

          /LOG K \(NEW\) GRID \(25-50-75-100/ && do {
	    $log_grid_flag = $TRUE;
	    last LGCASE;
	  };

	  $log_grid_flag && do {
            # Read in two lines (rows) of four columns each.
              $lg_list = $_;
              @lg_nums = split(' ', $lg_list);
              push @lg_stack, @lg_nums;
	      $log_grid_lines_read++;
              if ($log_grid_lines_read == 2) {$log_grid = $FALSE};
              last LGCASE;
          };

    } #LGCASE
    
    } # foreach
    
    if ($keytype eq 'gas') {
      $specie = $namefield1;
      $common_name = " ";
    }
    else {  # keytype is liquid
      $common_name = $namefield1;
      $specie = $namefield2;
    }

    return ($common_name, $specie, $keytype, $prtrval, $moleweight,
              $elcount, \%el_pairs, $reaction_species, \%rs_pairs, \@lg_stack);

}  
###

####
sub parse_aa {

# Designed to accept the contents of the auxiliary and aqueous
# sub-superblocks within the aqueous species superblock.
#
# In turn calls parse_aa_block, the headers comments of which show
# an example of a species block.

# Example non-basis aqueous block
# Note that the form is the same as the auxiliary specie block.
#Ag(CO3)2-3
#    date last revised =  19-feb-1991
# keys   = aqueous          active
#    charge =   -3.0
#*   azero  =    4.0
#*   mwt    =  227.88660 g/mol
#     3 chemical elements =
#      1.0000 Ag             2.0000 C              6.0000 O
#     4 species in reaction =
#    -1.0000  Ag(CO3)2-3                  -2.0000  h+
#     2.0000  HCO3-                        1.0000  Ag+
#*    LOG K (NEW) GRID (25-50-75-100/175-250-300-350)=Ag(CO3)2-3
#        18.2980   18.2960   18.3230   18.3750
#        18.6620   19.1590   19.6620   20.4970


# Example auxiliary subblock
#NH4+
#    date last revised =  26-jun-1990
# keys   = aux              active
#    charge =    1.0
#*   azero  =    3.0
#*   mwt    =   18.03850 g/mol
#     2 chemical elements =
#      4.0000 H              1.0000 N
#     5 species in reaction =
#    -1.0000  NH4+                        -2.0000  o2(g)
#     1.0000  h2o                          1.0000  NO3-
#     2.0000  h+
#*    LOG K (NEW) GRID (25-50-75-100/175-250-300-350)=NH4+
#        46.9390   42.0510   37.8410   34.1730
#        25.5220   19.1880   15.7510   12.5990



    my $inblock = $_[0];    # the block to be parsed
    my @asblock = ();

    my $endofblock = $FALSE;
    my $blocktopline = 0;
    my $specieblock;
    my $specieblocksize;
    my $blockbottomline;
    my $blockcnt = 0;

    tie my %comph, 'Tie::IxHash';
    %comph = ();

# Make a local copy of the input block
    foreach $inline (@{$inblock}) {
	push @asblock, $inline;
    }

    until ($endofblock == $TRUE) {

#  $blocktopline is modified each time by &capture_block.
#  $blockbottom is returned but isn't needed.

	($specieblock, $blocktopline, $blockbottomline) =
          capture_block('\+--------------','\+--------------',
          $blocktopline, \@asblock);

        my $specieblocksize= @{$specieblock};
        $specieblocksize= @{$specieblock};
        if ($specieblocksize > 1) {
         # Turn the species block into a hash table
         my ($sp_eqpt, $bt, $cv, $mw, $ec, $eh, $rc, $rh, $lg, $az, $sp_alt) = parse_aa_block($specieblock);



	 print STDERR "Undefined $mw\n" unless (defined($sp_eqpt));
 
        # Now add the particulars for this species into the hash
        #for all (sub)species.
	$sp = expand_charge_in_name($sp_eqpt);
	$sp = comma2paren($sp);
        $comph{$sp}{name} = $sp;
        $comph{$sp}{eqpt_name} = $sp_eqpt;
        $comph{$sp}{keytype} = $bt;
        $comph{$sp}{chargeval}=$cv;
        $comph{$sp}{moleweight}=$mw;
        $comph{$sp}{elcount}=$ec;
        $comph{$sp}{elements}=$eh;
        $comph{$sp}{reactionSpecies} = $rc;
        $comph{$sp}{reactionPairs} = $rh;
        $comph{$sp}{log_grid_vals} = $lg;
        $comph{$sp}{comments}=$specieblock;
	$comph{$sp}{this_azero}=$az;

        if (defined($sp_alt)) {        
	    $comph{$sp}{alt_eqpt_name}= expand_charge_in_name($sp_alt);
        }


	$blockcnt++;

        if ($blockbottomline >= $#asblock) {$endofblock = $TRUE;}
    }
    }

    return (\%comph, $blockcnt);

}

####
sub parse_aa_block {

    my $passedblock = $_[0];

    my $line;
    my %comph = ();
    my @comments = ();
    my $elflag = $FALSE;
    my @el_els = ();
    my %el_pairs = ();
    my %elh = ();
    my $el;
    my %rs_pairs = ();
    my @rs_els = ();
    my $rs_list;
    my @inblock=();
    my $i;
    my $rsflag=$FALSE;
    my $rslinecount=0;
    my $rslinestodo=0;
    my $log_grid_flag = $FALSE;
    my $lg_list;
    my @lg_nums;
    my @lg_stack = ();
    my $log_grid_lines_read = 0;
    my ($specie,$keytype,$chargeval,$moleweight,$elcount,$reaction_species);
    my $this_azero;

# Make a local copy of the input block
    foreach $passedline (@{$passedblock}) {
	push @inblock, $passedline;
    }

    foreach $line (@inblock) {

	$_ = $line;

#	print "$_";

      AACASE: {

           # We want the species name, which is alone on its own line
           # except for acetic acid. For now I'll ignore this singular
           # case. 7/2007 Below add a case for acetic acid and acetate
          /^([a-zA-Z0-9_+\-,()]*)\s+$/ && do {
	      $specie = $1;
              last AACASE;
          };

          /^(ACET[-\w]+,AQ)\s+([a-zA-Z0-9_+\-,()]*)\s+$/ && do {
              $specie = $1;
              $speciealt = $2;
              last AACASE;
          };
              

          /keys\s+=\s+(\w+)/ && do {
              $keytype = $1;
              last AACASE;
	  };
 
          /charge\s+=\s+(-{0,1}\d+.\d+)/ && do {
	      $chargeval = $1;
              last AACASE;
	  };

	  /(\d) chemical elements\s+=/ && do {
              $elcount = $1;
              $elflag = $TRUE;
              last AACASE;
	  };

          /mwt\s+=\s+(\d+.\d+)\s+g\/mol/ && do {
              $moleweight = $1;
              last AACASE;
          };

          /azero\s+=\s+(\d+.\d+)/ && do {
              $this_azero = $1;
              last AACASE;
	  };

          $elflag && do {
              $els_list = $_;
#              print "LIST:$els_list";
	      @el_els = split(' ', $els_list);
              for ($i=1; $i<=$#el_els; $i+=2) {$el_pairs{$el_els[$i]}=$el_els[$i-1];}
              $elflag = $FALSE;
              last AACASE;
          };

          /(\d+)\s+species in reaction/ && do {
              $reaction_species = $1;
              $rsflag = $TRUE;
              $rslinestodo = ceil($reaction_species/2);
              last AACASE;
          };

          $rsflag && do {
	    $rs_list = $_;
	    # On the flipping of $rsflag $rslinestodo lines should be
            # read.
	    @rs_els = split(' ', $rs_list);
            for ($i=1; $i<=$#rs_els; $i+=2) {
              $rs_pairs{expand_charge_in_name($rs_els[$i])}=$rs_els[$i-1];
            }
            $rslinecount++;
            if ($rslinecount == $rslinestodo) {$rsflag = $FALSE;}
	    last AACASE;
	  };

          /LOG K \(NEW\) GRID \(25-50-75-100/ && do {
	    $log_grid_flag = $TRUE;
	    last AACASE;
	  };

	  $log_grid_flag && do {
            # Read in two lines (rows) of four columns each.
              $lg_list = $_;
              @lg_nums = split(' ', $lg_list);
              push @lg_stack, @lg_nums;
	      $log_grid_lines_read++;
              if ($log_grid_lines_read == 2) {$log_grid = $FALSE};
              last AACASE;
          };

    } #AACASE
    
    } # foreach
    
 #print "\n=========== $specie ==============\n";
             
    return ($specie, $keytype, $chargeval, $moleweight, $elcount,
             \%el_pairs, $reaction_species, \%rs_pairs, \@lg_stack, $this_azero, $speciealt);

}

####
sub parse_basis {

# Designed to accept the contets of one of the 3 sub-superblocks
# within the aqueous species superblock, NOT the entire superblock.
#

# Example species subblock FROM THE BASIS TYPE. The other types
# (auxiliary basis and non-basis aqueous) have different forms.
#+--------------------------------------------------------------------
#h2o
#    date last revised =  13-jul-1990
# keys   = basis            active
#    charge =    0.0
#*   azero  =    3.0
#*   mwt    =   18.01528 g/mol
#     2 chemical elements =
#      2.0000 H              1.0000 O
#+--------------------------------------------------------------------

# Note that each species is delimited by a line with '+-----------'.
# Yank out each successive specie and process

    my $inblock = $_[0];    # the block to be parsed
    my @asblock = ();

    my $endofblock = $FALSE;
    my $blocktopline = 0;
    my $specieblock;
    my $specieblocksize;
    my $blockbottomline;
    my $blockcnt = 0;

    tie my %comph, 'Tie::IxHash';
    %comph = ();

# Make a local copy of the input block
    foreach $inline (@{$inblock}) {
	push @asblock, $inline;
    }

#    open OUTB, ">outbasis.txt";
#    print OUTB "@asblock";
#    close OUTB;

    until ($endofblock == $TRUE) {

#  $blocktopline is modified each time by &capture_block.
#  $blockbottom is returned but isn't needed.

	($specieblock, $blocktopline, $blockbottomline) =
          capture_block('\+--------------','\+--------------',
          $blocktopline, \@asblock);

#              print "\n$blocktopline:$#asblock\n\n";

        $specieblocksize= @{$specieblock};
       if ($specieblocksize > 1) {

        # Turn the species block into a hash table
        my ($sp_eqpt, $bt, $cv, $mw, $ec, $eh) = parse_basis_block($specieblock);
        
        # Now add the particulars for this species into the hash
        #for all (sub)species.
        $sp = expand_charge_in_name($sp_eqpt);
        $comph{$sp}{name} = $sp;
        $comph{$sp}{eqpt_name} = $sp_eqpt;
        $comph{$sp}{keytype} = $bt;
        $comph{$sp}{chargeval}=$cv;
        $comph{$sp}{moleweight} = $mw;
        $comph{$sp}{elcount}=$ec;
        $comph{$sp}{elements}=$eh;
        $comph{$sp}{comments}=$specieblock;

	$blockcnt++;

        if ($blockbottomline >= $#asblock) {$endofblock = $TRUE;}
    }
    }

    return (\%comph, $blockcnt);

}

####
sub parse_basis_block {

    my $passedblock = $_[0];

# Example species subblock
#
#h2o
#    date last revised =  13-jul-1990
# keys   = basis            active
#    charge =    0.0
#*   azero  =    3.0
#*   mwt    =   18.01528 g/mol
#     2 chemical elements =
#      2.0000 H              1.0000 O
#

    my $line;
    my %comph = ();
    my @comments = ();
    my $elflag = $FALSE;
    my @inblock=();
    my $i;
    my @el_els = ();
    my %el_pairs = ();

# Make a local copy of the input block
    foreach $passedline (@{$passedblock}) {
	push @inblock, $passedline;
    }

    foreach $line (@inblock) {

	$_ = $line;

#	print "$_";

      SPSCASE: {

           # We want the species name, which is alone on its own line
          /^([a-zA-Z0-9_+\-,()]*)\s+$/ && do {
	      $specie = $1;
              last SPSCASE;
          };

          /keys\s+=\s+(\w+)/ && do {
              $keytype = $1;
              last SPSCASE;
	  };
 
          /charge\s+=\s+(-{0,1}\d+.\d+)/ && do {
	      $chargeval = $1;
              last SPSCASE;
	  };

	  /(\d) chemical elements\s+=/ && do {
              $elcount = $1;
              $elflag = $TRUE;
              last SPSCASE;
	  };

          /mwt\s+=\s+(\d+.\d+)\s+g\/mol/ && do {
              $moleweight = $1;
              last SPSCASE;
          };

          $elflag && do {

              $els_list = $_;
	      @el_els = split(" ", $els_list);
              if ($#el_els != ($elcount*2-1)){warn "Element count and pair issue: count=$elcount;pairs=$#el_els\n";}
              for ($i=1; $i<=$#el_els; $i+=2) {$el_pairs{$el_els[$i]}=$el_els[$i-1];}
              $elflag = $FALSE;
              last SPSCASE;
          };
    } #SPSCASE
    
    } # foreach
    
    return ($specie, $keytype, $chargeval, $moleweight, $elcount, \%el_pairs);

} 


####
sub parse_solids {

#   'lsg' to indicate that this is for liquids, solids, and gases.
    my $inblock = $_[0];    # the block to be parsed
    my @lblock = ();

    my $endofblock = $FALSE;
    my $blocktopline = 0;
    my $specieblock;
    my $specieblocksize;
    my $blockbottomline;
    my $blockcnt = 0;

    tie my %comph, 'Tie::IxHash';
    %comph = ();

# Make a local copy of the input block
    foreach $inline (@{$inblock}) {
	push @lblock, $inline;
    }

    until ($endofblock == $TRUE) {

      ($specieblock, $blocktopline, $blockbottomline) =
          capture_block('\+--------------','\+--------------',
          $blocktopline, \@lblock);

        my $specieblocksize= @{$specieblock};
       if ($specieblocksize > 1) {
        # Turn the species block into a hash table
        my ($cn, $sp_eqpt, $bt, $cv, $mw, $ec, $eh, $rc, $rp, $lg) = parse_solids_block($specieblock);
        # Now add the particulars for this species into the hash
        #for all (sub)species.

        # The solids section is inconsistent with naming: usually
        # there's both a common name and a chemical compound string AKA $sp;
        # sometimes there's just a common name, e.g., AMORPHOUS-SILICA
        # and NICKEL. Sometimes a specie will have the same chemical
        # formula and is differentiated by the contents of the common name,
        # as in 'naalsi3o8' denoting all three of ALBITE; ALBITE,HIGH; and
        # ALBITE,LOW. These have the same mole weights, elemental makeup
        # (except for --, HIGH, and LOW), but the Log K values differ. So
        # we gotta keep all three.
        # So...
        #  It looks like the organizing element of the solids hash will
        # have to be the common name. Blech.

	$sp = expand_charge_in_name($sp_eqpt);
        $comph{$cn}{name} = $sp;
        $comph{$cn}{eqpt_name} = $sp_eqpt;
        $comph{$cn}{commonname} = $cn;
        $comph{$cn}{keytype} = $bt;
        $comph{$cn}{prtrval}=$cv;
        $comph{$cn}{moleweight}=$mw;
        $comph{$cn}{elcount}=$ec;
        $comph{$cn}{elements}=$eh;
        $comph{$cn}{reactionSpecies} = $rc;
        $comph{$cn}{reactionPairs} = $rp;
        $comph{$cn}{log_grid_vals} = $lg;
        $comph{$cn}{comments}=$specieblock;

	$blockcnt++;

        if ($blockbottomline >= $#lblock) {$endofblock = $TRUE;}
    }
    }

    return (\%comph, $blockcnt);

}

####
sub parse_solids_block {

# Note that unlike of most the aqueous & aux species the name lines
# contain a common name AND a chemical designation.
#
# Example solids specie block
#ACANTHITE               ag2s
#    date last revised =  22-jun-1999
# keys   = solid            active
#       v0prtr =    34.200 cm**3/mol  (source = supcrt92                )
#*      mwt    =   247.80240 g/mol
#     2 chemical elements =
#      2.0000 Ag             1.0000 Sred
#
#     4 species in data0 reaction
#    -1.0000  ACANTHITE                   -1.0000  h+
#     1.0000  HS-                          2.0000  Ag+
#*    LOG K (NEW) GRID (25-50-75-100/175-250-300-350)=ACANTHITE
#       -35.9710  -32.9200  -30.3250  -28.0940
#       -22.9910  -19.6290  -18.0960  -17.2320

    my $passedblock = $_[0];
    my $line;
    my @comments = ();
    my $elflag = $FALSE;
    my $els_list;
    my @el_els = ();
    my %el_pairs = ();
    my %rs_pairs = ();
    my @rs_els = ();
    my $rs_list;
    my @inblock=();
    my $i;
    my $rsflag=$FALSE;
    my $rslinecount=0;
    my $rslinestodo=0;
    my $log_grid_flag = $FALSE;
    my $lg_list;
    my @lg_nums;
    my @lg_stack = ();
    my $log_grid_lines_read = 0;
    my ($specie,$keytype,$prtrval,$moleweight,$elcount,$reaction_species);

# Make a local copy of the input block
    foreach $passedline (@{$passedblock}) {
	push @inblock, $passedline;
    }

    undef $namefield1;
    undef $namefield2;
    my $namefield1;
    my $namefield2;

    foreach $line (@inblock) {

	$_ = $line;

      SOLCASE: {

          # If there are two names, the first is common and the
          # second is the specie. If one, what it is depends on
          # what keytype (gas or solid) we're processing.
          /^([a-zA-Z0-9_+\-,()]*)\s+([a-zA-Z0-9_+\-,():]*)\s+$/ && do {
              $namefield1 = $1;
              $namefield2 = $2;
              last SOLCASE;
          };

          /keys\s+=\s+(\w+)/ && do {
              $keytype = $1;
              last SOLCASE;
	  };
 
          /v0prtr\s+=\s+(-{0,1}\d+.\d+) cm\*\*3\/mol/ && do {
	      $prtrval = $1;
              last SOLCASE;
	  };

	  /(\d) chemical elements\s+=/ && do {
              $elcount = $1;
              $elflag = $TRUE;
              last SOLCASE;
	  };

          /mwt\s+=\s+(\d+.\d+)\s+g\/mol/ && do {
              $moleweight = $1;
              last SOLCASE;
          };

          $elflag && do {
              $els_list = $_;
#              print "LIST:$els_list";
	      @el_els = split(' ', $els_list);
              for ($i=1; $i<=$#el_els; $i+=2) {$el_pairs{$el_els[$i]}=$el_els[$i-1];}
              $elflag = $FALSE;
              last SOLCASE;
          };

          /(\d+)\s+species in data0 reaction/ && do {
              $reaction_species = $1;
              $rsflag = $TRUE;
              $rslinestodo = ceil($reaction_species/2);
              last SOLCASE;
          };

          $rsflag && do {
	    $rs_list = $_;
	    # On the flipping of $rsflag $rslinestodo lines should be
            # read.
	    @rs_els = split(' ', $rs_list);
            for ($i=1; $i<=$#rs_els; $i+=2) {
              $rs_pairs{expand_charge_in_name($rs_els[$i])}=$rs_els[$i-1];
            }
            $rslinecount++;
            if ($rslinecount == $rslinestodo) {$rsflag = $FALSE;}
	    last SOLCASE;
	  };

          /LOG K \(NEW\) GRID \(25-50-75-100/ && do {
	    $log_grid_flag = $TRUE;
	    last SOLCASE;
	  };

	  $log_grid_flag && do {
            # Read in two lines (rows) of four columns each.
              $lg_list = $_;
              @lg_nums = split(' ', $lg_list);
              push @lg_stack, @lg_nums;
	      $log_grid_lines_read++;
              if ($log_grid_lines_read == 2) {$log_grid = $FALSE};
              last SOLCASE;
          };

    } #SOLCASE
    
    } # foreach
    
    # Inconsistencies in the name line (within solids and between solids
    # and gases) means some post-processing of cases. Gases namelines have
    # one field that should be assigned to the specie name. Solids
    # namelines might have one field or two (usually two). If two, the
    # first is a common name and the second is a specie name. If one, it's
    # a common name. This leads to hash table problems later on, because
    # the organizing item in the hash tables produced by this code is the
    # specie name, not the common name. Short term: use the common name
    # as the organizing item. Longer term: perhaps develop a equating
    # hash table to assign a specie name given the common name. Another
    # possibility is to capture a list of species for which only a common
    # name is given and to send the list back to the user for correction.

    # Another problem is that some species have the same chemical
    # symbols but different parameters and common names, e.g., ALBITE.
    # This also messes up the organization.
    if ($keytype eq 'gas') {
	$specie = $namefield1;
    }
    elsif ($keytype eq 'solid') {
	if (defined($namefield1) and defined($namefield2)) {
	  $common_name = $namefield1;
          $specie = $namefield2;
        }
        else {
          $common_name = $namefield1;
          $specie = "n/a";
        }
    }
     else {        # keytype is liquid
	$common_name = $namefield1;
        $specie = $namefield2;
     }
            
    return ($common_name, $specie, $keytype, $prtrval, $moleweight,
             $elecount, \%el_pairs, $reaction_species, \%rs_pairs, \@lg_stack);

}
###

sub upper_down {

# Converts a word in all uppercase to one thats in lower
# case with a capital first letter.

$inword = $_[0];

$newstr = substr $inword, 1;  # letter in position 2 to the end.
$newstr =~ tr/[A-Z]/[a-z]/;

$newstr2 = sprintf("%s%s", substr($inword, 0,1), $newstr);

return $newstr2;
}

sub comma2paren {

# Converts the EQPT name, which uses an ELEMENT, g or
# ELEMENT, aq convention to denote gaseous or aqueous
# states. GWB wants ELEMENT(g) and ELEMENT(aq).

    my $comm_str = $_[0];
    my $alt_str = "";
    # in case the instring isn't changed.
    my $new_str = sprintf("%s", $comm_str);   

    $_ = $comm_str;

    /(.+),(\w+)/ && do {
	my $specie = $1;
        my $g_or_aq = $2;

        $new_str = sprintf("%s(%s)",$specie, $g_or_aq);
#        print STDERR "$new_str ";
    };

    return($new_str);
}



sub filter_specie_from_reaction_pair {

    my $inh = $_[0];
    my %newlist;
    my @duplist;
    my $newcount;
    my $root1;
    my $root2;
    my $in_sp="";

  foreach $in_sp (keys %{$inh}) {
      %newlist=();
      %duplist = ();
      $newcount=0;
      $root1 = expand_charge_in_name($in_sp);

      print STDERR "Undefined key!\n" unless (defined($in_sp));
      
      ROOT1CASE: {
	$_ = $root1;

#        /[-+]/ && do {
#	    $root1 =~ s/[+-]//g;
#            print STDERR "1: $root1\n";
#            last ROOT1CASE;
#        };

        (/\(aq\)/ or /\(AQ\)/) && do {
            last ROOT1CASE;
        };

        /,aq/ && do {
            $root1 =~ s/,aq/(aq)/;
            last ROOT1CASE;
        };

        /,AQ/ && do {
            $root1 =~ s/,AQ/(AQ)/;
            last ROOT1CASE;
        };

        /,G/ && do {
            $root1 =~ s/,G/(G)/;
            last ROOT1CASE;
        };

        /g/ && do {
            $root1 =~ s/,g/(g)/;
            last ROOT1CASE;
        };

        /,G/ && do {
            $root1 =~ s/,G/(G)/;
            last ROOT1CASE;
        };

        (/\(g\)/ or /\(G\)/) && do {
            last ROOT1CASE;
        };

#        $flag=$TRUE;

        last ROOT1CASE;
         
    }; # ROOT1CASE

    while (($rp_name, $rp_num) = each %{$inh->{$in_sp}->{reactionPairs}}) {
	$root2 = expand_charge_in_name($rp_name);

      ROOT2CASE: {
	$_ = $root2;

#        /[-+]/ && do {
#	    $root2 =~ s/[+-]//g;
#            last ROOT2CASE;
#        };

        (/\(aq\)/ or /\(AQ\)/) && do {
            last ROOT2CASE;
        };

        /,aq/ && do {
            $root2 =~ s/,aq/(aq)/;
            last ROOT2CASE;
        };

        /,AQ/ && do {
            $root2 =~ s/,AQ/(AQ)/;
            last ROOT2CASE;
        };

        /,G/ && do {
            $root2 =~ s/,G/(G)/;
            last ROOT2CASE;
        };

        /g/ && do {
            $root2 =~ s/,g/(g)/;
            last ROOT2CASE;
        };

        /,G/ && do {
            $root2 =~ s/,G/(G)/;
            last ROOT2CASE;
        };

        (/\(g\)/ or /\(G\)/) && do {
            last ROOT2CASE;
        };

#        if ($flag==$TRUE){print STDERR "$in_sp $rp_name\n";}
        last ROOT2CASE;
         
    }; # ROOT2CASE

      if ($root1 eq $root2) {
#	  print STDERR "MATCH:$root1 $root2\n";
          $duplist{$rp_name}= $inh->{$in_sp}->{reactionPairs}->{$rp_name};
      } else {
#          print STDERR "$rp_name\n";
	  $newcount++;
          $newlist{$rp_name}= $inh->{$in_sp}->{reactionPairs}->{$rp_name};
	}

    }  # while

      $newref = \%newlist;
      $dupref = \%duplist;
      $inh->{$in_sp}->{reactionPairsNewList}= {%$newref};
      $inh->{$in_sp}->{reactionSpeciesNewCount}=$newcount;
      $inh->{$in_sp}->{reactionPairsDupList}={%$dupref};
    } # foreach

    return($inh);
}
    
