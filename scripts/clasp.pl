#!xPERLx -w
#
# Based on ~david/Surface_deformation/How_to_extract_cortical_surfaces.txt

use strict;

use Getopt::Tabular;
use MNI::Startup;
use MNI::FileUtilities;
use MNI::Spawn;

# Several of the sub-programs will look for static files using
# MNI::DataDir.  Set the directory here if the user does not have
# their own.
#
$ENV{MNI_DATAPATH} = 'xDATADIRx'
  unless defined $ENV{MNI_DATAPATH};


MNI::Spawn::RegisterPrograms
  ( [qw/  rm
     mask_cortical_white_matter
     surface_mask2 mask_volume
     extract_white_surface
     classify_correct
     make_asp_grid
     expand_from_white/ ] )
  or exit 1;


# --- set the help & usage strings ---
my $help = <<HELP;
The help string goes here.

HELP

my $usage = <<USAGE;
usage: $ProgramName [options] classified.mnc final.mnc
       $ProgramName -help to list options
USAGE

Getopt::Tabular::SetHelp( $help, $usage );

# --- process options ---
my( $masking_surface, $masked_input, $output_prefix );
my( $laplace, $cls_correct, $skelCSF );
my @options = 
  ( @DefaultArgs,     # from MNI::Startup
    ['-out', 'string', 1, \$output_prefix,
     'prefix for all output files'],
  );

GetOptions( \@options, \@ARGV ) 
  or exit 1;
die "$usage\n" unless @ARGV == 2 or @ARGV == 3;

# $input_volume should be a CLASSIFIED MINC volume
my $input_volume = shift;
my $t1_volume = shift;
my $objMask = shift;


# Masking surface and cortical-matter masked version of input.
# These are temporary files.
MNI::FileUtilities::check_output_path("${TmpDir}")
  or exit 1;
$masking_surface = "${TmpDir}input-mask.obj"
  unless defined $masking_surface;
$masked_input = "${TmpDir}input-masked"
  unless defined $masked_input;

# Output filename prefix
if ( !defined $output_prefix ) {
    my( $dir, $base, $ext ) = MNI::PathUtilities::split_path($input_volume);
    $output_prefix = $base;
}


#---------------------------------------------------------------------------
#  Step 1:   2 hours    Correcting a cls volume to classify putamen as WM
#                       Partial volume estimation
#---------------------------------------------------------------------------

# Partial volume estimation using pve2
$cls_correct = "${TmpDir}cls_correct.mnc";
$skelCSF = "${TmpDir}skel_CSF.mnc";
if( !defined $objMask ) {
  Spawn(["classify_correct", "-unmasked", $input_volume, $t1_volume, $cls_correct, $skelCSF]);
}
else {
  my $masked_input = "${TmpDir}masked_input.mnc";
  Spawn(["surface_mask2", $input_volume, $objMask, $masked_input]);
  Spawn(["classify_correct", $masked_input, $t1_volume, $cls_correct, $skelCSF]);
}
$input_volume = $cls_correct;


#---------------------------------------------------------------------------
#  Step 2:   1/2 hour   Creating a mask surface to chop of non-cortical white
#                       matter
#---------------------------------------------------------------------------

# Deforms from model file `white_matter_mask.obj' using new_fit_3d
Spawn( "mask_cortical_white_matter $input_volume  $masking_surface  2.5");


#---------------------------------------------------------------------------
#  Step 3:   2 minutes  Masking the classified volume with the masking_surface
#---------------------------------------------------------------------------

# conglomerate routine
Spawn("surface_mask2  $input_volume  $masking_surface  $masked_input");
#Spawn("rm -f $masking_surface");

# Set all voxels with value in range (1,2.5) to zero
# conglomerate routine
Spawn("mask_volume    $masked_input $masked_input $masked_input 1 2.5 0");


#---------------------------------------------------------------------------
#  Step 4:   20 hours   Shrink wrapping a sphere to the MASKED white matter
#                       matter,  (note that the second argument below is a
#                                 prefix)
#                       creating white_surface_{320,1280,5120,20480,81920}.obj
#                       you can delete the model files that are created:
#                            white_surface_m_{320,1280,5120,20480,81920}.obj
#---------------------------------------------------------------------------

# Deforms from model file `white_model_320.obj' using new_fit_3d
Spawn("extract_white_surface  $masked_input   ${output_prefix}_white  2.5");
#Spawn("rm -f $masked_input");


#---------------------------------------------------------------------------
#  Step 5:   1 hour     Create a Laplacian field from the WM surface to the
#                       outer boundary of gray matter
#---------------------------------------------------------------------------

$laplace = "${TmpDir}laplace.mnc";
Spawn(["make_asp_grid", $skelCSF, "${output_prefix}_white_81920.obj", $input_volume, $laplace]);


#---------------------------------------------------------------------------
#  Step 6:   20 hours   Expand a copy of the white surface out to the gray
#                       boundary.  Note that now we can delete the temporary
#                       masked white matter ($masked_input), and use the
#                       original classified volume for this step.
#---------------------------------------------------------------------------

# Model file is ellipsoid_${n_polygons}.obj.gz
# (n_polygons will be 81920)
#Spawn("expand_from_white  $input_volume ${output_prefix}_white_81920.obj"
#      . " ${output_prefix}_gray_81920.obj  1.5");
Spawn(["expand_from_white", $input_volume, "${output_prefix}_white_81920.obj", "${output_prefix}_gray_81920.obj", "1.5", $laplace]);