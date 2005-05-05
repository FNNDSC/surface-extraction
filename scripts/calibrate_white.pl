#!xPERLx -w

# Calibrate the white surface with a gradient field.
#
# Authors: June-sic Kim <luck3d@bic.mni.mcgill.ca> and
#
# May 2004

use strict;
use MNI::Startup;
use Getopt::Tabular;
use MNI::Spawn;
use MNI::DataDir;
use MNI::FileUtilities qw(check_output_dirs);


# ===== Global Variables =====
my ($usage, $help);
my ($mri, $cls, $white, $gradient, $gradient_masked);
my ($white_fix, $gradient_parm, $cls_tmp, $cls_gradient);
my ($cls_gradient_tmp, $cls_gm_gradient, $cls_wm_gradient);
my ($cls_gm, $cls_wm, $csf_skel, $cls_fix, $mask);
my ($mult_const, $max_num, $mean);
my ($filledImage, $mask_dil);
my ($dx, $dy, $dz, $kerndx, $kerndy, $kerndz);
my ($std, $gradient_weight, $gradient_threshold);

# ===== Argument Processing =====

$usage = "$ProgramName mri.mnc cls.mnc csf_skel.mnc white.obj output.obj\n";
$help = "Help still to be written";

my @leftOverArgs;
my @argTbl = 
    (
     @DefaultArgs,
     ,
     );
GetOptions(\@argTbl, \@ARGV, \@leftOverArgs) or die "\n";

$mri = shift @leftOverArgs or die $usage;
$cls = shift @leftOverArgs or die $usage;
$csf_skel = shift @leftOverArgs or die $usage;
$white = shift @leftOverArgs or die $usage;
$white_fix = shift @leftOverArgs or die $usage;

RegisterPrograms(["minccalc", "mincmath", "mincmorph",
                  "surface_mask2", "mincblur", "intensity_statistics",
                  "new_fit_3d", "mv", "rm", "mincresample", "dilate_volume"]);

if ($Clobber) {
    AddDefaultArgs("minccalc", ["-clobber"]);
}

# create necessary tmp directory
check_output_dirs($TmpDir);

$kerndx="${TmpDir}/dx.kern";
$kerndy="${TmpDir}/dy.kern";
$kerndz="${TmpDir}/dz.kern";
&CreateKernelFiles;
$dx="${TmpDir}/dx.mnc";
$dy="${TmpDir}/dy.mnc";
$dz="${TmpDir}/dz.mnc";

# compute gradient volume
sub gradient_volume{
  Spawn(["mincmorph","-clobber","-convolve", "-kernel", $kerndx, $_[0], $dx]);
  Spawn(["mincmorph","-clobber","-convolve", "-kernel", $kerndy, $_[0], $dy]);
  Spawn(["mincmorph","-clobber","-convolve", "-kernel", $kerndz, $_[0], $dz]);
  Spawn(["minccalc", "-clobber", "-expression", 'out=sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2])/2;', $dx, $dy, $dz, $_[1]]);
}

# Prepare cls file
$cls_fix = "${TmpDir}/cls.mnc";
Spawn(["mincresample", "-like", $mri, $cls, $cls_fix]);
$cls = "${TmpDir}/cls.mnc";

# Make a gradient mask
$mask = "${TmpDir}/mask.mnc";
$mask_dil = "${TmpDir}/mask_dil.mnc";
$filledImage = "${TmpDir}/filledImage.mnc";
Spawn(["minccalc", "-expression", 'out=1;', $cls, $filledImage]);
Spawn(["surface_mask2", "-binary", $filledImage, $white, $mask]);
Spawn(["dilate_volume", $mask, $mask_dil, "1", "26", "1"]);
Spawn(["mincresample", "-like", $cls, $mask_dil, "${TmpDir}/mask_tmp.mnc"]);
Spawn(["mv", "-f", "${TmpDir}/mask_tmp.mnc", $mask_dil]);

# Make gradient images of cls
$cls_fix = "${TmpDir}/cls_fix.mnc";
$cls_gm = "${TmpDir}/cls_gm.mnc";
$cls_gm_gradient = "${TmpDir}/cls_gm_gradient.mnc";
$cls_gradient = "${TmpDir}/cls_gradient.mnc";
$cls_gradient_tmp = "${TmpDir}/cls_gradient_tmp.mnc";
Spawn(["mincresample", "-like", $cls, $csf_skel, "${TmpDir}/csf_skel.mnc"]);
$csf_skel = "${TmpDir}/csf_skel.mnc";
Spawn(["minccalc", "-expression", 'if(A[1]>0){out=1;}else{out=A[0];}', $cls, $csf_skel, $cls_fix]);
Spawn(["minccalc", "-expression", 'if(A[0]>2){out=2;}else if(A[0]<1){out=1;}else{out=A[0];}', $cls_fix, $cls_gm]);
#Spawn(["make_gradient_volume", $cls_gm, $cls_gm_gradient, "1", "3"]);
#Spawn(["gradient_volume.sh", $cls_gm, $cls_gm_gradient]);
gradient_volume($cls_gm, $cls_gm_gradient);
Spawn(["minccalc", "-clobber", "-expression", 'if(A[0]>0.01){out=1;}else{out=0;}', $cls_gm_gradient, $cls_gradient_tmp]);
Spawn(["mincblur", "-clobber", "-fwhm", "1", $cls_gradient_tmp, "${TmpDir}/cls_gm_1"]);

# Make a gradient image of an MRI
$gradient = "${TmpDir}/gradient.mnc";
#Spawn(["make_gradient_volume", $mri, $gradient, "1", "3"]);
#Spawn(["gradient_volume.sh", $mri, $gradient]);
gradient_volume( $mri, $gradient );

# Composite a gradient-parameter image
$gradient_parm = "${TmpDir}/gradient_parm.mnc";
$gradient_masked = "${TmpDir}/gradient_masked.mnc";
Spawn(["minccalc", "-expression", 'if(A[0]>0.01 || A[2]==0){out=0;}else{out=A[1];}', "${TmpDir}/cls_gm_1_blur.mnc", $gradient, $mask_dil, $gradient_parm]);
#Spawn(["mincblur", "-clobber", "-fwhm", "1", $gradient_parm, "${TmpDir}/grad_1"]);
#Spawn(["mv", "-f", "${TmpDir}/grad_1_blur.mnc", $gradient_parm]);
Spawn(["surface_mask2", $gradient_parm, $white, $gradient_masked]);

# Mask the cls image with a gradient mask
Spawn(["minccalc", "-expression", 'if(A[1]>0){out=A[0];}else{out=0;}', $cls, $mask_dil, "${TmpDir}/cls_tmp.mnc"]);
Spawn(["mv", "-f", "${TmpDir}/cls_tmp.mnc", $cls]);

# Determine parameters for 'new_fit_3d'
#chomp( $max_num = qx(intensity_statistics ${gradient_masked} none | grep Max | cut -c 12-) );
chomp( $mean = qx(intensity_statistics ${gradient_masked} none | grep Mean | cut -c 12-) );
chomp( $std = qx(intensity_statistics ${gradient_masked} none | grep Std | cut -c 12-) );
$gradient_weight = 1e-2 / ($mean * $mean);
$gradient_threshold = $mean + $std;

Spawn(["new_fit_3d", "-surface", $white, $white_fix, "-stretch", "2", $white,
       "-.9", "0", "0", "0", "-boundary", "100", "1", $cls, "2.5", "-", "3", 
       "50", "0", "0", "1", "-self_intersect", "1", "0.75", "-self_intersect",
       "10", "0.667777777777778", "-self_intersect", "100", 
       "0.585555555555556", "-self_intersect", "1000", "0.503333333333333", 
       "-self_intersect", "10000", "0.421111111111111", "-self_intersect", 
       "100000", "0.338888888888889", "-self_intersect", "1000000", 
       "0.256666666666667", "-self_intersect", "10000000", "0.174444444444445",
       "-self_intersect", "100000000", "0.0922222222222222", "-self_intersect",
       "1e8", "0.01", "-print_deriv", "-step", "1", "-fitting", "300", "5", 
       "1e-10", "-ftol", "1e-10", "-stop", "0.03", "10", 
       "-gradient", $gradient_weight, $gradient_parm, "0", 
       $gradient_threshold, "0", $gradient_threshold, "1", "0", "1"],
       err_action => 'ignore');

#Create the kernel files used in taking the image derivative
sub CreateKernelFiles
{
 
  open DXFILE, "> $kerndx";
  print DXFILE "MNI Morphology Kernel File\n";
  print DXFILE "Kernel_Type = Normal_Kernel;\n";
  print DXFILE "Kernel =\n";
  print DXFILE "  -1.0  0.0  0.0  0.0  0.0     -0.5\n";
  print DXFILE "   1.0  0.0  0.0  0.0  0.0      0.5;\n";
  close DXFILE;


  open DYFILE, "> $kerndy";
  print DYFILE "MNI Morphology Kernel File\n";
  print DYFILE "Kernel_Type = Normal_Kernel;\n";
  print DYFILE "Kernel =\n";
  print DYFILE "   0.0 -1.0  0.0  0.0  0.0     -0.5\n";
  print DYFILE "   0.0  1.0  0.0  0.0  0.0      0.5;\n";
  close DYFILE;

  open DZFILE, "> $kerndz";
  print DZFILE "MNI Morphology Kernel File\n";
  print DZFILE "Kernel_Type = Normal_Kernel;\n";
  print DZFILE "Kernel =\n";
  print DZFILE "   0.0  0.0 -1.0  0.0  0.0     -0.5\n";
  print DZFILE "   0.0  0.0  1.0  0.0  0.0      0.5;\n";
  close DZFILE;

}
