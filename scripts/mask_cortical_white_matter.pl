#! xPERLx -w

# Part of the ASP cortex extraction suite. For a description of the
# algorithms used here see MacDonald et. al. "Automated 3-D Extraction
# of Inner and Outer Surfaces of Cerebral Cortex from MRI" in
# NeuroImage, 12, 340-356, 2000

# mask_cortical_white_matter takes a classified volume as an input,
# and masks out the skull

# Author: David MacDonald
# Last Modified: Sep 2001, Jason Lerch <jason@bic.mni.mcgill.ca>
#    * removed any hardcoded paths
#    * uses MNI::DataDir

require "xINCDIRx/deform_utils.pl";
use MNI::DataDir;

#    &print_out_script_info();

    $volume = shift;
    $output_file = shift;
    $isovalue = shift;
    $model = shift;
    $start_n = shift;
    $end_n = shift;
    $dont_copy = shift;

    if( ! defined($output_file) )
    {
        die "Usage: $0  input.mnc output.obj [isovalue] [model] [start] [end] [dont copy]\n";
    }

    if( !defined( $isovalue ) )
        { $isovalue = 2.5; }

#--------- initialize volumes

    $f = 5;
    $tmp_volume = "/tmp/volume_${$}.mnc";
    register_tmp_files( $tmp_volume );
    system_call( "mask_volume $volume $volume $tmp_volume 0 $isovalue 0");
    system_call( "box_filter_volume_nd $tmp_volume $tmp_volume $f $f $f" );

#---------

    $tmp_model = "/tmp/model_${$}.mnc";
    register_tmp_files( $tmp_model );

    if( !defined( $model ) )
    {
      $model_data_dir = MNI::DataDir::dir('CLASP');
      MNI::DataDir::check_data($model_data_dir, [qw(white_matter_mask.obj)]);
      $model = "${model_data_dir}/white_matter_mask.obj";
    }

    $threshold = 0;

    $fit = "new_fit_3d ";

    $self_dist2 = 0.01;
    $n_selfs = 9;
    $self_factor = 0.0;

    $stop_threshold = 3e-2;
    $stop_iters = 10;

    $n_per = 5;
    $tolerance = 1.0e-10;
    $f_tolerance = 1.0e-10;

    $iters_scale = 1.0;
    $break_scale = 1.0;
    $oo_scale = 1.0;
    $iters_override = 0;

    $stretch_scale = 1;
    $curvature_scale = 1;

    @schedule = (

#size    sw    cw  n_itr  inc   bw si   over   sw   self
#-----  ----  ---- -----  ----  -- --   ----  ----  ----
    5120, 1e4, 1e3,  500, 100,  1,  0.5,  1,   1e0, .5,
    5120, 1e3, 1e2, 1000, 100,  1,  0.5,  1,   1e0, .5,
    5120, 1e2, 1e1, 1000, 100,  1,  0.5,  1,   1e0, .5,
#    5120, 1e1, 1e0, 1000, 100,  1,  0.5,  1,   1e0, .5,
  );


    $sched_size =  10;
 
    $st = $schedule[0];

    if( ! $start_n )
        { $start_n = $st; }

    if( ! $end_n )
        { $end_n = @schedule[@schedule-$sched_size]; }

#------ loop over each schedule

    if( !$dont_copy && $start_n >= $schedule[0] )
    {
        system_call( "subdivide_polygons $model $output_file $st" );
    }

    $prev_n = $st;

    for( $i = 0;  $i < @schedule;  $i += $sched_size )
    {
        #--- get the 14 components of the deformation schedule entry

        ( $size, $sw, $cw, $n_iters, $iter_inc, $bw,
          $si_step, $oversample, $self_weight, $self_dist ) =
                     @schedule[$i..$i+$sched_size-1];

        $sw *= $stretch_scale;
        $cw *= $curvature_scale;
        $oversample = int( $oo_scale * $oversample + 0.5 );

        if( $iters_override > 0 )
            { $n_iters = $iters_override; }
        else
            { $n_iters = int( $n_iters * $iters_scale ); }

        $self_weight *= $self_factor;

        $self = get_self_intersect( $self_weight, $n_selfs, $self_dist,
                                    $self_dist2 );

        if( $size > $start_n || $size < $end_n )
            { $prev_n = $size;  next; }

        if( $size != $prev_n && (!$dont_copy || $size < $start_n) )
        {
            system_call( "subdivide_polygons $output_file $output_file $size" );
        }

        system_call( "subdivide_polygons $model $tmp_model $size" );

        $prev_n = $size;

        #--- if the schedule size is greater than the current number of
        #--- polygons in the deforming surface, subdivide the deforming surface

        print( "Fitting $size-sized polygons, " .
               "max $n_iters iters.\n" );

        $iter_inc *= $break_scale;
        if( $iter_inc <= 0 )  { $iter_inc = $n_iters; }

        $n_failed = 0;

        for( $iter = 0;  $iter < $n_iters;  $iter += $iter_inc )
        {
            system( "echo Step ${size}: $iter / $n_iters    $sw $cw" );

            $ni = $n_iters - $iter;
            if( $ni > $iter_inc )  { $ni = $iter_inc; }

            $b1 = "    -value $bw $tmp_volume 0 $threshold -1 0 0 $oversample ";

            $surf1_info = "-surface $output_file $output_file " .
                 " -stretch $sw $tmp_model 1 0 0 0" .
                 " -curvature $cw $tmp_model 1 0 0 0" .
                 " $b1 " .
                 " $self ";

            $command = "$fit $surf1_info  ".
                       " -print_deriv " .
                       " -step $si_step " .
                       " -fitting $ni $n_per $tolerance " .
                       " -ftol $f_tolerance " .
                       " -stop $stop_threshold $stop_iters ";

            $ret = system_call( "$command", 1 );

            if( $ret == 1 )
            {
                ++$n_failed;
                if( $n_failed == 2 )
                    { last; }
            }
            else
            {
                $n_failed = 0;
            }
        }
    }

    print( "Surface extraction finished.\n" );

    clean_up();
