// Downloaded on 2022-04-30 from 
// http://www.bic.mni.mcgill.ca/users/claude/inflate_to_sphere_implicit.c

/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/
/*
   inflate_to_sphere_implicit.c

   Inflate a .obj as a sphere. 
   - iterative smoothing of highest Gaussian curvature points first
   - alternating growth on convex/concave part of surface
   - implicit geometric smoothing on regions with small triangles
   - overall smoothing with equalizing triangle areas

   inflate_to_sphere_implicit in.obj out.obj [n_inflate [n_smooth]]

   Values: in.obj = input object file
           out.obj = output object file
           n_inflate = number of inflation iterations (default 2000)
           n_smooth = number of smoothing iterations (default 2000)

   Author: Claude Lepage, July 2019.

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <volume_io.h>
#include <bicpl.h>

Real Pi = 3.14159265358979323846;

int dbg = 0;

// Prototypes of functions in this file.

static void usage( char * );
static Status read_surface_obj( STRING, int *, Point *[],
                                Vector *[], int *, int *[], int *[], int **[] );
static Status get_surface_neighbours( polygons_struct *, int *[],
                                      int ** [] );
static void save_surface_obj( STRING, int, Point *, Vector *, int, int []);

static void compute_surface_normals( int, int, int *, Point *, Vector * );
static void compute_triangle_normal( Point, Point, Point, Real[3] );
static Real gaussian_curvature( int, int, int *, Point * );
void area_weights( int, int, int *, int *, Point *, Real * );
void implicit( int, int, Point *, Vector *, int *, int *, int **, Real );
Real smooth_iso( int n_points, int n_elems, int n_iters, Real relax,
                 Real aspect_ratio, int * connec,
                 Point * coords, int * n_ngh, int * ngh[] );
Real smooth_area( int n_points, int n_iters, Real relax, Real * ref_areas, 
                  Point * coords, int * n_ngh, int * ngh[] );

// Main program.

int main( int argc, char * argv[] ) {

  int      i, j, k, jj;

  int      n_points;           // number of grid points per object
  int      n_elems;            // number of triangles per object
  Point  * coords;             // coordinates
  Vector * normals;            // normal vectors
  int    * connec;             // connectivity
  int    * n_ngh = NULL;       // node neighbours (inverse connectivity)
  int   ** ngh = NULL;

  if( argc < 3 ) {
    usage( argv[0] );
    return( 1 );
  }

  // Read the surface file.
  if( read_surface_obj( argv[1], &n_points, &coords, &normals,
                        &n_elems, &connec, &n_ngh, &ngh ) != OK ) {
    return 1;
  }

  int n_iter = 2000;
  int n_smooth = 2000;
  if( argc >= 4 ) {
    n_iter = atoi( argv[3] );
    if( argc >= 5 ) {
      n_smooth = atoi( argv[4] );
    }
  }

  int * mask = (int *)malloc( n_points * sizeof( int ) );
  for( i = 0; i < n_points; i++ ) {
    mask[i] = 1;
  }
  Real * weight = (Real *)malloc( n_points * sizeof( Real ) );

  Real * mc = (Real *)malloc( n_points * sizeof( Real ) );
  for( i = 0; i < n_points; i++ ) {
    mc[i] = gaussian_curvature( i, n_ngh[i], ngh[i], coords );
  }

  // Compute reference vertex-based areas on the input surface.
  // The areas on the sphere will be matched to those areas so
  // that the sphere is area-preserving locally.

  Real * ref_areas = (Real *)malloc( n_points * sizeof( Real ) );

  Real total_area = 0.0;
  for( i = 0; i < n_points; i++ ) {
    ref_areas[i] = 0.0;
    for( k = 0; k < n_ngh[i]; k++ ) {
      int k1 = ngh[i][k];
      int k2 = ngh[i][(k+1)%n_ngh[i]];
      Real len01 = sqrt( ( coords[i].coords[0] - coords[k1].coords[0] ) *
                         ( coords[i].coords[0] - coords[k1].coords[0] ) +
                         ( coords[i].coords[1] - coords[k1].coords[1] ) *
                         ( coords[i].coords[1] - coords[k1].coords[1] ) +
                         ( coords[i].coords[2] - coords[k1].coords[2] ) *
                         ( coords[i].coords[2] - coords[k1].coords[2] ) );
      Real len02 = sqrt( ( coords[i].coords[0] - coords[k2].coords[0] ) *
                         ( coords[i].coords[0] - coords[k2].coords[0] ) +
                         ( coords[i].coords[1] - coords[k2].coords[1] ) *
                         ( coords[i].coords[1] - coords[k2].coords[1] ) +
                         ( coords[i].coords[2] - coords[k2].coords[2] ) *
                         ( coords[i].coords[2] - coords[k2].coords[2] ) );
      Real len12 = sqrt( ( coords[k1].coords[0] - coords[k2].coords[0] ) *
                         ( coords[k1].coords[0] - coords[k2].coords[0] ) +
                         ( coords[k1].coords[1] - coords[k2].coords[1] ) *
                         ( coords[k1].coords[1] - coords[k2].coords[1] ) +
                         ( coords[k1].coords[2] - coords[k2].coords[2] ) *
                         ( coords[k1].coords[2] - coords[k2].coords[2] ) );

      Real s = 0.5 * ( len01 + len02 + len12 );
      ref_areas[i] += sqrt( fabs( s * ( s - len01 ) * ( s - len02 ) *
                                  ( s - len12 ) ) + 1.0e-16 );
    }
    total_area += ref_areas[i];
  }

  // Scale areas of like surface to that of a unit sphere (4*pi).
  // Each triangle on the sphere surface covers the same percentage
  // of total area as on the brain surface. We are trying to minimize
  //    ( area[i] - areas_like[i] )**2
  // to have matching vertex-based normalized areas.

  Real factor = ( 4.0 * Pi ) / ( total_area / 3.0 );
  for( i = 0; i < n_points; i++ ) {
    ref_areas[i] *= factor;
  }

  Real relax = -0.75;  // this is SOR=1.75. Use relax=0.25 for no SOR.

  for( k = 1; k <= n_iter; k++ ) {

    // Re-center the surface and update surface normals (not every iteration).
    if( k % 10 == 0 ) {
      for( j = 0; j < 3; j++ ) {
        Real  sum = 0.0;
        for( i = 0; i < n_points; i++ ) {
          sum += coords[i].coords[j];
        }
        sum /= (Real)n_points;
        for( i = 0; i < n_points; i++ ) {
          coords[i].coords[j] -= sum;
        }
      }
      // this is to rescale in x, y, z so that the blob is sphere like

      Real axes[3] = { 0.0, 0.0, 0.0 };
      for( i = 0; i < n_points; i++ ) {
        for( j = 0; j < 3; j++ ) {
          axes[j] += fabs( coords[i].coords[j] );
        }
      }
      for( j = 0; j < 3; j++ ) {
        axes[j] /= (Real)n_points;
      }
      for( i = 0; i < n_points; i++ ) {
        for( j = 0; j < 3; j++ ) {
          coords[i].coords[j] /= axes[j];
        }
      }

      compute_surface_normals( n_elems, n_points, connec, coords, normals );
    }

    for( i = 0; i < n_points; i++ ) {
      int update_mc = 0;

      if( ( k % 2 == 0 && mc[i] < 0.0 ) || ( k % 2 == 1 && mc[i] > 0.0 ) ) {
        update_mc = 1;
      } else {
        // always smooth at a local max in abs(mc).
        int local_max = 1;
        for( j = 0; j < n_ngh[i]; j++ ) { 
          int vj = ngh[i][j];
          if( fabs( mc[vj] ) >= fabs( mc[i] ) ) local_max = 0;
        }
        update_mc = local_max;
      }
      if( update_mc ) {
        // geometric smoothing
        Real sum[3] = { 0.0, 0.0, 0.0 };
        for( j = 0; j < n_ngh[i]; j++ ) { 
          int vj = ngh[i][j];
          sum[0] += coords[vj].coords[0];
          sum[1] += coords[vj].coords[1];
          sum[2] += coords[vj].coords[2];
        }
        for( jj = 0; jj < 3; jj++ ) { 
          coords[i].coords[jj] += ( 1.0 - relax ) * 
                                  ( sum[jj] / (Real)n_ngh[i] - coords[i].coords[jj] );
        }
      }
    }

    // find area threshold based on standard deviation. We could use
    // the median, but this isn't so critical.
    if( k % 100 == 0 ) {
      area_weights( n_points, n_elems, connec, n_ngh, coords, weight );

      Real mean = 0.0;
      Real stdev = 0.0;
      for( i = 0; i < n_points; i++ ) {
        mean += weight[i];
        stdev += weight[i] * weight[i];
      }
      mean /= (Real)n_points;
      stdev = sqrt( stdev / (Real)n_points - mean * mean );
      Real area_thresh = mean - 2.0 * stdev;

      // printf( "Area: Mean = %g Stdev = %g Threshold = %g\n", 
      //         mean, stdev, area_thresh );

      int count = 0;
      for( i = 0; i < n_points; i++ ) {
        if( weight[i] < area_thresh ) {
          mask[i] = 1;
          count++;
        } else {
          mask[i] = 0;
        }
      }
      if( count > 0 ) {
        printf( "Implicit smoothing on %d vertices with small triangles\n", count );
        implicit( n_points, 50, coords, normals, mask, n_ngh, ngh, -0.75 );
      }
    }

    // smoothing to equalize area of triangles
    smooth_iso( n_points, n_elems, 10, 0.25, 0.02, connec, coords, n_ngh, ngh );

    for( i = 0; i < n_points; i++ ) {
      mc[i] = gaussian_curvature( i, n_ngh[i], ngh[i], coords );
    }

    int num_inverted = -1;
    if( k % 10 == 0 ) {

      Point * new_coords = (Point *)malloc( n_points * sizeof( Point ) );
      if( !new_coords ) {
        printf( "Error allocating memory for new_coords.\n" );
        exit( 1 );
      }

      for( i = 0; i < n_points; i++ ) {
        Real mag = 1.0 / sqrt( coords[i].coords[0] * coords[i].coords[0] +
                               coords[i].coords[1] * coords[i].coords[1] +
                               coords[i].coords[2] * coords[i].coords[2] );
        new_coords[i].coords[0] = coords[i].coords[0] / mag;
        new_coords[i].coords[1] = coords[i].coords[1] / mag;
        new_coords[i].coords[2] = coords[i].coords[2] / mag;
      } 

      num_inverted = 0;
      compute_surface_normals( n_elems, n_points, connec, new_coords, normals );
      for( i = 0; i < n_points; i++ ) {
        Real dot = new_coords[i].coords[0] * normals[i].coords[0] +
                   new_coords[i].coords[1] * normals[i].coords[1] +
                   new_coords[i].coords[2] * normals[i].coords[2];
        if( dot < 0.0 ) num_inverted++;
      }
      free( new_coords );

      printf( "Iter = %d  inverted = %d\n", k, num_inverted );

      compute_surface_normals( n_elems, n_points, connec, coords, normals );
    }

    if( num_inverted == 0 ) {
      break;
    }
  }
  free( mc );
  free( mask );
  free( weight );

  // Normalize the coordinates onto the unit sphere centered at zero
  // and perform smoothing to obtain equal area triangles.

  relax = -0.10;
  for( k = 1; k <= n_smooth; k++ ) {
    Real res = smooth_area( n_points, 10, relax, ref_areas,
                            coords, n_ngh, ngh );
    if( k % 10 == 0 ) printf( "Res(%d) = %g\n", k, res );

    for( i = 0; i < n_points; i++ ) {
      Real mag = 1.0 / sqrt( coords[i].coords[0] * coords[i].coords[0] +
                             coords[i].coords[1] * coords[i].coords[1] +
                             coords[i].coords[2] * coords[i].coords[2] );
      coords[i].coords[0] *= mag;
      coords[i].coords[1] *= mag;
      coords[i].coords[2] *= mag;
    } 
  } 

  // if n_iter=0 and n_smooth=0, this will simply normalize to a sphere.
  if( n_smooth == 0 ) {
    for( i = 0; i < n_points; i++ ) {
      Real mag = 1.0 / sqrt( coords[i].coords[0] * coords[i].coords[0] +
                             coords[i].coords[1] * coords[i].coords[1] +
                             coords[i].coords[2] * coords[i].coords[2] );
      coords[i].coords[0] *= mag;
      coords[i].coords[1] *= mag;
      coords[i].coords[2] *= mag;
    } 
  }

  free( ref_areas );

  compute_surface_normals( n_elems, n_points, connec, coords, normals );

  // test for inverted mesh

  int inverted = 0;

  for( i = 0; i < n_points; i++ ) {
    Real dot = coords[i].coords[0] * normals[i].coords[0] +
               coords[i].coords[1] * normals[i].coords[1] +
               coords[i].coords[2] * normals[i].coords[2];
    if( dot < 0.0 ) inverted++;
  }
  printf( "%d vertices with inverted triangles\n", inverted );
  
  save_surface_obj( argv[2], n_points, coords, normals, n_elems, connec );

  FREE( connec );
  FREE( coords );
  FREE( normals );
  FREE( n_ngh );
  FREE( ngh );

  return 0;
}


// Do smoothing on the coordinates to achieve equal area triangles.
// Taken from equalize_sphere.c.

Real smooth_iso( int n_points, int n_elems, int n_iters, Real relax,
                 Real aspect_ratio, int * connec,
                 Point * coords, int * n_ngh, int * ngh[] ) {

  int i, j, k, kk;
  Real res;
  Real * new_coords = (Real *)malloc( 3 * n_points * sizeof( Real ) );
  if( !new_coords ) {
    printf( "Error allocating memory for new_coords.\n" );
    exit( 1 );
  }
  Real * weight = (Real *)malloc( n_points * sizeof( Real ) );
  if( !weight ) {
    printf( "Error allocating memory for weight.\n" );
    exit( 1 );
  }
  Real * sumweight = (Real *)malloc( n_points * sizeof( Real ) );
  if( !sumweight ) {
    printf( "Error allocating memory for sumweight.\n" );
    exit( 1 );
  }

  for( kk = 1; kk <= n_iters; kk++ ) {
    for( i = 0; i < 3*n_points; i++ ) {
      new_coords[i] = 0.0;
    }

    area_weights( n_points, n_elems, connec, n_ngh, coords, weight );

    Real min_edge_threshold;
    Real total_area = 0.0;
    for( i = 0; i < n_points; i++ ) {
      total_area += weight[i];
    }
    min_edge_threshold = 2.0 * aspect_ratio *
                         sqrt( total_area / ( sqrt( 3.0 ) * n_points ) );

    for( i = 0; i < n_points; i++ ) {
      sumweight[i] = 0.0;
    }

    for( i = 0; i < n_elems; i++ ) {

      for( k = 0; k < 3; k++ ) {       // k is 3 verts of triangle
        int k0 = connec[3*i+k];
        int k1 = connec[3*i+(k+1)%3];
        int k2 = connec[3*i+(k+2)%3];

        // Stability condition:
        // Above term abs((w1-w0)/(w1+w0))<=1 for w1>=0, w0>=0.
        // The weight below acts as a fraction of movement along
        // the edge (x1,x0). So for a value of 10, we can a maximum
        // motion of (x1-x0)/10 or 10% of the edge length. For
        // stability, we must ensure that the maximum displacement
        // does not exceed the surrounding edge lengths. A value
        // of 2.0 is quite aggressive and gives fast stable convergence.

        if( weight[k1] - weight[k0] > 0.0 ) {
          Real factor = ( weight[k1] - weight[k0] ) /
                        ( weight[k1] + weight[k0] );
          for( j = 0; j < 3; j++ ) {     // j is x,y,z
            // this is a displacement vector from vertex k0
            new_coords[3*k0+j] += factor * ( coords[k1].coords[j] -
                                             coords[k0].coords[j] );
          }
          sumweight[k0] += 1.0;
        }
        if( weight[k2] - weight[k0] > 0.0 ) {
          Real factor = ( weight[k2] - weight[k0] ) /
                        ( weight[k2] + weight[k0] );
          for( j = 0; j < 3; j++ ) {     // j is x,y,z
            new_coords[3*k0+j] += factor * ( coords[k2].coords[j] -
                                             coords[k0].coords[j] );
          }
          sumweight[k0] += 1.0;
        }
      }
    }
    res = 0.0;
    for( i = 0; i < n_points; i++ ) {

      // Find maximum allowable distance of propagation of vertex k0
      // inside its connected triangles.

      Real mindist = 1.0e20;
      for( k = 0; k < n_ngh[i]; k++ ) {
        int k0 = i;
        int k1 = ngh[i][k];
        int k2 = ngh[i][(k+1)%n_ngh[i]];

        Real t = ( ( coords[k1].coords[0] - coords[k0].coords[0] ) *
                   ( coords[k1].coords[0] - coords[k2].coords[0] ) +
                   ( coords[k1].coords[1] - coords[k0].coords[1] ) *
                   ( coords[k1].coords[1] - coords[k2].coords[1] ) +
                   ( coords[k1].coords[2] - coords[k0].coords[2] ) *
                   ( coords[k1].coords[2] - coords[k2].coords[2] ) ) /
                 ( ( coords[k1].coords[0] - coords[k2].coords[0] ) *
                   ( coords[k1].coords[0] - coords[k2].coords[0] ) +
                   ( coords[k1].coords[1] - coords[k2].coords[1] ) *
                   ( coords[k1].coords[1] - coords[k2].coords[1] ) +
                   ( coords[k1].coords[2] - coords[k2].coords[2] ) *
                   ( coords[k1].coords[2] - coords[k2].coords[2] ) );
        if( t < 0.0 ) {
          t = 0.0;
        } else if( t > 1.0 ) {
          t = 1.0;
        }
        Real xt = ( 1.0 - t ) * coords[k1].coords[0] + t * coords[k2].coords[0];
        Real yt = ( 1.0 - t ) * coords[k1].coords[1] + t * coords[k2].coords[1];
        Real zt = ( 1.0 - t ) * coords[k1].coords[2] + t * coords[k2].coords[2];
        Real dist = ( xt - coords[k0].coords[0] ) * ( xt - coords[k0].coords[0] ) +
                    ( yt - coords[k0].coords[1] ) * ( yt - coords[k0].coords[1] ) +
                    ( zt - coords[k0].coords[2] ) * ( zt - coords[k0].coords[2] );
        if( dist < mindist ) {
          mindist = dist;
        }
      }
      mindist = 0.25 * sqrt( mindist );

      if( mindist > min_edge_threshold && sumweight[i] > 0.0 ) {
        Real mag = sqrt( new_coords[3*i+0] * new_coords[3*i+0] +
                         new_coords[3*i+1] * new_coords[3*i+1] +
                         new_coords[3*i+2] * new_coords[3*i+2] ) / sumweight[i];

        if( mag > mindist ) {
          sumweight[i] *= mag / mindist;
        }

        for( j = 0; j < 3; j++ ) {
          new_coords[3*i+j] /= sumweight[i];
          res += fabs( new_coords[3*i+j] );
          coords[i].coords[j] += ( 1.0 - relax ) * new_coords[3*i+j];
        }
      } else {
        // do regular isotropic smoothing on this vertex

        sumweight[i] = 0.0;
        for( j = 0; j < 3; j++ ) {
          new_coords[3*i+j] = 0.0;
        }
        for( k = 0; k < n_ngh[i]; k++ ) {
          for( j = 0; j < 3; j++ ) {
            new_coords[3*i+j] += weight[ngh[i][k]] *
                                 coords[ngh[i][k]].coords[j];
          }
          sumweight[i] += weight[ngh[i][k]];
        }
        for( j = 0; j < 3; j++ ) {
          coords[i].coords[j] = relax * coords[i].coords[j] +
                                ( 1.0 - relax ) * new_coords[3*i+j] / sumweight[i];
        }
      }
    }
    res /= (Real)(3.0 * n_points);
  }

  free( new_coords );
  free( weight );
  free( sumweight );
  return( res );
}

// Do smoothing on the coordinates to achieve matching distribution
// of vertex-based areas.

Real smooth_area( int n_points, int n_iters, Real relax, Real * ref_areas, 
                  Point * coords, int * n_ngh, int * ngh[] ) {

  int i, j, k, kk;
  Real res;

  Real * areas = (Real *)malloc( n_points * sizeof( Real ) );
  if( !areas ) {
    printf( "Error allocating memory for areas.\n" );
    exit( 1 );
  }

  for( kk = 1; kk <= n_iters; kk++ ) {

    res = 0.0;

    for( i = 0; i < n_points; i++ ) {

      areas[i] = 0.0;
      for( k = 0; k < n_ngh[i]; k++ ) {
        int k1 = ngh[i][k];
        int k2 = ngh[i][(k+1)%n_ngh[i]];

        Real len01 = sqrt( ( coords[i].coords[0] - coords[k1].coords[0] ) *
                           ( coords[i].coords[0] - coords[k1].coords[0] ) +
                           ( coords[i].coords[1] - coords[k1].coords[1] ) *
                           ( coords[i].coords[1] - coords[k1].coords[1] ) +
                           ( coords[i].coords[2] - coords[k1].coords[2] ) *
                           ( coords[i].coords[2] - coords[k1].coords[2] ) );
        Real len02 = sqrt( ( coords[i].coords[0] - coords[k2].coords[0] ) *
                           ( coords[i].coords[0] - coords[k2].coords[0] ) +
                           ( coords[i].coords[1] - coords[k2].coords[1] ) *
                           ( coords[i].coords[1] - coords[k2].coords[1] ) +
                           ( coords[i].coords[2] - coords[k2].coords[2] ) *
                           ( coords[i].coords[2] - coords[k2].coords[2] ) );
        Real len12 = sqrt( ( coords[k1].coords[0] - coords[k2].coords[0] ) *
                           ( coords[k1].coords[0] - coords[k2].coords[0] ) +
                           ( coords[k1].coords[1] - coords[k2].coords[1] ) *
                           ( coords[k1].coords[1] - coords[k2].coords[1] ) +
                           ( coords[k1].coords[2] - coords[k2].coords[2] ) *
                           ( coords[k1].coords[2] - coords[k2].coords[2] ) );

        Real s = 0.5 * ( len01 + len02 + len12 );
        areas[i] += sqrt( fabs( s * ( s - len01 ) * ( s - len02 ) *
                           ( s - len12 ) ) + 1.0e-16 );
      }
      res += ( ref_areas[i] - areas[i] ) * ( ref_areas[i] - areas[i] );
    }
    res = sqrt( res ) / (Real)n_points;

    // compute displacements
    Real new_coords[3];
    for( i = 0; i < n_points; i++ ) {

      int ii = ( kk % 2 == 1 ) ? i : n_points - i - 1;

      for( j = 0; j < 3; j++ ) {
        new_coords[j] = 0.0;
      }

      Real sum_weight = 0.0;
      for( k = 0; k < n_ngh[ii]; k++ ) {
        int k1 = ngh[ii][k];
        Real ww = areas[k1] / ref_areas[k1];
        sum_weight += ww;
        for( j = 0; j < 3; j++ ) {
          new_coords[j] += ww * ( coords[k1].coords[j] - coords[ii].coords[j] );
        }
      }
      for( j = 0; j < 3; j++ ) {
        coords[ii].coords[j] += ( 1.0 - relax ) * new_coords[j] / sum_weight;
      }
    }
  }

  free( areas );
  return( res );
}


// Do smoothing on the coordinates (simple averaging). No area weights.
// With relax < 0, we get the accelerated SOR scheme. 

void implicit( int n_points, int n_iters, Point * coords, 
               Vector * normals, 
               int * mask, int * n_ngh, int ** ngh, Real relax ) {

  int i, j, k, kk;

  // Real SOR = 1.75;  // higher values might be faster, but riskier.
  Real SOR = 1.0 - relax;  // want SOR = 1.75 for linear coefficients

  for( kk = 1; kk <= n_iters; kk++ ) {
    Real res = 0.0;
    for( i = 0; i < n_points; i++ ) {
      if( mask[i] > 0 && n_ngh[i] > 0 ) {
        Real old_coords[3];

        for( j = 0; j < 3; j++ ) {     // j is x,y,z
          old_coords[j] = coords[i].coords[j];
          coords[i].coords[j] = 0.0;
          for( k = 0; k < n_ngh[i]; k++ ) {
            coords[i].coords[j] += coords[ngh[i][k]].coords[j];
          }
          coords[i].coords[j] /= (Real)n_ngh[i];
          coords[i].coords[j] = old_coords[j] + SOR * ( coords[i].coords[j] -
                                                        old_coords[j] );
          res += fabs( coords[i].coords[j] - old_coords[j] );
        }
      }
    }
    res /= (Real)n_points;
    if( kk % 100 == 0 ) printf( "Iter %d  res %g\n", kk, res );
    // if( res < 1.0e-05 ) break;  // should be relative residual
  }
}

// Compute area weights for mesh smoothing. It returns the average
// area of triangles around each vertex.

void area_weights( int n_points, int n_elems, int * connec, int * n_ngh,
                   Point * coords, Real * weight ) {

  int i, j, k, kk;
  Real len[3];

  for( i = 0; i < n_points; i++ ) {
    weight[i] = 0.0;
  }
  for( i = 0; i < n_elems; i++ ) {
    for( j = 0; j < 3; j++ ) {
      int k0 = connec[3*i+j];
      int k1 = connec[3*i+(j+1)%3];
      len[j] = sqrt( ( coords[k0].coords[0] - coords[k1].coords[0] ) *
                     ( coords[k0].coords[0] - coords[k1].coords[0] ) +
                     ( coords[k0].coords[1] - coords[k1].coords[1] ) *
                     ( coords[k0].coords[1] - coords[k1].coords[1] ) +
                     ( coords[k0].coords[2] - coords[k1].coords[2] ) *
                     ( coords[k0].coords[2] - coords[k1].coords[2] ) );
    }
    Real s = 0.5 * ( len[0] + len[1] + len[2] );
    Real area = sqrt( fabs( s * ( s - len[0] ) * ( s - len[1] ) * ( s - len[2] ) ) + 1.0e-10 );
    for( j = 0; j < 3; j++ ) {
      weight[connec[3*i+j]] += area;
    }
  }

  // Average area of triangle around a vertex.
  for( i = 0; i < n_points; i++ ) {
    weight[i] /= (Real)n_ngh[i];
  }

}


// Recompute the surface normals at the nodes. Simply average the normal
// vector of the neighbouring faces at each node.
//
static void compute_surface_normals( int n_elems, int n_points, int * connec, 
                                     Point * coords, Vector * normals ) {

    int  i, j, v0, v1, v2;
    Real norm[3], mag;

    for( i = 0; i < n_points; i++ ) {
      normals[i].coords[0] = 0.0;
      normals[i].coords[1] = 0.0;
      normals[i].coords[2] = 0.0;
    }

    for( i = 0; i < n_elems; i++ ) {
      v0 = connec[3*i];
      v1 = connec[3*i+1];
      v2 = connec[3*i+2];
      compute_triangle_normal( coords[v0], coords[v1], coords[v2], norm );
      for( j = 0; j < 3; j++ ) {
        normals[v0].coords[j] += norm[j];
        normals[v1].coords[j] += norm[j];
        normals[v2].coords[j] += norm[j];
      }
    }

    for( i = 0; i < n_points; i++ ) {
      mag = sqrt( normals[i].coords[0] * normals[i].coords[0] +
                  normals[i].coords[1] * normals[i].coords[1] +
                  normals[i].coords[2] * normals[i].coords[2] );
      normals[i].coords[0] /= mag;
      normals[i].coords[1] /= mag;
      normals[i].coords[2] /= mag;
    }
}

// Compute the normal vector to a triangular face.

static void compute_triangle_normal( Point v1, Point v2, Point v3, Real norm[3] ) {

  Real a1x, a1y, a1z, a2x, a2y, a2z, mag;

  a1x = v2.coords[0] - v1.coords[0];
  a1y = v2.coords[1] - v1.coords[1];
  a1z = v2.coords[2] - v1.coords[2];

  a2x = v3.coords[0] - v1.coords[0];
  a2y = v3.coords[1] - v1.coords[1];
  a2z = v3.coords[2] - v1.coords[2];

  norm[0] = a1y * a2z - a1z * a2y;
  norm[1] = a1z * a2x - a1x * a2z;
  norm[2] = a1x * a2y - a1y * a2x;
  mag = sqrt( norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2] );
  norm[0] /= mag;
  norm[1] /= mag;
  norm[2] /= mag;
}


// -------------------------------------------------------------------
// Compute Gaussian curvature at a vertex.

static Real gaussian_curvature( int v0, int n_ngh, int * ngh, Point * coords ) {

  int i;

  Real area = 0.0, angle = 2.0 * Pi;
  Real prod[3];

  for( i = 0; i < n_ngh; i++ ) {
    int v1 = ngh[i];
    int v2 = ngh[(i+1)%n_ngh];

    // area of triangle

    Real a1x = coords[v1].coords[0] - coords[v0].coords[0];
    Real a1y = coords[v1].coords[1] - coords[v0].coords[1];
    Real a1z = coords[v1].coords[2] - coords[v0].coords[2];
    Real a2x = coords[v2].coords[0] - coords[v0].coords[0];
    Real a2y = coords[v2].coords[1] - coords[v0].coords[1];
    Real a2z = coords[v2].coords[2] - coords[v0].coords[2];

    prod[0] = a1y * a2z - a1z * a2y;
    prod[1] = a1z * a2x - a1x * a2z;
    prod[2] = a1x * a2y - a1y * a2x;
    area += sqrt( prod[0] * prod[0] + prod[1] * prod[1] + prod[2] * prod[2] );

    // angle with middle vertex

    Real mag = sqrt( ( a1x * a1x + a1y * a1y + a1z * a1z ) *
                     ( a2x * a2x + a2y * a2y + a2z * a2z ) );
    angle -= acos( ( a1x * a2x + a1y * a2y + a1z * a2z ) / mag ); 
  }

  // K = 3 * (2Pi - sum(angle)) / sum(area)
  // drop the constants 3 and 0.5 in area

  return( angle / area );
}


// -------------------------------------------------------------------
// Help message on how to use this module.
//
static void usage( char * executable_name ) {

  STRING  usage_format = "\
Description: Inflate a .obj as a sphere. \n\
   - iterative smoothing of highest Gaussian curvature points first\n\
   - alternating growth on convex/concave part of surface\n\
   - implicit geometric smoothing on regions with small triangles\n\
   - overall smoothing with equalizing triangle areas\n\
Usage: %s in.obj out.obj [n_inflate [n_smooth]]\n\
Values: in.obj = input object file\n\
        out.obj = output object file\n\
        n_inflate = number of inflation iterations (default 2000)\n\
        n_smooth = number of smoothing iterations (default 2000)\n\n\
Copyright Alan C. Evans\n\
Professor of Neurology\n\
McGill University\n\n";

  print_error( usage_format, executable_name );
}


// -------------------------------------------------------------------
// Load the cortical surface.
//
// filename: name of the .obj file
// n_points: the number of the vertices
// points: (x,y,z) coordinates
// normals: normal vectors
// n_elem: number of triangles
// connec: connectivity of triangles
// n_neighbours: number of vertices around each node
// neighbours: the set of ordered triangle consisting of the vertices
//
static Status read_surface_obj( STRING filename,
                                 int * n_points,
                                 Point * points[],
                                 Vector * normals[],
                                 int * n_elem,
                                 int * connec[],
                                 int * n_neighbours[],
                                 int ** neighbours[] ) {

  int               i, n_objects;
  object_struct  ** object_list;
  polygons_struct * surface;
  File_formats      format;
  STRING            expanded;

  expanded = expand_filename( filename );   // why?????

  int err = input_graphics_file( expanded, &format, &n_objects,
                                 &object_list );

  if( err != OK ) {
    print_error( "Error reading file %s\n", expanded );
    return( ERROR );
  }

  if( n_objects != 1 || 
      ( n_objects == 1 && get_object_type(object_list[0]) != POLYGONS ) ) {
    print_error( "Error in contents of file %s\n", expanded );
    return( ERROR );
  }

  delete_string( expanded );

  surface = get_polygons_ptr( object_list[0] );

  int ntri = 0, nquad = 0, unknown = 0;
  int start_ind = 0;
  for( i = 0; i < surface->n_items; i++ ) {
    int nn = surface->end_indices[i] - start_ind;
    start_ind = surface->end_indices[i];
    if( nn == 3 ) {
      ntri++;
    } else {
     if( nn == 4 ) {
       nquad++;
     } else {
       unknown++;
       printf( "face with %d nodes\n", nn );
     }
   }
  }
  printf( "%d triangles, %d quads, %d unknown faces in mesh\n", ntri, nquad, unknown );

  // Check if all polygons are triangles.

  if( 3 * surface->n_items != surface->end_indices[surface->n_items-1] ) {
    printf( "Error: Surface must contain only triangular polygons.\n" );
    delete_object_list( n_objects, object_list );
    return ERROR;
  }

  // Make a copy of the coordinates, the normals, and the
  // connectivity since delete_object_list will destroy them.

  *n_points = surface->n_points;
  *n_elem = surface->n_items;
  ALLOC( *points, surface->n_points );
  if( !(*points) ) {
    printf( "Error allocating memory for points.\n" );
    exit( 1 );
  }
  ALLOC( *normals, surface->n_points );
  if( !(*normals) ) {
    printf( "Error allocating memory for normals.\n" );
    exit( 1 );
  }

  for( i = 0; i < *n_points; i++ ) {
    (*points)[i].coords[0] = surface->points[i].coords[0];
    (*points)[i].coords[1] = surface->points[i].coords[1];
    (*points)[i].coords[2] = surface->points[i].coords[2];
    (*normals)[i].coords[0] = surface->normals[i].coords[0];
    (*normals)[i].coords[1] = surface->normals[i].coords[1];
    (*normals)[i].coords[2] = surface->normals[i].coords[2];
  }

  if( connec ) {
    ALLOC( *connec, 3*surface->n_items );
    if( !(*connec) ) {
      printf( "Error allocating memory for connec.\n" );
      exit( 1 );
    }
    for( i = 0; i < 3*surface->n_items; i++ ) {
      (*connec)[i] = surface->indices[i];
    }
  }

  if( n_neighbours && neighbours ) {
    get_surface_neighbours( surface, n_neighbours, neighbours );
  }

  delete_object_list( n_objects, object_list );

  return( OK );
}


// -------------------------------------------------------------------
// Construct the edges around each node. The edges are sorted to
// make an ordered closed loop.
//
private Status get_surface_neighbours( polygons_struct * surface,
                                       int * n_neighbours_return[],
                                       int ** neighbours_return[] ) {

  int    i, j, k, jj;
  int  * tri;
  int  * n_ngh;
  int ** ngh;
  int  * ngh_array;

  // Check if all polygons are triangles.

  if( 3 * surface->n_items != surface->end_indices[surface->n_items-1] ) {
    printf( "Surface must contain only triangular polygons.\n" );
    return ERROR;
  }

  // Check if the node numbering starts at 0 or 1.

  int min_idx, max_idx;

  min_idx = 100*surface->n_points;  // anything big
  max_idx = 0;                      // anything small

  for( i = 0; i < 3*surface->n_items; i++ ) {
    if( surface->indices[i] < min_idx ) min_idx = surface->indices[i];
    if( surface->indices[i] > max_idx ) max_idx = surface->indices[i];
  }

  // Shift numbering to start at zero, for array indexing. Note
  // that we don't care if surface->indices array is modified.

  if( min_idx != 0 ) {
    for( i = 0; i < 3*surface->n_items; i++ ) {
      surface->indices[i] -= min_idx;
    }
  }

  // Count number of triangles attached to each node.

  ALLOC( n_ngh, surface->n_points );
  if( !n_ngh ) {
    printf( "Error allocating memory for n_ngh.\n" );
    exit( 1 );
  }
  ALLOC( ngh, surface->n_points );
  if( !ngh ) {
    printf( "Error allocating memory for ngh.\n" );
    exit( 1 );
  }
  ALLOC( ngh_array, 3*surface->n_items );
  if( !ngh_array ) {
    printf( "Error allocating memory for ngh_array.\n" );
    exit( 1 );
  }

  for( i = 0; i < surface->n_points; i++ ) {
    n_ngh[i] = 0;
  }

  for( i = 0; i < 3*surface->n_items; i++ ) {
    n_ngh[surface->indices[i]]++;
    ngh_array[i] = -1;
  }

  int max_ngh = 0;
  int sum_ngh = 0;
  for( i = 0; i < surface->n_points; i++ ) {
    ngh[i] = &(ngh_array[sum_ngh]);
    sum_ngh += n_ngh[i];
    max_ngh = MAX( max_ngh, n_ngh[i] );
  }

  // At first, store the indices of the triangles in the neighbours.
  for( i = 0; i < surface->n_items; i++ ) {
    for( j = 0; j < 3; j++ ) {
      jj = surface->indices[3*i+j];
      for( k = 0; k < n_ngh[jj]; k++ ) {
        if( ngh[jj][k] == -1 ) {
          ngh[jj][k] = i;
          break;
        }
      }
    }
  }

  // Now create a sort closed loop of the node neighbours.
  // This is needed by the parametric=0 FEM algorithm.
  //
  //         1 ----- 2
  //          /\   /\
  //         /  \ /  \
  //       0 ----P---- 3
  //         \  / \  /
  //          \/   \/
  //         5 ----- 4
  //

  int * tmp;
  ALLOC( tmp, 2*max_ngh );
  if( !tmp ) {
    printf( "Error allocating memory for tmp.\n" );
    exit( 1 );
  }

  for( i = 0; i < surface->n_points; i++ ) {
    for( k = 0; k < n_ngh[i]; k++ ) {
      tri = &(surface->indices[3*ngh[i][k]]);
      for( j = 0; j < 3; j++ ) {
        if( tri[j] == i ) break;
      }
      tmp[2*k+0] = tri[(j+1)%3];
      tmp[2*k+1] = tri[(j+2)%3];
    }

    ngh[i][0] = tmp[0];
    ngh[i][1] = tmp[1];
    for( k = 2; k < n_ngh[i]; k++ ) {
      for( j = 1; j < n_ngh[i]; j++ ) {
        if( tmp[2*j] == ngh[i][k-1] || tmp[2*j+1] == ngh[i][k-1] ) {
          if( tmp[2*j] == ngh[i][k-1] ) {
            ngh[i][k] = tmp[2*j+1];
          } else {
            ngh[i][k] = tmp[2*j];
          }
          tmp[2*j] = -1;
          tmp[2*j+1] = -1;
          break;
        }
      }
    }
  }

  *n_neighbours_return = n_ngh;
  *neighbours_return = ngh;

  FREE( tmp );

  return OK;

}

 
// Save an .obj file.
 
static void save_surface_obj( STRING filename,
                              int n_points,
                              Point coords[],
                              Vector normals[],
                              int n_elems,
                              int connec[] ) {

  int               i, j;

  FILE * fp = fopen( filename, "w" );
  fprintf( fp, "P 0.3 0.3 0.4 10 1 %d\n", n_points );

  // print the coords
  for( j = 0; j < n_points; j++ ) {
    fprintf( fp, "%g %g %g\n", coords[j].coords[0], 
             coords[j].coords[1], coords[j].coords[2] );
    // fprintf( fp, "%g %g %g\n", normals[j].coords[0], 
    //          normals[j].coords[1], normals[j].coords[2] );
  }
  fprintf( fp, "\n" );

  // print the normals
  for( j = 0; j < n_points; j++ ) {
    fprintf( fp, "%g %g %g\n", normals[j].coords[0], 
             normals[j].coords[1], normals[j].coords[2] );
  }

  // The connectivity - part 1.
  fprintf( fp, "\n" );
  fprintf( fp, "%d\n", n_elems );
  fprintf( fp, "0 1 1 1 1\n\n" );

  for( i = 1; i <= n_elems; i++ ) {
    fprintf( fp, "%d ", 3*i );
    if( i%8 == 0 ) fprintf( fp, "\n" );
  }

  // The connectivity - part 2.

  int count = 0;
  for( j = 0; j < 3*n_elems; j++ ) {
    if( count%8 == 0 ) fprintf( fp, "\n" );
    fprintf( fp, "%d ", connec[j] );
    count++;
  }
  fclose( fp );
}
