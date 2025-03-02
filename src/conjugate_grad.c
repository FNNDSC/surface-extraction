/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include  <conjugate.h>

struct  conjugate_struct
{
    BOOLEAN  first_call;
    int      n_parameters;
    Real     *g;
    Real     *h;
};

public  conjugate   initialize_conjugate_gradient(
    int       n_parameters )
{
    conjugate  con;

    ALLOC( con, 1 );

    con->n_parameters = n_parameters;
    con->first_call = TRUE;

    if( n_parameters > 0 )
    {
        ALLOC( con->g, n_parameters );
        ALLOC( con->h, n_parameters );
    }

    return( con );
}

public  void   delete_conjugate_gradient(
    conjugate   con )
{
    if( con->n_parameters >= 0 )
    {
        FREE( con->g );
        FREE( con->h );
    }

    FREE( con );
}

public  BOOLEAN  get_conjugate_unit_direction(
    conjugate   con,
    Real        derivative[],
    Real        unit_dir[],
    BOOLEAN     *changed )
{
    int      p;
    Real     len, gg, dgg, gam;
    BOOLEAN  found;

    found = TRUE;

    if( con->first_call )
    {
        con->first_call = FALSE;

        for_less( p, 0, con->n_parameters )
        {
            con->g[p] = -derivative[p];
            con->h[p] = -derivative[p];
            unit_dir[p] = -derivative[p];
        }

        *changed = TRUE;
    }
    else
    {
        *changed = FALSE;

        gg = 0.0;
        dgg = 0.0;
        for_less( p, 0, con->n_parameters )
        {
            gg += con->g[p] * con->g[p];
            dgg += (derivative[p] + con->g[p]) * derivative[p];
/*
            dgg += derivative[p] * derivative[p];
*/
        }

        if( gg != 0.0 )
            gam = dgg / gg;
        else
        {
            gam = 0.0;
            found = FALSE;
        }

        for_less( p, 0, con->n_parameters )
        {
            con->g[p] = -derivative[p];
            con->h[p] = con->g[p] + gam * con->h[p];
            unit_dir[p] = con->h[p];
        }
    }

    len = 0.0;
    for_less( p, 0, con->n_parameters )
        len += unit_dir[p] * unit_dir[p];

    if( len == 0.0 )
        return( FALSE );

    len = sqrt( len );

    for_less( p, 0, con->n_parameters )
        unit_dir[p] /= len;

    return( found );
}
