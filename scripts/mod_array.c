#define _GNU_SOURCE

#include <alloca.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <sqlite3ext.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

SQLITE_EXTENSION_INIT1

/* Defines array access functions for SQLite3.

   The idea is to store the bulk data (tables 'expresses' and 'covers')
   in some compact arrays.  Originally, those were stored in BLOBs
   inside the database, but accessing them there is awkward at best.
   Now they are stored in an auxiliary file, and the 'samples' table
   stores pointers into that file.

   The auxiliary file is written by stringtie-db during import of
   samples.  From SQLite, it is only possible to read it.  To do so,
   first load this extension, then execute `attachArray(<path>)`.  In
   Python, it might look like this:

        con = sqlite3.connect(database_file)
        con.enable_load_extension(True)
        con.load_extension("mod_array.so")
        con.enable_load_extension(False)            # to prevent attacks from user provided SQL
        con.execute("select attachArray(?)", (array_file,))


   After that, the array file is memory mapped and can be accessed
   through four custom functions:

   - getTpm( off, num, tid )

     This returns the expression of a transcript in a sample.  'off' and
     'num' must be the values of 'expresses_off' and 'expresses_num' of
     the sample, tid is the rowid of a transcript.

   - getFpkm( off, num, tid )

     Same, but returns FPKM.

   - getCov( off, num, tid )

     Same, but returns COV.

   - getECov( off, num, eid )

     This returns the coverage of an exon by a sample.  'off' and
     'num' must be the values of 'covers_off' and 'covers_num' of
     the sample, eid is the rowid of an exon.
*/



static uint16_t to_f16( double x )
{
    if( x < 1e-6 ) return 0 ;
    if( x > 9993152 ) return 65535 ;
    return (log(x)/log(10) + 6) / 13 * 65536 + 0.5 ;
}

static float from_f16_( uint16_t y )
{
    return exp( ((y * 13.0 / 65536.0) - 6) * log(10) ) ;
}

static const char* array_base = 0 ;
static long long   array_size = 0 ;
static float       from_f16[65536] ;

static void read_external( sqlite3_context *context, int argc, sqlite3_value **argv )
{
    if( argc != 3 )
        sqlite3_result_error( context, "expected three arguments", -1 ) ;
    else if( sqlite3_value_type( argv[0] ) == SQLITE_NULL )
        sqlite3_result_null( context ) ;
    else if( sqlite3_value_type( argv[0] ) != SQLITE_INTEGER )
        sqlite3_result_error( context, "first argument must be an integer", -1 ) ;
    else if( sqlite3_value_type( argv[1] ) != SQLITE_INTEGER )
        sqlite3_result_error( context, "second argument must be an integer", -1 ) ;
    else if( sqlite3_value_type( argv[2] ) != SQLITE_INTEGER )
        sqlite3_result_error( context, "third argument must be an integer", -1 ) ;
    else {
        sqlite3_int64 offset = sqlite3_value_int64( argv[0] ) ;
        sqlite3_int64 length = sqlite3_value_int64( argv[1] ) ;
        int index  = sqlite3_value_int( argv[2] ) ;

        const uint64_t *bitvec = (const uint64_t*)(array_base + offset) ;

        int block = index >> 9 ;
        int word  = (index & 0x1ff) >> 6 ;
        int bit   = index & 0x3f ;

        if( index < 0 || index > length || !((bitvec[ 10*block + 2 + word ] >> bit) & 1) )
        {
            sqlite3_result_double( context, 0 ) ;
        }
        else
        {
            int total = bitvec[ 10*block ] ;
            int subtotal = (bitvec[ 10*block + 1 ] >> (63-9*word)) & 0x1ff ;
            int count = __builtin_popcountll( bitvec[ 10*block + 2 + word ] << (63-bit) ) ;

            int nwords = (length+64) >> 6 ;
            int nblocks = (nwords+7) >> 3 ;
            int jndex = total+subtotal+count-1 ;

            int stride = (uintptr_t)sqlite3_user_data(context) >> 4 ;
            int field  = (uintptr_t)sqlite3_user_data(context) & 0xf ;

            const uint16_t *floatvec = (const uint16_t*)( bitvec + 2*nblocks + nwords ) ;
            sqlite3_result_double( context, from_f16[ floatvec[ stride*jndex + field ] ] ) ;
        }
    }
}

static void attach_array( sqlite3_context *context, int argc, sqlite3_value **argv )
{
    if( argc != 1 )
        sqlite3_result_error( context, "expected one argument", -1 ) ;
    else if( sqlite3_value_type( argv[0] ) != SQLITE_TEXT )
        sqlite3_result_error( context, "first argument must be text", -1 ) ;

    if( array_base ) munmap( (void*)array_base, array_size ) ;

    const char *fn = sqlite3_value_text( argv[0] ) ;
    int fd = open( fn, O_RDONLY ) ;

    if( fd < 0 )
        sqlite3_result_error( context, strerror(errno), -1 ) ;
    else
    {
        struct stat the_stat ;
        fstat( fd, &the_stat ) ;
        array_size = the_stat.st_size ;

        array_base = mmap( 0, array_size, PROT_READ, MAP_SHARED, fd, 0 ) ;
        if( array_base == (void*)(-1) )
            sqlite3_result_error( context, strerror(errno), -1 ) ;
        else
        {
            char *msg ;
            int msg_len = asprintf( &msg, "Attached %s at %p, length %lld.", fn, array_base, array_size ) ;
            sqlite3_result_text( context, msg, msg_len, free ) ;
        }
        close( fd ) ;
    }
}

// standard deviation aggregation function, effectively from
// http://sqlite.org/contrib (extension-functions.c)

struct stdev_state
{
    double m ;
    double s ;
    long long k ;
} ;

static void stdev_step( sqlite3_context *context, int argc, sqlite3_value **argv )
{
    if( argc != 1 ) sqlite3_result_error( context, "expected two arguments", -1 ) ;

    if( sqlite3_value_type(argv[0]) != SQLITE_NULL )
    {
        double x = sqlite3_value_double(argv[0]) ;
        struct stdev_state *p = sqlite3_aggregate_context(context, sizeof(*p));

        double d = x - p->m ;
        p->k++ ;
        p->m += d / p->k ;
        p->s += d * (x - p->m) ;
    }
}

static void stdev_final( sqlite3_context *context )
{
    struct stdev_state *p = sqlite3_aggregate_context(context, 0);
    if( !p || p->k < 2 )
        sqlite3_result_null( context ) ;
    else
        sqlite3_result_double( context, sqrt( p->s / (p->k - 1) ) ) ;
}


int sqlite3_modarray_init( sqlite3 *db, char **pzErrMsg, const sqlite3_api_routines *pApi )
{
    SQLITE_EXTENSION_INIT2(pApi);

    // populate the lookup table
    for( int i = 0 ; i != 65536 ; ++i ) from_f16[i] = from_f16_( i ) ;

    return sqlite3_create_function( db, "attachArray", 1, 0, 0, attach_array, 0, 0 )
         | sqlite3_create_function( db, "stdev",   1, SQLITE_DETERMINISTIC | SQLITE_INNOCUOUS, 0, 0, stdev_step, stdev_final )
         | sqlite3_create_function( db, "getCov",  3, SQLITE_DETERMINISTIC, (void*)0x30, read_external, 0, 0 )
         | sqlite3_create_function( db, "getFpkm", 3, SQLITE_DETERMINISTIC, (void*)0x31, read_external, 0, 0 )
         | sqlite3_create_function( db, "getTpm",  3, SQLITE_DETERMINISTIC, (void*)0x32, read_external, 0, 0 )
         | sqlite3_create_function( db, "getECov", 3, SQLITE_DETERMINISTIC, (void*)0x10, read_external, 0, 0 ) ;
}

