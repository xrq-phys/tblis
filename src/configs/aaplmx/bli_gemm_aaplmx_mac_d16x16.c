/* 

   Kernel for the Apple matrix coprocessor.

*/

#include "blis.h"
#include "amx.h"
#include "amx_ext.h"
#include <stdlib.h>
#include <assert.h>

// Const memory for zeroizing as
//  zeroizing instruction is not found yet.
const uint8_t amx_zeros[64] = { 0 };

// BLIs Single-precision GEMM for Apple Matrix Coprocessor,
//   implemented with MACros, of size 32 * 32.
void bli_sgemm_aaplmx_mac_32x32
     (
       dim_t               k,
       float*     restrict alpha,
       float*     restrict a,
       float*     restrict b,
       float*     restrict beta,
       float*     restrict c0, inc_t rs_c0, inc_t cs_c0,
       auxinfo_t* restrict data,
       cntx_t*    restrict cntx
     )
{
    void* a_next = bli_auxinfo_next_a( data );
    void* b_next = bli_auxinfo_next_b( data );
    inc_t rs_c = rs_c0;
    inc_t cs_c = cs_c0;
    float* c = c0;
    float* c_= c0;

    // As current the RE work has not discovered any
    //  broadcasting-load instruction yet, use this
    //  halfway solution for alpha & beta.
    static float alphac[16] = { 0 };
    static float beta_c[16] = { 0 };

    // Isotropic kernel.
    if ( cs_c == 1 && rs_c != 1 )
        return bli_sgemm_aaplmx_mac_32x32
            (
              k, alpha, b, a,
              beta, c, cs_c, rs_c,
              data, cntx
            );

    // Direct support for generic strided storage is not available yet due to no
    //  instruction supporting scattered-storage is found til now.
    if ( rs_c != 1 )
    {
        c = ( float* )malloc( 32*32 * sizeof( float ) );
        c_= c;
        rs_c = 1;
        cs_c = 32;
        if ( *beta != 0.0 )
            for ( int j = 0; j < 32; ++j )
                for ( int i = 0; i < 32; ++i )
                    c_[ i * rs_c + j * cs_c ] = c0[ i * rs_c0 + j * cs_c0 ];
    }

    // Duplicate alpha & beta.
    if ( alphac[0] != *alpha )
        for ( int i = 0; i < 16; ++i )
            alphac[i] = *alpha;

    if ( beta_c[0] != *beta )
        for ( int i = 0; i < 16; ++i )
            beta_c[i] = *beta;

    // As the BLIS API has no self-finalization,
    //  AMX_START and AMX_STOP has to be included within
    //  kernels. They do not seem to take much time,
    //  either.
    AMX_START();

    // Zeroize Z.
    // AMX_START seems clearing already?
    AMX_MEM( LDY, amx_zeros, 0 );
    AMX_FMUL32_COMMON_REGALIGNED( 0, 0, 0 );
    AMX_FMUL32_COMMON_REGALIGNED( 0, 0, 1 );
    AMX_FMUL32_COMMON_REGALIGNED( 0, 0, 2 );
    AMX_FMUL32_COMMON_REGALIGNED( 0, 0, 3 );

#pragma nounroll
    for ( ; k >= 4; k -= 4 )
    {
        AMX_MEM( LDX, a + 32 * 0 + 0 , 0 ); // A column 0.
        AMX_MEM( LDX, a + 32 * 0 + 16, 1 );
        AMX_MEM( LDX, a + 32 * 1 + 0 , 2 ); // A column 1.
        AMX_MEM( LDX, a + 32 * 1 + 16, 3 );
        AMX_MEM( LDX, a + 32 * 2 + 0 , 4 ); // A column 2.
        AMX_MEM( LDX, a + 32 * 2 + 16, 5 );
        AMX_MEM( LDX, a + 32 * 3 + 0 , 6 ); // A column 3.
        AMX_MEM( LDX, a + 32 * 3 + 16, 7 );

        AMX_MEM( LDY, b + 32 * 0 + 0 , 0 ); // B row 0.
        AMX_MEM( LDY, b + 32 * 0 + 16, 1 );
        AMX_MEM( LDY, b + 32 * 1 + 0 , 2 ); // B row 1.
        AMX_MEM( LDY, b + 32 * 1 + 16, 3 );
        AMX_MEM( LDY, b + 32 * 2 + 0 , 4 ); // B row 2.
        AMX_MEM( LDY, b + 32 * 2 + 16, 5 );
        AMX_MEM( LDY, b + 32 * 3 + 0 , 6 ); // B row 3.
        AMX_MEM( LDY, b + 32 * 3 + 16, 7 );

        AMX_FMA32_COMMON_REGALIGNED( 0, 0, 0 ); // Block (0, 0)
        AMX_FMA32_COMMON_REGALIGNED( 1, 0, 1 ); // Block (1, 0)
        AMX_FMA32_COMMON_REGALIGNED( 0, 1, 2 ); // Block (0, 1)
        AMX_FMA32_COMMON_REGALIGNED( 1, 1, 3 ); // Block (1, 1)

        AMX_FMA32_COMMON_REGALIGNED( 2, 2, 0 );
        AMX_FMA32_COMMON_REGALIGNED( 3, 2, 1 );
        AMX_FMA32_COMMON_REGALIGNED( 2, 3, 2 );
        AMX_FMA32_COMMON_REGALIGNED( 3, 3, 3 );

        AMX_FMA32_COMMON_REGALIGNED( 4, 4, 0 );
        AMX_FMA32_COMMON_REGALIGNED( 5, 4, 1 );
        AMX_FMA32_COMMON_REGALIGNED( 4, 5, 2 );
        AMX_FMA32_COMMON_REGALIGNED( 5, 5, 3 );

        AMX_FMA32_COMMON_REGALIGNED( 6, 6, 0 );
        AMX_FMA32_COMMON_REGALIGNED( 7, 6, 1 );
        AMX_FMA32_COMMON_REGALIGNED( 6, 7, 2 );
        AMX_FMA32_COMMON_REGALIGNED( 7, 7, 3 );

        // Address forward.
        a += 4 * 32;
        b += 4 * 32;
    }
#pragma nounroll
    for ( ; k >= 1; k -= 1 )
    {
        AMX_MEM( LDX, a + 0 , 0 );
        AMX_MEM( LDX, a + 16, 1 );

        AMX_MEM( LDY, b + 0 , 0 );
        AMX_MEM( LDY, b + 16, 1 );

        AMX_FMA32_COMMON_REGALIGNED( 0, 0, 0 );
        AMX_FMA32_COMMON_REGALIGNED( 1, 0, 1 );
        AMX_FMA32_COMMON_REGALIGNED( 0, 1, 2 );
        AMX_FMA32_COMMON_REGALIGNED( 1, 1, 3 );

        a += 32;
        b += 32;
    }

    // Load alpha & beta.
    AMX_MEM( LDY, alphac, 0 );
    AMX_MEM( LDY, beta_c, 1 );

    // Multiply by alpha.
    if ( *alpha != 1.0 )
        for ( int i = 0; i < 16; ++i ) {
            AMX_EXTRX_REGALIGNED( i * 4 + 0, 4 );
            AMX_EXTRX_REGALIGNED( i * 4 + 1, 5 );
            AMX_EXTRX_REGALIGNED( i * 4 + 2, 6 );
            AMX_EXTRX_REGALIGNED( i * 4 + 3, 7 );

            AMX_FMUL32_SELCOL_REGALIGNED( i, 4, 0, 0 );
            AMX_FMUL32_SELCOL_REGALIGNED( i, 5, 0, 1 );
            AMX_FMUL32_SELCOL_REGALIGNED( i, 6, 0, 2 );
            AMX_FMUL32_SELCOL_REGALIGNED( i, 7, 0, 3 );
        }

    __asm__ volatile
    (
      "prfm PLDL2STRM, [%[a_next], 64*0] \n\t"
      "prfm PLDL2STRM, [%[a_next], 64*1] \n\t"
      "prfm PLDL2STRM, [%[a_next], 64*2] \n\t"
      "prfm PLDL2STRM, [%[a_next], 64*3] \n\t"
      "prfm PLDL2STRM, [%[b_next], 64*0] \n\t"
      "prfm PLDL2STRM, [%[b_next], 64*1] \n\t"
      "prfm PLDL2STRM, [%[b_next], 64*2] \n\t"
      "prfm PLDL2STRM, [%[b_next], 64*3] \n\t"
      :
      : [a_next] "r" (a_next),
        [b_next] "r" (b_next)
    );

    // Load and multiply by beta.
    // Write into Z registers.
    if ( *beta != 0.0 )
    {
        float *c_ldr = c;
        if ( rs_c == 1 )
        {
            // Load blocks (0, 0) and (1, 0).
            for (int i = 0; i < 16; ++i) {
                AMX_MEM( LDX, c_ldr + i * cs_c + 0 , 2 );
                AMX_MEM( LDX, c_ldr + i * cs_c + 16, 3 );

                AMX_FMA32_SELCOL_REGALIGNED( i, 2, 1, 0 );
                AMX_FMA32_SELCOL_REGALIGNED( i, 3, 1, 1 );
            }
            c_ldr += 16 * cs_c;

            // Load blocks (0, 1) and (1, 1).
            for (int i = 0; i < 16; ++i) {
                AMX_MEM( LDX, c_ldr + i * cs_c + 0 , 2 );
                AMX_MEM( LDX, c_ldr + i * cs_c + 16, 3 );

                AMX_FMA32_SELCOL_REGALIGNED( i, 2, 1, 2 );
                AMX_FMA32_SELCOL_REGALIGNED( i, 3, 1, 3 );
            }
        }
    }

    if ( rs_c == 1 )
    {
        // Store blocks (0, 0) and (1, 0).
        for (int i = 0; i < 16; ++i) {
            AMX_MEM( STZ, c + i * cs_c + 0 , i * 4 + 0 );
            AMX_MEM( STZ, c + i * cs_c + 16, i * 4 + 1 );
        }
        c += 16 * cs_c;

        // Store blocks (0, 1) and (1, 1).
        for (int i = 0; i < 16; ++i) {
            AMX_MEM( STZ, c + i * cs_c + 0 , i * 4 + 2 );
            AMX_MEM( STZ, c + i * cs_c + 16, i * 4 + 3 );
        }
    }

    AMX_STOP();

    // Compatibility layer.
    if ( rs_c0 != 1 )
    {
        for ( int j = 0; j < 32; ++j )
            for ( int i = 0; i < 32; ++i )
                c0[ i * rs_c0 + j * cs_c0 ] = c_[ i * rs_c + j * cs_c ];

        free(c_);
    }
}


void bli_dgemm_aaplmx_mac_16x16
     (
       dim_t               k,
       double*    restrict alpha,
       double*    restrict a,
       double*    restrict b,
       double*    restrict beta,
       double*    restrict c, inc_t rs_c, inc_t cs_c,
       auxinfo_t* restrict data,
       cntx_t*    restrict cntx
     )
{
    void* a_next = bli_auxinfo_next_a( data );
    void* b_next = bli_auxinfo_next_b( data );

    // As current the RE work has not discovered any
    //  broadcasting-load instruction yet, use this
    //  halfway solution for alpha & beta.
    static double alphac[8] = { 0 };
    static double beta_c[8] = { 0 };

    // Isotropic kernel.
    if ( cs_c == 1 && rs_c != 1 )
        return bli_dgemm_aaplmx_mac_16x16
            (
              k, alpha, b, a,
              beta, c, cs_c, rs_c,
              data, cntx
            );

    // TODO: Support generic strided storage.
    assert( rs_c == 1 );

    // Duplicate alpha & beta.
    if ( alphac[0] != *alpha )
        for ( int i = 0; i < 8; ++i )
            alphac[i] = *alpha;

    if ( beta_c[0] != *beta )
        for ( int i = 0; i < 8; ++i )
            beta_c[i] = *beta;

    // As the BLIS API has no self-finalization,
    //  AMX_START and AMX_STOP has to be included within
    //  kernels. They do not seem to take much time,
    //  either.
    AMX_START();

    // Zeroize Z.
    AMX_MEM( LDY, amx_zeros, 0 );
    AMX_FMUL64_COMMON_REGALIGNED( 0, 0, 0 );
    AMX_FMUL64_COMMON_REGALIGNED( 0, 0, 1 );
    AMX_FMUL64_COMMON_REGALIGNED( 0, 0, 2 );
    AMX_FMUL64_COMMON_REGALIGNED( 0, 0, 3 );

#pragma nounroll
    for ( ; k >= 4; k -= 4 )
    {
        // NOTE: Comments on this nanokernel is written in a column-major fasion:
        //
        // An AMX_FMA multiplies vector column x and vector row y^T into a matrix
        //  formed by vector columns of z:
        //
        // e.g.
        //   AMX_FMA64_COMMON_REGALIGNED( i, j, l ); where 0<=i,j<8, 0<=l<8
        //  takes the i-th vector from 8 registers for x (x_i),
        //  j-th vector from 8 registers for y (y_j) and performs:
        //   z_{8 * 0 + l} += x_i * y_i[0];
        //   z_{8 * 1 + l} += x_i * y_i[1];
        //   z_{8 * 2 + l} += x_i * y_i[2];
        //   z_{8 * 3 + l} += x_i * y_i[3];
        //   z_{8 * 4 + l} += x_i * y_i[4];
        //   z_{8 * 5 + l} += x_i * y_i[5];
        //   z_{8 * 6 + l} += x_i * y_i[6];
        //   z_{8 * 7 + l} += x_i * y_i[7];
        //  within 1 clock cycle.
        //
        // ---
        // Notation:
        //  v_i or v_{i}: the i-th register in the register file for v;
        //  v[k]: element k in some vector.

        AMX_MEM( LDX, a + 16 * 0 + 0, 0 ); // Load A column 0.
        AMX_MEM( LDX, a + 16 * 0 + 8, 1 );
        AMX_MEM( LDX, a + 16 * 1 + 0, 2 ); // Load A column 1.
        AMX_MEM( LDX, a + 16 * 1 + 8, 3 );
        AMX_MEM( LDX, a + 16 * 2 + 0, 4 ); // Load A column 2.
        AMX_MEM( LDX, a + 16 * 2 + 8, 5 );
        AMX_MEM( LDX, a + 16 * 3 + 0, 6 ); // Load A column 3.
        AMX_MEM( LDX, a + 16 * 3 + 8, 7 );

        AMX_MEM( LDY, b + 16 * 0 + 0, 0 ); // Load B row 0.
        AMX_MEM( LDY, b + 16 * 0 + 8, 1 );
        AMX_MEM( LDY, b + 16 * 1 + 0, 2 ); // Load B row 1.
        AMX_MEM( LDY, b + 16 * 1 + 8, 3 );
        AMX_MEM( LDY, b + 16 * 2 + 0, 4 ); // Load B row 2.
        AMX_MEM( LDY, b + 16 * 2 + 8, 5 );
        AMX_MEM( LDY, b + 16 * 3 + 0, 6 ); // Load B row 3.
        AMX_MEM( LDY, b + 16 * 3 + 8, 7 );

        AMX_FMA64_COMMON_REGALIGNED( 0, 0, 0 ); // C[ 0:15, 0:15] += A[ 0:15, 0] * B[ 0:15, 0]
        AMX_FMA64_COMMON_REGALIGNED( 1, 0, 1 ); // C[16:31, 0:15] += A[16:31, 0] * B[ 0:15, 0]
        AMX_FMA64_COMMON_REGALIGNED( 0, 1, 2 ); // C[ 0:15,16:31] += A[ 0:15, 0] * B[16:31, 0]
        AMX_FMA64_COMMON_REGALIGNED( 1, 1, 3 ); // C[16:31,16:31] += A[16:31, 0] * B[16:31, 0]

        AMX_FMA64_COMMON_REGALIGNED( 2, 2, 0 ); // C[ 0:15, 0:15] += A[ 0:15, 1] * B[ 0:15, 1]
        AMX_FMA64_COMMON_REGALIGNED( 3, 2, 1 ); // C[16:31, 0:15] += A[16:31, 1] * B[ 0:15, 1]
        AMX_FMA64_COMMON_REGALIGNED( 2, 3, 2 ); // C[ 0:15,16:31] += A[ 0:15, 1] * B[16:31, 1]
        AMX_FMA64_COMMON_REGALIGNED( 3, 3, 3 ); // C[16:31,16:31] += A[16:31, 1] * B[16:31, 1]

        AMX_FMA64_COMMON_REGALIGNED( 4, 4, 0 ); // C[ 0:15, 0:15] += A[ 0:15, 2] * B[ 0:15, 2]
        AMX_FMA64_COMMON_REGALIGNED( 5, 4, 1 ); // C[16:31, 0:15] += A[16:31, 2] * B[ 0:15, 2]
        AMX_FMA64_COMMON_REGALIGNED( 4, 5, 2 ); // C[ 0:15,16:31] += A[ 0:15, 2] * B[16:31, 2]
        AMX_FMA64_COMMON_REGALIGNED( 5, 5, 3 ); // C[16:31,16:31] += A[16:31, 2] * B[16:31, 2]

        AMX_FMA64_COMMON_REGALIGNED( 6, 6, 0 ); // C[ 0:15, 0:15] += A[ 0:15, 3] * B[ 0:15, 3]
        AMX_FMA64_COMMON_REGALIGNED( 7, 6, 1 ); // C[16:31, 0:15] += A[16:31, 3] * B[ 0:15, 3]
        AMX_FMA64_COMMON_REGALIGNED( 6, 7, 2 ); // C[ 0:15,16:31] += A[ 0:15, 3] * B[16:31, 3]
        AMX_FMA64_COMMON_REGALIGNED( 7, 7, 3 ); // C[16:31,16:31] += A[16:31, 3] * B[16:31, 3]

        // Address forward.
        a += 4 * 16;
        b += 4 * 16;
    }
#pragma nounroll
    for ( ; k >= 1; k -= 1 )
    {
        AMX_MEM( LDX, a + 0, 0 );
        AMX_MEM( LDX, a + 8, 1 );

        AMX_MEM( LDY, b + 0, 0 );
        AMX_MEM( LDY, b + 8, 1 );

        AMX_FMA64_COMMON_REGALIGNED( 0, 0, 0 );
        AMX_FMA64_COMMON_REGALIGNED( 1, 0, 1 );
        AMX_FMA64_COMMON_REGALIGNED( 0, 1, 2 );
        AMX_FMA64_COMMON_REGALIGNED( 1, 1, 3 );

        a += 16;
        b += 16;
    }

    // Load alpha & beta.
    AMX_MEM( LDY, alphac, 0 );
    AMX_MEM( LDY, beta_c, 1 );

    // Multiply by alpha.
    if ( *alpha != 1.0 )
        for ( int i = 0; i < 8; ++i ) {
            AMX_EXTRX_REGALIGNED( i * 8 + 0, 4 );
            AMX_EXTRX_REGALIGNED( i * 8 + 1, 5 );
            AMX_EXTRX_REGALIGNED( i * 8 + 2, 6 );
            AMX_EXTRX_REGALIGNED( i * 8 + 3, 7 );

            AMX_FMUL64_SELCOL_REGALIGNED( i, 4, 0, 0 );
            AMX_FMUL64_SELCOL_REGALIGNED( i, 5, 0, 1 );
            AMX_FMUL64_SELCOL_REGALIGNED( i, 6, 0, 2 );
            AMX_FMUL64_SELCOL_REGALIGNED( i, 7, 0, 3 );
        }

    __asm__ volatile
    (
      "prfm PLDL2STRM, [%[a_next], 64*0] \n\t"
      "prfm PLDL2STRM, [%[a_next], 64*1] \n\t"
      "prfm PLDL2STRM, [%[a_next], 64*2] \n\t"
      "prfm PLDL2STRM, [%[a_next], 64*3] \n\t"
      "prfm PLDL2STRM, [%[b_next], 64*0] \n\t"
      "prfm PLDL2STRM, [%[b_next], 64*1] \n\t"
      "prfm PLDL2STRM, [%[b_next], 64*2] \n\t"
      "prfm PLDL2STRM, [%[b_next], 64*3] \n\t"
      :
      : [a_next] "r" (a_next),
        [b_next] "r" (b_next)
    );

    // Load and multiply by beta.
    // Write into Z registers.
    if ( *beta != 0.0 )
    {
        double *c_ldr = c;
        if ( rs_c == 1 )
        {
            // Load blocks (0, 0) and (1, 0).
            for (int i = 0; i < 8; ++i) {
                AMX_MEM( LDX, c_ldr + i * cs_c + 0, 2 );
                AMX_MEM( LDX, c_ldr + i * cs_c + 8, 3 );

                AMX_FMA64_SELCOL_REGALIGNED( i, 2, 1, 0 );
                AMX_FMA64_SELCOL_REGALIGNED( i, 3, 1, 1 );
            }
            c_ldr += 8 * cs_c;

            // Load blocks (0, 1) and (1, 1).
            for (int i = 0; i < 8; ++i) {
                AMX_MEM( LDX, c_ldr + i * cs_c + 0, 2 );
                AMX_MEM( LDX, c_ldr + i * cs_c + 8, 3 );

                AMX_FMA64_SELCOL_REGALIGNED( i, 2, 1, 2 );
                AMX_FMA64_SELCOL_REGALIGNED( i, 3, 1, 3 );
            }
        }
    }

    if ( rs_c == 1 )
    {
        // Store blocks (0, 0) and (1, 0).
        for (int i = 0; i < 8; ++i) {
            AMX_MEM( STZ, c + i * cs_c + 0, i * 8 + 0 );
            AMX_MEM( STZ, c + i * cs_c + 8, i * 8 + 1 );
        }
        c += 8 * cs_c;

        // Store blocks (0, 1) and (1, 1).
        for (int i = 0; i < 8; ++i) {
            AMX_MEM( STZ, c + i * cs_c + 0, i * 8 + 2 );
            AMX_MEM( STZ, c + i * cs_c + 8, i * 8 + 3 );
        }
    }

    AMX_STOP();
}

