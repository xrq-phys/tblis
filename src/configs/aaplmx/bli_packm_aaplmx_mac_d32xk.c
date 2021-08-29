/*

   Kernel for the Apple matrix coprocessor.

*/

#include "blis.h"

#include "amx.h"
#include "amx_ext.h"


void bli_dpackm_aaplmx_mac_32xk_simp
     (
       conj_t              conja,
       dim_t               k0,
       double*    restrict kappa,
       double*    restrict a, inc_t inca, inc_t lda,
       double*    restrict p,             inc_t ldp
     )
{
    const int gs    = ( inca != 1 && lda != 1 );
    const int unitk = ( *kappa == 1.0 );

    // As current the RE work has not discovered any
    //  broadcasting-load instruction yet, use this
    //  halfway solution for kappa.
    static double kappa_c[8] = { 0 };
    if ( kappa_c[0] != *kappa )
        for ( int i = 0; i < 8; ++i )
            kappa_c[i] = *kappa;

    // Local copy.
    dim_t k = k0;

    // -------------------------------------------------------------------------

    if ( !gs && unitk )
    {
        AMX_START();
        AMX_MEM( LDX, kappa_c, 0 );
        AMX_MEM( LDY, kappa_c, 0 );

        if ( inca == 1 )
        {
            for ( ; k >= 2; k -= 2 )
            {
                AMX_MEM( LDX, a + lda * 0 + 0 , 0 );
                AMX_MEM( LDX, a + lda * 0 + 8 , 1 );
                AMX_MEM( LDX, a + lda * 0 + 16, 2 );
                AMX_MEM( LDX, a + lda * 0 + 24, 3 );
                AMX_MEM( LDX, a + lda * 1 + 0 , 4 );
                AMX_MEM( LDX, a + lda * 1 + 8 , 5 );
                AMX_MEM( LDX, a + lda * 1 + 16, 6 );
                AMX_MEM( LDX, a + lda * 1 + 24, 7 );

                if ( unitk )
                {
                    // Write as-is.
                    AMX_MEM( STX, p + ldp * 0 + 0 , 0 );
                    AMX_MEM( STX, p + ldp * 0 + 8 , 1 );
                    AMX_MEM( STX, p + ldp * 0 + 16, 2 );
                    AMX_MEM( STX, p + ldp * 0 + 24, 3 );
                    AMX_MEM( STX, p + ldp * 1 + 0 , 4 );
                    AMX_MEM( STX, p + ldp * 1 + 8 , 5 );
                    AMX_MEM( STX, p + ldp * 1 + 16, 6 );
                    AMX_MEM( STX, p + ldp * 1 + 24, 7 );
                }
                else // if ( !unitk )
                {
                    // Scale and write Z.
                    AMX_FMUL64_SELCOL_REGALIGNED( 0, 0, 0, 0 );
                    AMX_FMUL64_SELCOL_REGALIGNED( 1, 1, 0, 0 );
                    AMX_FMUL64_SELCOL_REGALIGNED( 2, 2, 0, 0 );
                    AMX_FMUL64_SELCOL_REGALIGNED( 3, 3, 0, 0 );
                    AMX_FMUL64_SELCOL_REGALIGNED( 4, 4, 0, 0 );
                    AMX_FMUL64_SELCOL_REGALIGNED( 5, 5, 0, 0 );
                    AMX_FMUL64_SELCOL_REGALIGNED( 6, 6, 0, 0 );
                    AMX_FMUL64_SELCOL_REGALIGNED( 7, 7, 0, 0 );

                    AMX_MEM( STZ, p + ldp * 0 + 0 , 0 * 8 );
                    AMX_MEM( STZ, p + ldp * 0 + 8 , 1 * 8 );
                    AMX_MEM( STZ, p + ldp * 0 + 16, 2 * 8 );
                    AMX_MEM( STZ, p + ldp * 0 + 24, 3 * 8 );
                    AMX_MEM( STZ, p + ldp * 1 + 0 , 4 * 8 );
                    AMX_MEM( STZ, p + ldp * 1 + 8 , 5 * 8 );
                    AMX_MEM( STZ, p + ldp * 1 + 16, 6 * 8 );
                    AMX_MEM( STZ, p + ldp * 1 + 24, 7 * 8 );
                }

                a += lda * 2;
                p += ldp * 2;
            }
            for ( ; k >= 1; k -= 1 )
            {
                AMX_MEM( LDX, a + 0 , 0 );
                AMX_MEM( LDX, a + 8 , 1 );
                AMX_MEM( LDX, a + 16, 2 );
                AMX_MEM( LDX, a + 24, 3 );

                if ( unitk )
                {
                    AMX_MEM( STX, p + 0 , 0 );
                    AMX_MEM( STX, p + 8 , 1 );
                    AMX_MEM( STX, p + 16, 2 );
                    AMX_MEM( STX, p + 24, 3 );
                }
                else // if ( !unitk )
                {
                    AMX_FMUL64_SELCOL_REGALIGNED( 0, 0, 0, 0 );
                    AMX_FMUL64_SELCOL_REGALIGNED( 1, 1, 0, 0 );
                    AMX_FMUL64_SELCOL_REGALIGNED( 2, 2, 0, 0 );
                    AMX_FMUL64_SELCOL_REGALIGNED( 3, 3, 0, 0 );

                    AMX_MEM( STZ, p + 0 , 0 * 8 );
                    AMX_MEM( STZ, p + 8 , 1 * 8 );
                    AMX_MEM( STZ, p + 16, 2 * 8 );
                    AMX_MEM( STZ, p + 24, 3 * 8 );
                }

                a += lda;
                p += ldp;
            }
        }
        else // if ( inca != 1 ) hence ( lda == 1 )
        {
            for ( ; k >= 8; k -= 8 )
            {
                // Load and transpose: block 1.
                AMX_MEM( LDY, a + inca * 0, 0 );
                AMX_MEM( LDY, a + inca * 1, 1 );
                AMX_MEM( LDY, a + inca * 2, 2 );
                AMX_MEM( LDY, a + inca * 3, 3 );
                AMX_MEM( LDY, a + inca * 4, 4 );
                AMX_MEM( LDY, a + inca * 5, 5 );
                AMX_MEM( LDY, a + inca * 6, 6 );
                AMX_MEM( LDY, a + inca * 7, 7 );
                AMX_FMUL64_SELROW_REGALIGNED( 0, 0, 0, 0 );
                AMX_FMUL64_SELROW_REGALIGNED( 1, 0, 1, 0 );
                AMX_FMUL64_SELROW_REGALIGNED( 2, 0, 2, 0 );
                AMX_FMUL64_SELROW_REGALIGNED( 3, 0, 3, 0 );
                AMX_FMUL64_SELROW_REGALIGNED( 4, 0, 4, 0 );
                AMX_FMUL64_SELROW_REGALIGNED( 5, 0, 5, 0 );
                AMX_FMUL64_SELROW_REGALIGNED( 6, 0, 6, 0 );
                AMX_FMUL64_SELROW_REGALIGNED( 7, 0, 7, 0 );

                // Load and transpose: block 2.
                AMX_MEM( LDY, a + inca * (8 + 0), 0 );
                AMX_MEM( LDY, a + inca * (8 + 1), 1 );
                AMX_MEM( LDY, a + inca * (8 + 2), 2 );
                AMX_MEM( LDY, a + inca * (8 + 3), 3 );
                AMX_MEM( LDY, a + inca * (8 + 4), 4 );
                AMX_MEM( LDY, a + inca * (8 + 5), 5 );
                AMX_MEM( LDY, a + inca * (8 + 6), 6 );
                AMX_MEM( LDY, a + inca * (8 + 7), 7 );
                AMX_FMUL64_SELROW_REGALIGNED( 0, 0, 0, 1 );
                AMX_FMUL64_SELROW_REGALIGNED( 1, 0, 1, 1 );
                AMX_FMUL64_SELROW_REGALIGNED( 2, 0, 2, 1 );
                AMX_FMUL64_SELROW_REGALIGNED( 3, 0, 3, 1 );
                AMX_FMUL64_SELROW_REGALIGNED( 4, 0, 4, 1 );
                AMX_FMUL64_SELROW_REGALIGNED( 5, 0, 5, 1 );
                AMX_FMUL64_SELROW_REGALIGNED( 6, 0, 6, 1 );
                AMX_FMUL64_SELROW_REGALIGNED( 7, 0, 7, 1 );

                // Load and transpose: block 3.
                AMX_MEM( LDY, a + inca * (16 + 0), 0 );
                AMX_MEM( LDY, a + inca * (16 + 1), 1 );
                AMX_MEM( LDY, a + inca * (16 + 2), 2 );
                AMX_MEM( LDY, a + inca * (16 + 3), 3 );
                AMX_MEM( LDY, a + inca * (16 + 4), 4 );
                AMX_MEM( LDY, a + inca * (16 + 5), 5 );
                AMX_MEM( LDY, a + inca * (16 + 6), 6 );
                AMX_MEM( LDY, a + inca * (16 + 7), 7 );
                AMX_FMUL64_SELROW_REGALIGNED( 0, 0, 0, 2 );
                AMX_FMUL64_SELROW_REGALIGNED( 1, 0, 1, 2 );
                AMX_FMUL64_SELROW_REGALIGNED( 2, 0, 2, 2 );
                AMX_FMUL64_SELROW_REGALIGNED( 3, 0, 3, 2 );
                AMX_FMUL64_SELROW_REGALIGNED( 4, 0, 4, 2 );
                AMX_FMUL64_SELROW_REGALIGNED( 5, 0, 5, 2 );
                AMX_FMUL64_SELROW_REGALIGNED( 6, 0, 6, 2 );
                AMX_FMUL64_SELROW_REGALIGNED( 7, 0, 7, 2 );

                // Load and transpose: block 4.
                AMX_MEM( LDY, a + inca * (24 + 0), 0 );
                AMX_MEM( LDY, a + inca * (24 + 1), 1 );
                AMX_MEM( LDY, a + inca * (24 + 2), 2 );
                AMX_MEM( LDY, a + inca * (24 + 3), 3 );
                AMX_MEM( LDY, a + inca * (24 + 4), 4 );
                AMX_MEM( LDY, a + inca * (24 + 5), 5 );
                AMX_MEM( LDY, a + inca * (24 + 6), 6 );
                AMX_MEM( LDY, a + inca * (24 + 7), 7 );
                AMX_FMUL64_SELROW_REGALIGNED( 0, 0, 0, 3 );
                AMX_FMUL64_SELROW_REGALIGNED( 1, 0, 1, 3 );
                AMX_FMUL64_SELROW_REGALIGNED( 2, 0, 2, 3 );
                AMX_FMUL64_SELROW_REGALIGNED( 3, 0, 3, 3 );
                AMX_FMUL64_SELROW_REGALIGNED( 4, 0, 4, 3 );
                AMX_FMUL64_SELROW_REGALIGNED( 5, 0, 5, 3 );
                AMX_FMUL64_SELROW_REGALIGNED( 6, 0, 6, 3 );
                AMX_FMUL64_SELROW_REGALIGNED( 7, 0, 7, 3 );

                for ( int i = 0; i < 8; ++i )
                {
                    AMX_MEM( STZ, p + ldp * i + 0 , i * 8 + 0 );
                    AMX_MEM( STZ, p + ldp * i + 8 , i * 8 + 1 );
                    AMX_MEM( STZ, p + ldp * i + 16, i * 8 + 2 );
                    AMX_MEM( STZ, p + ldp * i + 24, i * 8 + 3 );
                }
                a += 8;
                p += ldp * 8;
            }
            for ( ; k >= 1; k -= 1 )
            {
                // Plain C.
                if ( unitk )
                    for ( int i = 0; i < 32; ++i )
                        p[i] = *( a + inca * i );
                else
                    for ( int i = 0; i < 32; ++i )
                        p[i] = ( *kappa ) * *( a + inca * i );

                a += 1;
                p += ldp;
            }
        }

        AMX_STOP();
    }
    else
    {
        // Handled by caller.
    }
}

