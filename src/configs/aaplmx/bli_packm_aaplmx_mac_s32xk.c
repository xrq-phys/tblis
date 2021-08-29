/*

   Kernel for the Apple matrix coprocessor.

*/

#include "blis.h"

#include "amx.h"
#include "amx_ext.h"


void bli_spackm_aaplmx_mac_32xk_simp
     (
       conj_t              conja,
       dim_t               k0,
       float*     restrict kappa,
       float*     restrict a, inc_t inca, inc_t lda,
       float*     restrict p,             inc_t ldp
     )
{
    const int gs    = ( inca != 1 && lda != 1 );
    const int unitk = ( *kappa == 1.0 );

    // As current the RE work has not discovered any
    //  broadcasting-load instruction yet, use this
    //  halfway solution for kappa.
    static float kappa_c[16] = { 0 };
    if ( kappa_c[0] != *kappa )
        for ( int i = 0; i < 16; ++i )
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
            for ( ; k >= 4; k -= 4 )
            {
                AMX_MEM( LDX, a + lda * 0 + 0 , 0 );
                AMX_MEM( LDX, a + lda * 0 + 16, 1 );
                AMX_MEM( LDX, a + lda * 1 + 0 , 2 );
                AMX_MEM( LDX, a + lda * 1 + 16, 3 );
                AMX_MEM( LDX, a + lda * 2 + 0 , 4 );
                AMX_MEM( LDX, a + lda * 2 + 16, 5 );
                AMX_MEM( LDX, a + lda * 3 + 0 , 6 );
                AMX_MEM( LDX, a + lda * 3 + 16, 7 );

                if ( unitk )
                {
                    // Write as-is.
                    AMX_MEM( STX, p + ldp * 0 + 0 , 0 );
                    AMX_MEM( STX, p + ldp * 0 + 16, 1 );
                    AMX_MEM( STX, p + ldp * 1 + 0 , 2 );
                    AMX_MEM( STX, p + ldp * 1 + 16, 3 );
                    AMX_MEM( STX, p + ldp * 2 + 0 , 4 );
                    AMX_MEM( STX, p + ldp * 2 + 16, 5 );
                    AMX_MEM( STX, p + ldp * 3 + 0 , 6 );
                    AMX_MEM( STX, p + ldp * 3 + 16, 7 );
                }
                else // if ( !unitk )
                {
                    // Scale and write Z.
                    AMX_FMUL32_SELCOL_REGALIGNED( 0, 0, 0, 0 );
                    AMX_FMUL32_SELCOL_REGALIGNED( 1, 1, 0, 0 );
                    AMX_FMUL32_SELCOL_REGALIGNED( 2, 2, 0, 0 );
                    AMX_FMUL32_SELCOL_REGALIGNED( 3, 3, 0, 0 );
                    AMX_FMUL32_SELCOL_REGALIGNED( 4, 4, 0, 0 );
                    AMX_FMUL32_SELCOL_REGALIGNED( 5, 5, 0, 0 );
                    AMX_FMUL32_SELCOL_REGALIGNED( 6, 6, 0, 0 );
                    AMX_FMUL32_SELCOL_REGALIGNED( 7, 7, 0, 0 );

                    AMX_MEM( STZ, p + ldp * 0 + 0 , 0 * 4 );
                    AMX_MEM( STZ, p + ldp * 0 + 16, 1 * 4 );
                    AMX_MEM( STZ, p + ldp * 1 + 0 , 2 * 4 );
                    AMX_MEM( STZ, p + ldp * 1 + 16, 3 * 4 );
                    AMX_MEM( STZ, p + ldp * 2 + 0 , 4 * 4 );
                    AMX_MEM( STZ, p + ldp * 2 + 16, 5 * 4 );
                    AMX_MEM( STZ, p + ldp * 3 + 0 , 6 * 4 );
                    AMX_MEM( STZ, p + ldp * 3 + 16, 7 * 4 );
                }

                a += lda * 4;
                p += ldp * 4;
            }
            for ( ; k >= 1; k -= 1 )
            {
                AMX_MEM( LDX, a + 0 , 0 );
                AMX_MEM( LDX, a + 16, 1 );

                if ( unitk )
                {
                    AMX_MEM( STX, p + 0 , 0 );
                    AMX_MEM( STX, p + 16, 1 );
                }
                else // if ( !unitk )
                {
                    AMX_FMUL32_SELCOL_REGALIGNED( 0, 0, 0, 0 );
                    AMX_FMUL32_SELCOL_REGALIGNED( 1, 1, 0, 0 );

                    AMX_MEM( STZ, p + 0 , 0 * 4 );
                    AMX_MEM( STZ, p + 16, 1 * 4 );
                }

                a += lda;
                p += ldp;
            }
        }
        else // if ( inca != 1 ) hence ( lda == 1 )
        {
            for ( ; k >= 16; k -= 16 )
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
                AMX_FMUL32_SELROW_REGALIGNED( 0, 0, 0, 0 );
                AMX_FMUL32_SELROW_REGALIGNED( 1, 0, 1, 0 );
                AMX_FMUL32_SELROW_REGALIGNED( 2, 0, 2, 0 );
                AMX_FMUL32_SELROW_REGALIGNED( 3, 0, 3, 0 );
                AMX_FMUL32_SELROW_REGALIGNED( 4, 0, 4, 0 );
                AMX_FMUL32_SELROW_REGALIGNED( 5, 0, 5, 0 );
                AMX_FMUL32_SELROW_REGALIGNED( 6, 0, 6, 0 );
                AMX_FMUL32_SELROW_REGALIGNED( 7, 0, 7, 0 );

                AMX_MEM( LDY, a + inca * 8 , 0 );
                AMX_MEM( LDY, a + inca * 9 , 1 );
                AMX_MEM( LDY, a + inca * 10, 2 );
                AMX_MEM( LDY, a + inca * 11, 3 );
                AMX_MEM( LDY, a + inca * 12, 4 );
                AMX_MEM( LDY, a + inca * 13, 5 );
                AMX_MEM( LDY, a + inca * 14, 6 );
                AMX_MEM( LDY, a + inca * 15, 7 );
                AMX_FMUL32_SELROW_REGALIGNED( 8 , 0, 0, 0 );
                AMX_FMUL32_SELROW_REGALIGNED( 9 , 0, 1, 0 );
                AMX_FMUL32_SELROW_REGALIGNED( 10, 0, 2, 0 );
                AMX_FMUL32_SELROW_REGALIGNED( 11, 0, 3, 0 );
                AMX_FMUL32_SELROW_REGALIGNED( 12, 0, 4, 0 );
                AMX_FMUL32_SELROW_REGALIGNED( 13, 0, 5, 0 );
                AMX_FMUL32_SELROW_REGALIGNED( 14, 0, 6, 0 );
                AMX_FMUL32_SELROW_REGALIGNED( 15, 0, 7, 0 );

                // Load and transpose: block 2.
                AMX_MEM( LDY, a + inca * (16 + 0), 0 );
                AMX_MEM( LDY, a + inca * (16 + 1), 1 );
                AMX_MEM( LDY, a + inca * (16 + 2), 2 );
                AMX_MEM( LDY, a + inca * (16 + 3), 3 );
                AMX_MEM( LDY, a + inca * (16 + 4), 4 );
                AMX_MEM( LDY, a + inca * (16 + 5), 5 );
                AMX_MEM( LDY, a + inca * (16 + 6), 6 );
                AMX_MEM( LDY, a + inca * (16 + 7), 7 );
                AMX_FMUL32_SELROW_REGALIGNED( 0, 0, 0, 1 );
                AMX_FMUL32_SELROW_REGALIGNED( 1, 0, 1, 1 );
                AMX_FMUL32_SELROW_REGALIGNED( 2, 0, 2, 1 );
                AMX_FMUL32_SELROW_REGALIGNED( 3, 0, 3, 1 );
                AMX_FMUL32_SELROW_REGALIGNED( 4, 0, 4, 1 );
                AMX_FMUL32_SELROW_REGALIGNED( 5, 0, 5, 1 );
                AMX_FMUL32_SELROW_REGALIGNED( 6, 0, 6, 1 );
                AMX_FMUL32_SELROW_REGALIGNED( 7, 0, 7, 1 );

                AMX_MEM( LDY, a + inca * (16 + 8 ), 0 );
                AMX_MEM( LDY, a + inca * (16 + 9 ), 1 );
                AMX_MEM( LDY, a + inca * (16 + 10), 2 );
                AMX_MEM( LDY, a + inca * (16 + 11), 3 );
                AMX_MEM( LDY, a + inca * (16 + 12), 4 );
                AMX_MEM( LDY, a + inca * (16 + 13), 5 );
                AMX_MEM( LDY, a + inca * (16 + 14), 6 );
                AMX_MEM( LDY, a + inca * (16 + 15), 7 );
                AMX_FMUL32_SELROW_REGALIGNED( 8 , 0, 0, 1 );
                AMX_FMUL32_SELROW_REGALIGNED( 9 , 0, 1, 1 );
                AMX_FMUL32_SELROW_REGALIGNED( 10, 0, 2, 1 );
                AMX_FMUL32_SELROW_REGALIGNED( 11, 0, 3, 1 );
                AMX_FMUL32_SELROW_REGALIGNED( 12, 0, 4, 1 );
                AMX_FMUL32_SELROW_REGALIGNED( 13, 0, 5, 1 );
                AMX_FMUL32_SELROW_REGALIGNED( 14, 0, 6, 1 );
                AMX_FMUL32_SELROW_REGALIGNED( 15, 0, 7, 1 );

                for ( int i = 0; i < 16; ++i )
                {
                    AMX_MEM( STZ, p + ldp * i + 0 , i * 4 + 0 );
                    AMX_MEM( STZ, p + ldp * i + 16, i * 4 + 1 );
                }
                a += 16;
                p += ldp * 16;
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

