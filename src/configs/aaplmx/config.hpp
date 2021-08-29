#ifndef _TBLIS_CONFIGS_AAPLMX_CONFIG_HPP_
#define _TBLIS_CONFIGS_AAPLMX_CONFIG_HPP_

#include "configs/config_builder.hpp"

EXTERN_BLIS_GEMM_UKR(bli_sgemm_aaplmx_mac_32x32);
EXTERN_BLIS_GEMM_UKR(bli_dgemm_aaplmx_mac_32x16);

namespace tblis
{

EXTERN_PACK_NN_UKR(float, aaplmx_spackm_asm_32xk_sidem);
EXTERN_PACK_NN_UKR(float, aaplmx_spackm_asm_32xk_siden);
EXTERN_PACK_NN_UKR(double, aaplmx_dpackm_asm_32xk);
EXTERN_PACK_NN_UKR(double, aaplmx_dpackm_asm_16xk);

extern int aaplmx_check();

TBLIS_BEGIN_CONFIG(aaplmx)

    TBLIS_CONFIG_GEMM_MR(  32,   32, _, _)
    TBLIS_CONFIG_GEMM_NR(  32,   16, _, _)
    TBLIS_CONFIG_GEMM_KR(   4,    4, _, _)
    TBLIS_CONFIG_GEMM_MC( 512,  512, _, _)
    TBLIS_CONFIG_GEMM_NC(4096, 4096, _, _)
    TBLIS_CONFIG_GEMM_KC(2048, 2048, _, _)

    // TBLIS_CONFIG_M_THREAD_RATIO(_,3,_,_)
    // TBLIS_CONFIG_N_THREAD_RATIO(_,2,_,_)
    // TBLIS_CONFIG_MR_MAX_THREAD(_,1,_,_)
    // TBLIS_CONFIG_NR_MAX_THREAD(_,4,_,_)

    TBLIS_CONFIG_GEMM_WRAP_UKR(bli_sgemm_aaplmx_mac_32x32,
                               bli_dgemm_aaplmx_mac_32x16,
                               _,
                               _)
    TBLIS_CONFIG_PACK_NN_MR_UKR(aaplmx_spackm_asm_32xk_sidem, aaplmx_dpackm_asm_32xk, _, _)
    TBLIS_CONFIG_PACK_NN_NR_UKR(aaplmx_spackm_asm_32xk_siden, aaplmx_dpackm_asm_16xk, _, _)

    TBLIS_CONFIG_GEMM_ROW_MAJOR(false, false, _, _)
    TBLIS_CONFIG_GEMM_FLIP_UKR(true, true, _, _)

    TBLIS_CONFIG_CHECK(aaplmx_check)

TBLIS_END_CONFIG

}

#endif
