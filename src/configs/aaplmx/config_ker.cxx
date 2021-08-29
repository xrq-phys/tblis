#include "config.hpp"

#include "util/cpuid.hpp"

#define PACKM_PARAMS(ctype) \
     ( \
       conj_t       conja, \
       dim_t        n_,    \
       const ctype* kappa, \
       const ctype* a, inc_t inca_, inc_t lda_, \
       ctype*       p,              inc_t ldp_  \
     )
extern "C" void bli_dpackm_aaplmx_mac_16xk_simp PACKM_PARAMS(double);
extern "C" void bli_dpackm_aaplmx_mac_32xk_simp PACKM_PARAMS(double);
extern "C" void bli_spackm_aaplmx_mac_32xk_simp PACKM_PARAMS(float);

namespace tblis
{

void aaplmx_spackm_asm_32xk_sidem(len_type m, len_type k,
                            const void* alpha, bool conj,
                            const void* p_a, stride_type rs_a, stride_type cs_a,
                            const void* p_d, stride_type inc_d,
                            const void* p_e, stride_type inc_e,
                            void* p_ap)
{
    int gs    = rs_a != 1 && cs_a != 1;
    int unitk = *((float *)alpha) == float(1.0);
    if (m == 32 && !gs && unitk)
    {
        bli_spackm_aaplmx_mac_32xk_simp(conj ? BLIS_CONJUGATE : BLIS_NO_CONJUGATE, k,
                                        reinterpret_cast<const float*>(alpha),
                                        reinterpret_cast<const float*>(p_a), rs_a, cs_a,
                                        reinterpret_cast<float*>(p_ap), 32);
    }
    else
    {
        pack_nn_ukr_def<aaplmx_config, float, matrix_constants::MAT_A>
            (m, k, alpha, conj, p_a, rs_a, cs_a, p_d, inc_d, p_e, inc_e, p_ap);
    }
}

void aaplmx_spackm_asm_32xk_siden(len_type m, len_type k,
                            const void* alpha, bool conj,
                            const void* p_a, stride_type rs_a, stride_type cs_a,
                            const void* p_d, stride_type inc_d,
                            const void* p_e, stride_type inc_e,
                            void* p_ap)
{
    int gs    = rs_a != 1 && cs_a != 1;
    int unitk = *((float *)alpha) == float(1.0);
    if (m == 32 && !gs && unitk)
    {
        bli_spackm_aaplmx_mac_32xk_simp(conj ? BLIS_CONJUGATE : BLIS_NO_CONJUGATE, k,
                                        reinterpret_cast<const float*>(alpha),
                                        reinterpret_cast<const float*>(p_a), rs_a, cs_a,
                                        reinterpret_cast<float*>(p_ap), 32);
    }
    else
    {
        pack_nn_ukr_def<aaplmx_config, float, matrix_constants::MAT_B>
            (m, k, alpha, conj, p_a, rs_a, cs_a, p_d, inc_d, p_e, inc_e, p_ap);
    }
}

void aaplmx_dpackm_asm_32xk(len_type m, len_type k,
                            const void* alpha, bool conj,
                            const void* p_a, stride_type rs_a, stride_type cs_a,
                            const void* p_d, stride_type inc_d,
                            const void* p_e, stride_type inc_e,
                            void* p_ap)
{
    int gs    = rs_a != 1 && cs_a != 1;
    int unitk = *((double *)alpha) == double(1.0);
    if (m == 32 && !gs && unitk)
    {
        bli_dpackm_aaplmx_mac_32xk_simp(conj ? BLIS_CONJUGATE : BLIS_NO_CONJUGATE, k,
                                        reinterpret_cast<const double*>(alpha),
                                        reinterpret_cast<const double*>(p_a), rs_a, cs_a,
                                        reinterpret_cast<double*>(p_ap), 32);
    }
    else
    {
        pack_nn_ukr_def<aaplmx_config, double, matrix_constants::MAT_A>
            (m, k, alpha, conj, p_a, rs_a, cs_a, p_d, inc_d, p_e, inc_e, p_ap);
    }
}

void aaplmx_dpackm_asm_16xk(len_type m, len_type k,
                            const void* alpha, bool conj,
                            const void* p_a, stride_type rs_a, stride_type cs_a,
                            const void* p_d, stride_type inc_d,
                            const void* p_e, stride_type inc_e,
                            void* p_ap)
{
    int gs    = rs_a != 1 && cs_a != 1;
    int unitk = *((double *)alpha) == double(1.0);
    if (m == 16 && !gs && unitk)
    {
        bli_dpackm_aaplmx_mac_16xk_simp(conj ? BLIS_CONJUGATE : BLIS_NO_CONJUGATE, k,
                                        reinterpret_cast<const double*>(alpha),
                                        reinterpret_cast<const double*>(p_a), rs_a, cs_a,
                                        reinterpret_cast<double*>(p_ap), 16);
    }
    else
    {
        pack_nn_ukr_def<aaplmx_config, double, matrix_constants::MAT_B>
            (m, k, alpha, conj, p_a, rs_a, cs_a, p_d, inc_d, p_e, inc_e, p_ap);
    }
}

TBLIS_CONFIG_INSTANTIATE(aaplmx);

}
