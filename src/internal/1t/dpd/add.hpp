#ifndef _TBLIS_INTERNAL_1T_DPD_ADD_HPP_
#define _TBLIS_INTERNAL_1T_DPD_ADD_HPP_

#include "util/thread.h"
#include "util/basic_types.h"
#include "configs/configs.hpp"

namespace tblis
{
namespace internal
{

void add(type_t type, const communicator& comm, const config& cfg,
         const scalar& alpha, bool conj_A, const dpd_varray_view<char>& A,
         const dim_vector& idx_A,
         const dim_vector& idx_A_AB,
         const scalar&  beta, bool conj_B, const dpd_varray_view<char>& B,
         const dim_vector& idx_B,
         const dim_vector& idx_B_AB);

}
}

#endif
