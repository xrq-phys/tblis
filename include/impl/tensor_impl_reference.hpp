#ifndef _TENSOR_IMPL_TENSOR_IMPL_REFERENCE_HPP_
#define _TENSOR_IMPL_TENSOR_IMPL_REFERENCE_HPP_

#include "impl/tensor_impl.hpp"

namespace tensor
{
namespace impl
{

template <typename T>
int tensor_mult_reference(T alpha, const Tensor<T>& A, const std::string& idx_A,
                                   const Tensor<T>& B, const std::string& idx_B,
                          T  beta,       Tensor<T>& C, const std::string& idx_C);

template <typename T>
int tensor_contract_reference(T alpha, const Tensor<T>& A, const std::string& idx_A,
                                       const Tensor<T>& B, const std::string& idx_B,
                              T  beta,       Tensor<T>& C, const std::string& idx_C);

template <typename T>
int tensor_weight_reference(T alpha, const Tensor<T>& A, const std::string& idx_A,
                                     const Tensor<T>& B, const std::string& idx_B,
                            T  beta,       Tensor<T>& C, const std::string& idx_C);

template <typename T>
int tensor_outer_prod_reference(T alpha, const Tensor<T>& A, const std::string& idx_A,
                                         const Tensor<T>& B, const std::string& idx_B,
                                T  beta,       Tensor<T>& C, const std::string& idx_C);

template <typename T>
int tensor_sum_reference(T alpha, const Tensor<T>& A, const std::string& idx_A,
                         T  beta,       Tensor<T>& B, const std::string& idx_B);

template <typename T>
int tensor_trace_reference(T alpha, const Tensor<T>& A, const std::string& idx_A,
                           T  beta,       Tensor<T>& B, const std::string& idx_B);

template <typename T>
int tensor_replicate_reference(T alpha, const Tensor<T>& A, const std::string& idx_A,
                               T  beta,       Tensor<T>& B, const std::string& idx_B);

template <typename T>
int tensor_transpose_reference(T alpha, const Tensor<T>& A, const std::string& idx_A,
                               T  beta,       Tensor<T>& B, const std::string& idx_B);

template <typename T>
int tensor_dot_reference(const Tensor<T>& A, const std::string& idx_A,
                         const Tensor<T>& B, const std::string& idx_B, T& val);

template <typename T>
int tensor_scale_reference(T alpha, Tensor<T>& A, const std::string& idx_A);

template <typename T>
int tensor_reduce_reference(reduce_t op, const Tensor<T>& A, const std::string& idx_A, T& val, dim_t& idx);

}
}

#endif
