/*

   BLIS
   An object-based framework for developing high-performance BLAS-like
   libraries.

   Copyright (C) 2014, The University of Texas at Austin
   Copyright (C) 2020, Dept. Physics, The University of Tokyo

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:
    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    - Neither the name(s) of the copyright holder(s) nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


*/
#include "blis.h"

// Arm SVE composite instructions.
#include "armsve_asm_macros.h"

// 2vx10 microkernels.
#include "armsve_asm_2vx10.h"

#define PASTEMAC_(ch,op) bli_ ## ch  ## op
#define PASTEMAC(ch,op)  PASTEMAC_(ch,op)

#ifdef _A64FX
#define MEMTAGGING(a,a_next,b,b_next,c) \
  a      = (uint64_t)a      | ((uint64_t)1 << 56); \
  a_next = (uint64_t)a_next | ((uint64_t)1 << 56); \
  b      = (uint64_t)b      | ((uint64_t)2 << 56); \
  b_next = (uint64_t)b_next | ((uint64_t)2 << 56); \
  c      = (uint64_t)c      | ((uint64_t)3 << 56);
#else
#define MEMTAGGING(a,a_next,b,b_next,c)
#endif

#define SVE_TEMPLATED_KERNEL_2VX10(tchar,ctype) \
void PASTEMAC(tchar,gemm_armsve_asm_2vx10_unindexed) \
     ( \
       dim_t               k0,    \
       ctype*     restrict alpha, \
       ctype*     restrict a,     \
       ctype*     restrict b,     \
       ctype*     restrict beta,  \
       ctype*     restrict c, inc_t rs_c0, inc_t cs_c0, \
       auxinfo_t* restrict data,  \
       cntx_t*    restrict cntx   \
     ) \
{ \
  void* a_next = bli_auxinfo_next_a( data ); \
  void* b_next = bli_auxinfo_next_b( data ); \
\
  MEMTAGGING(a,a_next,b,b_next,c) \
\
  /* Typecast local copies of integers in case dim_t and inc_t are a
     different size than is expected by load instructions.
   */ \
  uint64_t k_mker = k0 / 4; \
  uint64_t k_left = k0 % 4; \
  uint64_t rs_c   = rs_c0;  \
  uint64_t cs_c   = cs_c0;  \
\
  __asm__ volatile ( \
" ldr             x0, %[a]                        \n\t" \
" ldr             x1, %[b]                        \n\t" \
" mov             x2, xzr                         \n\t" \
 INC "            x2, ALL, MUL #2                 \n\t" /* Column-skip of A. */ \
" mov             x3, #10                         \n\t" /* Row-skip of B. */ \
"                                                 \n\t" \
" ldr             x5, %[c]                        \n\t" \
" ldr             x6, %[rs_c]                     \n\t" /* Row-skip of C. */ \
" ldr             x7, %[cs_c]                     \n\t" /* Column-skip of C. */ \
"                                                 \n\t" \
" mov             x8, #" SZ "                     \n\t" /* Multiply some address skips by sizeof(double). */ \
" madd            x2, x8, x2, xzr                 \n\t" /* cs_a */ \
" madd            x3, x8, x3, xzr                 \n\t" /* rs_b */ \
" madd            x7, x8, x7, xzr                 \n\t" /* cs_c */ \
" ptrue           p0.b                            \n\t" \
"                                                 \n\t" \
" ldr             x4, %[k_mker]                   \n\t" /* Number of loops. */ \
" ldr             x8, %[k_left]                   \n\t" \
"                                                 \n\t" \
 DT "LOAD_ABC:                                    \n\t" \
" cmp             x4, #0                          \n\t" /* Don't preload if no microkernel there. */ \
" b.eq            " DT "END_CCOL_PRFM             \n\t" \
 LD1R "           z20." DT ", p0/z, [x1]          \n\t" /* Load 8/10 of first B row. */ \
 LD1R "           z21." DT ", p0/z, [x1, " SZ "*1]\n\t" \
 LD1R "           z22." DT ", p0/z, [x1, " SZ "*2]\n\t" \
 LD1R "           z23." DT ", p0/z, [x1, " SZ "*3]\n\t" \
 LD1R "           z24." DT ", p0/z, [x1, " SZ "*4]\n\t" \
 LD1R "           z25." DT ", p0/z, [x1, " SZ "*5]\n\t" \
 LD1R "           z26." DT ", p0/z, [x1, " SZ "*6]\n\t" \
 LD1R "           z27." DT ", p0/z, [x1, " SZ "*7]\n\t" \
"                                                 \n\t" \
GEMM_ACOL_CONTIGUOUS_LOAD(z28,z29,p0,p0,x0) \
"                                                 \n\t" \
 DT "CCOL_PRFM:                                   \n\t" \
" cmp             x6, #1                          \n\t" \
" b.ne            " DT "END_CCOL_PRFM             \n\t" /* Do not prefetch for generic C storage. */ \
" mov             x16, x5                         \n\t" \
" prfm            PLDL1KEEP, [x16]                \n\t" \
" add             x16, x16, x7                    \n\t" \
" prfm            PLDL1KEEP, [x16]                \n\t" \
" add             x16, x16, x7                    \n\t" \
" prfm            PLDL1KEEP, [x16]                \n\t" \
" add             x16, x16, x7                    \n\t" \
" prfm            PLDL1KEEP, [x16]                \n\t" \
" add             x16, x16, x7                    \n\t" \
" prfm            PLDL1KEEP, [x16]                \n\t" \
" add             x16, x16, x7                    \n\t" \
" prfm            PLDL1KEEP, [x16]                \n\t" \
" add             x16, x16, x7                    \n\t" \
" prfm            PLDL1KEEP, [x16]                \n\t" \
" add             x16, x16, x7                    \n\t" \
" prfm            PLDL1KEEP, [x16]                \n\t" \
" add             x16, x16, x7                    \n\t" \
" prfm            PLDL1KEEP, [x16]                \n\t" \
" add             x16, x16, x7                    \n\t" \
" prfm            PLDL1KEEP, [x16]                \n\t" \
 DT "END_CCOL_PRFM:                               \n\t" \
"                                                 \n\t" \
CLEAR_COL20(z0,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15,z16,z17,z18,z19) \
"                                                 \n\t" \
" cmp             x4, #0                          \n\t" /* If no 4-microkernel can be applied. */ \
" b.eq            " DT "K_LEFT_LOOP               \n\t" \
"                                                 \n\t" \
 DT "K_MKER_LOOP:                                 \n\t" \
"                                                 \n\t" \
" add             x0, x0, x2                      \n\t" /* Forward A's address to the next column. */ \
GEMM_ACOL_CONTIGUOUS_LOAD(z30,z31,p0,p0,x0) \
GEMM_2VX10_MKER_LOOP_PLAIN_C_1(z0,z2,z4,z6,z8,z10,z12,z14,z16,z18,z1,z3,z5,z7,z9,z11,z13,z15,z17,z19,p0,z28,z29,z20,z21,z22,z23,z24,z25,z26,z27,x1,x3) \
"                                                 \n\t" \
" add             x0, x0, x2                      \n\t" /* Forward A's address to the next column. */ \
GEMM_ACOL_CONTIGUOUS_LOAD(z28,z29,p0,p0,x0) \
GEMM_2VX10_MKER_LOOP_PLAIN_C_2(z0,z2,z4,z6,z8,z10,z12,z14,z16,z18,z1,z3,z5,z7,z9,z11,z13,z15,z17,z19,p0,z30,z31,z20,z21,z22,z23,z24,z25,z26,z27,x1,x3) \
"                                                 \n\t" \
" add             x0, x0, x2                      \n\t" /* Forward A's address to the next column. */ \
GEMM_ACOL_CONTIGUOUS_LOAD(z30,z31,p0,p0,x0) \
GEMM_2VX10_MKER_LOOP_PLAIN_C_3(z0,z2,z4,z6,z8,z10,z12,z14,z16,z18,z1,z3,z5,z7,z9,z11,z13,z15,z17,z19,p0,z28,z29,z20,z21,z22,z23,z24,z25,z26,z27,x1,x3) \
"                                                 \n\t" \
" subs            x4, x4, #1                      \n\t" /* Decrease counter before final replica. */ \
" b.eq            " DT "FIN_MKER_LOOP             \n\t" /* Branch early to avoid reading excess mem. */ \
"                                                 \n\t" \
" add             x0, x0, x2                      \n\t" /* Forward A's address to the next column. */ \
GEMM_ACOL_CONTIGUOUS_LOAD(z28,z29,p0,p0,x0) \
GEMM_2VX10_MKER_LOOP_PLAIN_C_4(z0,z2,z4,z6,z8,z10,z12,z14,z16,z18,z1,z3,z5,z7,z9,z11,z13,z15,z17,z19,p0,z30,z31,z20,z21,z22,z23,z24,z25,z26,z27,x1,x3) \
" b               " DT "K_MKER_LOOP               \n\t" \
"                                                 \n\t" \
 DT "FIN_MKER_LOOP:                               \n\t" \
GEMM_2VX10_MKER_LOOP_PLAIN_C_4_RESIDUAL(z0,z2,z4,z6,z8,z10,z12,z14,z16,z18,z1,z3,z5,z7,z9,z11,z13,z15,z17,z19,p0,z30,z31,z20,z21,z22,z23,z24,z25,z26,z27,x1,x3) \
" add             x0, x0, x2                      \n\t" /* Forward A to fill the blank. */ \
"                                                 \n\t" \
 DT "K_LEFT_LOOP:                                 \n\t" \
" cmp             x8, #0                          \n\t" /* End of execution. */ \
" b.eq            " DT "WRITE_MEM_PREP            \n\t" \
"                                                 \n\t" \
GEMM_ACOL_CONTIGUOUS_LOAD(z30,z31,p0,p0,x0) \
 LD1R "           z20." DT ", p0/z, [x1]          \n\t" /* Load 8/10 of first B row. */ \
 LD1R "           z21." DT ", p0/z, [x1, " SZ "*1]\n\t" \
 LD1R "           z22." DT ", p0/z, [x1, " SZ "*2]\n\t" \
 LD1R "           z23." DT ", p0/z, [x1, " SZ "*3]\n\t" \
 LD1R "           z24." DT ", p0/z, [x1, " SZ "*4]\n\t" \
 LD1R "           z25." DT ", p0/z, [x1, " SZ "*5]\n\t" \
 LD1R "           z26." DT ", p0/z, [x1, " SZ "*6]\n\t" \
 LD1R "           z27." DT ", p0/z, [x1, " SZ "*7]\n\t" \
 LD1R "           z28." DT ", p0/z, [x1, " SZ "*8]\n\t" \
 LD1R "           z29." DT ", p0/z, [x1, " SZ "*9]\n\t" \
GEMM_FMLA2(z0,z1,p0,z30,z31,z20) \
GEMM_FMLA2(z2,z3,p0,z30,z31,z21) \
GEMM_FMLA2(z4,z5,p0,z30,z31,z22) \
GEMM_FMLA2(z6,z7,p0,z30,z31,z23) \
GEMM_FMLA2(z8,z9,p0,z30,z31,z24) \
GEMM_FMLA2(z10,z11,p0,z30,z31,z25) \
GEMM_FMLA2(z12,z13,p0,z30,z31,z26) \
GEMM_FMLA2(z14,z15,p0,z30,z31,z27) \
GEMM_FMLA2(z16,z17,p0,z30,z31,z28) \
GEMM_FMLA2(z18,z19,p0,z30,z31,z29) \
" add             x0, x0, x2                      \n\t" /* Forward A. */ \
" add             x1, x1, x3                      \n\t" /* Forward B. */ \
" sub             x8, x8, #1                      \n\t" \
" b               " DT "K_LEFT_LOOP               \n\t" /* Next column / row. */ \
"                                                 \n\t" \
 DT "WRITE_MEM_PREP:                              \n\t" \
"                                                 \n\t" \
" ldr             x4, %[alpha]                    \n\t" /* Load alpha & beta (address). */ \
" ldr             x8, %[beta]                     \n\t" \
" ldr             " HF "4, [x4]                   \n\t" /* Load alpha & beta (value). */ \
" ldr             " HF "8, [x8]                   \n\t" \
" dup             z30." DT ", " HF "4             \n\t" /* Broadcast alpha & beta into vectors. */ \
" dup             z31." DT ", " HF "8             \n\t" \
" fmov            " DT "28, #1.0                  \n\t" /* Prepare FP 1.0. */ \
" fmov            " HF "16, " DT "28              \n\t" \
"                                                 \n\t" \
 DT "PREFETCH_ABNEXT:                             \n\t" \
" ldr             x0, %[a_next]                   \n\t" \
" ldr             x1, %[b_next]                   \n\t" \
" prfm            PLDL1STRM, [x0]                 \n\t" \
" prfm            PLDL1STRM, [x0, 256*1]          \n\t" \
/* " prfm            PLDL2KEEP, [x0, 256*2]          \n\t"
   " prfm            PLDL2KEEP, [x0, 256*3]          \n\t"
   " prfm            PLDL2KEEP, [x0, 256*4]          \n\t"
   " prfm            PLDL2KEEP, [x0, 256*5]          \n\t"
   " prfm            PLDL2KEEP, [x0, 256*6]          \n\t"
   " prfm            PLDL2KEEP, [x0, 256*7]          \n\t"
   " prfm            PLDL2KEEP, [x0, 256*8]          \n\t"
   " prfm            PLDL2KEEP, [x0, 256*9]          \n\t"
   " prfm            PLDL2KEEP, [x0, 256*10]         \n\t"
   " prfm            PLDL2KEEP, [x0, 256*11]         \n\t"
   " prfm            PLDL2KEEP, [x0, 256*12]         \n\t"
   " prfm            PLDL2KEEP, [x0, 256*13]         \n\t"
   " prfm            PLDL2KEEP, [x0, 256*14]         \n\t"
   " prfm            PLDL2KEEP, [x0, 256*15]         \n\t"
   */ \
" prfm            PLDL1STRM, [x1]                 \n\t" \
" prfm            PLDL1STRM, [x1, 256*1]          \n\t" \
/* " prfm            PLDL2KEEP, [x1, 256*2]          \n\t"
   " prfm            PLDL2KEEP, [x1, 256*3]          \n\t"
   " prfm            PLDL2KEEP, [x1, 256*4]          \n\t"
   " prfm            PLDL2KEEP, [x1, 256*5]          \n\t"
   " prfm            PLDL2KEEP, [x1, 256*6]          \n\t"
   " prfm            PLDL2KEEP, [x1, 256*7]          \n\t"
   " prfm            PLDL2KEEP, [x1, 256*8]          \n\t"
   " prfm            PLDL2KEEP, [x1, 256*9]          \n\t"
   */ \
"                                                 \n\t" \
" mov             x9, x5                          \n\t" /* C address for loading. */ \
"                                                 \n\t" /* C address for storing is x5 itself. */ \
" cmp             x6, #1                          \n\t" /* Preload first half of C for contiguous case. */ \
" b.ne            " DT "WRITE_MEM                 \n\t" \
GEMM_C_LOAD_UKER_C(z20,z22,z24,z26,z28,z21,z23,z25,z27,z29,p0,p0,x9,x7) \
"                                                 \n\t" \
 DT "WRITE_MEM:                                   \n\t" \
"                                                 \n\t" \
" cmp             " HF "16, " HF "4               \n\t" \
" b.eq            " DT "UNIT_ALPHA                \n\t" \
"                                                 \n\t" \
SCALE_COL20(z0,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15,z16,z17,z18,z19,z30) \
"                                                 \n\t" \
 DT "UNIT_ALPHA:                                  \n\t" \
" cmp             x6, #1                          \n\t" \
" b.ne            " DT "WRITE_MEM_G               \n\t" \
"                                                 \n\t" \
 DT "WRITE_MEM_C:                                 \n\t" /* Available scratch: Z[20-30]. */ \
"                                                 \n\t" /* Here used scratch: Z[20-29]. */ \
/* First half of C is already loaded in this case. */ \
/* GEMM_C_FMAD_LOAD_UKER_C(z20,z22,z24,z26,z28,z21,z23,z25,z27,z29,p0,p0,z0,z2,z4,z6,z8,z1,z3,z5,z7,z9,z31,x9,x7) */ \
" fcmp            " DT "31, #0.0                  \n\t" /* Skip loading for *beta = 0.0. */ \
" b.eq            " DT "BETA_ZERO_C               \n\t" \
GEMM_C_FMLA_UKER(z0,z2,z4,z6,z8,z1,z3,z5,z7,z9,p0,z20,z22,z24,z26,z28,z21,z23,z25,z27,z29,z31) \
GEMM_C_LOAD_UKER_C(z20,z22,z24,z26,z28,z21,z23,z25,z27,z29,p0,p0,x9,x7) \
GEMM_C_FMLA_UKER(z10,z12,z14,z16,z18,z11,z13,z15,z17,z19,p0,z20,z22,z24,z26,z28,z21,z23,z25,z27,z29,z31) \
"                                                 \n\t" \
 DT "BETA_ZERO_C:                                 \n\t" \
GEMM_C_STORE_UKER_C(z0,z2,z4,z6,z8,z1,z3,z5,z7,z9,p0,p0,x5,x7) \
GEMM_C_STORE_UKER_C(z10,z12,z14,z16,z18,z11,z13,z15,z17,z19,p0,p0,x5,x7) \
" b               " DT "END_WRITE_MEM             \n\t" \
"                                                 \n\t" \
 DT "WRITE_MEM_G:                                 \n\t" /* Available scratch: Z[20-30]. */ \
"                                                 \n\t" /* Here used scratch: Z[20-30] - Z30 as index. */ \
" mov             x8, xzr                         \n\t" \
" incb            x8                              \n\t" \
" madd            x8, x8, x6, xzr                 \n\t" /* C-column's logical 1-vector skip. */ \
" index           z30." DT ", " HF "zr, " HF "6   \n\t" /* Skips passed to index is not multiplied by 8. */ \
"                                                 \n\t" \
" fcmp            " DT "31, #0.0                  \n\t" /* Skip loading for *beta = 0.0. */ \
" b.eq            " DT "BETA_ZERO_G               \n\t" \
GEMM_C_LOAD_UKER_G(z20,z22,z24,z26,z28,z21,z23,z25,z27,z29,z30,p0,p0,x9,x7,x8,x16) \
GEMM_C_FMLA_UKER(z0,z2,z4,z6,z8,z1,z3,z5,z7,z9,p0,z20,z22,z24,z26,z28,z21,z23,z25,z27,z29,z31) \
GEMM_C_LOAD_UKER_G(z20,z22,z24,z26,z28,z21,z23,z25,z27,z29,z30,p0,p0,x9,x7,x8,x16) \
GEMM_C_FMLA_UKER(z10,z12,z14,z16,z18,z11,z13,z15,z17,z19,p0,z20,z22,z24,z26,z28,z21,z23,z25,z27,z29,z31) \
"                                                 \n\t" \
 DT  "BETA_ZERO_G:                                \n\t" \
GEMM_C_STORE_UKER_G(z0,z2,z4,z6,z8,z1,z3,z5,z7,z9,z30,p0,p0,x5,x7,x8,x16) \
GEMM_C_STORE_UKER_G(z10,z12,z14,z16,z18,z11,z13,z15,z17,z19,z30,p0,p0,x5,x7,x8,x16) \
"                                                 \n\t" \
 DT "END_WRITE_MEM:                               \n\t" \
" b               " DT "END_EXEC                  \n\t" \
"                                                 \n\t" \
 DT "END_ERROR:                                   \n\t" \
" mov             x0, #1                          \n\t" /* Return error. */ \
 DT "END_EXEC:                                    \n\t" \
" mov             x0, #0                          \n\t" /* Return normal. */ \
: \
: [a]      "m" (a),      \
  [b]      "m" (b),      \
  [c]      "m" (c),      \
  [rs_c]   "m" (rs_c),   \
  [cs_c]   "m" (cs_c),   \
  [k_mker] "m" (k_mker), \
  [k_left] "m" (k_left), \
  [alpha]  "m" (alpha),  \
  [beta]   "m" (beta),   \
  [a_next] "m" (a_next), \
  [b_next] "m" (b_next)  \
: "x0","x1","x2","x3","x4","x5","x6","x7","x8",  \
  "x9","x16",                                    \
  "z0","z1","z2","z3","z4","z5","z6","z7",       \
  "z8","z9","z10","z11","z12","z13","z14","z15", \
  "z16","z17","z18","z19", \
  "z20","z21","z22","z23", \
  "z24","z25","z26","z27", \
  "z28","z29","z30","z31"  \
   ); \
}

// Instantiate (double).
#define HF    "x"
#define DT    "d"
#define LD1   "ld1d"
#define ST1   "st1d"
#define LD1R  "ld1rd"
#define PRFG  "prfd"
#define INC   "incd"
#define SZ    "8"
#define OFFS  "lsl #3"
SVE_TEMPLATED_KERNEL_2VX10(d,double)
#undef HF
#undef DT
#undef LD1
#undef ST1
#undef LD1R
#undef PRFG
#undef INC
#undef SZ
#undef OFFS

// Instantiate (float).
#define HF    "w"
#define DT    "s"
#define LD1   "ld1w"
#define ST1   "st1w"
#define LD1R  "ld1rw"
#define PRFG  "prfw"
#define INC   "incw"
#define SZ    "4"
#define OFFS  "uxtw #2"
SVE_TEMPLATED_KERNEL_2VX10(s,float)

