/* 

   Additional macros for AMX instructions

   TODO: Clean up. Merge duplicates.

*/
#define AMX_MEM(OP, ADDR, REG) \
  AMX_## OP( (uint64_t)(ADDR) | ((uint64_t)(REG) << 56) )

#define AMX_EXTRX_REGALIGNED(ZREG, XREG) \
  AMX_EXTRX( ((uint64_t)(ZREG) << 20) | ((uint64_t)(XREG) << 16) )

#define AMX_EXTRY64_REGALIGNED(ZREG, YREG) \
  AMX_EXTRY( ((uint64_t)(ZREG) << 20) | ((uint64_t)(YREG) << 6) )

#define AMX_FMA32_COMMON(PADROWS, PADCOLS, ZREGS) \
  AMX_FMA32( (uint64_t)(PADCOLS) | ((uint64_t)(PADROWS) << 10) | ((uint64_t)(ZREGS)<< 20) )

#define AMX_FMA32_COMMON_REGALIGNED(XREG, YREG, ZREGS) \
  AMX_FMA32_COMMON( ((uint64_t)(XREG) << 6), ((uint64_t)(YREG) << 6), (ZREGS) )

#define AMX_FMUL32_COMMON(PADROWS, PADCOLS, ZREGS) \
  AMX_FMA32(  (uint64_t)(PADCOLS) | ((uint64_t)(PADROWS) << 10) | ((uint64_t)(ZREGS)<< 20) \
           | ((uint64_t) 1 << 27) )

#define AMX_FMUL32_COMMON_REGALIGNED(XREG, YREG, ZREGS) \
  AMX_FMUL32_COMMON( ((uint64_t)(XREG) << 6), ((uint64_t)(YREG) << 6), (ZREGS) )

#define AMX_FMA32_SEL(SELCOLS, SELROWS, PADROWS, PADCOLS, ZREGS) \
  AMX_FMA32(  (uint64_t)(PADCOLS) \
           | ((uint64_t)(PADROWS) << 10) | ((uint64_t)(ZREGS) << 20) \
           | ((uint64_t)(SELCOLS) << 32) | ((uint64_t)(SELROWS) << 41) )

#define AMX_FMA32_SELCOL_REGALIGNED(COLIDX, XREG, YREG, ZREGS) \
  AMX_FMA32_SEL( ((uint64_t)0x20 | COLIDX), 0, \
                 ((uint64_t)(XREG) << 6), ((uint64_t)(YREG) << 6), ZREGS )

#define AMX_FMUL32_SEL(SELCOLS, SELROWS, PADROWS, PADCOLS, ZREGS) \
  AMX_FMA32(  (uint64_t)(PADCOLS) \
           | ((uint64_t)(PADROWS) << 10) | ((uint64_t)(ZREGS) << 20) \
           | ((uint64_t)(SELCOLS) << 32) | ((uint64_t)(SELROWS) << 41) \
           | ((uint64_t) 1 << 27) )

#define AMX_FMUL32_SELCOL_REGALIGNED(COLIDX, XREG, YREG, ZREGS) \
  AMX_FMUL32_SEL( ((uint64_t)0x20 | COLIDX), 0, \
                  ((uint64_t)(XREG) << 6), ((uint64_t)(YREG) << 6), ZREGS )

#define AMX_FMUL32_SELROW_REGALIGNED(ROWIDX, XREG, YREG, ZREGS) \
  AMX_FMUL32_SEL( 0, ((uint64_t)0x20 | ROWIDX), \
                  ((uint64_t)(XREG) << 6), ((uint64_t)(YREG) << 6), ZREGS )

#define AMX_FMA64_COMMON(PADROWS, PADCOLS, ZREGS) \
  AMX_FMA64( (uint64_t)(PADCOLS) | ((uint64_t)(PADROWS) << 10) | ((uint64_t)(ZREGS)<< 20) )

#define AMX_FMA64_COMMON_REGALIGNED(XREG, YREG, ZREGS) \
  AMX_FMA64_COMMON( ((uint64_t)XREG << 6), ((uint64_t)YREG << 6), ZREGS )

#define AMX_FMUL64_COMMON(PADROWS, PADCOLS, ZREGS) \
  AMX_FMA64(  (uint64_t)(PADCOLS) | ((uint64_t)(PADROWS) << 10) | ((uint64_t)(ZREGS)<< 20) \
           | ((uint64_t) 1 << 27))

#define AMX_FMUL64_COMMON_REGALIGNED(XREG, YREG, ZREGS) \
  AMX_FMUL64_COMMON( ((uint64_t)(XREG) << 6), ((uint64_t)(YREG) << 6), (ZREGS) )

#define AMX_FMA64_SEL(SELCOLS, SELROWS, PADROWS, PADCOLS, ZREGS) \
  AMX_FMA64(  (uint64_t)(PADCOLS) \
           | ((uint64_t)(PADROWS) << 10) | ((uint64_t)(ZREGS)<< 20) \
           | ((uint64_t)(SELCOLS) << 32) | ((uint64_t)(SELROWS) << 41) )

#define AMX_FMA64_SELCOL_REGALIGNED(COLIDX, XREG, YREG, ZREGS) \
  AMX_FMA64_SEL( ((uint64_t)0x20 | COLIDX), 0, \
                 ((uint64_t)(XREG) << 6), ((uint64_t)(YREG) << 6), ZREGS )

#define AMX_FMA64_SELROW_REGALIGNED(ROWIDX, XREG, YREG, ZREGS) \
  AMX_FMA64_SEL( 0, ((uint64_t)0x20 | ROWIDX), \
                 ((uint64_t)(XREG) << 6), ((uint64_t)(YREG) << 6), ZREGS )

#define AMX_FMUL64_SEL(SELCOLS, SELROWS, PADROWS, PADCOLS, ZREGS) \
  AMX_FMA64(  (uint64_t)(PADCOLS) \
           | ((uint64_t)(PADROWS) << 10) | ((uint64_t)(ZREGS)<< 20) \
           | ((uint64_t)(SELCOLS) << 32) | ((uint64_t)(SELROWS) << 41) \
           | ((uint64_t) 1 << 27) )

#define AMX_FMUL64_SELCOL_REGALIGNED(COLIDX, XREG, YREG, ZREGS) \
  AMX_FMUL64_SEL( ((uint64_t)0x20 | COLIDX), 0, \
                  ((uint64_t)(XREG) << 6), ((uint64_t)(YREG) << 6), ZREGS )

#define AMX_FMUL64_SELROW_REGALIGNED(ROWIDX, XREG, YREG, ZREGS) \
  AMX_FMUL64_SEL( 0, ((uint64_t)0x20 | ROWIDX), \
                  ((uint64_t)(XREG) << 6), ((uint64_t)(YREG) << 6), ZREGS )


