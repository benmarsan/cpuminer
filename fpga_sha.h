#ifndef FPGA_SHA_H
#define FPGA_SHA_H

#include <inttypes.h>

extern void fpga_setup();

extern int fpga_scanhash_sha256d(int thr_id, uint32_t *pdata,
                            const uint32_t *ptarget, uint32_t max_nonce, unsigned long *hashes_done);

#endif //FPGA_SHA_H
