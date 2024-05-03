/**
 *
 */

#include <inttypes.h>
#include <string.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/mman.h>
#include <stdint.h>


#include "fpga_sha.h"
#include "miner.h"

#define MAP_SIZE 4096UL
#define MAP_MASK (MAP_SIZE - 1)
#define MINE_CORE 0x80000000 // Add the correct physical address here

/* Elementary functions used by SHA256 */
#define Ch(x, y, z)     ((x & (y ^ z)) ^ z)
#define Maj(x, y, z)    ((x & (y | z)) | (y & z))
#define ROTR(x, n)      ((x >> n) | (x << (32 - n)))
#define S0(x)           (ROTR(x, 2) ^ ROTR(x, 13) ^ ROTR(x, 22))
#define S1(x)           (ROTR(x, 6) ^ ROTR(x, 11) ^ ROTR(x, 25))
#define s0(x)           (ROTR(x, 7) ^ ROTR(x, 18) ^ (x >> 3))
#define s1(x)           (ROTR(x, 17) ^ ROTR(x, 19) ^ (x >> 10))

/* SHA256 round function */
#define RND(a, b, c, d, e, f, g, h, k) \
	do { \
		t0 = h + S1(e) + Ch(e, f, g) + k; \
		t1 = S0(a) + Maj(a, b, c); \
		d += t0; \
		h  = t0 + t1; \
	} while (0)

/* Adjusted round function for rotating state */
#define RNDr(S, W, i) \
	RND(S[(64 - i) % 8], S[(65 - i) % 8], \
	    S[(66 - i) % 8], S[(67 - i) % 8], \
	    S[(68 - i) % 8], S[(69 - i) % 8], \
	    S[(70 - i) % 8], S[(71 - i) % 8], \
	    W[i] + sha256_k[i])

static const uint32_t sha256d_hash1[16] = {
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x80000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000100
};

static const uint32_t sha256_h[8] = {
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
        0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
};

static const uint32_t sha256_k[64] = {
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
        0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
        0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
        0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
        0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
        0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
        0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
        0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
};

static inline void sha256d_preextend(uint32_t *W)
{
    W[16] = s1(W[14]) + W[ 9] + s0(W[ 1]) + W[ 0];
    W[17] = s1(W[15]) + W[10] + s0(W[ 2]) + W[ 1];
    W[18] = s1(W[16]) + W[11]             + W[ 2];
    W[19] = s1(W[17]) + W[12] + s0(W[ 4]);
    W[20] =             W[13] + s0(W[ 5]) + W[ 4];
    W[21] =             W[14] + s0(W[ 6]) + W[ 5];
    W[22] =             W[15] + s0(W[ 7]) + W[ 6];
    W[23] =             W[16] + s0(W[ 8]) + W[ 7];
    W[24] =             W[17] + s0(W[ 9]) + W[ 8];
    W[25] =                     s0(W[10]) + W[ 9];
    W[26] =                     s0(W[11]) + W[10];
    W[27] =                     s0(W[12]) + W[11];
    W[28] =                     s0(W[13]) + W[12];
    W[29] =                     s0(W[14]) + W[13];
    W[30] =                     s0(W[15]) + W[14];
    W[31] =                     s0(W[16]) + W[15];
}

static inline void sha256d_prehash(uint32_t *S, const uint32_t *W)
{
    uint32_t t0, t1;
    RNDr(S, W, 0);
    RNDr(S, W, 1);
    RNDr(S, W, 2);
}

static inline void sha256d_ms(uint32_t *hash, uint32_t *W,
                              const uint32_t *midstate, const uint32_t *prehash)
{
    uint32_t S[64];
    uint32_t t0, t1;
    int i;

    S[18] = W[18];
    S[19] = W[19];
    S[20] = W[20];
    S[22] = W[22];
    S[23] = W[23];
    S[24] = W[24];
    S[30] = W[30];
    S[31] = W[31];

    W[18] += s0(W[3]);
    W[19] += W[3];
    W[20] += s1(W[18]);
    W[21]  = s1(W[19]);
    W[22] += s1(W[20]);
    W[23] += s1(W[21]);
    W[24] += s1(W[22]);
    W[25]  = s1(W[23]) + W[18];
    W[26]  = s1(W[24]) + W[19];
    W[27]  = s1(W[25]) + W[20];
    W[28]  = s1(W[26]) + W[21];
    W[29]  = s1(W[27]) + W[22];
    W[30] += s1(W[28]) + W[23];
    W[31] += s1(W[29]) + W[24];
    for (i = 32; i < 64; i += 2) {
        W[i]   = s1(W[i - 2]) + W[i - 7] + s0(W[i - 15]) + W[i - 16];
        W[i+1] = s1(W[i - 1]) + W[i - 6] + s0(W[i - 14]) + W[i - 15];
    }

    memcpy(S, prehash, 32);

    RNDr(S, W,  3);
    RNDr(S, W,  4);
    RNDr(S, W,  5);
    RNDr(S, W,  6);
    RNDr(S, W,  7);
    RNDr(S, W,  8);
    RNDr(S, W,  9);
    RNDr(S, W, 10);
    RNDr(S, W, 11);
    RNDr(S, W, 12);
    RNDr(S, W, 13);
    RNDr(S, W, 14);
    RNDr(S, W, 15);
    RNDr(S, W, 16);
    RNDr(S, W, 17);
    RNDr(S, W, 18);
    RNDr(S, W, 19);
    RNDr(S, W, 20);
    RNDr(S, W, 21);
    RNDr(S, W, 22);
    RNDr(S, W, 23);
    RNDr(S, W, 24);
    RNDr(S, W, 25);
    RNDr(S, W, 26);
    RNDr(S, W, 27);
    RNDr(S, W, 28);
    RNDr(S, W, 29);
    RNDr(S, W, 30);
    RNDr(S, W, 31);
    RNDr(S, W, 32);
    RNDr(S, W, 33);
    RNDr(S, W, 34);
    RNDr(S, W, 35);
    RNDr(S, W, 36);
    RNDr(S, W, 37);
    RNDr(S, W, 38);
    RNDr(S, W, 39);
    RNDr(S, W, 40);
    RNDr(S, W, 41);
    RNDr(S, W, 42);
    RNDr(S, W, 43);
    RNDr(S, W, 44);
    RNDr(S, W, 45);
    RNDr(S, W, 46);
    RNDr(S, W, 47);
    RNDr(S, W, 48);
    RNDr(S, W, 49);
    RNDr(S, W, 50);
    RNDr(S, W, 51);
    RNDr(S, W, 52);
    RNDr(S, W, 53);
    RNDr(S, W, 54);
    RNDr(S, W, 55);
    RNDr(S, W, 56);
    RNDr(S, W, 57);
    RNDr(S, W, 58);
    RNDr(S, W, 59);
    RNDr(S, W, 60);
    RNDr(S, W, 61);
    RNDr(S, W, 62);
    RNDr(S, W, 63);

    for (i = 0; i < 8; i++)
        S[i] += midstate[i];

    W[18] = S[18];
    W[19] = S[19];
    W[20] = S[20];
    W[22] = S[22];
    W[23] = S[23];
    W[24] = S[24];
    W[30] = S[30];
    W[31] = S[31];

    memcpy(S + 8, sha256d_hash1 + 8, 32);
    S[16] = s1(sha256d_hash1[14]) + sha256d_hash1[ 9] + s0(S[ 1]) + S[ 0];
    S[17] = s1(sha256d_hash1[15]) + sha256d_hash1[10] + s0(S[ 2]) + S[ 1];
    S[18] = s1(S[16]) + sha256d_hash1[11] + s0(S[ 3]) + S[ 2];
    S[19] = s1(S[17]) + sha256d_hash1[12] + s0(S[ 4]) + S[ 3];
    S[20] = s1(S[18]) + sha256d_hash1[13] + s0(S[ 5]) + S[ 4];
    S[21] = s1(S[19]) + sha256d_hash1[14] + s0(S[ 6]) + S[ 5];
    S[22] = s1(S[20]) + sha256d_hash1[15] + s0(S[ 7]) + S[ 6];
    S[23] = s1(S[21]) + S[16] + s0(sha256d_hash1[ 8]) + S[ 7];
    S[24] = s1(S[22]) + S[17] + s0(sha256d_hash1[ 9]) + sha256d_hash1[ 8];
    S[25] = s1(S[23]) + S[18] + s0(sha256d_hash1[10]) + sha256d_hash1[ 9];
    S[26] = s1(S[24]) + S[19] + s0(sha256d_hash1[11]) + sha256d_hash1[10];
    S[27] = s1(S[25]) + S[20] + s0(sha256d_hash1[12]) + sha256d_hash1[11];
    S[28] = s1(S[26]) + S[21] + s0(sha256d_hash1[13]) + sha256d_hash1[12];
    S[29] = s1(S[27]) + S[22] + s0(sha256d_hash1[14]) + sha256d_hash1[13];
    S[30] = s1(S[28]) + S[23] + s0(sha256d_hash1[15]) + sha256d_hash1[14];
    S[31] = s1(S[29]) + S[24] + s0(S[16])             + sha256d_hash1[15];
    for (i = 32; i < 60; i += 2) {
        S[i]   = s1(S[i - 2]) + S[i - 7] + s0(S[i - 15]) + S[i - 16];
        S[i+1] = s1(S[i - 1]) + S[i - 6] + s0(S[i - 14]) + S[i - 15];
    }
    S[60] = s1(S[58]) + S[53] + s0(S[45]) + S[44];

    sha256_init(hash);

    RNDr(hash, S,  0);
    RNDr(hash, S,  1);
    RNDr(hash, S,  2);
    RNDr(hash, S,  3);
    RNDr(hash, S,  4);
    RNDr(hash, S,  5);
    RNDr(hash, S,  6);
    RNDr(hash, S,  7);
    RNDr(hash, S,  8);
    RNDr(hash, S,  9);
    RNDr(hash, S, 10);
    RNDr(hash, S, 11);
    RNDr(hash, S, 12);
    RNDr(hash, S, 13);
    RNDr(hash, S, 14);
    RNDr(hash, S, 15);
    RNDr(hash, S, 16);
    RNDr(hash, S, 17);
    RNDr(hash, S, 18);
    RNDr(hash, S, 19);
    RNDr(hash, S, 20);
    RNDr(hash, S, 21);
    RNDr(hash, S, 22);
    RNDr(hash, S, 23);
    RNDr(hash, S, 24);
    RNDr(hash, S, 25);
    RNDr(hash, S, 26);
    RNDr(hash, S, 27);
    RNDr(hash, S, 28);
    RNDr(hash, S, 29);
    RNDr(hash, S, 30);
    RNDr(hash, S, 31);
    RNDr(hash, S, 32);
    RNDr(hash, S, 33);
    RNDr(hash, S, 34);
    RNDr(hash, S, 35);
    RNDr(hash, S, 36);
    RNDr(hash, S, 37);
    RNDr(hash, S, 38);
    RNDr(hash, S, 39);
    RNDr(hash, S, 40);
    RNDr(hash, S, 41);
    RNDr(hash, S, 42);
    RNDr(hash, S, 43);
    RNDr(hash, S, 44);
    RNDr(hash, S, 45);
    RNDr(hash, S, 46);
    RNDr(hash, S, 47);
    RNDr(hash, S, 48);
    RNDr(hash, S, 49);
    RNDr(hash, S, 50);
    RNDr(hash, S, 51);
    RNDr(hash, S, 52);
    RNDr(hash, S, 53);
    RNDr(hash, S, 54);
    RNDr(hash, S, 55);
    RNDr(hash, S, 56);

    hash[2] += hash[6] + S1(hash[3]) + Ch(hash[3], hash[4], hash[5])
               + S[57] + sha256_k[57];
    hash[1] += hash[5] + S1(hash[2]) + Ch(hash[2], hash[3], hash[4])
               + S[58] + sha256_k[58];
    hash[0] += hash[4] + S1(hash[1]) + Ch(hash[1], hash[2], hash[3])
               + S[59] + sha256_k[59];
    hash[7] += hash[3] + S1(hash[0]) + Ch(hash[0], hash[1], hash[2])
               + S[60] + sha256_k[60]
               + sha256_h[7];
}

static void sha256d_80_swap(uint32_t *hash, const uint32_t *data)
{
    uint32_t S[16];
    int i;

    sha256_init(S);
    sha256_transform(S, data, 0);
    sha256_transform(S, data + 16, 0);
    memcpy(S + 8, sha256d_hash1 + 8, 32);
    sha256_init(hash);
    sha256_transform(hash, S, 0);
    for (i = 0; i < 8; i++)
        hash[i] = swab32(hash[i]);
}

// Memory mapped device registers
volatile uint32_t* MappedBase;
volatile uint32_t* Control;
volatile uint32_t* Header;
volatile uint32_t* Nonce;
volatile uint32_t* Difficulty;

void fpga_setup()
{
    int fd = open("/dev/mem", O_RDWR|O_SYNC);
    if(fd == -1) {
        perror("Error opening /dev/mem");
        return;
    }

    MappedBase = (volatile uint32_t*) mmap(NULL, MAP_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, fd, MINE_CORE);
    if (MappedBase == MAP_FAILED) {
        perror("Error mapping memory");
        close(fd);
        return;
    }

    Control = MappedBase;
    Header  = MappedBase + 1;
    Nonce   = MappedBase + 21;
    Difficulty = MappedBase + 22;

    *Control = 0;
    // Reset device
    *Control |= (1 << 4); // Set bit 4
    *Control &= ~(1 << 4); // Reset bit 4
}

int fpga_scanhash_sha256d(int thr_id, uint32_t *pdata, const uint32_t *ptarget,
                     uint32_t max_nonce, unsigned long *hashes_done)
{
    uint32_t data[64] __attribute__((aligned(128)));
    uint32_t hash[8] __attribute__((aligned(32)));
    uint32_t midstate[8] __attribute__((aligned(32)));
    uint32_t prehash[8] __attribute__((aligned(32)));
    uint32_t n = pdata[19] - 1;
    const uint32_t first_nonce = pdata[19];
    const uint32_t Htarg = ptarget[7];

    printf("Input header data (bytes): ");
    for(int i = 0; i < 80; i++) {
        printf("%02x ", ((uint8_t *)pdata)[i]);
    }
    printf("\r\n");

    printf("Input header data (BEHex): ");
    for(int i = 19; i >= 0; i--) {
        printf("%08x", pdata[i]);
    }
    printf("\r\n");

    printf("Input target: ");
    for(int i = 0; i < 8*4; i++) {
        printf("%02x ", ((uint8_t *)ptarget)[i]);
    }
    printf("\r\n");

    printf("Input target (BEHex): ");
    for(int i = 7; i >= 0; i--) {
        printf("%08x", ptarget[i]); // SWAP BYTES?
    }
    printf("\r\n");

    uint32_t target_be[8];
    char target_str[65];

    for (int i = 0; i < 8; i++) {
        be32enc(target_be + i, ptarget[7 - i]);
    }
    bin2hex(target_str, (unsigned char *)target_be, 32);
    printf("fpga_sha target: %s\r\n", target_str);

    *Control = 0;
    // Reset device
    *Control |= (1 << 4); // Set bit 4
    *Control &= ~(1 << 4); // Reset bit 4

    // Copy expanded difficulty to registers
    printf("Setting difficulty target to: 0x");
    for (int i = 0; i < 8; i++) {
        //printf("%i, %i\n", expanded_difficulty[i], i);
        Difficulty[i] = ptarget[7-i];
        printf("%08X", ptarget[7-i]);
    }
    printf("\r\n");

    // Handle block header
    for (size_t i = 0; i < 80/4; i++) {
        // Header in from pdata is already byte-swapped, we don't need to swap again
        Header[i] = pdata[i];
    }

    // set bit 0 to start mining
    *Control |= (1 << 0);

    do {
        // Poll to check if accelerator has finished
        if (*Control & (1 << 8)) {
            uint32_t nonce_from_fpga = *Nonce;
            printf("Found valid nonce: 0x%08x\r\n", nonce_from_fpga);
            *hashes_done = (nonce_from_fpga % (286331153)) * 15;

            pdata[19] = nonce_from_fpga;
            printf("Nonce from fpga: 0x%08x\r\n", nonce_from_fpga);

            printf("Resulting header bytes: ");
            for(int i = 0; i < 80; i++) {
                printf("%02x ", ((uint8_t *)pdata)[i]);
            }
            printf("\r\n");

            printf("maybe good header (bytes): ");
            for(int i = 0; i < 80; i++) {
                printf("%02x ", ((uint8_t *)pdata)[i]);
            }
            printf("\r\n");

            // Check if hash is indeed valid
            memcpy(data, pdata + 16, 64);
            sha256d_preextend(data);
            sha256_init(midstate);
            sha256_transform(midstate, pdata, 0);
            memcpy(prehash, midstate, 32);
            sha256d_prehash(prehash, pdata + 16);
            data[3] = nonce_from_fpga;
            sha256d_ms(hash, data, midstate, prehash);
            pdata[19] = data[3];
            sha256d_80_swap(hash, pdata);
            if (fulltest(hash, ptarget)) {
                printf("Hash from FPGA meets target!\r\n");
                return 1;
            }
            printf("Error! Accelerator reported false positive!\r\n");
            return 0;
        } else if (*Control & (1 << 12)) {
            printf("Exhausted possible nonce values\r\n");
            *hashes_done = 4294967295; // 2^32-1
            return 0;
        }

        usleep(1000);
    } while (!work_restart[thr_id].restart);

    // Stop mining
    *Control &= ~(1 << 0);

    *hashes_done = (*Nonce % (286331153)) * 15;
    pdata[19] = *Nonce;
    return 0;
}
