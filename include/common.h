#include <cstddef>

#include <gmpxx.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <memory>
#include <nfl.hpp>

#include "bench.h"
#include <sys/random.h>

using namespace std;

#ifndef COMMON_H
#define COMMON_H

/* Parameter v in the commitment  scheme (laximum l1-norm of challs). */
#define NONZERO     36
/* Security level to attain. */
#define LEVEL       128
/* The \infty-norm bound of certain elements. */
#define BETA 	    1
/* Width k of the commitment matrix. */
#define WIDTH 	    4
/* Height n of the commitment matrix. */
#define HEIGHT 	    1
/* Dimension of the committed messages. */
#ifndef SIZE
#define SIZE        2
#endif
/* Small modulus. */
#define PRIMEP      3
/* Degree of the irreducible polynomial. */
#define DEGREE      4096
/* Sigma for the commitment gaussian distribution. */
#define SIGMA_C     (1u << 12)
/* Parties that run the distributed decryption protocol. */
#define PARTIES     5
/* Sigma for the OR proof. */
#define SIGMA_O     (1555)
/* Bound for the OR proof. */
#define BOUND_O     (18662)

#if USERS == 1024
    /* Large modulus. */
    #define PRIMEQ      "316912650057057350374175870977"
    /* Parameters for encoding messages modulo q. */
    #define DELTA       562949953421312
    #define DELTA_R     (DELTA-139266)
    /* Sigma for the boundness proof. */
    #define SIGMA_B     (11906898151735774419091456.0)
    /* Norm bound for boundness proof. */
    #define BOUND_B     "144379862274704472782405632"
    /* Bound for Distributed Decryption. */
    #define BOUND_D     "1155038898197635782259245056"
#elif USERS == 10000
    /* Large modulus. */
    #define PRIMEQ      "2658455991569831745807614120560865281"
    /* Parameter for encoding messages modulo q. */
    #define DELTA       2251799813685248
    #define DELTA_R     (DELTA-565252)
    /* Sigma for the boundness proof. */
    #define SIGMA_B     (355504921678021018857242624.0)
    /* Norm bound for boundness proof. */
    #define BOUND_B     "11649185273545392745914126303232"
    /* Bound for Distributed Decryption. */
    #define BOUND_D     "818235810536127948211699224848367616"
#endif

namespace params {
    using poly_p = nfl::poly_from_modulus<uint32_t, DEGREE, 30>;
    using poly_q = nfl::poly_from_modulus<uint64_t, DEGREE, 124>;
    using poly_2q = nfl::poly_from_modulus<uint64_t, 2*DEGREE, 124>;
}

/*============================================================================*/
/* Type definitions                                                           */
/*============================================================================*/

/* Class that represents a commitment key pair. */
class comkey_t {
    public:
       params::poly_q A1[HEIGHT][WIDTH - HEIGHT];
       params::poly_q A2[SIZE][WIDTH];
};

/* Class that represents a commitment in CRT representation. */
class commit_t {
    public:
      params::poly_q c1;
      vector<params::poly_q> c2;
};

/* Class that represents a BGV key pair. */
class bgvkey_t {
    public:
       params::poly_q a;
       params::poly_q b;
};

class bgvenc_t {
    public:
       params::poly_q u;
       params::poly_q v;
};

/* Class that represents a GHL key pair. */
class ghlkey_t {
    public:
       params::poly_2q a;
       params::poly_2q b;
};

class ghlenc_t {
    public:
       params::poly_2q u;
       params::poly_2q v;
};

#include "falcon.h"

typedef struct {
	unsigned logn;
	shake256_context rng;
	uint8_t *tmp;
	size_t tmp_len;
	uint8_t *pk;
	uint8_t *sk;
	uint8_t *esk;
	uint8_t *sig;
	size_t sig_len;
	uint8_t *sigct;
	size_t sigct_len;
} falcon_t;

#include "util.hpp"

void bdlop_sample_rand(vector<params::poly_q>& r);
void bdlop_sample_chal(params::poly_q& f);
bool bdlop_test_norm(params::poly_q r, uint64_t sigma_sqr);
void bdlop_commit(commit_t& com, vector<params::poly_q> m, comkey_t& key, vector<params::poly_q> r);
int  bdlop_open(commit_t& com, vector<params::poly_q> m, comkey_t& key, vector<params::poly_q> r, params::poly_q& f);
void bdlop_keygen(comkey_t& key);

void bgv_sample_message(params::poly_p& r);
void bgv_sample_short(params::poly_q& r);
void bgv_keygen(bgvkey_t& pk, params::poly_q& sk);
void bgv_encrypt(bgvenc_t &c, bgvkey_t& pk, params::poly_p& m);
void bgv_decrypt(params::poly_p& m, bgvenc_t& c, params::poly_q & sk);

void ghl_sample_message(params::poly_q & m);
void ghl_encode(params::poly_2q & _m, params::poly_q & m);
void ghl_decode(params::poly_q & _m, params::poly_2q & m);
void ghl_keygen(ghlkey_t & pk, params::poly_2q & sk);
void ghl_keyshare(params::poly_2q s[], size_t shares, params::poly_2q & sk);
void ghl_encrypt(ghlenc_t & c, ghlkey_t & pk, params::poly_2q & m);
void ghl_encrypt(ghlenc_t & c, ghlkey_t & pk, params::poly_2q & m, params::poly_2q & r);
void ghl_decrypt(params::poly_2q & m, ghlenc_t & c, params::poly_2q & sk);
void ghl_add(ghlenc_t & c, ghlenc_t & d, ghlenc_t & e);
void ghl_distdec(params::poly_2q & tj, ghlenc_t & c, params::poly_2q & sj);
void ghl_comb(params::poly_2q & m, ghlenc_t & c, params::poly_2q t[], size_t shares);

#endif /* COMMON_H */
