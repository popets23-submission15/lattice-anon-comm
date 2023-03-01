#include "blake3.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "common.h"
#include "test.h"
#include "bench.h"
#include "assert.h"
#include "sample_z_small.h"

void pior_prehash(uint8_t hash[BLAKE3_OUT_LEN],
		params::poly_q A0[SIZE][SIZE + 1], params::poly_q A1[SIZE],
		params::poly_q t1) {
	blake3_hasher hasher;
	blake3_hasher_init(&hasher);

	for (size_t i = 0; i < SIZE; i++) {
		for (size_t j = 0; j < SIZE + 1; j++) {
			blake3_hasher_update(&hasher, A0[i][j].data(), 16 * DEGREE);
		}
		blake3_hasher_update(&hasher, A1[i].data(), 16 * DEGREE);
	}
	blake3_hasher_update(&hasher, t1.data(), 16 * DEGREE);

	blake3_hasher_finalize(&hasher, hash, BLAKE3_OUT_LEN);
}

void pior_hash(uint8_t hash[BLAKE3_OUT_LEN], uint8_t prehash[BLAKE3_OUT_LEN],
		params::poly_q w0[SIZE], params::poly_q w1, params::poly_q t0[SIZE]) {
	blake3_hasher hasher;

	blake3_hasher_init(&hasher);

	for (size_t i = 0; i < SIZE; i++) {
		blake3_hasher_update(&hasher, w0[i].data(), 16 * DEGREE);
		blake3_hasher_update(&hasher, t0[i].data(), 16 * DEGREE);
	}
	blake3_hasher_update(&hasher, w1.data(), 16 * DEGREE);
	blake3_hasher_update(&hasher, prehash, BLAKE3_OUT_LEN);

	blake3_hasher_finalize(&hasher, hash, BLAKE3_OUT_LEN);
}

/**
 * Test if the l2-norm is within bounds (4 * sigma * sqrt(N)).
 *
 * @param[in] r 			- the polynomial to compute the l2-norm.
 * @return the computed norm.
 */
bool pior_test_norm(params::poly_q r, uint64_t sigma_sqr) {
	array < mpz_t, params::poly_q::degree > coeffs;
	mpz_t norm, qDivBy2, tmp;

	/// Constructors
	mpz_inits(norm, qDivBy2, tmp, nullptr);
	for (size_t i = 0; i < params::poly_q::degree; i++) {
		mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
	}

	r.poly2mpz(coeffs);
	mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);
	mpz_set_ui(norm, 0);
	for (size_t i = 0; i < params::poly_q::degree; i++) {
		util::center(coeffs[i], coeffs[i],
				params::poly_q::moduli_product(), qDivBy2);
		mpz_mul(tmp, coeffs[i], coeffs[i]);
		mpz_add(norm, norm, tmp);
	}

	// Compare to (2 * sigma * sqrt(2N))^2 = 8 * sigma^2 * N.
	uint64_t bound = 2 * sigma_sqr * params::poly_q::degree;
	int result = mpz_cmp_ui(norm, bound) <= 0;

	mpz_clears(norm, qDivBy2, tmp, nullptr);
	for (size_t i = 0; i < params::poly_q::degree; i++) {
		mpz_clear(coeffs[i]);
	}

	return result;
}

// Sample a challenge.
void pior_sample_chall(params::poly_q & f) {
	f = nfl::ZO_dist();
	f.ntt_pow_phi();
}

void pior_prove(uint8_t h[2][BLAKE3_OUT_LEN], uint8_t prehash[BLAKE3_OUT_LEN],
		params::poly_q w0[SIZE], params::poly_q & w1,
		params::poly_q z0[SIZE + 1], params::poly_q z1[SIZE],
		params::poly_q s0[SIZE + 1], params::poly_q s1[SIZE],
		params::poly_q A0[SIZE][SIZE + 1], params::poly_q A1[SIZE],
		params::poly_q t0[SIZE], params::poly_q t1, int b) {
	std::array < mpz_t, params::poly_q::degree > coeffs;
	params::poly_q c[2], y0[SIZE + 1], y1[SIZE];

	for (size_t i = 0; i < params::poly_q::degree; i++) {
		mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
	}

	if (b == 0) {
		/* Prover samples z_1-b,y_b from Gaussian. */
		for (size_t i = 0; i < SIZE; i++) {
			for (size_t k = 0; k < params::poly_q::degree; k++) {
				int64_t coeff = sample_z(0.0, SIGMA_O);
				mpz_set_si(coeffs[k], coeff);
			}
			z1[i].mpz2poly(coeffs);
			z1[i].ntt_pow_phi();
			for (size_t k = 0; k < params::poly_q::degree; k++) {
				int64_t coeff = sample_z(0.0, SIGMA_O);
				mpz_set_si(coeffs[k], coeff);
			}
			y0[i].mpz2poly(coeffs);
			y0[i].ntt_pow_phi();
		}
		for (size_t k = 0; k < params::poly_q::degree; k++) {
			int64_t coeff = sample_z(0.0, SIGMA_O);
			mpz_set_si(coeffs[k], coeff);
		}
		y0[SIZE].mpz2poly(coeffs);
		y0[SIZE].ntt_pow_phi();
	} else {
		/* Prover samples z_1-b,y_b from Gaussian. */
		for (size_t i = 0; i < SIZE; i++) {
			for (size_t k = 0; k < params::poly_q::degree; k++) {
				int64_t coeff = sample_z(0.0, SIGMA_O);
				mpz_set_si(coeffs[k], coeff);
			}
			z0[i].mpz2poly(coeffs);
			z0[i].ntt_pow_phi();
			for (size_t k = 0; k < params::poly_q::degree; k++) {
				int64_t coeff = sample_z(0.0, SIGMA_O);
				mpz_set_si(coeffs[k], coeff);
			}
			y1[i].mpz2poly(coeffs);
			y1[i].ntt_pow_phi();
		}
		for (size_t k = 0; k < params::poly_q::degree; k++) {
			int64_t coeff = sample_z(0.0, SIGMA_O);
			mpz_set_si(coeffs[k], coeff);
		}
		z0[SIZE].mpz2poly(coeffs);
		z0[SIZE].ntt_pow_phi();
	}

	nfl::fastrandombytes(h[1 - b], BLAKE3_OUT_LEN);
	nfl::fastrandombytes_seed(h[1 - b]);
	pior_sample_chall(c[1 - b]);

	if (b == 0) {
		w1 = A1[0] * z1[0] + A1[1] * z1[1] - c[1] * t1;
		w0[0] = A0[0][0] * y0[0] + A0[0][1] * y0[1] + A0[0][2] * y0[2];
		w0[1] = A0[1][0] * y0[0] + A0[1][1] * y0[1] + A0[1][2] * y0[2];
	} else {
		w0[0] = A0[0][0] * z0[0] + A0[0][1] * z0[1] + A0[0][2] * z0[2] -
				c[0] * t0[0];
		w0[1] = A0[1][0] * z0[0] + A0[1][1] * z0[1] + A0[1][2] * z0[2] -
				c[0] * t0[1];
		w1 = A1[0] * y1[0] + A1[1] * y1[1];
	}

	pior_hash(h[b], prehash, w0, w1, t0);
	for (size_t i = 0; i < BLAKE3_OUT_LEN; i++) {
		h[b][i] ^= h[1 - b][i];
	}
	nfl::fastrandombytes_seed(h[b]);
	pior_sample_chall(c[b]);

	if (b == 0) {
		z0[0] = c[b] * s0[0] + y0[0];
		z0[1] = c[b] * s0[1] + y0[1];
		z0[2] = c[b] * s0[2] + y0[2];
	} else {
		z1[0] = c[b] * s1[0] + y1[0];
		z1[1] = c[b] * s1[1] + y1[1];
	}

	for (size_t i = 0; i < params::poly_q::degree; i++) {
		mpz_clear(coeffs[i]);
	}
}

int pior_verify(uint8_t h[2][BLAKE3_OUT_LEN], uint8_t prehash[BLAKE3_OUT_LEN],
		params::poly_q w0[SIZE], params::poly_q w1,
		params::poly_q z0[SIZE + 1], params::poly_q z1[SIZE],
		params::poly_q A0[SIZE][SIZE + 1], params::poly_q A1[SIZE],
		params::poly_q t0[SIZE], params::poly_q t1) {
	int result = 0;
	uint8_t hash[BLAKE3_OUT_LEN];
	params::poly_q c[2], y0[2], y1;

	pior_hash(hash, prehash, w0, w1, t0);
	for (size_t i = 0; i < BLAKE3_OUT_LEN; i++) {
		hash[i] ^= h[0][i];
		hash[i] ^= h[1][i];
	}

	nfl::fastrandombytes_seed(h[0]);
	pior_sample_chall(c[0]);
	nfl::fastrandombytes_seed(h[1]);
	pior_sample_chall(c[1]);

	y0[0] = A0[0][0] * z0[0] + A0[0][1] * z0[1] + A0[0][2] * z0[2] -
			c[0] * t0[0];
	y0[1] = A0[1][0] * z0[0] + A0[1][1] * z0[1] + A0[1][2] * z0[2] -
			c[0] * t0[1];
	y1 = A1[0] * z1[0] + A1[1] * z1[1] - c[1] * t1;

	result = ((w1 - y1) == 0);
	result &= ((w0[0] - y0[0]) == 0);
	result &= ((w0[1] - y0[1]) == 0);

	for (int i = 0; i < SIZE + 1; i++) {
		z0[i].invntt_pow_invphi();
		result &= pior_test_norm(z0[i], 4 * SIGMA_O * SIGMA_O);
	}
	z1[0].invntt_pow_invphi();
	z1[1].invntt_pow_invphi();
	result &= pior_test_norm(z1[0], 4 * SIGMA_O * SIGMA_O);
	result &= pior_test_norm(z1[1], 4 * SIGMA_O * SIGMA_O);

	return result;
}

#ifdef MAIN
static void test() {
	params::poly_q A0[SIZE][SIZE + 1], A1[SIZE], s0[SIZE + 1], s1[SIZE];
	params::poly_q z0[SIZE + 1], z1[SIZE], t0[SIZE], t1, w0[SIZE], w1;
	uint8_t ph[BLAKE3_OUT_LEN];
	uint8_t h[2][BLAKE3_OUT_LEN];

	for (int i = 0; i < SIZE + 1; i++) {
		A0[0][i] = nfl::uniform();
		A0[1][i] = nfl::uniform();
		s0[i] = nfl::ZO_dist();
		s0[i].ntt_pow_phi();
	}
	for (int i = 0; i < SIZE; i++) {
		A1[i] = nfl::uniform();
		s1[i] = nfl::ZO_dist();
		s1[i].ntt_pow_phi();
	}

	t0[0] = A0[0][0] * s0[0] + A0[0][1] * s0[1] + A0[0][2] * s0[2];
	t0[1] = A0[1][0] * s0[0] + A0[1][1] * s0[1] + A0[1][2] * s0[2];
	t1 = A1[0] * s1[0] + A1[1] * s1[1];

	pior_prehash(ph, A0, A1, t1);

	TEST_ONCE("OR proof is consistent") {
		pior_prove(h, ph, w0, w1, z0, z1, s0, s1, A0, A1, t0, t1, 0);
		TEST_ASSERT(pior_verify(h, ph, w0, w1, z0, z1, A0, A1, t0, t1) == 1,
				end);
		t1 = nfl::uniform();
		pior_prove(h, ph, w0, w1, z0, z1, s0, s1, A0, A1, t0, t1, 0);
		TEST_ASSERT(pior_verify(h, ph, w0, w1, z0, z1, A0, A1, t0, t1) == 1,
				end);
		t0[0] = nfl::uniform();
		t0[0] = nfl::uniform();
		t1 = A1[0] * s1[0] + A1[1] * s1[1];
		pior_prove(h, ph, w0, w1, z0, z1, s0, s1, A0, A1, t0, t1, 1);
		TEST_ASSERT(pior_verify(h, ph, w0, w1, z0, z1, A0, A1, t0, t1) == 1,
				end);
		t1 = nfl::uniform();
		pior_prove(h, ph, w0, w1, z0, z1, s0, s1, A0, A1, t0, t1, 1);
		TEST_ASSERT(pior_verify(h, ph, w0, w1, z0, z1, A0, A1, t0, t1) == 0,
				end);
		pior_prove(h, ph, w0, w1, z0, z1, s0, s1, A0, A1, t0, t1, 1);
		TEST_ASSERT(pior_verify(h, ph, w0, w1, z0, z1, A0, A1, t0, t1) == 0,
				end);
	} TEST_END;

  end:
	return;
}

static void bench() {
	params::poly_q A0[SIZE][SIZE + 1], A1[SIZE], s0[SIZE + 1], s1[SIZE];
	params::poly_q z0[SIZE + 1], z1[SIZE], t0[SIZE], t1, w0[SIZE], w1;
	uint8_t ph[BLAKE3_OUT_LEN];
	uint8_t h[2][BLAKE3_OUT_LEN];
	blake3_hasher hasher;

	for (int i = 0; i < SIZE + 1; i++) {
		A0[0][i] = nfl::uniform();
		A0[1][i] = nfl::uniform();
		s0[i] = nfl::ZO_dist();
		s0[i].ntt_pow_phi();
	}
	for (int i = 0; i < SIZE; i++) {
		A1[i] = nfl::uniform();
		s1[i] = nfl::ZO_dist();
		s1[i].ntt_pow_phi();
	}

	t0[0] = A0[0][0] * s0[0] + A0[0][1] * s0[1] + A0[0][2] * s0[2];
	t0[1] = A0[1][0] * s0[0] + A0[1][1] * s0[1] + A0[1][2] * s0[2];
	t1 = A1[0] * s1[0] + A1[1] * s1[1];
	pior_prehash(ph, A0, A1, t1);

	blake3_hasher_init(&hasher);
	BENCH_BEGIN("blake3 hashing (65536 bytes)") {
		BENCH_ADD(blake3_hasher_update(&hasher, t1.data(), 16 * DEGREE);
				blake3_hasher_finalize(&hasher, h[0], BLAKE3_OUT_LEN));
	} BENCH_END;

	BENCH_BEGIN("pior_prove") {
		BENCH_ADD(pior_prove(h, ph, w0, w1, z0, z1, s0, s1, A0, A1, t0, t1, 0));
	} BENCH_END;

	BENCH_BEGIN("pior_verify") {
		BENCH_ADD(pior_verify(h, ph, w0, w1, z0, z1, A0, A1, t0, t1));
	} BENCH_END;
}

int main(int argc, char *argv[]) {
	printf("\n** Tests for lattice-based OR proof:\n\n");
	test();

	printf("\n** Benchmarks for lattice-based OR proof:\n\n");
	bench();
	printf("\nMultiply prover by 3 due to rejection sampling.\n");
}
#endif
