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

void pibdn_hash(uint8_t hash[BLAKE3_OUT_LEN], params::poly_q w[SIZE],
		params::poly_q A[SIZE][WIDTH], params::poly_q t[SIZE]) {
	blake3_hasher hasher;

	blake3_hasher_init(&hasher);

	for (size_t i = 0; i < SIZE; i++) {
		blake3_hasher_update(&hasher, w[i].data(), 16 * DEGREE);
		for (size_t j = 0; j < WIDTH; j++) {
			blake3_hasher_update(&hasher, A[i][j].data(), 16 * DEGREE);
		}
		blake3_hasher_update(&hasher, t[i].data(), 16 * DEGREE);
	}

	blake3_hasher_finalize(&hasher, hash, BLAKE3_OUT_LEN);
}

/**
 * Test if the l2-norm is within bounds (4 * sigma * sqrt(N)).
 *
 * @param[in] r 			- the polynomial to compute the l2-norm.
 * @return the computed norm.
 */
bool pibdn_test_norm1(params::poly_q r, uint64_t sigma_sqr) {
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

bool pibdn_test_norm2(params::poly_q r) {
	array < mpz_t, params::poly_q::degree > coeffs;
	mpz_t norm, qDivBy2, tmp, q;

	/// Constructors
	mpz_inits(norm, qDivBy2, tmp, q, nullptr);
	for (size_t i = 0; i < params::poly_q::degree; i++) {
		mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
	}

	r.poly2mpz(coeffs);
	mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);
	mpz_set_ui(norm, 0);
	mpz_set_str(q, PRIMEQ, 10);
	for (size_t i = 0; i < params::poly_q::degree; i++) {
		mpz_mod(coeffs[i], coeffs[i], q);
		util::center(coeffs[i], coeffs[i],
				params::poly_q::moduli_product(), qDivBy2);
		mpz_mul(tmp, coeffs[i], coeffs[i]);
		mpz_add(norm, norm, tmp);
	}

	mpz_set_str(tmp, BOUND_B, 10);
	mpz_mul(tmp, tmp, tmp);
	int result = mpz_cmp(norm, tmp) <= 0;

	mpz_clears(norm, qDivBy2, tmp, q, nullptr);
	for (size_t i = 0; i < params::poly_q::degree; i++) {
		mpz_clear(coeffs[i]);
	}

	return !result;
}

// Sample a challenge.
void pibdn_sample_chall(params::poly_q & f) {
	f = nfl::ZO_dist();
	f.ntt_pow_phi();
}

void pibdn_sample_large(params::poly_q & s) {
	std::array < mpz_t, DEGREE > coeffs;
	mpz_t qDivBy2, bound;

	mpz_init(qDivBy2);
	mpz_init(bound);
	for (size_t i = 0; i < params::poly_q::degree; i++) {
		mpz_init2(coeffs[i], params::poly_q::bits_in_moduli_product() << 2);
	}
	mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);
	mpz_set_str(bound, BOUND_D, 10);

	s = nfl::uniform();
	s.poly2mpz(coeffs);
	for (size_t i = 0; i < params::poly_q::degree; i++) {
		util::center(coeffs[i], coeffs[i], params::poly_q::moduli_product(),
				qDivBy2);
		mpz_mod(coeffs[i], coeffs[i], bound);
	}
	s.mpz2poly(coeffs);
	s.ntt_pow_phi();
}

void pibdn_prove(params::poly_q w[WIDTH], params::poly_q z[WIDTH],
		params::poly_q A[SIZE][WIDTH], params::poly_q s[WIDTH],
		params::poly_q t[SIZE]) {
	std::array < mpz_t, params::poly_q::degree > coeffs;
	params::poly_q c, y[WIDTH];
	uint8_t hash[BLAKE3_OUT_LEN];

	for (size_t i = 0; i < params::poly_q::degree; i++) {
		mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
	}

	/* Prover samples z_1-b,y_b from Gaussian. */
	for (size_t i = 0; i < WIDTH - 1; i++) {
		for (size_t k = 0; k < params::poly_q::degree; k++) {
			int64_t coeff = sample_z(0.0, SIGMA_O);
			mpz_set_si(coeffs[k], coeff);
		}
		y[i].mpz2poly(coeffs);
		y[i].ntt_pow_phi();
	}
	for (size_t k = 0; k < params::poly_q::degree; k++) {
		int64_t coeff = sample_z((__float128) 0.0, (__float128) SIGMA_B);
		mpz_set_si(coeffs[k], coeff);
	}
	y[WIDTH - 1].mpz2poly(coeffs);
	y[WIDTH - 1].ntt_pow_phi();

	w[0] = y[0] + A[0][1] * y[1] + A[0][2] * y[2];
	w[1] = y[1] + A[1][2] * y[2] + y[3];

	pibdn_hash(hash, w, A, t);
	nfl::fastrandombytes_seed(hash);
	pibdn_sample_chall(c);

	for (size_t i = 0; i < WIDTH; i++) {
		z[i] = c * s[i] + y[i];
	}

	for (size_t i = 0; i < params::poly_q::degree; i++) {
		mpz_clear(coeffs[i]);
	}
}

int pibdn_verify(params::poly_q w[SIZE], params::poly_q z[WIDTH],
		params::poly_q A[SIZE][WIDTH], params::poly_q t[SIZE]) {
	int result = 0;
	uint8_t hash[BLAKE3_OUT_LEN];
	params::poly_q c, y[SIZE];

	pibdn_hash(hash, w, A, t);
	nfl::fastrandombytes_seed(hash);
	pibdn_sample_chall(c);

	y[0] = z[0] + A[0][1] * z[1] + A[0][2] * z[2] - c * t[0];
	y[1] = z[1] + A[1][2] * z[2] + z[3] - c * t[1];

	result = ((w[0] - y[0]) == 0);
	result &= ((w[1] - y[1]) == 0);

	for (int i = 0; i < WIDTH - 1; i++) {
		z[i].invntt_pow_invphi();
		result &= pibdn_test_norm1(z[i], 4 * SIGMA_O * SIGMA_O);
	}
	z[WIDTH - 1].invntt_pow_invphi();
	result &= pibdn_test_norm2(z[WIDTH - 1]);

	return result;
}

#ifdef MAIN
static void test() {
	params::poly_q one, A[SIZE][WIDTH], s[WIDTH], t[SIZE], z[WIDTH], w[SIZE];

	one = 1;
	one.ntt_pow_phi();

	A[1][0] = A[0][3] = 0;
	A[0][0] = A[1][1] = A[1][3] = one;
	A[0][1] = nfl::uniform();
	A[0][2] = nfl::uniform();
	A[1][2] = nfl::uniform();
	for (int i = 0; i < WIDTH - 1; i++) {
		s[i] = nfl::ZO_dist();
		s[i].ntt_pow_phi();
	}
	pibdn_sample_large(s[WIDTH - 1]);

	t[0] = s[0] + A[0][1] * s[1] + A[0][2] * s[2];
	t[1] = s[1] + A[1][2] * s[2] + s[3];

	TEST_ONCE("boundedness proof is consistent") {
		pibdn_prove(w, z, A, s, t);
		TEST_ASSERT(pibdn_verify(w, z, A, t) == 1, end);
	} TEST_END;

  end:
	return;
}

static void bench() {
	params::poly_q one, A[SIZE][WIDTH], s[WIDTH], t[SIZE], z[WIDTH], w[SIZE];

	one = 1;
	one.ntt_pow_phi();

	A[1][0] = A[0][3] = 0;
	A[0][0] = A[1][1] = A[1][3] = one;
	A[0][1] = nfl::uniform();
	A[0][2] = nfl::uniform();
	A[1][2] = nfl::uniform();
	for (int i = 0; i < WIDTH - 1; i++) {
		s[i] = nfl::ZO_dist();
		s[i].ntt_pow_phi();
	}
	pibdn_sample_large(s[WIDTH - 1]);

	t[0] = s[0] + A[0][1] * s[1] + A[0][2] * s[2];
	t[1] = s[1] + A[1][2] * s[2] + s[3];

	BENCH_BEGIN("pibdn_prove") {
		BENCH_ADD(pibdn_prove(w, z, A, s, t));
	} BENCH_END;

	BENCH_BEGIN("pibdn_verify") {
		pibdn_prove(w, z, A, s, t);
		BENCH_ADD(pibdn_verify(w, z, A, t));
	} BENCH_END;
}

int main(int argc, char *argv[]) {
	printf("\n** Tests for lattice-based boundness proof:\n\n");
	test();

	printf("\n** Benchmarks for lattice-based boundness proof:\n\n");
	bench();
	printf("\nMultiply prover by 3 due to rejection sampling.\n");
}
#endif
