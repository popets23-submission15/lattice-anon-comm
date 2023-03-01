#include <math.h>
#include <stdlib.h>

#include <flint/flint.h>
#include <flint/fmpz_mod_poly.h>

#include "common.h"
#include "test.h"
#include "bench.h"
#include "assert.h"
#include "blake3.h"
#include "sample_z_small.h"

#define MSGS 2

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

void linear_hash(params::poly_q & beta, comkey_t & key, commit_t x,
		commit_t y, params::poly_q alpha[2], params::poly_q & u,
		params::poly_q t, params::poly_q _t) {
	uint8_t hash[BLAKE3_OUT_LEN];
	blake3_hasher hasher;

	blake3_hasher_init(&hasher);

	/* Hash public key. */
	for (size_t i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH - HEIGHT; j++) {
			blake3_hasher_update(&hasher, key.A1[i][j].data(), 16 * DEGREE);
		}
		for (size_t j = 0; j < WIDTH; j++) {
			blake3_hasher_update(&hasher, key.A2[i][j].data(), 16 * DEGREE);
		}
	}

	/* Hash alpha, beta from linear relation. */
	for (size_t i = 0; i < 2; i++) {
		blake3_hasher_update(&hasher, alpha[i].data(), 16 * DEGREE);
	}

	blake3_hasher_update(&hasher, x.c1.data(), 16 * DEGREE);
	blake3_hasher_update(&hasher, y.c1.data(), 16 * DEGREE);
	for (size_t i = 0; i < x.c2.size(); i++) {
		blake3_hasher_update(&hasher, x.c2.data(), 16 * DEGREE);
		blake3_hasher_update(&hasher, y.c2.data(), 16 * DEGREE);
	}

	blake3_hasher_update(&hasher, u.data(), 16 * DEGREE);
	blake3_hasher_update(&hasher, t.data(), 16 * DEGREE);
	blake3_hasher_update(&hasher, _t.data(), 16 * DEGREE);

	blake3_hasher_finalize(&hasher, hash, BLAKE3_OUT_LEN);

	/* Sample challenge from RNG seeded with hash. */
	nfl::fastrandombytes_seed(hash);
	bdlop_sample_chal(beta);
	nfl::fastrandombytes_reseed();
}

static void poly_inverse(params::poly_q & inv, params::poly_q p) {
	std::array < mpz_t, params::poly_q::degree > coeffs;
	fmpz_t q;
	fmpz_mod_poly_t poly, irred;
	fmpz_mod_ctx_t ctx_q;

	fmpz_init(q);
	for (size_t i = 0; i < params::poly_q::degree; i++) {
		mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
	}

	fmpz_set_mpz(q, params::poly_q::moduli_product());
	fmpz_mod_ctx_init(ctx_q, q);
	fmpz_mod_poly_init(poly, ctx_q);
	fmpz_mod_poly_init(irred, ctx_q);

	p.poly2mpz(coeffs);
	fmpz_mod_poly_set_coeff_ui(irred, params::poly_q::degree, 1, ctx_q);
	fmpz_mod_poly_set_coeff_ui(irred, 0, 1, ctx_q);

	for (size_t i = 0; i < params::poly_q::degree; i++) {
		fmpz_mod_poly_set_coeff_mpz(poly, i, coeffs[i], ctx_q);
	}
	fmpz_mod_poly_invmod(poly, poly, irred, ctx_q);

	for (size_t i = 0; i < params::poly_q::degree; i++) {
		fmpz_mod_poly_get_coeff_mpz(coeffs[i], poly, i, ctx_q);
	}

	inv.mpz2poly(coeffs);

	fmpz_clear(q);
	for (size_t i = 0; i < params::poly_q::degree; i++) {
		mpz_clear(coeffs[i]);
	}
}

static void simul_inverse(params::poly_q inv[MSGS], params::poly_q m[MSGS]) {
	params::poly_q u, t[MSGS];
	inv[0] = m[0];
	t[0] = m[0];

	for (size_t i = 1; i < MSGS; i++) {
		t[i] = m[i];
		inv[i] = inv[i - 1] * m[i];
	}

	u = inv[MSGS - 1];
	u.invntt_pow_invphi();
	poly_inverse(u, u);
	u.ntt_pow_phi();

	for (size_t i = MSGS - 1; i > 0; i--) {
		inv[i] = u * inv[i - 1];
		u = u * t[i];
	}
	inv[0] = u;
}

static void linear_prove(params::poly_q y[WIDTH], params::poly_q _y[WIDTH],
		params::poly_q & t, params::poly_q & _t, params::poly_q & u,
		commit_t x, commit_t _x, params::poly_q alpha[2],
		comkey_t & key, vector < params::poly_q > r) {
	params::poly_q beta;
	std::array < mpz_t, params::poly_q::degree > coeffs;
	mpz_t qDivBy2;

	mpz_init(qDivBy2);
	for (size_t i = 0; i < params::poly_q::degree; i++) {
		mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
	}
	mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);

	/* Prover samples y,y' from Gaussian. */
	for (int i = 0; i < WIDTH; i++) {
		for (size_t k = 0; k < params::poly_q::degree; k++) {
			int64_t coeff = sample_z(0.0, SIGMA_C);
			mpz_set_si(coeffs[k], coeff);
		}
		y[i].mpz2poly(coeffs);
		y[i].ntt_pow_phi();
		for (size_t k = 0; k < params::poly_q::degree; k++) {
			int64_t coeff = sample_z(0.0, SIGMA_C);
			mpz_set_si(coeffs[k], coeff);
		}
		_y[i].mpz2poly(coeffs);
		_y[i].ntt_pow_phi();
	}

	t = y[0];
	_t = _y[0];
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH - HEIGHT; j++) {
			t = t + key.A1[i][j] * y[j + HEIGHT];
			_t = _t + key.A1[i][j] * _y[j + HEIGHT];
		}
	}

	u = 0;
	for (int i = 0; i < WIDTH; i++) {
		u = u + alpha[0] * (key.A2[0][i] * y[i]) - (key.A2[0][i] * _y[i]);
	}

	/* Sample challenge. */
	linear_hash(beta, key, x, _x, alpha, u, t, _t);

	/* Prover */
	for (int i = 0; i < WIDTH; i++) {
		y[i] = y[i] + beta * r[i];
		_y[i] = _y[i] + beta * r[i];
	}

	for (size_t i = 0; i < params::poly_q::degree; i++) {
		mpz_clear(coeffs[i]);
	}
	mpz_clear(qDivBy2);
}

static int linear_verify(commit_t x, commit_t _x, params::poly_q z[WIDTH],
		params::poly_q _z[WIDTH], params::poly_q t, params::poly_q _t,
		params::poly_q u, comkey_t & key, params::poly_q alpha[2], int l) {
	params::poly_q beta, v, _v, tmp, zero = 0;
	int result = 1;

	/* Sample challenge. */
	linear_hash(beta, key, x, _x, alpha, u, t, _t);

	/* Verifier checks norm, reconstruct from NTT representation. */
	for (int i = 0; i < WIDTH; i++) {
		v = z[i];
		v.invntt_pow_invphi();
		result &= bdlop_test_norm(v, 4 * SIGMA_C * SIGMA_C);
		v = _z[i];
		v.invntt_pow_invphi();
		result &= bdlop_test_norm(v, 4 * SIGMA_C * SIGMA_C);
	}

	/* Verifier computes A1z and A1z'. */
	v = z[0];
	_v = _z[0];
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH - HEIGHT; j++) {
			v = v + key.A1[i][j] * z[j + HEIGHT];
			_v = _v + key.A1[i][j] * _z[j + HEIGHT];
		}
	}

	tmp = t + beta * x.c1 - v;
	tmp.invntt_pow_invphi();
	result &= (tmp == zero);
	tmp = _t + beta * _x.c1 - _v;
	tmp.invntt_pow_invphi();
	result &= (tmp == zero);

	v = 0;
	for (int i = 0; i < WIDTH; i++) {
		v = v + alpha[0] * (key.A2[0][i] * z[i]) - (key.A2[0][i] * _z[i]);
	}
	t = (alpha[0] * x.c2[0] + alpha[1] - _x.c2[0]) * beta + u;

	t.invntt_pow_invphi();
	v.invntt_pow_invphi();

	result &= ((t - v) == 0);
	return result;
}

static int run(commit_t com[MSGS], vector < vector < params::poly_q >> m,
		vector < vector < params::poly_q >> _m, comkey_t & key, comkey_t & _key,
		vector < params::poly_q > r) {
	int result = 1;
	params::poly_q ms[MSGS], _ms[MSGS], inv[MSGS];
	commit_t d[MSGS], cs[MSGS];
	vector < params::poly_q > t0(1);
	params::poly_q one, t1, rho[SIZE], theta[MSGS], s[MSGS];
	params::poly_q y[WIDTH], _y[WIDTH], t, _t, u;
	size_t messages = m.size();

	rho[0] = 1;
	rho[0].ntt_pow_phi();
	for (size_t j = 1; j < SIZE; j++) {
		rho[j] = nfl::uniform();
	}
	for (size_t i = 0; i < messages; i++) {
		ms[i] = m[i][0];
		ms[i].ntt_pow_phi();
		_ms[i] = _m[i][0];
		_ms[i].ntt_pow_phi();
		cs[i].c1 = com[i].c1;
		cs[i].c2.resize(1);
		cs[i].c2[0] = com[i].c2[0] * rho[0];
		for (size_t j = 1; j < SIZE; j++) {
			cs[i].c2[0] = cs[i].c2[0] + com[i].c2[j] * rho[j];
			t1 = m[i][j];
			t1.ntt_pow_phi();
			ms[i] = ms[i] + t1 * rho[j];
			t1 = _m[i][j];
			t1.ntt_pow_phi();
			_ms[i] = _ms[i] + t1 * rho[j];
		}
	}

	for (size_t i = 0; i < HEIGHT; i++) {
		for (size_t j = HEIGHT; j < WIDTH; j++) {
			_key.A1[i][j - HEIGHT] = key.A1[i][j - HEIGHT];
		}
	}
	_key.A2[0][0] = 0;
	_key.A2[0][1] = rho[0];
	for (size_t j = 2; j < WIDTH; j++) {
		_key.A2[0][j] = key.A2[0][j];
	}
	for (size_t i = 1; i < SIZE; i++) {
		for (size_t j = 2; j < WIDTH; j++) {
			_key.A2[0][j] = _key.A2[0][j] + rho[i] * key.A2[i][j];
		}
	}

	/* Prover samples theta_i and computes commitments D_i. */
	theta[0] = nfl::ZO_dist();
	theta[0].ntt_pow_phi();
	t0[0] = theta[0] * _ms[0];
	t0[0].invntt_pow_invphi();
	bdlop_commit(d[0], t0, _key, r);
	t0[0] = theta[messages - 2] * ms[messages - 1];
	t0[0].invntt_pow_invphi();

	//Check relationship here
	params::poly_q beta = nfl::uniform();
	simul_inverse(inv, _ms);
	s[0] = theta[0] * _ms[0] - beta * ms[0];
	s[0] = s[0] * inv[0];

	/* Now run \Prod_LIN instances, one for each commitment. */
	t0[0] = s[0] * _ms[0];

	params::poly_q alpha[2] = { beta, t0[0] };
	linear_prove(y, _y, t, _t, u, cs[0], d[0], alpha, _key, r);
	result = linear_verify(cs[0], d[0], y, _y, t, _t, u, _key, alpha, 0);

	BENCH_SMALL("linear proof", linear_prove(y, _y, t, _t, u, cs[0], d[0],
					alpha, _key, r));
	BENCH_SMALL("linear verify", linear_verify(cs[0], d[0], y, _y, t, _t, u,
					_key, alpha, 0));

	return result;
}

#ifdef MAIN
static void test() {
	comkey_t key, _key;
	commit_t com[MSGS];
	vector < vector < params::poly_q >> m(MSGS), _m(MSGS);
	vector < params::poly_q > r(WIDTH);

	/* Generate commitment key-> */
	bdlop_keygen(key);
	bdlop_sample_rand(r);
	for (int i = 0; i < MSGS; i++) {
		m[i].resize(SIZE);
		for (int j = 0; j < SIZE; j++) {
			m[i][j] = nfl::ZO_dist();
		}
		bdlop_commit(com[i], m[i], key, r);
	}

	/* Prover shuffles messages (only a circular shift for simplicity). */
	for (int i = 0; i < MSGS; i++) {
		_m[i].resize(SIZE);
		for (int j = 0; j < SIZE; j++) {
			_m[i][j] = m[i][j];
		}
	}

	TEST_ONCE("linear proof is consistent") {
		TEST_ASSERT(run(com, m, _m, key, _key, r) == 1, end);
	} TEST_END;

  end:
	return;
}

int main(int argc, char *argv[]) {
	printf("\n** Tests and benchmarks for lattice-based linear proof:\n\n");
	test();

	printf("\nMultiply prover by 3 due to rejection sampling.\n");
}
#endif
