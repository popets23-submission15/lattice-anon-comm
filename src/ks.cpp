#include "test.h"
#include "bench.h"
#include "common.h"
#include "blake3.h"
#include "falcon.h"
#include <sys/random.h>

#define FALCON 9

static uint8_t *
xmalloc(size_t len)
{
	uint8_t *buf;

	if (len == 0) {
		return NULL;
	}
	buf = (uint8_t *)malloc(len);
	if (buf == NULL) {
		fprintf(stderr, "memory allocation error\n");
		exit(EXIT_FAILURE);
	}
	return buf;
}

static void
xfree(uint8_t *buf)
{
	if (buf != NULL) {
		free(buf);
	}
}

static inline size_t
maxsz(size_t a, size_t b)
{
	return a > b ? a : b;
}

void ks_comh(uint8_t *hash, params::poly_2q &m) {
	blake3_hasher hasher;

	blake3_hasher_init(&hasher);
    blake3_hasher_update(&hasher, m.data(), 16 * DEGREE);
    blake3_hasher_finalize(&hasher, hash, BLAKE3_OUT_LEN);
}

void ks_comh(uint8_t *hash, ghlenc_t &c) {
	blake3_hasher hasher;

	blake3_hasher_init(&hasher);
    blake3_hasher_update(&hasher, c.u.data(), 16 * DEGREE);
    blake3_hasher_update(&hasher, c.v.data(), 16 * DEGREE);
    blake3_hasher_finalize(&hasher, hash, BLAKE3_OUT_LEN);
}

void ks_perm(size_t *p, size_t n) {
	size_t i, j, k;

	for (i = 0; i < n; i++) {
		getrandom((uint8_t *)&k, sizeof(size_t), 0);
		j = k % (i+1);
		p[i] = p[j];
		p[j] = i;
	}
}

static void ks_client_submit(ghlenc_t *ctx, falcon_t *sig, uint8_t h[USERS][BLAKE3_OUT_LEN],
		bgvkey_t *cpk, params::poly_q *csk, ghlkey_t &spk, falcon_t *fkp) {
    params::poly_2q _m;

    for (size_t i = 0; i < USERS; i++) {
        params::poly_2q r = nfl::ZO_dist();
        bgv_keygen(cpk[i], csk[i]);
        ghl_encode(_m, cpk[i].b);
        ks_comh(h[i], r);
        ghl_encrypt(ctx[i], spk, _m, r);
		falcon_keygen_make(&fkp->rng, fkp->logn, fkp->sk,
			FALCON_PRIVKEY_SIZE(fkp->logn), fkp->pk,
			FALCON_PUBKEY_SIZE(fkp->logn), fkp->tmp, fkp->tmp_len);
		uint8_t hash[BLAKE3_OUT_LEN];
		ks_comh(hash, ctx[i]);
		falcon_sign_dyn(&sig->rng,
			sig->sig, &sig->sig_len, FALCON_SIG_COMPRESSED,
			sig->sk, FALCON_PRIVKEY_SIZE(sig->logn),
			hash, BLAKE3_OUT_LEN, sig->tmp, sig->tmp_len);
    }
}

static void ks_server_submit(falcon_t *fkp) {
	for (size_t i = 0; i < PARTIES; i++) {
		falcon_keygen_make(&fkp->rng, fkp->logn, fkp->sk,
			FALCON_PRIVKEY_SIZE(fkp->logn), fkp->pk,
			FALCON_PUBKEY_SIZE(fkp->logn), fkp->tmp, fkp->tmp_len);
	}
}

static void ks_shuffle(ghlenc_t *out, falcon_t *sig, uint8_t h[BLAKE3_OUT_LEN], uint8_t hs[USERS][BLAKE3_OUT_LEN], ghlenc_t *in, ghlkey_t &spk, falcon_t fkp) {
    params::poly_2q _m, rho = nfl::ZO_dist();
    params::poly_q zero = 0;
    ghlenc_t enc_zero;
    size_t shuffle[USERS];

    ks_comh(h, rho);
    ks_perm(shuffle, USERS);
    ghl_encode(_m, zero);
    ghl_encrypt(enc_zero, spk, _m, rho);

    for (size_t i = 0; i < USERS; i++) {
		ks_comh(hs[i], in[i]);
		falcon_verify(
			sig->sig, sig->sig_len, FALCON_SIG_COMPRESSED,
			sig->pk, FALCON_PUBKEY_SIZE(sig->logn),
			hs[i], BLAKE3_OUT_LEN, sig->tmp, sig->tmp_len);
        ghl_add(out[i], in[i], enc_zero);
        ks_comh(hs[i], out[i]);
		falcon_sign_dyn(&sig->rng,
			sig->sig, &sig->sig_len, FALCON_SIG_COMPRESSED,
			sig->sk, FALCON_PRIVKEY_SIZE(sig->logn),
			hs[i], BLAKE3_OUT_LEN, sig->tmp, sig->tmp_len);
    }
}

static void ks_distdec(params::poly_2q *t, uint8_t hs[USERS][BLAKE3_OUT_LEN], ghlenc_t *ctx, params::poly_2q &ssk) {
    for (size_t i = 0; i < USERS; i++) {
        ghl_distdec(t[i], ctx[i], ssk);
        ks_comh(hs[i], t[i]);
    }
}

static void test() {
	falcon_t client_fkp[USERS], server_fkp[PARTIES], sig[USERS];
    ghlenc_t ctx[USERS], out[USERS];
    ghlkey_t spk;
    bgvkey_t cpk[USERS];
    params::poly_q csk[USERS];
    params::poly_2q sk, ssk[PARTIES], acc, t[PARTIES][USERS], _t[PARTIES], m, _m;
    uint8_t h[BLAKE3_OUT_LEN], hs[USERS][BLAKE3_OUT_LEN];
	size_t len;

    TEST_ONCE("key scheduling setup is correct") {
        ghl_keygen(spk, sk);
        ghl_keyshare(ssk, PARTIES, sk);
        acc = ssk[0];
        for (size_t j = 1; j < PARTIES; j++) {
            acc = acc + ssk[j];
        }
        TEST_ASSERT(sk - acc == 0, end);

		for (size_t i = 0; i < USERS; i++) {
			sig[i].logn = FALCON;
			shake256_init_prng_from_system(&sig[i].rng);
			len = FALCON_TMPSIZE_KEYGEN(FALCON);
			len = maxsz(len, FALCON_TMPSIZE_SIGNDYN(FALCON));
			len = maxsz(len, FALCON_TMPSIZE_SIGNTREE(FALCON));
			len = maxsz(len, FALCON_TMPSIZE_EXPANDPRIV(FALCON));
			len = maxsz(len, FALCON_TMPSIZE_VERIFY(FALCON));
			sig[i].tmp = xmalloc(len);
			sig[i].tmp_len = len;
			sig[i].pk = xmalloc(FALCON_PUBKEY_SIZE(FALCON));
			sig[i].sk = xmalloc(FALCON_PRIVKEY_SIZE(FALCON));
			sig[i].esk = xmalloc(FALCON_EXPANDEDKEY_SIZE(FALCON));
			sig[i].sig = xmalloc(FALCON_SIG_COMPRESSED_MAXSIZE(FALCON));
			sig[i].sig_len = 0;
			sig[i].sigct = xmalloc(FALCON_SIG_CT_SIZE(FALCON));
			sig[i].sigct_len = 0;

			client_fkp[i].logn = FALCON;
			shake256_init_prng_from_system(&client_fkp[i].rng);
			client_fkp[i].tmp = xmalloc(len);
			client_fkp[i].tmp_len = len;
			client_fkp[i].pk = xmalloc(FALCON_PUBKEY_SIZE(FALCON));
			client_fkp[i].sk = xmalloc(FALCON_PRIVKEY_SIZE(FALCON));
			client_fkp[i].esk = xmalloc(FALCON_EXPANDEDKEY_SIZE(FALCON));
			client_fkp[i].sig = xmalloc(FALCON_SIG_COMPRESSED_MAXSIZE(FALCON));
			client_fkp[i].sig_len = 0;
			client_fkp[i].sigct = xmalloc(FALCON_SIG_CT_SIZE(FALCON));
			client_fkp[i].sigct_len = 0;
		}

		for (size_t i = 0; i < PARTIES; i++) {
			server_fkp[i].logn = FALCON;
			shake256_init_prng_from_system(&server_fkp[i].rng);
			server_fkp[i].tmp = xmalloc(len);
			server_fkp[i].tmp_len = len;
			server_fkp[i].pk = xmalloc(FALCON_PUBKEY_SIZE(FALCON));
			server_fkp[i].sk = xmalloc(FALCON_PRIVKEY_SIZE(FALCON));
			server_fkp[i].esk = xmalloc(FALCON_EXPANDEDKEY_SIZE(FALCON));
			server_fkp[i].sig = xmalloc(FALCON_SIG_COMPRESSED_MAXSIZE(FALCON));
			server_fkp[i].sig_len = 0;
			server_fkp[i].sigct = xmalloc(FALCON_SIG_CT_SIZE(FALCON));
			server_fkp[i].sigct_len = 0;
		}

        ks_client_submit(ctx, sig, hs, cpk, csk, spk, client_fkp);
		ks_server_submit(server_fkp);
    } TEST_END;

    TEST_ONCE("key scheduling is correct") {
        for (size_t j = 0; j < PARTIES; j++) {
			ks_shuffle(out, sig, h, hs, ctx, spk, server_fkp[j]);
			for (size_t i = 0; i < USERS; i++) {
				ctx[i] = out[i];
			}
			ks_distdec(t[j], hs, ctx, ssk[j]);
		}
        for (size_t i = 0; i < USERS; i++) {
            for (size_t j = 0; j < PARTIES; j++) {
                _t[j] = t[j][i];
            }
            ghl_comb(_m, out[i], _t, PARTIES);
            ghl_decrypt(m, ctx[i], sk);
            TEST_ASSERT(m - _m == 1, end);
        }
    } TEST_END;

    end:
      return;
}

static void bench() {
	falcon_t client_fkp[USERS], server_fkp[PARTIES], sig[USERS];
    ghlenc_t ctx[USERS], out[USERS];
    ghlkey_t spk;
    bgvkey_t cpk[USERS];
    params::poly_q csk[USERS], m;
    params::poly_2q sk, ssk[PARTIES], acc, t[USERS], _m;
    uint8_t h[BLAKE3_OUT_LEN], hs[USERS][BLAKE3_OUT_LEN];
	size_t len;

    ghl_keygen(spk, sk);
    ghl_keyshare(ssk, PARTIES, sk);
    ghl_sample_message(m);

    BENCH_BEGIN("ks_comh") {
		BENCH_ADD(ks_comh(h, _m));
	} BENCH_END;

    ghl_encode(_m, m);
    ghl_encrypt(ctx[0], spk, _m);
    BENCH_BEGIN("ks_comh") {
		BENCH_ADD(ks_comh(h, ctx[0]));
	} BENCH_END;

	for (size_t i = 0; i < USERS; i++) {
		sig[i].logn = FALCON;
		shake256_init_prng_from_system(&sig[i].rng);
		len = FALCON_TMPSIZE_KEYGEN(FALCON);
		len = maxsz(len, FALCON_TMPSIZE_SIGNDYN(FALCON));
		len = maxsz(len, FALCON_TMPSIZE_SIGNTREE(FALCON));
		len = maxsz(len, FALCON_TMPSIZE_EXPANDPRIV(FALCON));
		len = maxsz(len, FALCON_TMPSIZE_VERIFY(FALCON));
		sig[i].tmp = xmalloc(len);
		sig[i].tmp_len = len;
		sig[i].pk = xmalloc(FALCON_PUBKEY_SIZE(FALCON));
		sig[i].sk = xmalloc(FALCON_PRIVKEY_SIZE(FALCON));
		sig[i].esk = xmalloc(FALCON_EXPANDEDKEY_SIZE(FALCON));
		sig[i].sig = xmalloc(FALCON_SIG_COMPRESSED_MAXSIZE(FALCON));
		sig[i].sig_len = 0;
		sig[i].sigct = xmalloc(FALCON_SIG_CT_SIZE(FALCON));
		sig[i].sigct_len = 0;

		client_fkp[i].logn = FALCON;
		shake256_init_prng_from_system(&client_fkp[i].rng);
		client_fkp[i].tmp = xmalloc(len);
		client_fkp[i].tmp_len = len;
		client_fkp[i].pk = xmalloc(FALCON_PUBKEY_SIZE(FALCON));
		client_fkp[i].sk = xmalloc(FALCON_PRIVKEY_SIZE(FALCON));
		client_fkp[i].esk = xmalloc(FALCON_EXPANDEDKEY_SIZE(FALCON));
		client_fkp[i].sig = xmalloc(FALCON_SIG_COMPRESSED_MAXSIZE(FALCON));
		client_fkp[i].sig_len = 0;
		client_fkp[i].sigct = xmalloc(FALCON_SIG_CT_SIZE(FALCON));
		client_fkp[i].sigct_len = 0;
	}

	for (size_t i = 0; i < PARTIES; i++) {
		server_fkp[i].logn = FALCON;
		shake256_init_prng_from_system(&server_fkp[i].rng);
		server_fkp[i].tmp = xmalloc(len);
		server_fkp[i].tmp_len = len;
		server_fkp[i].pk = xmalloc(FALCON_PUBKEY_SIZE(FALCON));
		server_fkp[i].sk = xmalloc(FALCON_PRIVKEY_SIZE(FALCON));
		server_fkp[i].esk = xmalloc(FALCON_EXPANDEDKEY_SIZE(FALCON));
		server_fkp[i].sig = xmalloc(FALCON_SIG_COMPRESSED_MAXSIZE(FALCON));
		server_fkp[i].sig_len = 0;
		server_fkp[i].sigct = xmalloc(FALCON_SIG_CT_SIZE(FALCON));
		server_fkp[i].sigct_len = 0;
	}

    BENCH_SMALL("ks_client_submit", ks_client_submit(ctx, sig, hs, cpk, csk, spk, client_fkp));

	BENCH_SMALL("ks_server_submit", ks_server_submit(server_fkp));

    BENCH_SMALL("ks_shuffle", ks_shuffle(out, sig, h, hs, ctx, spk, server_fkp[0]));

    BENCH_SMALL("ks_distdec", ks_distdec(t, hs, out, ssk[0]));
}

#ifdef MAIN
int main(int argc, char *argv[]) {
	printf("\n** Tests and benchmarks for lattice-based key schedule:\n\n");
	test();

    printf("\n** Benchmarks for lattice-based key scheduling:\n\n");
    bench();
}
#endif
