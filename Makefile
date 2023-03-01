CC = clang
CPP = g++
CFLAGS = -O3 -march=native -mtune=native -Wall -ggdb -I NFLlib/include/ -I NFLlib/include/nfl -I NFLlib/include/nfl/prng -I include -I src/falcon -DNFL_OPTIMIZED=ON -DNTT_AVX2
INCLUDES = include/bench.h include/cpucycles.h
BENCH = src/bench.c src/cpucycles.c
TEST = src/test.c
BLAKE3 = src/blake3/blake3.c src/blake3/blake3_dispatch.c src/blake3/blake3_portable.c src/blake3/blake3_sse2_x86-64_unix.S src/blake3/blake3_sse41_x86-64_unix.S src/blake3/blake3_avx2_x86-64_unix.S src/blake3/blake3_avx512_x86-64_unix.S
LIBS = deps/libnfllib_static.a -lgmp -lmpfr -L deps/ -lflint -lquadmath deps/falcon.a

all: deps/falcon.a bdlop1 bgv1 ghl1 pior1 pibdn1 linear1 ks1 bdlop2 bgv2 ghl2 pior2 pibdn2 linear2 ks2

bdlop1: src/bdlop.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -DUSERS=1024 -c src/bgv.cpp -o bgv1.o
	${CPP} ${CFLAGS} -DMAIN src/bdlop.cpp bgv1.o ${TEST} ${BENCH} -o bdlop1 ${LIBS}

bdlop2: src/bdlop.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -DUSERS=10000 -c src/bgv.cpp -o bgv2.o
	${CPP} ${CFLAGS} -DMAIN src/bdlop.cpp bgv2.o ${TEST} ${BENCH} -o bdlop2 ${LIBS}

bgv1: src/bgv.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -DUSERS=1024 -DMAIN src/bgv.cpp ${TEST} ${BENCH} -o bgv1 ${LIBS}

bgv2: src/bgv.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -DUSERS=10000 -DMAIN src/bgv.cpp ${TEST} ${BENCH} -o bgv2 ${LIBS}

ghl1: src/ghl.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -DUSERS=1024 -DMAIN src/ghl.cpp ${TEST} ${BENCH} -o ghl1 ${LIBS}

ghl2: src/ghl.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -DUSERS=10000 -DMAIN src/ghl.cpp ${TEST} ${BENCH} -o ghl2 ${LIBS}

pior1: src/pior.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -c src/sample_z_small.c -o sample_z_small.o
	${CPP} ${CFLAGS} -DUSERS=1024 -DMAIN src/pior.cpp ${BLAKE3} sample_z_small.o ${TEST} ${BENCH} -o pior1 ${LIBS}

pior2: src/pior.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -c src/sample_z_small.c -o sample_z_small.o
	${CPP} ${CFLAGS} -DUSERS=10000 -DMAIN src/pior.cpp ${BLAKE3} sample_z_small.o ${TEST} ${BENCH} -o pior2 ${LIBS}

pibdn1: src/pibdn.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -c src/sample_z_small.c -o sample_z_small.o
	${CPP} ${CFLAGS} -c src/sample_z_large.c -o sample_z_large.o
	${CPP} ${CFLAGS} -DUSERS=1024 -DMAIN src/pibdn.cpp ${BLAKE3} sample_z_small.o sample_z_large.o ${TEST} ${BENCH} -o pibdn1 ${LIBS}

pibdn2: src/pibdn.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -c src/sample_z_small.c -o sample_z_small.o
	${CPP} ${CFLAGS} -c src/sample_z_large.c -o sample_z_large.o
	${CPP} ${CFLAGS} -DUSERS=10000 -DMAIN src/pibdn.cpp ${BLAKE3} sample_z_small.o sample_z_large.o ${TEST} ${BENCH} -o pibdn2 ${LIBS}

linear1: src/linear.cpp src/linear.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -c src/sample_z_small.c -o sample_z_small.o
	${CPP} ${CFLAGS} -c src/bdlop.cpp -o bdlop.o
	${CPP} ${CFLAGS} -DUSERS=1024 -DMAIN src/linear.cpp ${BLAKE3} sample_z_small.o bdlop.o ${TEST} ${BENCH} -o linear1 ${LIBS}

linear2: src/linear.cpp src/linear.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -c src/sample_z_small.c -o sample_z_small.o
	${CPP} ${CFLAGS} -c src/bdlop.cpp -o bdlop.o
	${CPP} ${CFLAGS} -DUSERS=10000 -DMAIN src/linear.cpp ${BLAKE3} sample_z_small.o bdlop.o ${TEST} ${BENCH} -o linear2 ${LIBS}

ks1: src/ks.cpp src/ghl.cpp
	${CPP} ${CFLAGS} -DUSERS=1024 -c src/bgv.cpp -o bgv.o
	${CPP} ${CFLAGS} -DUSERS=1024 -c src/ghl.cpp -o ghl.o
	${CPP} ${CFLAGS} -DUSERS=1024 -DMAIN src/ks.cpp ${BLAKE3} ghl.o bgv.o ${TEST} ${BENCH} -o ks1 ${LIBS}

ks2: src/ks.cpp src/ghl.cpp
	${CPP} ${CFLAGS} -DUSERS=10000 -c src/bgv.cpp -o bgv.o
	${CPP} ${CFLAGS} -DUSERS=10000 -c src/ghl.cpp -o ghl.o
	${CPP} ${CFLAGS} -DUSERS=10000 -DMAIN src/ks.cpp ${BLAKE3} ghl.o bgv.o ${TEST} ${BENCH} -o ks2 ${LIBS}

deps/falcon.a: src/falcon/*.c
	${CC} -Wall -Wextra -Wshadow -Wundef -O3  -c -o src/falcon/codec.o src/falcon/codec.c
	${CC} -Wall -Wextra -Wshadow -Wundef -O3  -c -o src/falcon/common.o src/falcon/common.c
	${CC} -Wall -Wextra -Wshadow -Wundef -O3  -c -o src/falcon/falcon.o src/falcon/falcon.c
	${CC} -Wall -Wextra -Wshadow -Wundef -O3  -c -o src/falcon/fft.o src/falcon/fft.c
	${CC} -Wall -Wextra -Wshadow -Wundef -O3  -c -o src/falcon/fpr.o src/falcon/fpr.c
	${CC} -Wall -Wextra -Wshadow -Wundef -O3  -c -o src/falcon/keygen.o src/falcon/keygen.c
	${CC} -Wall -Wextra -Wshadow -Wundef -O3  -c -o src/falcon/rng.o src/falcon/rng.c
	${CC} -Wall -Wextra -Wshadow -Wundef -O3  -c -o src/falcon/shake.o src/falcon/shake.c
	${CC} -Wall -Wextra -Wshadow -Wundef -O3  -c -o src/falcon/sign.o src/falcon/sign.c
	${CC} -Wall -Wextra -Wshadow -Wundef -O3  -c -o src/falcon/vrfy.o src/falcon/vrfy.c
	ar -rv deps/falcon.a src/falcon/*.o

clean:
	rm *.o bdlop1 bgv1 ghl1 pior1 pibdn1 linear1 bdlop2 bgv2 ghl2 pior2 pibdn2 linear2 ks1 ks2
