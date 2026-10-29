# Gosset-KEM — build, test and benchmark
# make            -> libraries + test programs
# make test       -> full test suite (NIST ACVP, Gosset-KEM KATs, reference == optimized, masked == unmasked,
#                    X25519, cross-platform digest)
# make bench      -> benchmarks and constant-time (dudect-style) tests
CC      ?= gcc
CFLAGS  ?= -O3 -mavx2 -Wall -Wno-unused-function
CFLAGS  += -Isrc -Isrc/masked
LDLIBS  = -lm
B       = build

CORE    = src/mlkem.c src/fips202.c src/e8fast.c src/e8codec.c
FULL    = $(CORE) src/dem.c src/x25519.c src/hybrid.c
REF     = src/ref/mlkem_ref.c src/fips202.c src/e8fast.c src/e8codec.c
MASKED  = src/masked/masked_decaps.c src/masked/masked_core.c src/masked/masked_keccak.c src/masked/masked_e8bs.c \
          src/fips202.c src/e8fast.c src/e8codec.c

all: $(B)/libgosset.so $(B)/libgosset_ref.so $(B)/libgosset_masked.so $(B)/xplat_kat $(B)/demo

$(B):
	mkdir -p $(B)
$(B)/libgosset.so: $(FULL) | $(B)
	$(CC) $(CFLAGS) -shared -fPIC -o $@ $(FULL) $(LDLIBS)
$(B)/libgosset_ref.so: $(REF) | $(B)
	$(CC) $(CFLAGS) -shared -fPIC -o $@ $(REF) $(LDLIBS)
$(B)/libgosset_masked.so: $(MASKED) | $(B)
	$(CC) $(CFLAGS) -shared -fPIC -o $@ $(MASKED) $(LDLIBS)
$(B)/xplat_kat: tests/xplat_kat.c $(FULL) | $(B)
	$(CC) $(CFLAGS) -o $@ tests/xplat_kat.c $(FULL) $(LDLIBS)
$(B)/demo: bench/demo.c $(CORE) src/dem.c | $(B)
	$(CC) $(CFLAGS) -o $@ bench/demo.c $(CORE) src/dem.c $(LDLIBS)

test: all
	python3 tests/run_tests.py

BENCH = $(B)/bench_kem $(B)/bench_maxsec $(B)/bench_ct $(B)/bench_fast $(B)/bench_hybrid $(B)/bench_mdecaps
bench: $(BENCH)
	$(B)/bench_kem && $(B)/bench_maxsec && $(B)/bench_ct && $(B)/bench_fast && $(B)/bench_hybrid
$(B)/bench_kem: bench/bench_kem.c $(CORE) | $(B)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)
$(B)/bench_maxsec: bench/bench_maxsec.c $(CORE) | $(B)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)
$(B)/bench_ct: bench/bench_ct.c src/e8fast.c src/e8codec.c | $(B)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)
$(B)/bench_fast: bench/bench_fast.c src/e8fast.c src/e8codec.c | $(B)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)
$(B)/bench_hybrid: bench/bench_hybrid.c $(FULL) | $(B)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)
$(B)/bench_mdecaps: bench/bench_mdecaps.c $(MASKED) | $(B)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -rf $(B)
.PHONY: all test bench clean
