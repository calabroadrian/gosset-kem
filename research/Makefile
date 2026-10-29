# Builds the shared libraries the research scripts load (run `make` here, then run each script from its own folder,
# e.g. `cd v13 && python3 failure_boost.py`). The production code lives in ../src; these are the historical copies.
CC ?= gcc
F = -O3 -mavx2 -shared -fPIC
LIBS = v4/libe8codec.so v4/libe8fast_chk.so v5/libleech.so v7/libmlkem.so v7/libmlkem_opt.so \
       v9/libx25519.so v9/libhybrid.so v10/libe8masked.so v11/libmdecaps.so
all: $(LIBS)
v4/libe8codec.so: v4/e8codec.c
	$(CC) $(F) -o $@ $^
v4/libe8fast_chk.so: v4/e8fast.c v4/e8codec.c
	$(CC) $(F) -o $@ $^
v5/libleech.so: v5/leech.c
	$(CC) $(F) -o $@ $^
v7/libmlkem.so: v7/mlkem.c v7/fips202.c v4/e8fast.c v4/e8codec.c
	$(CC) $(F) -o $@ $^
v7/libmlkem_opt.so: v7/mlkem_opt.c v7/fips202.c v4/e8fast.c v4/e8codec.c
	$(CC) $(F) -o $@ $^
v9/libx25519.so: v9/x25519.c
	$(CC) $(F) -o $@ $^
v9/libhybrid.so: v9/hybrid.c v9/x25519.c v7/mlkem_opt.c v7/fips202.c v4/e8fast.c v4/e8codec.c
	$(CC) $(F) -o $@ $^
v10/libe8masked.so: v10/e8_masked.c v4/e8fast.c v4/e8codec.c
	$(CC) $(F) -o $@ $^ -lm
v11/libmdecaps.so: v11/masked_decaps.c v11/masked_core.c v11/masked_keccak.c v11/masked_e8bs.c v7/fips202.c v4/e8fast.c v4/e8codec.c
	$(CC) $(F) -o $@ $^ -lm
clean:
	rm -f $(LIBS)
