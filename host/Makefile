# Copyright (c) Open Enclave SDK contributors.
# Licensed under the MIT License.

include ../config.mk

CFLAGS=$(shell pkg-config oehost-$(COMPILER) --cflags)
LDFLAGS=$(shell pkg-config oehost-$(COMPILER) --libs)
INCDIR=$(shell pkg-config oehost-$(COMPILER) --variable=includedir)

build:
	@ echo "Compilers used: $(CC), $(CXX)"
	oeedger8r ../fast_pca.edl --untrusted \
		--search-path $(INCDIR) \
		--search-path $(INCDIR)/openenclave/edl/sgx
	$(CXX) -g -c $(CFLAGS) host.cpp myEigenFunctions.cpp
	$(CC) -g -c $(CFLAGS) fast_pca_u.c
	$(CXX) -o fast_pca_host fast_pca_u.o host.o myEigenFunctions.o $(LDFLAGS)

clean:
	rm -f fast_pca_host host.o fast_pca_u.o fast_pca_u.c fast_pca_u.h fast_pca_args.h
