# Copyright (c) Open Enclave SDK contributors.
# Licensed under the MIT License.

.PHONY: all build clean run simulate

OE_CRYPTO_LIB := mbedtls
export OE_CRYPTO_LIB

sim: 
	make build;
	make run

all: build

build:
	$(MAKE) -C enclave
	$(MAKE) -C host

clean:
	$(MAKE) -C enclave clean
	$(MAKE) -C host clean
	rm coffeeAI_eigenvalues.csv
	rm coffeeAI_principal_components.csv

run:
	host/fast_pca_host ./enclave/enclave.signed

simulate:
	host/fast_pca_host ./enclave/enclave.signed --simulate

diff_coffee:
	diff ../fast-pca/coffeeAI_eigenvalues.csv coffeeAI_eigenvalues.csv
	diff ../fast-pca/coffeeAI_principal_components.csv coffeeAI_principal_components.csv

diff_milk:
	diff ../fast-pca/milk_eigenvalues.csv milk_eigenvalues.csv
	diff ../fast-pca/milk_principal_components.csv milk_principal_components.csv