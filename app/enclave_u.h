#ifndef ENCLAVE_U_H__
#define ENCLAVE_U_H__

#include <stdint.h>
#include <wchar.h>
#include <stddef.h>
#include <string.h>
#include "sgx_edger8r.h" /* for sgx_status_t etc. */

#include "myEigenFunctions.h"

#include <stdlib.h> /* for size_t */

#define SGX_CAST(type, item) ((type)(item))

#ifdef __cplusplus
extern "C" {
#endif

#ifndef _struct_foo_t
#define _struct_foo_t
typedef struct struct_foo_t {
	uint32_t struct_foo_0;
	uint64_t struct_foo_1;
} struct_foo_t;
#endif

typedef enum enum_foo_t {
	ENUM_FOO_0 = 0,
	ENUM_FOO_1 = 1,
} enum_foo_t;

#ifndef _union_foo_t
#define _union_foo_t
typedef union union_foo_t {
	uint32_t union_foo_0;
	uint32_t union_foo_1;
	uint64_t union_foo_3;
} union_foo_t;
#endif

#ifndef OCALL_POINTER_USER_CHECK_DEFINED__
#define OCALL_POINTER_USER_CHECK_DEFINED__
void SGX_UBRIDGE(SGX_NOCONVENTION, ocall_pointer_user_check, (int* val));
#endif
#ifndef OCALL_POINTER_IN_DEFINED__
#define OCALL_POINTER_IN_DEFINED__
void SGX_UBRIDGE(SGX_NOCONVENTION, ocall_pointer_in, (int* val));
#endif
#ifndef OCALL_POINTER_OUT_DEFINED__
#define OCALL_POINTER_OUT_DEFINED__
void SGX_UBRIDGE(SGX_NOCONVENTION, ocall_pointer_out, (int* val));
#endif
#ifndef OCALL_POINTER_IN_OUT_DEFINED__
#define OCALL_POINTER_IN_OUT_DEFINED__
void SGX_UBRIDGE(SGX_NOCONVENTION, ocall_pointer_in_out, (int* val));
#endif
#ifndef OCALL_FUNCTION_ALLOW_DEFINED__
#define OCALL_FUNCTION_ALLOW_DEFINED__
void SGX_UBRIDGE(SGX_NOCONVENTION, ocall_function_allow, (void));
#endif

sgx_status_t ecall_DotProduct_aa(sgx_enclave_id_t eid, double* retval, double** A, double** B, int n, int m, size_t len);
sgx_status_t ecall_DotProduct_vv(sgx_enclave_id_t eid, double* retval, double* A, double* B, int n, size_t len);
sgx_status_t ecall_DotProduct_av(sgx_enclave_id_t eid, double** retval, double** A, double* v, int n, size_t len1, size_t len2);
sgx_status_t ecall_CenterMatrix(sgx_enclave_id_t eid, double*** retval, double** A, int n, int m, size_t len);
sgx_status_t ecall_CovarianceMatrix(sgx_enclave_id_t eid, double*** retval, double** A, int n, int m, size_t len);
sgx_status_t ecall_power_method(sgx_enclave_id_t eid, Eigenpair* retval, double** A, double* v, int n, double tol, size_t len1, size_t len2);
sgx_status_t ecall_type_char(sgx_enclave_id_t eid, char val);
sgx_status_t ecall_type_int(sgx_enclave_id_t eid, int val);
sgx_status_t ecall_type_float(sgx_enclave_id_t eid, float val);
sgx_status_t ecall_type_double(sgx_enclave_id_t eid, double val);
sgx_status_t ecall_type_size_t(sgx_enclave_id_t eid, size_t val);
sgx_status_t ecall_type_wchar_t(sgx_enclave_id_t eid, wchar_t val);
sgx_status_t ecall_type_struct(sgx_enclave_id_t eid, struct struct_foo_t val);
sgx_status_t ecall_type_enum_union(sgx_enclave_id_t eid, enum enum_foo_t val1, union union_foo_t* val2);
sgx_status_t ecall_pointer_user_check(sgx_enclave_id_t eid, size_t* retval, void* val, size_t sz);
sgx_status_t ecall_pointer_in(sgx_enclave_id_t eid, int* val);
sgx_status_t ecall_pointer_out(sgx_enclave_id_t eid, int* val);
sgx_status_t ecall_pointer_in_out(sgx_enclave_id_t eid, int* val);
sgx_status_t ecall_pointer_string(sgx_enclave_id_t eid, char* str);
sgx_status_t ecall_pointer_string_const(sgx_enclave_id_t eid, const char* str);
sgx_status_t ecall_pointer_size(sgx_enclave_id_t eid, void* ptr, size_t len);
sgx_status_t ecall_pointer_count(sgx_enclave_id_t eid, int* arr, size_t cnt);
sgx_status_t ecall_pointer_isptr_readonly(sgx_enclave_id_t eid, buffer_t buf, size_t len);
sgx_status_t ocall_pointer_attr(sgx_enclave_id_t eid);
sgx_status_t ecall_array_user_check(sgx_enclave_id_t eid, int arr[4]);
sgx_status_t ecall_array_in(sgx_enclave_id_t eid, int arr[4]);
sgx_status_t ecall_array_out(sgx_enclave_id_t eid, int arr[4]);
sgx_status_t ecall_array_in_out(sgx_enclave_id_t eid, int arr[4]);
sgx_status_t ecall_array_isary(sgx_enclave_id_t eid, array_t arr);
sgx_status_t ecall_function_public(sgx_enclave_id_t eid);
sgx_status_t ecall_function_private(sgx_enclave_id_t eid, int* retval);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
