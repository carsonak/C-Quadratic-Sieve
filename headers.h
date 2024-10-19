#ifndef FAC_HEADERS
#define FAC_HEADERS

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <assert.h>
#include <unistd.h>
#include <sys/time.h>

// Quadratic sieve integers.
typedef uint32_t qs_sm; // small size, like a factor base prime number (32-bit)
typedef uint64_t qs_md; // medium size, like a factor base prime number squared (64-bit)
typedef int64_t qs_tmp; // signed type to perform intermediates computations

// 64-bit factorization integers.
typedef unsigned long long int u64;

// 64-bit factorization structure

typedef struct {
	qs_md prime;
	int power;
} fac64_row;

// Factorization manager's structure

typedef struct {
	struct {
		int verbose;
		qs_md force;
		qs_md help;
		struct {
			qs_md custom;
			qs_md seed;
		} rand;
		qs_md generate[3];
		qs_md timeout;
		char *input_file;
		char *output_file;
		char output_format;
		qs_md qs_multiplier;
		qs_md qs_base_size;
		qs_md qs_sieve;
		qs_md qs_error_bits;
		qs_md qs_threshold;
		qs_md qs_alloc_mb;
		qs_md qs_large_prime;
		qs_md qs_poly;
		qs_md qs_laziness;
	} params;
	qs_md duration_ms ;
	FILE *in;
	FILE *out;
	int code;
	struct {
		size_t row_idx;
		size_t total_rows;
		size_t max_factors;
		size_t max_digits;
		size_t max_bits;
	} scale;
	struct {
		cint num;
		cint tmp[10];
		cint_sheet *sheet;
		qs_md trial_start;
		qs_md duration_ms;
		int power;
		char *input_string;
		char *output_string;
		struct {
			cint num;
			int power;
			int prime;
		} *res;
		struct {
			void *base;
			void *now;
		} mem;
	} session;
} state;

// Quadratic sieve structures

struct qs_relation {
	qs_sm id; // definitive relations have a non-zero id.
	cint *X;
	struct {
		qs_sm *data;
		qs_sm length;
	} Y;
	union {
		struct {
			qs_sm *data;
			qs_sm length;
		} Z;
		struct qs_relation *next;
	} axis;
	// axis :
	// - "Z" field is used by definitive relations.
	// - "next" is used by data that wait to be paired, it uses a linked list instead of a "Z" field.
};

typedef struct {

	// computation sheet
	cint_sheet *sheet;

	// numbers that are updated
	struct {
		cint N;
		cint FACTOR;
		cint X;
		cint KEY;
		cint VALUE;
		cint CYCLE;
		cint TEMP[5];
		cint MY[5]; // available for developers
	} vars;

	// polynomial vars
	struct {
		cint A;
		cint B;
		cint C;
		cint D;
		qs_sm d_bits;
		qs_sm offset;
		struct {
			qs_sm x_1;
			qs_sm half ;
			qs_sm x_2 ;
			qs_sm x_3 ;
		} span;
		qs_sm gray_max;
		qs_sm curves;
	} poly;

	// constants
	struct {
		cint kN;
		cint ONE;
		cint TOO_LARGE_PRIME;
		cint MULTIPLIER;
		cint M_HALF;
	} constants;

	// system
	struct {
		qs_sm bytes_allocated;
		void *base;
		void *now;
	} mem;

	// parameters and miscellaneous vars
	struct{
		qs_md start;
		qs_md end;
		qs_md now;
		qs_sm tick;
		char remaining[127];
	} time;
	qs_md seed;
	qs_md adjustor;
	qs_sm multiplier;
	qs_sm n_bits;
	qs_sm kn_bits;
	struct {
		uint8_t **positions[2];
		uint8_t *sieve;
		uint8_t *flags;
		qs_sm length;
		qs_sm length_half;
		qs_sm cache_size;
	} m;
	qs_sm iterative_list[5];
	qs_sm error_bits;
	struct {
		qs_sm value;
	} threshold;
	qs_sm sieve_again_perms;

	// useful data sharing same length
	struct {
		struct {
			qs_sm num;
			qs_sm size;
			qs_sm sqrt_kN_mod_prime;
			qs_sm root[2];
		} *data;
		qs_sm largest;
		qs_sm length;
	} base;

	// useful data sharing same length
	struct {
		qs_sm *A_indexes;
		struct {
			cint B_terms;
			qs_sm *A_inv_double_value_B_terms;
			qs_sm A_over_prime_mod_prime;
			qs_sm prime_index;
			qs_md prime_squared;
		} *data;
		struct {
			qs_sm defined;
			qs_sm subtract_one;
			qs_sm double_value;
		} values;
	} s;

	qs_sm *buffer[2]; // proportional to "length of factor base" (medium or large)

	// uniqueness trees : [ relations, cycle finder, divisors of N, ]
	struct avl_manager uniqueness[3];

	// data collection made by algorithm
	struct {
		struct qs_relation **data;
		struct {
			qs_sm prev;
			qs_sm now;
			qs_sm peak;
			qs_sm needs;
			qs_sm capacity;
			qs_sm by_partial;
		} length;
		qs_md too_large_prime;
	} relations;

	// pointers to the divisors of N are kept in a flat array
	struct {
		qs_sm processing_index;
		qs_sm total_primes;
		qs_sm length;
		cint **data;
	} divisors;

	// Lanczos has its own struct
	struct {
		qs_sm safe_length;
		struct {
			struct qs_relation *relation;
			qs_sm y_length;
		} *snapshot;
	} lanczos;

	state *state;

} qs_sheet;

// Factorization manager and file i/o utilities.
static void print_help(const char *);
static inline int cli_param_match(const char *, const char *, const char *);
static inline qs_md get_num(char *);
static int read_key_value(char **, state *);
static int read_flags(char **, state *);
static void generate_input_file(state *state);
static void output_csv(state *, int, int);
static void output_default(state *, int, int);
static void output_json_compact(state *, int, int);
static void output_json_pretty_print(state *, int, int);
static inline void display_progress(const char *, double);
static inline void output(state *);
static int validate_input_file(state *state);
static size_t prepare_file_descriptors(state *);
static int validate_string_number(const char *, state *);

// Misc utilities (including the factorization manager).
static inline void simple_inline_cint(cint *, size_t, void **);
static inline void simple_dup_cint(cint *, const cint *, void **);
static inline void simple_int_to_cint(cint *, qs_md);
static inline qs_md simple_cint_to_int(const cint *);
static struct avl_node *avl_cint_inserter(void *, const void *);
static inline void *mem_aligned(void *);
static qs_md get_time_ms();
static void prepare_sessions(state *state);
static void erase_session(state *state);
static void clear_sessions(state *state);
static void manager_add_factor(state *state, cint *num, int pow, int is_prime);
static inline void manager_add_simple_factor(state *state, qs_md num, int pow, int is_prime);
static void factorization_64_bits(state *);
static int factorization_trial_division(state *, int, int);
static int factorization_any_root_checker(state *, const cint *, cint *, cint *);
static int factorization_perfect_power_checker(state *, int);
static int factorization_prime_number_checker(state *, int);
static int factorization_give_up(state *, int);
static void factor(state *);
static void process_single(state *);
static void process_multi(int, char **, state *);

// Math functions and system utilities.
static qs_md xor_random(qs_md *);
static inline qs_md xor_rand(qs_md *, qs_md, qs_md);
static int is_prime_4669913(qs_sm);
static double log_computation(double);
static inline qs_sm multiplication_modulo(qs_md, qs_md, qs_sm);
static inline qs_sm power_modulo(qs_md, qs_md, qs_sm);
static int kronecker_symbol(qs_sm, qs_sm);
static qs_sm tonelli_shanks(qs_sm, qs_sm);
static qs_sm modular_inverse(qs_sm, qs_sm);

// 64-bit factorization.
static u64 mul_mod(u64, u64, u64);
static u64 pow_mod(u64, u64, u64);
static int is_prime_64_bits(u64);
static u64 pollard_rho(u64, uint64_t *);
static u64 nth_root(u64, u64);
static inline u64 square_extraction(u64 *, int *);
static void fac_64_worker(state *, u64, fac64_row *);

// Quadratic sieve.
static int factorization_quadratic_sieve(state *, int);
static int inner_continuation_condition(qs_sheet *);
static int outer_continuation_condition(qs_sheet *);
static void qs_parametrize(qs_sheet *);
static qs_sm linear_param_resolution(const double [][2], qs_sm);
static void preparation_part_1(qs_sheet *, state *, int);
static void preparation_part_2(qs_sheet *);
static void preparation_part_3(qs_sheet *);
static void preparation_part_3_default_proposition(qs_sheet *qs, qs_sm *, size_t);
static void preparation_part_3_alternative_proposition(qs_sheet *qs, qs_sm *, size_t);
static void preparation_part_4(qs_sheet *);
static void preparation_part_5(qs_sheet *);
static void preparation_part_6(qs_sheet *);
static void get_started_iteration(qs_sheet *);
static void iteration_part_1(qs_sheet *, const cint *, cint *);
static void iteration_part_2(qs_sheet *, const cint *, cint *);
static void iteration_part_3(qs_sheet *, const cint *, const cint *);
static qs_sm iteration_part_4(const qs_sheet *, qs_sm, qs_sm **, cint *);
static void iteration_part_5(qs_sheet *, const cint *, const cint *);
static void iteration_part_6(qs_sheet *, const cint *, const cint *, const cint *, cint *);
static void iteration_part_7(qs_sheet *, qs_sm, const qs_sm *restrict);
static void iteration_part_8(qs_sheet *, qs_sm, const qs_sm *);
static int qs_register_divisor(qs_sheet *qs);
static void register_relations(qs_sheet *, const cint *, const cint *, const cint *);
static void register_regular_relation(qs_sheet *, const cint *, const qs_sm *const restrict[4]);
static void register_partial_relation(qs_sheet *, const cint *, const cint *, const qs_sm *const restrict[4]);
static void finalization_part_1(qs_sheet *, const uint64_t *restrict);
static int finalization_part_2(qs_sheet *);

// Linear algebra with Block Lanczos algorithm.
static void lanczos_mul_MxN_Nx64(const qs_sheet *, const uint64_t *, qs_sm, uint64_t *);
static void lanczos_mul_trans_MxN_Nx64(const qs_sheet *, const uint64_t *, uint64_t *);
static void lanczos_mul_64xN_Nx64(const qs_sheet *, const uint64_t *, const uint64_t *, uint64_t *, uint64_t *);
static uint64_t lanczos_find_non_singular_sub(const uint64_t *, const uint64_t *, uint64_t *, uint64_t, uint64_t *);
static void lanczos_mul_Nx64_64x64_acc(qs_sheet *, const uint64_t *, const uint64_t *, uint64_t *, uint64_t *);
static void lanczos_mul_64x64_64x64(const uint64_t *, const uint64_t *, uint64_t *);
static void lanczos_transpose_vector(qs_sheet *, const uint64_t *, uint64_t **);
static void lanczos_combine_cols(qs_sheet *, uint64_t *, uint64_t *, uint64_t *, uint64_t *);
static void lanczos_build_array(qs_sheet *, uint64_t **, size_t, size_t);
static uint64_t *lanczos_block_worker(qs_sheet *);
static void lanczos_reduce_matrix(qs_sheet *);
static uint64_t *lanczos_block(qs_sheet *);

// Verbose Level 0: Just need a factorization of N, no debug messages.
// Verbose Level 1: Also print the factorization progress in percentage.
// Verbose Level 2: Also show critical issues, should usually show nothing.
// Verbose Level 3: Also show detailed information, for in-depth analysis.
// Verbose Level 4: Also show intermediate mathematical values (divisors of N).

#define DEBUG_CRITICAL(fmt, ...) debug_print(qs->state, 1, fmt, __VA_ARGS__)
#define DEBUG_NORMAL(fmt, ...) debug_print(qs->state, 2, fmt, __VA_ARGS__)
#define DEBUG_MATH(fmt, ...) debug_print(qs->state, 3, fmt, __VA_ARGS__)

#endif //FAC_HEADERS