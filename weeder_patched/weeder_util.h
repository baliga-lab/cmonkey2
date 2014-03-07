#pragma once
#ifndef __WEEDER_UTIL_H__
#define __WEEDER_UTIL_H__

/*
 * Functions shared by locator.c, weederTFBS.c, adviser.c and weederlauncher.c
 * Refactored and fixed:
 * - hamming_distance()/hamming() in the original adviser.c does not really compute the
 *   reverse hamming distance
 * - revpattern() in the original code leaks memory
 * - generic, polymorphic mergesort on linked lists
 * - flag type variables are now declared as BOOL
 * - replaced all pow() calls with POW_4_N or pow_4(n), since we always use base 4
 *   and most of the time even know the exponent
 */

/* File names used by the weeder programs that can be customized here */
#define FORMAT_MIX_FILE       "%s.mix"
#define FORMAT_WEE_FILE       "%s.wee"
#define FORMAT_HTML_FILE      "%s.html"

/* format params are (inputfile, locatepattern) */
#define FORMAT_LOCATOR_WEE_FILE       "%s.%s.wee"
#define FORMAT_LOCATOR_HTML_FILE      "%s.%s.html"

#define FORMAT_8MER_FREQ_FILE "./FreqFiles/%s.8.freq"
#define FORMAT_6MER_FREQ_FILE "./FreqFiles/%s.6.freq"

typedef int BOOL;
#define TRUE  1
#define FALSE 0

#ifndef min
#define min(a, b) ((a < b) ? a : b)
#endif /* min */

/*
 * Zero-cost power calculation: When we know base and exponent, we can just
 * use the result at compile time.
 */
#define POW_4_0  1
#define POW_4_1  4
#define POW_4_2  16
#define POW_4_3  64
#define POW_4_4  256
#define POW_4_5  1024
#define POW_4_6  4096
#define POW_4_7  16384
#define POW_4_8  65536
#define POW_4_9  262144
#define POW_4_10 1048576
#define POW_4_11 4194304
#define POW_4_12 16777216
#define POW_4_13 67108864
#define POW_4_14 268435456

/*
 * Fast integer power function to the base 4, weeder uses it quite a bit.
 * Only calculates up to 4^14, which should be sufficient.
 */
extern int pow_4(int exponent);
extern int chartoint(char c);
extern char revcomp(char c);
extern char inttochar(int num);
extern const char *reverse_pattern(const char *pattern, char *result);
extern BOOL hamming_distance_exceeds_threshold(const char *str1, const char *str2, int threshold);
extern BOOL inside(const char *str1, const char *str2, BOOL checkreverse);
extern BOOL overlap(const char *str1, const char *str2, BOOL checkreverse);
extern int hamming_distance(const char *str1, const char *str2, BOOL checkreverse);
extern BOOL isvalid(char c);

extern const char *string_toupper(const char *str, char *result);
extern void check_file_exists(const char *filename, const char *errormsg);

#define TAIL_SCORE -111111.0
#define EPS 0.0001
#define MAX_PATTERN_LENGTH 15

/*
 * A generic merge sort for linked lists as used in weeder. It it based on a common
 * list node layout. Clients of the merge sort functions need to define their node
 * structures using the same memory layout as base_node and cast the respective
 * node pointers into base_node pointers.
 * See weederTFBS.c and adviser.c for details.
 */
typedef struct _base_node
{
  struct _base_node *next;
  double score;
  int    counter;
  char   pat[MAX_PATTERN_LENGTH];
} BaseNode;

void print_node(BaseNode *node);
void print_list(BaseNode *head);
BaseNode *mergesort_linkedlist(BaseNode *head, BaseNode *tail, int length);

/* Define an order of bases that works across the utilities */
#define DNA_A 0
#define DNA_C 1
#define DNA_G 2
#define DNA_T 3
#define ALPHABET_SIZE 4
extern char ALPHABET[ALPHABET_SIZE];

extern const char *HTML_COLORS[ALPHABET_SIZE];

#endif /* __WEEDER_UTIL_H__ */
