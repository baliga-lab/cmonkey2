#include "weeder_util.h"
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <memory.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

BOOL file_exists(const char *filename);
BOOL is_tail_node(BaseNode *node);

char ALPHABET[ALPHABET_SIZE]           = "ACGT";
const char *HTML_COLORS[ALPHABET_SIZE] = {"green", "blue", "gold", "red"};

int chartoint(char c)
{
  switch (tolower(c)) {
    case 'a': return DNA_A;
    case 'c': return DNA_C;
    case 'g': return DNA_G;
    case 't': return DNA_T;
    case '$': return 4;
    default: return -1;
  }
}

char revcomp(char c)
{
  char ch = toupper(c);

  switch (ch) {
  case 'A': return 'T';
  case 'C': return 'G';
  case 'G': return 'C';
  case 'T': return 'A';
  default: return ch;
  }
}

char inttochar(int num)
{
  if (num == DNA_A) return 'A';
  if (num == DNA_C) return 'C';
  if (num == DNA_G) return 'G';
  if (num == DNA_T) return 'T';
  /*
    WARNING: the original code does not include a default case.
    We return a '?' here which indicates something went wrong.
    Better than returning an undefined value.
   */
  fprintf(stderr, "*WARNING*: inttochar() received unsupported number: %d\n", num);
  return '?';
}

const char *reverse_pattern(const char *pattern, char *result)
{
  int len;
  int i;

  len = strlen(pattern);
  bzero(result, len + 1);

  for (i = 0; i < len; i++)
    result[i] = toupper(revcomp(pattern[len - i - 1]));
  return result;
}

BOOL hamming_distance_exceeds_threshold(const char *str1, const char *str2, int threshold)
{
  int i, len;
  int hamming_dist = 0;

  len = strlen(str1);

  for (i = 0; i < len; i++) {
    if (str1[i] != str2[i]) hamming_dist++;
    if (hamming_dist > threshold) break;
  }
  return hamming_dist > threshold ? 1 : 0;
}

BOOL inside(const char *str1, const char *str2, BOOL checkreverse) {
  int len1, len2, i, j;
  BOOL result = FALSE, match;

  len1 = strlen(str1);
  len2 = strlen(str2);

  if ((len2 - len1) != 2) return 0;

  for (i = 0; i < len2 - len1 + 1; i++) {
    match = TRUE;

    for (j = i; j < i + len1; j++) {
      if (str1[j - i] != str2[j]) {
        match = FALSE;
        break;
      }
    }
    if (match) result = TRUE;
  }
  if (checkreverse) {
    char *str1_reversed = malloc(sizeof(char) * (len1 + 1));
    BOOL rev_result;

    reverse_pattern(str1, str1_reversed);
    rev_result = inside(str1_reversed, str2, FALSE);
    free(str1_reversed);
    if (rev_result) result = TRUE;
  }
  return result;
}

BOOL overlap(const char *a, const char *b, BOOL checkreverse)
{
  int l;
  int i;
  int len;
  BOOL overlapping, result = FALSE;

  len = strlen(a);
  overlapping = TRUE;

  for (l = 1; l < 3; l++) {
    for (i = 0; i < len - l; i++) {
      if (a[i + l] != b[i]) {
	      overlapping = FALSE;
	      break;
	    }
    }
    if (overlapping) result = TRUE;
    overlapping = TRUE;

    for (i = 0; i < len - l; i++) {
      if (a[i] != b[i + l]) {
	      overlapping = FALSE;
	      break;
	    }
    }
    if (overlapping) result = TRUE;
  }

  if (checkreverse) {
    BOOL reverse_result;
    char *str1_reversed = malloc(sizeof(char) * (len + 1));
    reverse_pattern(a, str1_reversed);  
    reverse_result = overlap(str1_reversed, b, 0);
    if (reverse_result) result = TRUE;
    free(str1_reversed);
  }
  return result;
}

#define HAMMING_MAX 9999

int hamming_distance(const char *str1, const char *str2, BOOL checkrev)
{
  int dist_forward = 0, dist_reverse = HAMMING_MAX, i;
  int strlen1 = strlen(str1);

  if (strlen1 != strlen(str2)) return HAMMING_MAX;
  if (strcmp(str1, str2) == 0) return HAMMING_MAX; /* Note: Why equal strings return HAMMING_MAX ? */

  for (i = 0; i < strlen1; i++)
    if (str1[i] != str2[i]) dist_forward++;

  if (!checkrev) return dist_forward;
  if (checkrev) {
    char *str1_reversed = malloc(sizeof(char) * (strlen1 + 1));
    dist_reverse = 0;
    reverse_pattern(str1, str1_reversed);

    for (i = 0; i < strlen1; i++)
      if (str1_reversed[i] != str2[i]) dist_reverse++;
    free(str1_reversed);
  }
  return dist_reverse < dist_forward ? dist_reverse : dist_forward;
}

BOOL isvalid(char c)
{
  return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

extern const char *string_toupper(const char *str, char *result)
{
  int n = strlen(str);
  int i;
  for (i = 0; i < n; i++) {
    result[i] = toupper(str[i]);
  }
  result[i] = 0;
  return result;
}

BaseNode *merge_linkedlist(BaseNode *list1, BaseNode *list2,
                           BaseNode *tail);

BaseNode *mergesort_linkedlist(BaseNode *list_head, BaseNode *list_tail, int length)
{
  int i;
  int first, second;
  double mid;
  BaseNode *left, *right;
  BaseNode *result, *tmp;

  if (length == 1) {
    return list_head;
  } else {
    left  = list_head;
    right = list_head;

    mid = (double) (length) / 2;

    if (mid > floor(mid)) {
      first = (int) floor(mid);
      second = first + 1;
    } else {
      first = (int) mid;
      second = (int) mid;
    }

    /* Split the list */
    for (i = 0; i < first - 1; i++) right = right->next;
    tmp = right;
    right = right->next;
    tmp->next = list_tail;

    left   = mergesort_linkedlist(left, list_tail, first);
    right  = mergesort_linkedlist(right, list_tail, second);
    result = merge_linkedlist(left, right, list_tail);
    return result;
  }
}

/* Comparison function */
static BOOL node_greater_equals(BaseNode *node1, BaseNode *node2)
{
  return (node1->score > node2->score) ||
    ((node1->score == node2->score) && (node1->counter >= node2->counter));
}

BOOL is_tail_node(BaseNode *node)
{
  return abs(node->score - TAIL_SCORE) < EPS;
}

BaseNode *merge_linkedlist(BaseNode *list1, BaseNode *list2, BaseNode *tail)
{
  BaseNode *tmp, *first;
  int tail_node_reached = 0;

  if (node_greater_equals(list1, list2)) {
    first = list1;
    list1 = list1->next;
    tmp = first;
    tmp->next = tail;
  } else {
    first = list2;
    list2 = list2->next;
    tmp = first;
    tmp->next = tail;
  }

  while (!tail_node_reached) {
    if (node_greater_equals(list1, list2)) {

      if (is_tail_node(list1)) tail_node_reached = 1;
      else {
	      tmp->next = list1;
	      list1 = list1->next;
	      tmp = tmp->next;
	      tmp->next = tail;
	    }
    } else {
      tmp->next = list2;
      list2 = list2->next;
      tmp = tmp->next;
      tmp->next = tail;
    }
  }
  return first;
}

void print_node(BaseNode *node)
{
  const char *label = node->pat ? node->pat : "";
  fprintf(stderr, "list node [name: '%s', score = %f, counter = %d]\n", label,
          node->score, node->counter);
}

void print_list(BaseNode *head)
{
  if (head != NULL) {
    print_node(head);
    print_list(head->next);
  }
}

BOOL file_exists(const char *filename)
{
  return access(filename, R_OK) == 0 ? TRUE : FALSE;
}

void check_file_exists(const char *filename, const char *errormsg)
{
  if (!file_exists(filename)) {
    fprintf(stderr, errormsg, filename);
    exit(1);
  }
}

static int powers_base4[] = {
  POW_4_0, POW_4_1, POW_4_2,  POW_4_3,  POW_4_4,  POW_4_5,  POW_4_6,
  POW_4_7, POW_4_8, POW_4_9,  POW_4_10, POW_4_11, POW_4_12, POW_4_13,
  POW_4_14
};

int pow_4(int exponent)
{
  return powers_base4[exponent];
}
