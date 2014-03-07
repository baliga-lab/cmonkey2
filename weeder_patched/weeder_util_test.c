#include "CuTest.h"
#include "weeder_util.h"
#include <stdio.h>
#include <string.h>

/* A learning test for constant definition */
void TestAlphabet(CuTest *tc)
{
  CuAssertIntEquals(tc, 'A', ALPHABET[0]);
  CuAssertIntEquals(tc, 'C', ALPHABET[1]);
  CuAssertIntEquals(tc, 'G', ALPHABET[2]);
  CuAssertIntEquals(tc, 'T', ALPHABET[3]);
}

void TestReversePattern(CuTest *tc)
{
  const char *input = "AGTC";
  char result[5];
  CuAssertStrEquals(tc, "GACT", reverse_pattern(input, result));
}

void TestMergeSortLengthOne(CuTest *tc)
{
  BaseNode node0;
  CuAssertPtrEquals(tc, &node0, mergesort_linkedlist(&node0, &node0, 1));
}

BaseNode *init_node(BaseNode *node,
                    const char *label, double score,
                    int counter, BaseNode *next)
{
  /* node->pat     = (char *) label; */
  strcpy(node->pat, label);
  node->score   = score;
  node->counter = counter;
  node->next    = next;
  return node;
}

void TestMergeSortLength2(CuTest *tc)
{
  BaseNode node0, node1, tail;
  BaseNode *sorted;

  init_node(&node0, "node0", 1.0, 1, &node1);
  init_node(&node1, "node1", 2.0, 1, &tail);
  init_node(&tail, "tail", TAIL_SCORE, 0, NULL);
  sorted = mergesort_linkedlist(&node0, &tail, 2);

  CuAssertPtrEquals(tc, &node1, sorted);
  CuAssertPtrEquals(tc, &node0, sorted->next);
  CuAssertPtrEquals(tc, &tail, sorted->next->next);
}

void TestMergeSortLength3(CuTest *tc)
{
  BaseNode node0, node1, node2, tail;
  BaseNode *sorted;

  init_node(&node0, "node0", 1.0, 1, &node1);
  init_node(&node1, "node1", 2.0, 1, &node2);
  init_node(&node2, "node2", 3.0, 1, &tail);
  init_node(&tail, "tail", TAIL_SCORE, 0, NULL);

  sorted = mergesort_linkedlist(&node0, &tail, 3);

  CuAssertPtrEquals(tc, &node2, sorted);
  CuAssertPtrEquals(tc, &node1, sorted->next);
  CuAssertPtrEquals(tc, &node0, sorted->next->next);
  CuAssertPtrEquals(tc, &tail, sorted->next->next->next);
}

typedef struct _ext_node {
  struct _ext_node *next;
  double score;
  int    counter;
  char   pat[MAX_PATTERN_LENGTH];

  int attr1;
  int attr2;
} ExtNode;

ExtNode *init_ext_node(ExtNode *node,
                               const char *label, double score,
                               int counter, ExtNode *next,
                               int attr1, int attr2)
{
  init_node((BaseNode *) node, label, score, counter, (BaseNode *) next);
  node->attr1 = attr1;
  node->attr2 = attr2;
  return node;
}

/*
  A unit test that demonstrates how to use merge sort on nodes that extend on
  base_node's layout.
*/
void TestMergeSortExtendedNode(CuTest *tc)
{
  ExtNode node0, node1, node2, tail;
  ExtNode *sorted;

  init_ext_node(&node0, "node0", 1.0, 1, &node1, 42, 43);
  init_ext_node(&node1, "node1", 2.0, 1, &node2, 45, 46);
  init_ext_node(&node2, "node2", 3.0, 1, &tail,  50, 51);
  init_ext_node(&tail, "tail", TAIL_SCORE, 0, NULL, 0, 0);

  sorted = (ExtNode *) mergesort_linkedlist((BaseNode *) &node0,
                                            (BaseNode *) &tail, 3);

  CuAssertPtrEquals(tc, &node2, sorted);
  CuAssertPtrEquals(tc, &node1, sorted->next);
  CuAssertPtrEquals(tc, &node0, sorted->next->next);
  CuAssertPtrEquals(tc, &tail, sorted->next->next->next);
}

void TestStringToUpper(CuTest *tc)
{
  const char *input    = "SomeString";
  char buffer[strlen(input) + 1];
  const char *result = string_toupper(input, buffer);
  CuAssertStrEquals(tc, "SOMESTRING", result);
}

void TestStringToUpperInputIsOutput(CuTest *tc)
{
  const char *input    = "SomeString";
  char buffer[strlen(input) + 1];
  strcpy(buffer, input);
  string_toupper(buffer, buffer);
  CuAssertStrEquals(tc, "SOMESTRING", buffer);
}

/*
  This is a learning test, to demonstrate that we can simply replace several
  array allocations with literals, which is cleaner and leaves less memory holes.
*/
void TestColorArray(CuTest *tc)
{
  const char *color[] = {"green", "blue", "red", "gold"};

  CuAssertIntEquals(tc, 4 * sizeof(const char *), sizeof(color));
  CuAssertStrEquals(tc, "green", color[0]);
  CuAssertStrEquals(tc, "blue",  color[1]);
  CuAssertStrEquals(tc, "red",   color[2]);
  CuAssertStrEquals(tc, "gold",  color[3]);
}

void Testpow_4(CuTest *tc)
{
  CuAssertIntEquals(tc, 1,  pow_4(0));
  CuAssertIntEquals(tc, 4,  pow_4(1));
  CuAssertIntEquals(tc, 16, pow_4(2));
  CuAssertIntEquals(tc, 64, pow_4(3));
}

CuSuite *WeederUtilGetSuite()
{
  CuSuite *suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, TestAlphabet);
  SUITE_ADD_TEST(suite, TestReversePattern);
  SUITE_ADD_TEST(suite, TestStringToUpper);
  SUITE_ADD_TEST(suite, TestStringToUpperInputIsOutput);
  SUITE_ADD_TEST(suite, TestColorArray);

  SUITE_ADD_TEST(suite, TestMergeSortLengthOne);
  SUITE_ADD_TEST(suite, TestMergeSortLength2);
  SUITE_ADD_TEST(suite, TestMergeSortLength3);
  SUITE_ADD_TEST(suite, TestMergeSortExtendedNode);

  SUITE_ADD_TEST(suite, Testpow_4);
  return suite;
}

int main(int argc, char **argv)
{
  CuString *output = CuStringNew();
  CuSuite *suite = WeederUtilGetSuite();
  CuSuiteRun(suite);
  CuSuiteSummary(suite, output);
  CuSuiteDetails(suite, output);
  printf("%s\n", output->buffer);
  return 0;
}
