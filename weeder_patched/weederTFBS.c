/***************************************/
/*      Weeder 1.4      
        see LICENSE.txt file for
        terms and conditions of use
*/
#define OUTFILE_MAX_NAME_LENGTH 800

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <math.h>
#include <signal.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>

#include "weeder_util.h"

/* n(permutations) = |{A, C, G, T}| ^ k, k in { 6, 8, 10} */
#define NUM_6MER_PERMUTATIONS  POW_4_6
#define NUM_8MER_PERMUTATIONS  POW_4_8
#define NUM_10MER_PERMUTATIONS POW_4_10

#define DEFAULT_MAX_RESULTS    10
#define SEQUENCE_NAME_LENGTH   20
typedef struct _pattern
{
  /* base_node elements */
  struct _pattern *next;
  double score;
  int    counter;
  char   pat[MAX_PATTERN_LENGTH];

  /* weederTFBS specific elements */
  double sig;
  unsigned short int searchcount;
} Pattern;

typedef struct
{
  int *besthere, *inbestseq;
  char **sequence, **name;
  BOOL **used;
} WorkBuffers;

/* Module-global variables */
/* head and tail node do not change, so we do not need to malloc() them. */
static Pattern head_node, tail_node;

static int  extfile = 0;
static int *sequence_lengths;
static double sixmers[NUM_6MER_PERMUTATIONS];
static double eightmers[NUM_8MER_PERMUTATIONS];
static double tenmers[NUM_10MER_PERMUTATIONS];

static int num_sequences_total = 0, sum_all_seq_lengths = 0;
static double repeat_percentage = -1;

static time_t starttime;

/* Options derived from command line arguments */
static struct
{
  char infile[128];
  char organism[80];
  BOOL verbose, reversed, onlybestflag;
  int min_repeat_percentage, max_mismatches, motif_length, max_results;
} options;

/* Module  functions */
void init_options(void);
Pattern *mergesrt(Pattern *list, int length);
void   readfreqs(void);
double compfreqs(char *pat, char *pos, int err);
double exactfreqs(char *pat, int errs);
double newtwmers(char *pat, char *pos, int err);
double tentotwelve(char *pat);
double exactfreqswithcheck(char *pat, int errs, char *check);
double compfreqswithcheck(char *pat, char *pos, int err, char *check);
double newtwmerswithcheck(char *inpat, char *pos, int err, char *check);
void   determine_num_sequences_and_lengths(char *sequence_buffer);
void   init_workbuffers(WorkBuffers *);
void   read_sequence_data(WorkBuffers *, char *sequence_buffer);
double compsix(char *pat, char *pos, int err);
int pattolong(const char *pat);

/* Error checking */
BOOL checkmotifanderror(int motif_length, int errs);
int parse_option(int argc, char *argv[], int option_index);
void print_usage(const char *progname);
BOOL has_option_value(int argc, char *argv[], int option_index);
void print_option_error(const char *progname);
void check_condition(BOOL condition, const char *errormsg, int exitcode);

/* Reporting */
void report_best_results(FILE *wee_fp, FILE *html_fp, FILE *mix_fp, Pattern *seekpat);
void report_time_taken(FILE *wee_fp, FILE *html_fp);

int main(int argc, char *argv[])
{
  int uu, start, end, maah, lastocc, bestsofar,
    bestpos, lastinseq = -1, howmanyseq = 0, howmany = 0, inbest = 0,
    resultlist_length = 0, localerror, seq_index, base_index, arg_index;
  BOOL flagover; /* Note: this flag is calculated but never read */
  double patfreqs[5], myscore, *besterr;
  char revseek[15], sequence_buffer[MAXSEQ];
  Pattern *seekpat;
  WorkBuffers buffers;

  char weefile[OUTFILE_MAX_NAME_LENGTH]    = "";
  char mixfile[OUTFILE_MAX_NAME_LENGTH]    = "";
  char htmlfile[OUTFILE_MAX_NAME_LENGTH] = "";
  FILE *wee_fp, *mix_fp, *html_fp;
 
  if (argc < 4) {
    print_usage(argv[0]);
    exit(1);
  }

  init_options();
  arg_index = 1;
  while (arg_index < argc) {
    arg_index = parse_option(argc, argv, arg_index);
  }
  check_condition(extfile != 0, "\nMissing input filename\n", 1);
  check_file_exists(options.infile, "\nNo such file : %s\n");
  check_condition(options.motif_length != 0, "\nIllegal motif length\n", 1);
  check_condition(options.max_mismatches >= 0, "\nMissing error number -e\n", 1);
  check_condition(repeat_percentage >= 0, "\nUndefined repeat percentage -R\n", 1);

  if (!checkmotifanderror(options.motif_length, options.max_mismatches)) {
    fprintf(stderr, "\nIllegal combination of length %d and mismatches %d\n",
             options.motif_length, options.max_mismatches);
      exit(1);
  }

  readfreqs();
  determine_num_sequences_and_lengths(sequence_buffer);

  /*
    now that we know how many sequences we have, we can allocate
    the buffers and re-read the input file
   */
  init_workbuffers(&buffers);
  read_sequence_data(&buffers, sequence_buffer);

  besterr = (double *) malloc((options.max_mismatches + 1) * sizeof(double));
  repeat_percentage = (repeat_percentage * num_sequences_total / 100.0);
  options.min_repeat_percentage = ceil(repeat_percentage);

  sprintf(mixfile,  FORMAT_MIX_FILE,  options.infile);
  sprintf(weefile,  FORMAT_WEE_FILE,  options.infile);
  sprintf(htmlfile, FORMAT_HTML_FILE, options.infile);

  wee_fp  = fopen(weefile, "a");
  mix_fp  = fopen(mixfile, "a");
  html_fp = fopen(htmlfile, "a");

  fprintf(stderr,
          "\nSearching for motifs of length %d with %d mutations in file %s.....\n",
          options.motif_length, options.max_mismatches, options.infile);
  fprintf(wee_fp, "\nSearching for motifs of length %d with %d mutations.....\n",
          options.motif_length, options.max_mismatches);
  fprintf(html_fp,
          "<br>Searching for motifs of length %d with %d mutations.....<br><br>",
          options.motif_length, options.max_mismatches);

  if (options.verbose)
    fprintf(stderr, "\nMinimum sequence repeats : %d out of %d\n",
            options.min_repeat_percentage,
            num_sequences_total);

  fprintf(wee_fp,  "(Actual commandline is : \n");
  fprintf(html_fp, "(Actual commandline is : <br>");

  for (arg_index = 0; arg_index < argc; arg_index++) {
    fprintf(wee_fp,  "%s ", argv[arg_index]);
    fprintf(html_fp, "%s ", argv[arg_index]);  
  }

  fprintf(wee_fp,  ")\n");
  fprintf(html_fp, ")<br>");

  for (seq_index = 0; seq_index < num_sequences_total; seq_index++) {
    for (base_index = 0; base_index < sequence_lengths[seq_index]; base_index++) {
      buffers.sequence[seq_index][base_index] =
        toupper(buffers.sequence[seq_index][base_index]);
      
      if (!isvalid(buffers.sequence[seq_index][base_index])) {
        int used_base_index;

        if (base_index - options.motif_length + 1 < 0) start = 0;
        else start = base_index - options.motif_length + 1;
        end = base_index;

        for (used_base_index = start; used_base_index <= end; used_base_index++)
          buffers.used[seq_index][used_base_index] = TRUE;
      }
    }
  }
  starttime = clock() / CLOCKS_PER_SEC;
  tail_node.score = TAIL_SCORE;

  head_node.next = &tail_node;

  seq_index = 0;
  while (seq_index < num_sequences_total) {
    if (options.verbose)
      fprintf(stderr, "Processing sequence %s\n", buffers.name[seq_index]);

    for (base_index = 0; base_index < sequence_lengths[seq_index] - options.motif_length;
         base_index++) {

      if (options.verbose && ((base_index % 100) == 0))
        fprintf(stderr, "Position %d\n", base_index);

      if (!buffers.used[seq_index][base_index]) {
        int motif_base_index;

	      buffers.used[seq_index][base_index] = TRUE;

	      seekpat = (Pattern *) malloc(sizeof(Pattern));
	      bzero(seekpat->pat, options.motif_length + 1);
	      seekpat->counter = seq_index;
	      seekpat->searchcount = base_index;
	      seekpat->sig = 1.0;

	      resultlist_length++;

	      for (motif_base_index = 0; motif_base_index < options.motif_length;
             motif_base_index++) {
          seekpat->pat[motif_base_index] =
            toupper(buffers.sequence[seq_index][base_index + motif_base_index]);
        }

	      seekpat->next = head_node.next;
	      head_node.next = seekpat;
        
	      seekpat = &head_node;
	      myscore = 0.0;

	      howmany = 0;
	      howmanyseq = 0;

	      if (options.reversed) {
          int seek_pattern_length_plus_one;
          bzero(revseek, 15);

          seek_pattern_length_plus_one = strlen(seekpat->next->pat) + 1;

          for (uu = 0; uu <= seek_pattern_length_plus_one - 2; uu++) {
            revseek[uu] = revcomp(seekpat->next->pat[seek_pattern_length_plus_one - 2 - uu]);
          }
        }

        /* QUI PER SINGOLO STRAND */
	      if (!options.reversed) {

          for (uu = 0; uu <= options.max_mismatches; uu++) {
            besterr[uu] = -1.0;
          }

          lastinseq = -1;

          for (uu = 0; uu < num_sequences_total; uu++) {
            buffers.inbestseq[uu] = -1;
            buffers.besthere[uu] = -1;
            inbest = 0;
            lastocc = -1;
            bestsofar = 1000;
            bestpos = -1;
            inbest = 0;
            localerror = options.max_mismatches;
            maah = sequence_lengths[uu] - strlen(seekpat->next->pat) + 1;

            if ((maah > 0) &&
                (howmanyseq + num_sequences_total - uu >= options.min_repeat_percentage)) {
              int j, num_mismatches;

              for (j = 0; j < maah; j++) {
                num_mismatches = 0;

                for (motif_base_index = 0; motif_base_index < options.motif_length;
                     motif_base_index++) {
                  if (toupper(buffers.sequence[uu][j + motif_base_index]) !=
                      seekpat->next->pat[motif_base_index]) {
                    num_mismatches++;
                  }
                  if (num_mismatches > localerror) break;
                }
                if (num_mismatches <= localerror) {
                  if (num_mismatches == 0) buffers.used[uu][j] = TRUE;
          
                  howmany++;
                  lastocc = j;

                  if ((num_mismatches == bestsofar) && (lastinseq == uu)) {
                    inbest++;
                  }
                  if (lastinseq != uu) {
                    lastinseq = uu;
                    howmanyseq++;
                    bestpos = j;
                    bestsofar = num_mismatches;
                    inbest = 1;
                  }

                  if (num_mismatches < bestsofar) {
                    bestpos = j;
                    bestsofar = num_mismatches;
                    inbest = 1;

                    if (options.onlybestflag) localerror = bestsofar;
                  }
                  num_mismatches = 0;
                }
              }
            }

            if (bestpos > -1) {
              if (besterr[bestsofar] == -1) besterr[bestsofar] = 1;
              buffers.besthere[uu] = bestsofar;
              buffers.inbestseq[uu] = inbest;
            }
          }

          if (howmanyseq >= options.min_repeat_percentage) {
            patfreqs[0] = exactfreqs(seekpat->next->pat, 0);

            for (uu = 1; uu <= options.max_mismatches; uu++) {
              if (options.motif_length != 12)
                patfreqs[uu] = exactfreqs(seekpat->next->pat, uu);
              else
                patfreqs[uu] = exactfreqs(seekpat->next->pat, uu) + patfreqs[uu - 1];
            }

            seekpat->next->score = 0.0;
            myscore = 0.0;

            for (uu = 0; uu < num_sequences_total; uu++) {
              double tmpfloat;
              if (buffers.besthere[uu] >= 0) {
                tmpfloat =
                  log((double) buffers.inbestseq[uu]) -
                  (log(patfreqs[buffers.besthere[uu]]) +
                   log((double) sequence_lengths[uu]));
              } else {
                tmpfloat = 0;
              }
              myscore += tmpfloat;
            }
            seekpat->next->score = myscore;

            if (!options.onlybestflag) 
              seekpat->next->score += log((double) howmany) -
                (log(patfreqs[options.max_mismatches]) + log((double) sum_all_seq_lengths));
          } else {
            seekpat->next->score = -1.0;
          }
        } else {
          for (uu = 0; uu <= options.max_mismatches; uu++) {
            besterr[uu] = -1.0;
          }
          howmany = 0;
          howmanyseq = 0;
          lastinseq = -1;
          myscore = 0.0;

          for (uu = 0; uu < num_sequences_total; uu++) {
            buffers.inbestseq[uu] = -1;
            buffers.besthere[uu] = -1;
            bestsofar = 1000;
            bestpos = -1;
            maah = sequence_lengths[uu] - options.motif_length + 1;
            localerror = options.max_mismatches;

            if ((maah > 0) &&
                (howmanyseq + num_sequences_total - uu >= options.min_repeat_percentage)) {
              int j, num_mismatches, num_mismatches_rev;
              for (j = 0; j < maah; j++) {
                num_mismatches     = 0;
                num_mismatches_rev = 0;

                for (motif_base_index = 0; motif_base_index < options.motif_length;
                     motif_base_index++) {
                  if (buffers.sequence[uu][j + motif_base_index] !=
                      seekpat->next->pat[motif_base_index]) {
                    num_mismatches++;
                  }
                  if (buffers.sequence[uu][j + motif_base_index] !=
                      revseek[motif_base_index]) {
                    num_mismatches_rev++;
                  }
                  if ((num_mismatches > localerror) && (num_mismatches_rev > localerror))
                    break;
                }

                if ((num_mismatches <= localerror) &&
                    (num_mismatches <= num_mismatches_rev)) {
                  howmany++;

                  if (num_mismatches == 0) buffers.used[uu][j] = TRUE;
                  flagover = FALSE;

                  if ((j - options.motif_length + 1 > lastocc) || (lastocc == -1)) {
                    flagover = TRUE;
                  }

                  if ((num_mismatches == bestsofar) && (lastinseq == uu))
                    inbest++;

                  if (lastinseq != uu) {
                    inbest = 1;
                    lastinseq = uu;
                    howmanyseq++;
                    bestpos = j;
                    bestsofar = num_mismatches;
                  }

                  if (num_mismatches < bestsofar) {
                    inbest = 1;
                    bestpos = j;
                    bestsofar = num_mismatches;

                    if (options.onlybestflag) localerror = bestsofar;
                  }
                }

                if ((num_mismatches_rev <= localerror) &&
                    (num_mismatches_rev < num_mismatches)) {
                  if (num_mismatches_rev == 0) buffers.used[uu][j] = TRUE;

                  howmany++;
                  flagover = FALSE;

                  if ((j - options.motif_length + 1 > lastocc) || (lastocc == -1)) {
                    flagover = TRUE;
                  }

                  if ((num_mismatches_rev == bestsofar) && (lastinseq == uu))
                    inbest++;

                  if (lastinseq != uu) {
                    inbest = 1;
                    lastinseq = uu;
                    howmanyseq++;
                    bestpos = j;
                    bestsofar = num_mismatches_rev;
                  }

                  if (num_mismatches_rev < bestsofar) {
                    inbest = 1;
                    bestpos = j;
                    bestsofar = num_mismatches_rev;
                    if (options.onlybestflag) localerror = bestsofar;
                  }
                  num_mismatches = 0;
                }
              }
            }
            if (bestpos > -1) {
              buffers.besthere[uu] = bestsofar;
              buffers.inbestseq[uu] = inbest;
            }
          }

          if (howmanyseq >= options.min_repeat_percentage) {
            patfreqs[0] =
              exactfreqs(seekpat->next->pat, 0) +
              exactfreqswithcheck(revseek, 0,
                                  seekpat->next->pat);
            
            for (uu = 1; uu <= options.max_mismatches; uu++) {
              if (options.motif_length != 12)
                patfreqs[uu] =
                  exactfreqs(seekpat->next->pat, uu) +
                  exactfreqswithcheck(revseek, uu, seekpat->next->pat);
              else
                patfreqs[uu] =
                  exactfreqs(seekpat->next->pat, uu) +
                  exactfreqswithcheck(revseek, 0, seekpat->next->pat) +
                  patfreqs[uu - 1];
            }

            seekpat->next->score = 0.0;
            myscore = 0.0;

            for (uu = 0; uu < num_sequences_total; uu++) {
              double tmpfloat;

              if (buffers.besthere[uu] > -1)
                tmpfloat = log((double) buffers.inbestseq[uu]) -
                  (log(patfreqs[buffers.besthere[uu]]) + log(sequence_lengths[uu]));
              else
                tmpfloat = 0.0;

              myscore += tmpfloat;
            }
            seekpat->next->score = myscore;

            if (!options.onlybestflag)
              seekpat->next->score +=
                log((double) howmany) - (log(patfreqs[options.max_mismatches]) +
                                         log((double) sum_all_seq_lengths));
          } else {
            seekpat->next->score = -1;
          }
        }
	    }
    }
    seq_index++;
  }

  if (head_node.next != &tail_node)
    head_node.next = mergesrt(head_node.next, resultlist_length);

  report_best_results(wee_fp, html_fp, mix_fp, &head_node);
  report_time_taken(wee_fp, html_fp);

  fflush(wee_fp);
  fflush(mix_fp);
  fflush(html_fp);
  fclose(wee_fp);
  fclose(mix_fp);
  fclose(html_fp);

  fprintf(stderr, "\nOutput written in files %s and %s\n", weefile, htmlfile);
  exit(0);
}

void report_best_results(FILE *wee_fp, FILE *html_fp, FILE *mix_fp, Pattern *seekpat)
{
  int i = 0;
  while ((seekpat->next->score > 0) && (i < options.max_results)) {
    i++;

    fprintf(stderr,  "%d) %s %.2f\n",    i,
            seekpat->next->pat, seekpat->next->score/num_sequences_total);
    fprintf(wee_fp,  "%d) %s %.2f \n",   i,
            seekpat->next->pat, seekpat->next->score/num_sequences_total);
    fprintf(html_fp, "%d) %s %.2f <br>", i,
            seekpat->next->pat, seekpat->next->score/num_sequences_total);
    fprintf(mix_fp,  "%d) %s %.2f %d\n", i,
            seekpat->next->pat, seekpat->next->score/num_sequences_total,
            options.max_mismatches);
    seekpat = seekpat->next;
  }
}

void report_time_taken(FILE *wee_fp, FILE *html_fp)
{
  time_t endtime = clock() / CLOCKS_PER_SEC;
  fprintf(wee_fp,  "\nElapsed time : %d min. %d sec.\n\n",
          (int)((endtime - starttime) / 60), (int) ((endtime - starttime) % 60));
  fprintf(html_fp, "<br>Elapsed time : %d min. %d sec.<br><br>",
          (int)((endtime - starttime) / 60), (int)((endtime - starttime) % 60));
  fprintf(stderr,  "Elapsed time : %d min. %d sec.\n",
          (int)((endtime - starttime) / 60), (int)((endtime - starttime) % 60));
}

void readfreqs(void)
{
  int pre, post, backpre, backpost, i, j, k, n, tmpint, sum, deccount,
    i1, i2, i3, i4, i5, i6, i7, i8;
  double backfreqpre, backfreqpost, decfreq, prefreq, postfreq,
    totalcountsix, totalcounteight, checksum = 0, decsum = 0;
  char freqfname[40] = "", decamer[10];
  FILE *freq_fp;

  if (options.verbose) fprintf(stderr, "\nComputing frequencies...\n");

  if (!strcmp(options.organism, "")) {
    fprintf(stderr, "\nUndefined organism (must add -O)\n");
    exit(0);
  }
  string_toupper(options.organism, options.organism);

  /* Read 8-mer freq file */
  sprintf(freqfname, FORMAT_8MER_FREQ_FILE, options.organism);
  freq_fp = fopen(freqfname, "r");

  if (freq_fp == NULL) {
    fprintf(stderr,
            "\nCannot find frequency file %s\nWrong organism code, or the FreqFiles directory is not correctly positioned\n",
            freqfname);
    exit(0);
  }

  i = 0;
  n = fscanf(freq_fp, "%*s %d\n", &tmpint);
  sum = 0;

  while (n != EOF) {
    sum += tmpint;
    eightmers[i] = (double) tmpint;
    i++;

    n = fscanf(freq_fp, "%*s %d\n", &tmpint);
  }
  fclose(freq_fp);
  totalcounteight = sum;

  checksum = 0;

  for (i = 0; i < NUM_8MER_PERMUTATIONS; i++) {
    eightmers[i] /= (double) totalcounteight;
    checksum += eightmers[i];
  }

  /* Read 6-mer freq file  */
  sprintf(freqfname, FORMAT_6MER_FREQ_FILE, options.organism);
  freq_fp = fopen(freqfname, "r");

  if (freq_fp == NULL) {
    fprintf(stderr, "\nCannot find frequency file %s\n", freqfname);
    exit(0);
  }

  i = 0;
  sum = 0;
  n = fscanf(freq_fp, "%*s %d\n", &tmpint);

  while (n != EOF) {
    sum += tmpint;
    sixmers[i] = (double) tmpint;
    i++;
    n = fscanf(freq_fp, "%*s %d\n", &tmpint);
  }

  totalcountsix = sum;

  for (i = 0; i < NUM_6MER_PERMUTATIONS; i++) sixmers[i] /= (double) totalcountsix;

  if (options.motif_length >= 10) {
    int eightmer;

    for (i = 0; i < ALPHABET_SIZE; i++) {
      decamer[0] = inttochar(i);
      eightmer = 0;

      for (i1 = 0; i1 < ALPHABET_SIZE; i1++) {
	      decamer[1] = inttochar (i1);
	      eightmer += POW_4_7 * i1;

	      for (i2 = 0; i2 < ALPHABET_SIZE; i2++) {
          decamer[2] = inttochar(i2);
          eightmer += POW_4_6 * i2;
          for (i3 = 0; i3 < ALPHABET_SIZE; i3++) {
            decamer[3] = inttochar(i3);
            eightmer += POW_4_5 * i3;
            for (i4 = 0; i4 < ALPHABET_SIZE; i4++) {
              decamer[4] = inttochar(i4);
              eightmer += POW_4_4 * i4;
              for (i5 = 0; i5 < ALPHABET_SIZE; i5++) {
                decamer[5] = inttochar(i5);
                eightmer += POW_4_3 * i5;
                for (i6 = 0; i6 < ALPHABET_SIZE; i6++) {
                  decamer[6] = inttochar(i6);
                  eightmer += POW_4_2 * i6;
                  for (i7 = 0; i7 < ALPHABET_SIZE; i7++) {
                    decamer[7] = inttochar(i7);
                    eightmer += POW_4_1 * i7;

                    for (i8 = 0; i8 < ALPHABET_SIZE; i8++) {
                      decamer[8] = inttochar(i8);
                      pre = eightmer / 4 + POW_4_7 * i;
                      backpre = 0;
                      backfreqpre = 0.0;

                      for (k = 0; k < ALPHABET_SIZE; k++) {
                        backpre = eightmer / 4 + POW_4_7 * k;
                        backfreqpre += eightmers[backpre];
                      }

                      prefreq = eightmers[pre] / backfreqpre;
                      eightmer += POW_4_0 * i8;

                      for (j = 0; j < ALPHABET_SIZE; j++) {
                        decfreq = prefreq * eightmers[eightmer];
                        decamer[9] = inttochar(j);

                        post = (eightmer - (POW_4_7 * i1)) * 4;
                        backpost = (eightmer - (POW_4_7 * i1)) * 4;
                        post += j;
                        backfreqpost = 0;

                        for (k = 0; k < ALPHABET_SIZE; k++) {
                          backpost = backpost + k;
                          backfreqpost += eightmers[backpost];
                          backpost -= k;
                        }
                        
                        decfreq = decfreq * (eightmers[post] / backfreqpost);
                        postfreq = (eightmers[post] / backfreqpost);
                        decsum += decfreq;
                        deccount = eightmer * 4 + (i * POW_4_9) + j;
                        tenmers[deccount] = decfreq;
                      }
                      eightmer -= POW_4_0 * i8;
                    }
                    eightmer -= POW_4_1 * i7;
                  }
                  eightmer -= POW_4_2 * i6;
                }
                eightmer -= POW_4_3 * i5;
              }
              eightmer -= POW_4_4 * i4;
            }
            eightmer -= POW_4_5 * i3;
          }
          eightmer -= POW_4_6 * i2;
        }
	      eightmer -= POW_4_7 * i1;
	    }
    }
  }
}

double compsix(char *pat, char *pos, int err)
{
  int i, jj, nn, totapp;
  double result = 0;
  int appro[POW_4_1 + 1]; /* err == 1 */

  nn = strlen(pat) - 1;
  totapp = pow_4(err);

  for (i = 0; i <= totapp; i++) appro[i] = 0;

  for (jj = 0; jj <= nn; jj++) {
    if (pos[nn - jj] != '*') {
      for (i = 0; i < totapp; i++)
        appro[i] += chartoint(pat[nn - jj]) * pow_4(jj);
    } else {
      for (i = 0; i < totapp; i++) appro[i] += (i % 4) * pow_4(jj);
    }
  }
  for (i = 0; i < totapp; i++) result += sixmers[appro[i]];
  return result;
}

double compfreqs(char *pat, char *pos, int err)
{
  int i, jj, nn, totapp, contast = 0;
  double result = 0;
  int appro[POW_4_3 + 1]; /* 1 <= err <= 3 */

  nn = strlen(pat) - 1;
  totapp = pow_4(err);

  for (i = 0; i <= totapp; i++) appro[i] = 0;

  for (jj = 0; jj <= nn; jj++) {
    if (pos[nn - jj] != '*') {
      for (i = 0; i < totapp; i++)
        appro[i] += chartoint(pat[nn - jj]) * pow_4(jj);
    } else {
      for (i = 0; i < totapp; i++)
        appro[i] += ((i / pow_4(contast)) % 4) * pow_4(jj);

      contast++;
    }
  }

  if (nn == 7) {
    for (i = 0; i < totapp; i++) {
      result += eightmers[appro[i]];
    }
  }
  if (nn == 9) {
    for (i = 0; i < totapp; i++) {
      result += tenmers[appro[i]];
    }
  }
  return result;
}

double exactfreqs(char *pat, int errs)
{
  int i, j, k, jj, kk, patno = 0;
  double result = 0.0;
  char ast[9], longast[13];
  char empty[9] = "        ";
  char longempty[13] = "            ";

  bzero(ast, 9);
  bzero(longast, 13);

  if (strlen(pat) == 6) {
    for (jj = 0; jj < 6; jj++)
      patno += chartoint(tolower(pat[5 - jj])) * pow_4(jj);

    if (errs == 0) return sixmers[patno];
    if (errs == 1) {
      for (j = 0; j < 6; j++) {
	      strcpy (ast, empty);
	      ast[j] = '*';

	      result += compsix(pat, ast, 1);
	    }
      return result;
    }
  }

  if (strlen(pat) == 8) {
    for (jj = 0; jj < 8; jj++)
      patno += chartoint(tolower(pat[7 - jj])) * pow_4(jj);

      if (errs == 0) return eightmers[patno];
      if (errs == 1) {
        for (j = 0; j < 8; j++) {
          strcpy(ast, empty);
          ast[j] = '*';
          result += compfreqs(pat, ast, 1);
        }
        return result;
      }

      if (errs == 2) {
        for (j = 0; j < 7; j++) {
          strcpy(ast, empty);
          ast[j] = '*';

          for (k = j + 1; k < 8; k++) {
            ast[k] = '*';
            if (k > j + 1) ast[k - 1] = ' ';

            result += compfreqs(pat, ast, 2);
          }
        }
        return result;
      }
  }

  if (strlen(pat) == 10) {
    if (errs == 0) {
      patno = 0;
        
      for (jj = 0; jj < 10; jj++)
        patno += chartoint(tolower(pat[9 - jj])) * pow_4(jj);
      result = tenmers[patno];

      return result;
    }
    if (errs == 1) {
      for (j = 0; j < 10; j++) {
        strcpy (longast, longempty);
        longast[j] = '*';
          
        result += compfreqs(pat, longast, 1);
      }
      return result;
    }
    if (errs == 2) {
      for (j = 0; j < 9; j++) {
        strcpy(longast, longempty);
        longast[j] = '*';

        for (k = j + 1; k < 10; k++) {
          longast[k] = '*';

          if (k > j + 1) longast[k - 1] = ' ';
          result += compfreqs(pat, longast, 2);
        }
      }
      return result;
    }
    if (errs == 3) {
      for (j = 0; j < 8; j++) {
        strcpy (longast, longempty);
        longast[j] = '*';

        for (k = j + 1; k < 9; k++) {
          longast[k] = '*';

          for (i = k + 1; i < 10; i++) {
            longast[i] = '*';

            result += compfreqs(pat, longast, 3);
            longast[i] = ' ';
          }
          longast[k] = ' ';
        }
      }
      return result;
    }
  }
  if (strlen(pat) == 12) {
    if (errs == 0) {
      strcpy(longast, longempty);
      result += newtwmers(pat, longast, 0);
      return result;
    }
    if (errs == 1) {
      for (j = 0; j < 12; j++) {
	      strcpy(longast, longempty);
	      longast[j] = '*';
	      result += newtwmers(pat, longast, 1);
	    }
      return result;
    }
    if (errs == 2) {
      for (j = 0; j < 11; j++) {
	      strcpy(longast, longempty);
	      longast[j] = '*';

	      for (k = j + 1; k < 12; k++) {
          longast[k] = '*';
          if (k > j + 1) longast[k - 1] = ' ';
          result += newtwmers(pat, longast, 2);
          longast[k] = ' ';
        }
	    }
      return result;
    }
    if (errs == 3) {
      for (j = 0; j < 10; j++) {
	      strcpy(longast, longempty);
	      longast[j] = '*';

	      for (k = j + 1; k < 11; k++) {
          longast[k] = '*';

          for (i = k + 1; i < 12; i++) {
            longast[i] = '*';
            result += newtwmers(pat, longast, 3);
            longast[i] = ' ';
          }
          longast[k] = ' ';
        }
	    }
      return result;
    }
    if (errs == 4) {
      for (j = 0; j < 9; j++) {
	      strcpy(longast, longempty);
	      longast[j] = '*';

	      for (k = j + 1; k < 10; k++) {
          longast[k] = '*';

          for (i = k + 1; i < 11; i++) {
            longast[i] = '*';

            for (kk = i + 1; kk < 12; kk++) {
              longast[kk] = '*';

              result += newtwmers(pat, longast, 4);
              longast[kk] = ' ';
            }
            longast[i] = ' ';
          }
          longast[k] = ' ';
        }
	    }
      return result;
    }
  }
  return 5859.0;
}

int pattolong(const char *pat)
{
  int i, result = 0, len;

  len = strlen(pat) - 1;
  for (i = 0; i <= len; i++) result += chartoint(pat[len - i]) * pow_4(i);

  return result;
}

double tentotwelve(char *pat) {
  int i, len, prelong, postlong, primo;
  char mid[11], post[11], pre[11], back[11];
  double startfreq = 0, frontfreq = 0, endfreq = 0, norm = 0;

  bzero(post, 11);
  bzero(pre, 11);
  bzero(back, 11);
  bzero(mid, 11);

  len = strlen(pat);

  for (i = 0; i < 10; i++) {
    mid[i]  = pat[i + 1];
    post[i] = pat[i + 2];
    pre[i]  = pat[i];
  }
  strcpy(back, pre);
  prelong = pattolong(pre);
  frontfreq = tenmers[prelong];
  norm = frontfreq;
  primo = chartoint(pre[0]);

  /* ATTENZIONE */
  for (i = 0; i < ALPHABET_SIZE; i++) {
    back[0] = ALPHABET[i];
    prelong -= (primo * POW_4_9);

    if (ALPHABET[i] != pre[0]) {
      primo = i;
      prelong += (primo * POW_4_9);
      norm += tenmers[prelong];
    }
  }

  frontfreq /= norm;
  strcpy(back, post);
  postlong = pattolong(post);
  endfreq = tenmers[postlong];
  primo = chartoint(post[9]);
  norm = endfreq;

  for (i = 0; i < ALPHABET_SIZE; i++) {
    back[9] = ALPHABET[i];
    postlong -= primo;

    if (ALPHABET[i] != post[9]) {
      primo = i;
      postlong += primo;
      norm += tenmers[postlong];
    }
  }
  endfreq /= norm;
  startfreq = tenmers[pattolong(mid)];
  return startfreq * frontfreq * endfreq;
}

double newtwmers(char *inpat, char *pos, int err)
{
  int i, j, jj, k, len, paterr, asts[4];
  char pre[9], post[9], pat[13];
  double result = 0; 

  len = strlen (inpat);
  strcpy(pat, inpat);
  j = 0;

  for (i = 0; i < len; i++) {
    if (pos[i] == '*') asts[j++] = i;
  }
  bzero(post, 9);
  bzero(pre, 9);

  if (err == 0) return tentotwelve(pat);

  if (err == 1) {
    for (i = 0; i < ALPHABET_SIZE; i++) {
      paterr = 0;
      pat[asts[0]] = ALPHABET[i];

      if (pat[asts[0]] != inpat[asts[0]]) result += tentotwelve(pat);
    }
    return result;
  }

  if (err == 2) {
    for (i = 0; i < ALPHABET_SIZE; i++) {
      pat[asts[0]] = ALPHABET[i];

      if (pat[asts[0]] != inpat[asts[0]]) {
	      for (j = 0; j < ALPHABET_SIZE; j++) {
          pat[asts[1]] = ALPHABET[j];

          if (pat[asts[1]] != inpat[asts[1]]) result += tentotwelve(pat);
        }
	    }
    }
    return result;
  }

  if (err == 3) {
    result += tentotwelve(inpat);

    for (i = 0; i < ALPHABET_SIZE; i++) {
      pat[asts[0]] = ALPHABET[i];

      if (pat[asts[0]] != inpat[asts[0]]) {
	      for (j = 0; j < ALPHABET_SIZE; j++) {
          pat[asts[1]] = ALPHABET[j];

          if (pat[asts[1]] != inpat[asts[1]]) {
            for (k = 0; k < ALPHABET_SIZE; k++) {
              pat[asts[2]] = ALPHABET[k];

              if (pat[asts[1]] != inpat[asts[1]]) result += tentotwelve(pat);
            }
          }
        }
	    }
    }
    return result;
  }
  if (err == 4) {
    result += tentotwelve(inpat);

    for (i = 0; i < ALPHABET_SIZE; i++) {
      pat[asts[0]] = ALPHABET[i];

      if (pat[asts[0]] != inpat[asts[0]]) {
	      for (j = 0; j < ALPHABET_SIZE; j++) {
          pat[asts[1]] = ALPHABET[j];

          if (pat[asts[1]] != inpat[asts[1]]) {
            for (k = 0; k < ALPHABET_SIZE; k++) {
              pat[asts[2]] = ALPHABET[k];

              if (pat[asts[2]] != inpat[asts[2]]) {
                for (jj = 0; jj < ALPHABET_SIZE; jj++) {
                  pat[asts[3]] = ALPHABET[jj];

                  if (pat[asts[3]] != inpat[asts[3]])
                    result += tentotwelve(pat);
                }
              }
            }
          }
        }
	    }
    }
    return result;
  }
  /* WARNING - the original code does not specify handling for this !! */
  fprintf(stderr, "*WARNING*: newtwmers() - unexpected return.\n");
  return result;
}

double exactfreqswithcheck(char *pat, int errs, char *check)
{
  int i, j, k, jj, kk, patno = 0;
  char partial[13];
  char ast[9], longast[13];
  char empty[9] = "        ";
  char longempty[13] = "            ";
  double result = 0.0;

  bzero(ast, 9);
  bzero(longast, 13);
  bzero(partial, 13);

  if (strlen(pat) == 6) {
    for (jj = 0; jj < 6; jj++)
      patno += chartoint(tolower(pat[5 - jj])) * pow_4(jj);

    if (errs == 0) {
      if (strcmp(pat, check) != 0) return sixmers[patno];
      else return 0.0;
    }

    if (errs == 1) {
      for (j = 0; j < 6; j++) {
	      strcpy(ast, empty);
	      ast[j] = '*';
	      result += compsix(pat, ast, 1);
	    }
      return result;
    }
  }

  if (strlen(pat) == 8) {
    for (jj = 0; jj < 8; jj++)
      patno += chartoint(tolower(pat[7 - jj])) * pow_4(jj);

    if (errs == 0) {
      if (strcmp(pat, check) != 0) return eightmers[patno];
      else return 0.0;
    }

    if (errs == 1) {
      for (j = 0; j < 8; j++) {
	      strcpy (ast, empty);
	      ast[j] = '*';
	      result += compfreqswithcheck(pat, ast, 1, check);
	    }
      return result;
    }

    if (errs == 2) {
      for (j = 0; j < 7; j++) {
	      strcpy (ast, empty);
	      ast[j] = '*';

	      for (k = j + 1; k < 8; k++) {
          ast[k] = '*';
          result += compfreqswithcheck(pat, ast, 2, check);
          ast[k] = ' ';
        }
	    }
      return result;
    }
  }

  if (strlen(pat) == 10) {
    if (errs == 0) {
      patno = 0;

      for (jj = 0; jj < 10; jj++)
        patno += chartoint(tolower(pat[9 - jj])) * pow_4(jj);
      result = tenmers[patno];
      return result;
    }
    if (errs == 1) {
      for (j = 0; j < 10; j++) {
	      strcpy(longast, longempty);
	      longast[j] = '*';
	      result += compfreqswithcheck(pat, longast, 1, check);
	    }
      return result;
    }
    if (errs == 2) {
      for (j = 0; j < 9; j++) {
	      strcpy(longast, longempty);
	      longast[j] = '*';

	      for (k = j + 1; k < 10; k++) {
          longast[k] = '*';
          result += compfreqswithcheck(pat, longast, 2, check);
          longast[k] = ' ';
        }
	    }
      return result;
    }
    if (errs == 3) {
      for (j = 0; j < 8; j++) {
	      strcpy(longast, longempty);
	      longast[j] = '*';

	      for (k = j + 1; k < 9; k++) {
          longast[k] = '*';

          for (i = k + 1; i < 10; i++) {
            longast[i] = '*';
            result += compfreqswithcheck(pat, longast, 3, check);
            longast[i] = ' ';
          }
          longast[k] = ' ';
        }
	    }
      return result;
    }
  }
  if (strlen(pat) == 12) {
    if (errs == 0) {
      strcpy(longast, longempty);
      if (hamming_distance_exceeds_threshold(pat, check, 0))
        result += newtwmers (pat, longast, 0);
      else result = 0;
      return result;
    }
    if (errs == 1) {
      for (j = 0; j < 12; j++) {
	      strcpy (longast, longempty);
	      longast[j] = '*';
	      result += newtwmerswithcheck(pat, longast, 1, check);
	    }
      return result;
    }
    if (errs == 2) {
      for (j = 0; j < 11; j++) {
	      strcpy (longast, longempty);
	      longast[j] = '*';

	      for (k = j + 1; k < 12; k++) {
          longast[k] = '*';
          result += newtwmerswithcheck(pat, longast, 2, check);
          longast[k] = ' ';
        }
	    }
      return result;
    }
    if (errs == 3) {
      for (j = 0; j < 10; j++) {
	      strcpy(longast, longempty);
	      longast[j] = '*';

	      for (k = j + 1; k < 11; k++) {
          longast[k] = '*';

          for (i = k + 1; i < 12; i++) {
            longast[i] = '*';
            result += newtwmerswithcheck(pat, longast, 3, check);
            longast[i] = ' ';
          }
          longast[k] = ' ';
        }
	    }
      return result;
    }
    if (errs == 4) {
      for (j = 0; j < 9; j++) {
        strcpy(longast, longempty);
        longast[j] = '*';

        for (k = j + 1; k < 10; k++) {
          longast[k] = '*';

          for (i = k + 1; i < 11; i++) {
            longast[i] = '*';
            for (kk = i + 1; kk < 12; kk++) {
              longast[kk] = '*';

              result += newtwmerswithcheck(pat, longast, 4, check);
              longast[kk] = ' ';
            }
            longast[i] = ' ';
          }
          longast[k] = ' ';
        }
      }
      return result;
    }
  }
  return 5859.0;
}

double compfreqswithcheck(char *pat, char *pos, int err, char *check)
{
  int i, jj, nn, totapp, contast = 0;
  int appro[POW_4_3 + 1]; /* 1 <= err <= 3 */
  char partial[POW_4_3 + 1][15]; /* max chars per row = 15 */
  double result = 0;

  nn = strlen(pat) - 1; /* nn == |pat| - 1, |pat| == |revseek| == 14 => 0 <= nn <= 13 */
  totapp = pow_4(err);

  for (i = 0; i <= totapp; i++) bzero(partial[i], nn + 2);

  for (i = 0; i <= totapp; i++) appro[i] = 0;
  for (i = 0; i <= totapp; i++) strcpy(partial[i], pat);

  for (jj = 0; jj <= nn; jj++) {
    if (pos[nn - jj] != '*') {
      for (i = 0; i < totapp; i++)
        appro[i] += chartoint(pat[nn - jj]) * pow_4(jj);
    } else {
      for (i = 0; i < totapp; i++) {
        appro[i] += ((i / pow_4(contast)) % 4) * pow_4(jj);
        partial[i][jj] = inttochar(((i / pow_4(contast)) % 4));
      }
      contast++;
    }
  }

  if (nn == 7) {
    for (i = 0; i < totapp; i++) {
      if (hamming_distance_exceeds_threshold(partial[i], check, err))
        result += eightmers[appro[i]];
    }
  }
  if (nn == 9) {
    for (i = 0; i < totapp; i++) {
      if (hamming_distance_exceeds_threshold(partial[i], check, err))
        result += tenmers[appro[i]];
    }
  }
  return result;
}

double newtwmerswithcheck(char *inpat, char *pos, int err, char *check)
{
  int i, j, jj, k, len, asts[4];
  char pre[9], post[9], pat[13];
  double result = 0;

  len = strlen(inpat);
  strcpy(pat, inpat);
  j = 0;

  for (i = 0; i < len; i++) {
    if (pos[i] == '*') asts[j++] = i;
  }

  bzero(post, 9);
  bzero(pre, 9);

  if (err == 0) {
    return (hamming_distance_exceeds_threshold(pat, check, err)) ?
      tentotwelve(pat) : 0;
  }

  if (err == 1) {
    for (i = 0; i < ALPHABET_SIZE; i++) {
      pat[asts[0]] = ALPHABET[i];

      if (pat[asts[0]] != inpat[asts[0]]) {
	      if (hamming_distance_exceeds_threshold(pat, check, err))
          result += tentotwelve(pat);
	    }
    }
    return result;
  }

  if (err == 2) {
    result += tentotwelve (inpat);
    for (i = 0; i < ALPHABET_SIZE; i++) {
      pat[asts[0]] = ALPHABET[i];
      
      if (pat[asts[0]] != inpat[asts[0]]) {
        for (j = 0; j < ALPHABET_SIZE; j++) {
          pat[asts[1]] = ALPHABET[j];

          if (pat[asts[1]] != inpat[asts[1]]) {
            if (hamming_distance_exceeds_threshold(pat, check, err)
                && strcmp(pat, inpat))
              result += tentotwelve(pat);
          }
        }
      }
    }
    return result;
  }

  if (err == 3) {
    result += tentotwelve(inpat);

    for (i = 0; i < ALPHABET_SIZE; i++) {
      pat[asts[0]] = ALPHABET[i];

      if (pat[asts[0]] != inpat[asts[0]]) {
	      for (j = 0; j < ALPHABET_SIZE; j++) {
          pat[asts[1]] = ALPHABET[j];

          if (pat[asts[1]] != inpat[asts[1]]) {
            for (k = 0; k < ALPHABET_SIZE; k++) {
              pat[asts[2]] = ALPHABET[k];

              if (pat[asts[2]] != inpat[asts[2]]) {
                if (hamming_distance_exceeds_threshold(pat, check, err)
                    && strcmp(pat, inpat))
                  result += tentotwelve(pat);
              }
            }
          }
        }
	    }
    }
    return result;
  }
  if (err == 4) {
    result += tentotwelve(inpat);
    for (i = 0; i < ALPHABET_SIZE; i++) {
      pat[asts[0]] = ALPHABET[i];
      
      if (pat[asts[0]] != inpat[asts[0]]) {
	      for (j = 0; j < ALPHABET_SIZE; j++) {
          pat[asts[1]] = ALPHABET[j];

          if (pat[asts[1]] != inpat[asts[1]]) {
            for (k = 0; k < ALPHABET_SIZE; k++) {
              pat[asts[2]] = ALPHABET[k];

              if (pat[asts[2]] != inpat[asts[2]]) {
                for (jj = 0; jj < ALPHABET_SIZE; jj++) {
                  pat[asts[3]] = ALPHABET[jj];

                  if (pat[asts[3]] != inpat[asts[3]]) {
                    if (hamming_distance_exceeds_threshold(pat, check, err)
                        && strcmp(pat, inpat))
                      result += tentotwelve(pat);
                  }
                }
              }
            }
          }
        }
	    }
    }
    return result;
  }

  /* WARNING: The original code does not include a default return */
  fprintf(stderr, "*WARNING*: newtwmerswithcheck() - unexpected return\n");
  return result;
}

BOOL checkmotifanderror(int motif_length, int errs) {
  switch (motif_length) {
  case 6:  return errs > 1 ? 0 : 1;
  case 8:  return errs > 3 ? 0 : 1;
  case 10: return errs > 4 ? 0 : 1;
  case 12: return errs > 4 ? 0 : 1;
  default: return FALSE;
  }
}

void print_usage(const char *progname)
{
  fprintf(stderr,
          "\nUsage : %s -f inputfilename -O speciescode -W motif_width -e error -R repeatspercentage <options>\n"
          "\nspeciescode: two-letter species code (e.g. HS, MM, RN, SC, and so on)\n"
          "-W : motif length\n"
          "-e : number of mutations allowed\n"
          "-R : sequences percentage on which motif has to appear\n",
          progname);
  fputs("\noptions (not required):\n"
        "-S: process both strands of input sequences (default: single strand)\n"
        "-M: multiple motif occurrences in each sequence (default: expect zero or one occurrences per sequence)\n"
        "-T <number>: report top <number> motifs\n"
        "-V: verbose mode\n", stderr);
}

void print_option_error(const char *progname)
{
  fprintf(stderr,
          "\nSomething wrong with the commandline\n"
          "\nUsage : %s -f inputfilename -O speciescode -W motif_width -e error -R repeatspercentage <options>\n"
          "\nspeciescode: two-letter species code (e.g. HS, MM, RN, SC, and so on)\n"
          "-W : motif length\n"
          "-e : number of mutations allowed\n"
          "-R : sequences percentage on which motif has to appear\n",
          progname);
  fputs("\noptions (not required):\n"
        "-S: process both strands of input sequences (default: single strand)\n"
        "-M: multiple motif occurrences in each sequence (default: expect zero or one occurrences per sequence)\n"
        "-T <number>: report top <number> motifs\n"
        "-V: verbose mode\n", stderr);
}

BOOL has_option_value(int argc, char *argv[], int option_index)
{
  return (option_index + 1) < argc && argv[option_index + 1][0] != '-';
}

int parse_option(int argc, char *argv[], int option_index)
{
  switch (argv[option_index][1]) {
  case 'O':
    check_condition(has_option_value(argc, argv, option_index),
                    "\nProblem in organism definition\n", 1);
    strcpy(options.organism, argv[option_index + 1]);
    option_index += 2;
    break;
  case 'W':
    check_condition(has_option_value(argc, argv, option_index),
                    "\nProblem in motif length definition\n", 1);
    options.motif_length = atoi(argv[option_index + 1]);
    check_condition(options.motif_length == 6 || options.motif_length == 8 ||
                    options.motif_length == 10 || options.motif_length == 12,
                    "\nAccepted motif lengths are even values between 6 and 12\n", 1);
    option_index += 2;
    break;
  case 'T':
    check_condition(has_option_value(argc, argv, option_index),
                    "\nProblem in number of reported motifs definition\n", 1);
    options.max_results = atoi(argv[option_index + 1]);
    option_index += 2;
    break;
  case 'e':
    check_condition(has_option_value(argc, argv, option_index),
                    "\nError in defining mismatch number\n", 1);
    options.max_mismatches = atoi(argv[option_index + 1]);
    check_condition(options.max_mismatches >= 0, "\nError in defining mismatch number\n", 1);
    option_index += 2;
    break;
	case 'R':
    check_condition(has_option_value(argc, argv, option_index),
                    "\nError in defining repeat percentage -R\n", 1);
    options.min_repeat_percentage = atoi(argv[option_index + 1]);

    check_condition(options.min_repeat_percentage >= 1 &&
                    options.min_repeat_percentage <= 100,
                    "\nError in defining repeat percentage -R\n", 1);
    option_index += 2;
    repeat_percentage = options.min_repeat_percentage;
	  break;
	case 'f':
	  if (!strcmp ("-f", argv[option_index])) {
      check_condition(option_index + 1 < argc,
                      "\nDid not understand filename\n", 1);
      strcpy(options.infile, argv[option_index + 1]);
      extfile = 1;
      option_index += 2;
    }
	  break;
	case 'V':
	  options.verbose = TRUE;
	  option_index++;
	  break;
	case 'S':
	  options.reversed = TRUE;
	  option_index++;
	  break;
	case 'N':
	  option_index++;
	  break;
	case 'M':
	  options.onlybestflag = FALSE;
	  option_index++;
	  break;
	default:
    print_option_error(argv[0]);
	  exit(1);
	  break;
  }
  return option_index;
}

Pattern *mergesrt(Pattern *list, int length)
{
  return (Pattern *) mergesort_linkedlist((BaseNode *) list,
                                          (BaseNode *) &tail_node,
                                          length);
}

void check_condition(BOOL condition, const char *errormsg, int exitcode)
{
  if (!condition) {
    fputs(errormsg, stderr);
    exit(exitcode);
  }
}

void determine_num_sequences_and_lengths(char *sequence_buffer)
{
  FILE *input_fp;
  int seq_index = 0, nchar = 0;

  num_sequences_total = 0;
  sequence_lengths = (int *) malloc(MAXSEQ * sizeof(int));
  bzero(sequence_buffer, MAXSEQ);

  input_fp = fopen(options.infile, "r");
  nchar = fscanf(input_fp, "\n%[^\n]", sequence_buffer);

  check_condition(sequence_buffer[0] == '>',
                  "\nInput file not in FASTA format\nExiting program\n", 0);

  while (nchar != EOF) {
    int partiallen = 0;

    num_sequences_total++;
    nchar = fscanf(input_fp, "\n%[^\n]", sequence_buffer);

    while (sequence_buffer[0] != '>' && nchar != EOF) {
      partiallen += strlen(sequence_buffer);
      bzero(sequence_buffer, MAXSEQ);
      nchar = fscanf(input_fp, "\n%[^\n]", sequence_buffer);
    }
    sequence_lengths[seq_index++] = partiallen;
  }

  /* Fine conteggio lunghezza sequenze */
  fclose(input_fp);

  if (num_sequences_total > MAXSEQ) {
    fprintf(stderr, "\nSorry! At most %d sequences can be input!\n", MAXSEQ);
    exit(1);
  }
}

void init_workbuffers(WorkBuffers *buffers)
{
  int seq_index, base_index;

  buffers->besthere  = (int *) malloc(num_sequences_total * sizeof(int));
  buffers->inbestseq = (int *) malloc(num_sequences_total * sizeof(int));
  buffers->sequence  = (char **) malloc((num_sequences_total + 1) * sizeof(char *));
  buffers->used      = (BOOL **) malloc((num_sequences_total + 1) * sizeof(BOOL *));
  buffers->name      = (char **) malloc((num_sequences_total + 1) * sizeof(char *));

  for (seq_index = 0; seq_index < num_sequences_total + 1; seq_index++) {
    buffers->name[seq_index]     = (char *) malloc(SEQUENCE_NAME_LENGTH * sizeof(char));
    buffers->sequence[seq_index] = (char *) malloc((sequence_lengths[seq_index] + 1 +
                                                    20 + options.motif_length) * sizeof(char));
    buffers->used[seq_index]     = (BOOL *)
      malloc(sequence_lengths[seq_index] * sizeof(BOOL));
  }
  for (seq_index = 0; seq_index < num_sequences_total; seq_index++) {
    for (base_index = 0; base_index < sequence_lengths[seq_index]; base_index++)
      buffers->used[seq_index][base_index] = FALSE;
  }
}

void read_sequence_data(WorkBuffers *buffers, char *sequence_buffer)
{
  int seq_index = 0, nchar;
  FILE *input_fp = fopen(options.infile, "r");
  bzero(sequence_buffer, MAXSEQ);
  nchar = fscanf(input_fp, "\n%[^\n]", sequence_buffer);
  strncpy(buffers->sequence[seq_index], sequence_buffer, 20);

  while (seq_index < num_sequences_total) {
    int seqlen;
    if (sequence_buffer[0] == '>') {
      strncpy(buffers->name[seq_index], sequence_buffer, SEQUENCE_NAME_LENGTH - 1);

      bzero(buffers->sequence[seq_index], sequence_lengths[seq_index] + 20 + 1);
      bzero(sequence_buffer, MAXSEQ);
      nchar = fscanf(input_fp, "\n%[^\n]", sequence_buffer);

      while (sequence_buffer[0] != '>' && nchar != EOF) {
	      strcat(buffers->sequence[seq_index], sequence_buffer);
	      bzero(sequence_buffer, MAXSEQ);
	      nchar = fscanf(input_fp, "\n%[^\n]", sequence_buffer);
	    }
    }

    seqlen = strlen(buffers->sequence[seq_index]);

    if (options.verbose)
      fprintf(stderr, "Sequence %s is %d letters long\n", buffers->name[seq_index], seqlen);
    sequence_lengths[seq_index] = strlen(buffers->sequence[seq_index]);
    sum_all_seq_lengths += sequence_lengths[seq_index];
    seq_index++;
  }
}

void init_options(void)
{
  options.infile[0]   = '\0';
  options.organism[0] = '\0';
  options.verbose = options.reversed = FALSE;
  options.onlybestflag = TRUE;
  options.min_repeat_percentage = 2;
  options.max_mismatches = -1;
  options.motif_length = 0;
  options.max_results = DEFAULT_MAX_RESULTS;
}
