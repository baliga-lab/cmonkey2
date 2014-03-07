/***************************************/
/*      Weeder 1.0
   see LICENSE.txt file for
  terms and conditions of use
*/
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

#define SEQUENCE_NAME_LENGTH_SIMPLE   40
#define SEQUENCE_NAME_LENGTH_MAX (SEQUENCE_NAME_LENGTH_SIMPLE * 2)

typedef struct _brother
{
  char pat[MAX_PATTERN_LENGTH];
  struct _brother *nextbro;
} Brother;

typedef struct _pattern
{
  /* base node */
  struct _pattern *next;
  double score;
  int    counter;
  char   pat[MAX_PATTERN_LENGTH];

  /* adviser specific members */
  BOOL    shadowed;
  int     rank;
  BOOL    overlap;
  int     hamming;
  BOOL    inside;
  BOOL    outside;
  int     hor;
  int     error;
  int     vert;
  double  max;
  double  min;
  int     **matrix;
  int     inshuffled;
  Brother *nextbro;
} Pattern;

typedef struct
{
  char **sequence, **name;
} WorkBuffers;

/* head_nodes and tail_node do not change, no need to malloc() them. */
static Pattern head_nodes[4], tail_node;
static int *sequence_lengths, *nuclei, *revnuclei;
static int uu, pp, qq, winlen = 8, MAXERR = 0;
static BOOL reversed = FALSE;

void init_resultlists(void);
void init_workbuffers(WorkBuffers *buffers, int num_sequences);
Pattern *new_pattern(const char *buffer);
Brother *new_brother(const char *pat, Brother *next);
Pattern *mergesrt(Pattern *list, int length);
void output_results(WorkBuffers *buffers, int num_sequences,
                    const char *weefile, const char *htmlfile, BOOL motif_found);

Pattern *new_pattern(const char *buffer)
{
  int i, j;
  Pattern *result = (Pattern *) malloc(sizeof(Pattern));
  bzero(result->pat, MAX_PATTERN_LENGTH);

  result->hor        = 0;
  result->vert       = 0;
  result->overlap    = FALSE;
  result->hamming    = 0;
  result->inshuffled = 0;
  result->inside     = FALSE;
  result->outside    = FALSE;
  result->shadowed   = FALSE;
  result->nextbro    = NULL;

  sscanf(buffer, "%d) %s %lf %d",
         &result->rank, result->pat, &result->score, &result->error);  
  result->matrix = (int **) malloc(MAX_PATTERN_LENGTH * sizeof(int *));

  for (i = 0; i < MAX_PATTERN_LENGTH; i++) {
    result->matrix[i] = (int *) malloc(ALPHABET_SIZE * sizeof(int));
    for (j = 0; j < ALPHABET_SIZE; j++) result->matrix[i][j] = 0;
  }
  return result;
}

Brother *new_brother(const char *pat, Brother *next)
{
  Brother *result = (Brother *) malloc(sizeof(Brother));
  strcpy(result->pat, pat);
  result->nextbro = next;
  return result;
}

void init_resultlists(void)
{
  int i;
  tail_node.score = TAIL_SCORE;

  for (i = 0; i < 4; i++) {
    head_nodes[i].next = &tail_node;
    head_nodes[i].counter = 0;
  }

}

void init_workbuffers(WorkBuffers *buffers, int num_sequences)
{
  int i;
  /* Reserve enough space to read all sequences (including reversed versions if necessary) */
  if (!reversed) buffers->sequence = (char **) malloc((num_sequences + 1) * sizeof(char *));
  else buffers->sequence = (char **) malloc((2 * num_sequences + 1) * sizeof(char *));

  /* Reserve space for the names */
  if (!reversed) buffers->name = (char **) malloc((num_sequences + 1) * sizeof(char *));
  else buffers->name = (char **) malloc((2 * num_sequences + 1) * sizeof(char *));

  for (i = 0; i < num_sequences + 1; i++) {
    buffers->name[i] = (char *) malloc(SEQUENCE_NAME_LENGTH_MAX * sizeof(char));

    if (!reversed) buffers->sequence[i] = (char *) malloc((sequence_lengths[i] + 1 + winlen) * sizeof(char));
    else buffers->sequence[i] = (char *) malloc(2 * (sequence_lengths[i] + 1 + winlen) * sizeof(char));
  }

  if (reversed) {
    for (i = 0; i < num_sequences; i++) {
      buffers->sequence[i + num_sequences] = (char *) malloc((sequence_lengths[i] + 1 + winlen) * sizeof(char));
      buffers->name[i + num_sequences] = (char *) malloc(SEQUENCE_NAME_LENGTH_MAX * sizeof (char));
    }
  }
}

int main(int argc, char *argv[])
{
  WorkBuffers buffers;
  int nchar, len, partiallen = 0, i, l, num_sequences;
  BOOL flagjump, motif_found = FALSE;

  char buffer[MAXSEQ];
  char fastafile[400] = "stdin";
  char mixfile[400] = "";

  Pattern *newpat, *tmppat;

  char weefile[200] = "";
  char htmlfile[200] = "";
  FILE *fastafile_fp, *mixfile_fp;

  init_resultlists();
  sequence_lengths = (int *) malloc((MAXSEQ * 2) * sizeof(int));

  /* Conteggio lunghezza sequenze */
  sprintf(fastafile, "%s", argv[1]);
  if ((argc > 2) && !strcmp(argv[2], "S")) reversed = TRUE;

  fastafile_fp = fopen(fastafile, "r");

  if (fastafile_fp == NULL) {
    fprintf(stderr, "\nError in opening file %s\n", fastafile);
    exit(0);
  }

  partiallen = 0;
  l = 0;
  nchar = 0;
  num_sequences = 0;

  bzero(buffer, MAXSEQ);
  nchar = fscanf(fastafile_fp, "\n%[^\n]", buffer);

  if (buffer[0] != '>') {
    fprintf(stderr, "\nInput file not in FASTA format\nExiting program\n");
    exit(0);
  }

  /* read sequence lengths */
  while (nchar != EOF) {
    num_sequences++;
    nchar = fscanf(fastafile_fp, "\n%[^\n]", buffer);

    while ((buffer[0] != '>') && (nchar != EOF)) {
      partiallen += strlen(buffer);
      bzero(buffer, MAXSEQ);
      nchar = fscanf(fastafile_fp, "\n%[^\n]", buffer);
    }
    sequence_lengths[l++] = partiallen;
    partiallen = 0;
  }

  /* Fine conteggio lunghezza sequenze */
  fclose(fastafile_fp);

  sprintf(mixfile,  FORMAT_MIX_FILE,  fastafile);
  sprintf(weefile,  FORMAT_WEE_FILE,  fastafile);
  sprintf(htmlfile, FORMAT_HTML_FILE, fastafile);

  fprintf(stderr,
          "\nRunning adviser on file %s\nWriting output on files %s and %s\n",
          mixfile, weefile, htmlfile);

  /* re-read file with known sequence lengths */
  fastafile_fp = fopen(fastafile, "r");

  init_workbuffers(&buffers, num_sequences);

  bzero(buffer, MAXSEQ);
  l = 0;
  nchar = fscanf(fastafile_fp, "\n%[^\n]", buffers.sequence[l]);

ERE:
  if (buffers.sequence[l][0] == '>') {
    strncpy(buffers.name[l], buffers.sequence[l], SEQUENCE_NAME_LENGTH_SIMPLE - 1);

    if (reversed) {
      strcpy(buffers.name[l + num_sequences], "REVCOMP ");
      strcat(buffers.name[l + num_sequences], buffers.name[l]);
    }

    bzero(buffers.sequence[l], sequence_lengths[l]);
    bzero(buffer, MAXSEQ);
    nchar = fscanf(fastafile_fp, "\n%[^\n]", buffer);

    while ((buffer[0] != '>') && (nchar != EOF)) {
      strcat(buffers.sequence[l], buffer);
      bzero(buffer, MAXSEQ);
      nchar = fscanf(fastafile_fp, "\n%[^\n]", buffer);
    }
  }

  bzero(buffers.sequence[l + 1], sequence_lengths[l + 1]);

  if (buffer[0] == '>') strcat(buffers.sequence[l + 1], buffer);

  for (i = 0; i < sequence_lengths[l]; i++) {
    buffers.sequence[l][i] = tolower(buffers.sequence[l][i]);
  }

  sequence_lengths[l] = strlen(buffers.sequence[l]);
  l++;

  if (l < num_sequences) goto ERE; /* should be replaced by do .. while */

  if (reversed) {
    for (l = 0; l < num_sequences; l++) {
      for (i = sequence_lengths[l] - 1; i >= 0; i--) {
        buffers.sequence[l + num_sequences][sequence_lengths[l] - i - 1] = revcomp(buffers.sequence[l][i]);
        buffers.sequence[l][sequence_lengths[l]+sequence_lengths[l]-i-1] =  revcomp(buffers.sequence[l][i]);
      }
      sequence_lengths[l + num_sequences] = sequence_lengths[l];
    }
  }
  fclose(fastafile_fp);

  /* Note that the mixfile is never closed !!*/
  mixfile_fp = fopen(mixfile, "r");

  nchar = fscanf(mixfile_fp, "\n%[^\n]", buffer);

  while (nchar != EOF) {
    flagjump    = FALSE;
    motif_found = TRUE;

    newpat = new_pattern(buffer);

    i = (strlen(newpat->pat) - 6) / 2;
    tmppat = head_nodes[i].next;

    while ((tmppat->score > -1111.0) && !flagjump) {
      if (strcmp(tmppat->pat, newpat->pat) == 0) {
        flagjump = TRUE;

        if (tmppat->rank == newpat->rank) {
          if (tmppat->score < newpat->score) {
            tmppat->score = newpat->score;
            tmppat->error = newpat->error;
          }
        }
	      if (tmppat->rank < newpat->rank) {
          tmppat->score = newpat->score;
          tmppat->error = newpat->error;
        }
	    }
      tmppat = tmppat->next;
    }

    if (!flagjump) {
      newpat->next = head_nodes[i].next;
      head_nodes[i].next = newpat;
      head_nodes[i].counter++;
    } else {
      free(newpat);
    }
    nchar = fscanf(mixfile_fp, "\n%[^\n]", buffer);  
  }

  for (i = 0; i < 4; i++) {
    if (head_nodes[i].counter > 0) {
      head_nodes[i].next = mergesrt(head_nodes[i].next, head_nodes[i].counter);
    }
  }

  output_results(&buffers, num_sequences, weefile, htmlfile, motif_found);
  fprintf(stderr, "\nJob completed\n");
  exit(0);
}

void output_results(WorkBuffers *buffers, int num_sequences,
                    const char *weefile, const char *htmlfile,
                    BOOL motif_found)
{
  FILE *wee_fp, *html_fp;
  int **newmatrix;
  double percentage, expected, tmpflo, mintmp, maxtmp, newmin, newmax;
  BOOL flagnope, has_advice = TRUE;
  int i, j, k, nn, jj, bestpos, maah, rankcount, patlimit = 4, howmany = 0, lastinseq = -1;
  Pattern *pointer[4], *seekpat;
  Pattern *tmppat;
  char entry[4][6];

  for (i = 0; i < 4; i++) bzero(entry[i], 6);
  wee_fp = fopen(weefile, "a");
  html_fp = fopen(htmlfile, "a");

  if (!motif_found) {
    fprintf(wee_fp, "\n\nSorry, no interesting motif found\n\n");
    fprintf(html_fp,
            "<span style='font-size:12.0pt;font-family:courier,garamond,arial,sans-serif;color:green;'>");
    fprintf(html_fp,
            "<br><br><b>Sorry, no interesting motif found</b><br><br>\n");
    fprintf(html_fp, "</span>");
    fclose(wee_fp);
    fclose(html_fp);
    exit(22);
  }

  fprintf(wee_fp, "\nYour sequences:\n\n");
  fprintf(html_fp,
          "\n<br>Your sequences:<br><br>\n\n"
          "<span style='font-size:10.0pt;font-family:courier,garamond,arial,sans-serif;color:black;'>");

  for (j = 0; j < num_sequences; j++) {
    fprintf(wee_fp,  "Sequence %d : %s\n", j + 1, buffers->name[j]);
    fprintf(html_fp, "Sequence %d : %s<br>\n", j + 1, buffers->name[j]);

  }
  fprintf(wee_fp, "\n\n**** MY ADVICE ****\n\n");
  fprintf(html_fp,
          "</span>"
          "<span style='font-size:12.0pt;font-family:courier,garamond,arial,sans-serif;color:green;'>"
          "<br><br><big><b>**** MY ADVICE ****</b></big><br><br>\n"
          "</span>");

  for (j = 0; j < 4; j++) {
    tmppat = head_nodes[j].next;

    while (tmppat->score > -111.0) {
      for (i = j; i < 4; i++) {
	      seekpat = head_nodes[i].next;

	      while (seekpat->score > -111.0) {
          if (strcmp (seekpat->pat, tmppat->pat) != 0) {
            if (overlap(seekpat->pat, tmppat->pat, reversed)
                && (tmppat->rank < seekpat->rank)) {

              tmppat->vert = 1;
              seekpat->vert = 1;
              seekpat->overlap = TRUE;
              tmppat->overlap  = TRUE;
              seekpat->inshuffled++;
              tmppat->inshuffled++;

              seekpat->nextbro = new_brother(tmppat->pat, seekpat->nextbro);
              tmppat->nextbro  = new_brother(seekpat->pat, tmppat->nextbro);
            }

            if ((hamming_distance(seekpat->pat, tmppat->pat, reversed) <= seekpat->error) &&
                (tmppat->rank < seekpat->rank)) {

              tmppat->hamming = 1;
              tmppat->vert = 1;
              tmppat->inshuffled++;

              seekpat->hamming = 1;
              seekpat->vert = 1;
              seekpat->inshuffled++;

              seekpat->nextbro = new_brother(tmppat->pat, seekpat->nextbro);
              tmppat->nextbro  = new_brother(seekpat->pat, tmppat->nextbro);
            }
            /* BUG: There is something wrong here !!!! */
            if (strlen(seekpat->pat) != strlen(tmppat->pat)) {
              if (inside(tmppat->pat, seekpat->pat, reversed)) {
                tmppat->hor      = 1;
                seekpat->hor     = 1;
                tmppat->inside   = TRUE;
                seekpat->outside = TRUE;
                seekpat->inshuffled++;
                tmppat->inshuffled++;

                seekpat->nextbro = new_brother(tmppat->pat, seekpat->nextbro);
                tmppat->nextbro  = new_brother(seekpat->pat, tmppat->nextbro);
              }
            }
          }
          seekpat = seekpat->next;
        }
	    }
      tmppat = tmppat->next;
    }
  }

  /* Start exact search */
  nuclei = (int *) malloc(winlen * (sizeof(int)));
  revnuclei = (int *) malloc(winlen * (sizeof(int)));

  for (i = 0; i < 4; i++) {
    seekpat = &(head_nodes[i]);

    while (seekpat->next->score > -111.0) {
      MAXERR = seekpat->next->error;

      howmany = 0;
      lastinseq = -1;
      expected = 0;

      if ((((seekpat->next->hor > 0) || (seekpat->next->vert > 2)) &&
           ((seekpat->next->vert > 0))) ||
          ((seekpat->next->rank == 1) && (strlen(seekpat->next->pat) > 6))) {
	      for (uu = 0; uu < num_sequences; uu++) {
          bestpos = -1;

          if (!reversed) maah = sequence_lengths[uu] - strlen(seekpat->next->pat) + 1;
          else maah = (2 * sequence_lengths[uu]) - strlen(seekpat->next->pat) + 1;

          for (j = 0; j < maah; j++) {
            k = 0;
            flagnope = FALSE;

            for (pp = 0; pp < strlen(seekpat->next->pat); pp++)
              if (chartoint(tolower(buffers->sequence[uu][j + pp])) < 0) flagnope = TRUE;

            if (!flagnope) {
              for (pp = 0; pp < strlen(seekpat->next->pat); pp++) {
                nuclei[pp] = chartoint(tolower(buffers->sequence[uu][j + pp]));

                if (tolower(buffers->sequence[uu][j + pp]) !=
                    tolower(seekpat->next->pat[pp])) {
                  k++;
                }
              }
            } else k = MAXERR + 1;

            if ((k <= MAXERR) &&
                (!reversed ||
                 ((j < sequence_lengths[uu] - strlen(seekpat->next->pat) + 1) || 
                  (j >= sequence_lengths[uu])))) {
              howmany++;

              nn = strlen(seekpat->next->pat) - 1;
              pp = 0;

              for (jj = 0; jj <= nn; jj++)
                pp += nuclei[nn - jj] * pow_4(jj);

              for (jj = 0; jj <= nn; jj++) {
                seekpat->next->matrix[jj][nuclei[jj]]++;
              }

              if (lastinseq != uu) {
                lastinseq = uu;
                bestpos = j;
              }
              k = 0;
            }
          }
        }
	      seekpat->next->min = 0;
	      seekpat->next->max = 0;

	      for (pp = 0; pp < strlen(seekpat->next->pat); pp++) {
          maxtmp = 0;
          mintmp = 1;

          for (qq = 0; qq < ALPHABET_SIZE; qq++) {

            if (((double) seekpat->next->matrix[pp][qq] / howmany) > maxtmp)
              maxtmp = (double) seekpat->next->matrix[pp][qq] / howmany;

            if (((double) seekpat->next->matrix[pp][qq] / howmany) < mintmp)
              mintmp = (double) seekpat->next->matrix[pp][qq] / howmany;
          }
          seekpat->next->min += log(mintmp + 0.001);
          seekpat->next->max += log(maxtmp);
        }
	      seekpat->next->counter = howmany;
	    }
      seekpat = seekpat->next;
    }
  }
  /* write results */  
  for (i = 0; i < 4; i++) pointer[i] = &(head_nodes[i]);

  if (head_nodes[2].counter <= 0) patlimit = 2;

  for (rankcount = 0; rankcount < head_nodes[0].counter / 2; rankcount++) {
    for (i = 0; i < 4; i++) {
      if ((rankcount == 0) && (i == 0)) {

	      fprintf(html_fp,
                "<span style='font-size:12.0pt;font-family:courier,garamond,arial,sans-serif;color:green'>"
                "<br><b> *** Interesting motifs (highest-ranking) seem to be : </b><br></span>");

	      fprintf(wee_fp, "\n *** Interesting motifs (highest-ranking) seem to be : \n");
	      has_advice = FALSE;
	    }
      if ((rankcount == 1) && (i == 0)) {
	      fprintf(html_fp,
                "<span style='font-size:12.0pt;font-family:courier,garamond,arial,sans-serif;color:green'>"
                "<br><b> *** Interesting motifs (not highest-ranking) can also be : </b><br></span>");

	      fprintf(wee_fp,
                "\n *** Interesting motifs (not highest-ranking) can also be : \n");
	      has_advice = FALSE;
	    }

      if (pointer[i]->next->score > 0) {
	      seekpat = pointer[i];
	      pointer[i] = pointer[i]->next;

	      /* while (seekpat->next->score > -111.0) */
	      {
          if ((((seekpat->next->hor > 0) ||
                (seekpat->next->vert > 2)) &&
               ((seekpat->next->vert > 0))) ||
              ((seekpat->next->rank == 1) &&
               (strlen(seekpat->next->pat) > 6))) {

            newmatrix = (int **) malloc(MAX_PATTERN_LENGTH * sizeof(int *));

            for (j = 0; j < MAX_PATTERN_LENGTH; j++)
              newmatrix[j] = (int *) malloc(ALPHABET_SIZE * sizeof(int));

            fprintf(wee_fp, "\n\n%s\n", seekpat->next->pat);
            if (reversed) {
              char *revpat = malloc(sizeof(char) *
                                    (strlen(seekpat->next->pat) + 1));
              fprintf(wee_fp, "%s\n ",
                      reverse_pattern(seekpat->next->pat, revpat));
              free(revpat);
            }

            fprintf(html_fp,
                    "<span style='font-size:14.0pt;font-family:courier,garamond,arial,sans-serif;color:red'>");
            fprintf(html_fp, "<br><br><b> %s </b><br>",
                    seekpat->next->pat);
            if (reversed) {
              char *revpat = malloc(sizeof(char) *
                                    (strlen(seekpat->next->pat) + 1));
              fprintf(html_fp, "<b>%s</b><br>",
                      reverse_pattern(seekpat->next->pat, revpat));
              free(revpat);
            }

            fprintf(wee_fp, "\n%d redundant motifs found:\n", seekpat->next->inshuffled);

            fprintf(html_fp, "</span>");
            fprintf(html_fp, "<br>%d redundant motifs found:<br>", seekpat->next->inshuffled);

            while (seekpat->next->nextbro != NULL) {
              fprintf(wee_fp,  "%s - ", seekpat->next->nextbro->pat);
              fprintf(html_fp, "%s - ", seekpat->next->nextbro->pat);

              seekpat->next->nextbro = seekpat->next->nextbro->nextbro;
            }
            fprintf(wee_fp, "\n\nBest occurrences: \nSeq  St  oligo  pos  match\n");
            fprintf(html_fp,
                    "<br><br>Best occurrences (match percentage): <br>"
                    "Seq  St  oligo  pos  match\n<br>");
            MAXERR = seekpat->next->error;
            
            for (uu = 0; uu < MAX_PATTERN_LENGTH; uu++)
              for (j = 0; j < ALPHABET_SIZE; j++) newmatrix[uu][j] = 0;

            howmany = 0;
            lastinseq = -1;
            expected = 0;

            for (uu = 0; uu < num_sequences; uu++) {
              bestpos = -1;
              if (!reversed)
                maah = (sequence_lengths[uu]) - strlen(seekpat->next->pat) + 1;
              else
                maah = (2*sequence_lengths[uu]) - strlen(seekpat->next->pat) + 1;
                       

              for (j = 0; j < maah; j++) {
                k = 0;
                flagnope = FALSE;

                for (pp = 0; pp < strlen(seekpat->next->pat); pp++)
                  if (chartoint(tolower(buffers->sequence[uu][j + pp])) < 0)
                    flagnope = TRUE;

                if (!flagnope) {
                  for (pp = 0; pp < strlen(seekpat->next->pat); pp++) {
                    nuclei[pp] = chartoint(tolower(buffers->sequence[uu][j + pp]));
                  }
                }
                else k = MAXERR + 1;

                if ((k <= MAXERR) &&
                    (!reversed ||
                     ((j < sequence_lengths[uu] - strlen(seekpat->next->pat)+1) ||
                      (j >= sequence_lengths[uu])))) {

                  nn = strlen(seekpat->next->pat) - 1;
                  tmpflo = 0;

                  for (jj = 0; jj <= nn; jj++) {
                    tmpflo +=
                      log((double) seekpat->next->matrix[jj][nuclei[jj]]) -
                      log((double) seekpat->next->counter);
                  }

                  percentage =
                    (100.0 * (tmpflo - seekpat->next->min) /
                     (seekpat->next->max - seekpat->next->min));

                  if (percentage > 85) {
                    howmany++;

                    if (!reversed) {
                      fprintf(wee_fp,  "%3d + ", uu + 1);
                      fprintf(html_fp, "%3d +<b> ", uu + 1);
                    } else {
                      if (j < sequence_lengths[uu]) {
                        fprintf(wee_fp,  "%3d + ", uu + 1);
                        fprintf(html_fp, "%3d +<b> ", uu + 1);
                      } else {
                        fprintf(wee_fp,  "%3d - ", uu + 1);
                        fprintf(html_fp, "%3d -<b> ", uu + 1);
                      }                                         
                    }
                    if (percentage < 90) {
                      fprintf(wee_fp,  "[");
                      fprintf(html_fp, "[");
                    } else {
                      fprintf(wee_fp,  " ");
                      fprintf(html_fp, ".");
                    }

                    for (jj = 0; jj <= nn; jj++) {

                      fprintf(wee_fp,  "%c", ALPHABET[nuclei[jj]]);
                      fprintf(html_fp, "<span style='color:%s'>%c</span>",
                              HTML_COLORS[nuclei[jj]], ALPHABET[nuclei[jj]]);

                      if (percentage >= 90) newmatrix[jj][nuclei[jj]]++;
                    }

                    if (percentage < 90) {
                      fprintf(wee_fp,  "]");
                      fprintf(html_fp, "]");
                    } else {
                      fprintf(wee_fp,  " ");
                      fprintf(html_fp, ".");
                    }
                    if (reversed) {

                      if (j > sequence_lengths[uu]) {
                        fprintf(wee_fp, " %4d ", (int)
                                (2 * sequence_lengths[uu] - j - strlen(seekpat->next->pat) + 1));
                        fprintf(html_fp, "</b> %4d ", (int)
                                (2 * sequence_lengths[uu] - j - strlen(seekpat->next->pat) + 1));
                      } else {
                        fprintf(wee_fp,  " %4d ", j + 1);
                        fprintf(html_fp, "</b> %4d ", j + 1);
                      }
                    } else {
                      fprintf(wee_fp,  " %4d ", j + 1);
                      fprintf(html_fp, "</b> %4d ", j + 1);
                    }
                    fprintf(wee_fp, " (%.2f) ", percentage);
                    fprintf(wee_fp, "\n");

                    fprintf(html_fp, " (%.2f) ", percentage);
                    fprintf(html_fp, "<br>\n");
                  }
                  k = 0;
                }
              }
            }
            newmin = 0;
            newmax = 0;

            for (pp = 0; pp < strlen(seekpat->next->pat); pp++) {
              maxtmp = 0;
              mintmp = 1;

              for (qq = 0; qq < ALPHABET_SIZE; qq++) {

                if (((double) newmatrix[pp][qq] / howmany) > maxtmp)
                  maxtmp = (double) newmatrix[pp][qq] / howmany;

                if (((double) newmatrix[pp][qq] / howmany) < mintmp)
                  mintmp = newmatrix[pp][qq] / howmany;
              }

              newmin += log(mintmp + 0.001);
              newmax += log(maxtmp);
            }
            fprintf(wee_fp,
                    "\n\t                   Frequency Matrix               \n"
                     "\t      All Occurrences  \t\t      Best Occurrences  \n"
                     "\t    A     C     G     T\t\t    A     C     G     T\n");

            fprintf(html_fp,
                    "<br><TABLE WIDTH = 250>"
                    "<CAPTION><b>Frequency Matrix</b></CAPTION>"
                    "<TR><TH><TH colspan=\"5\">  All Occs <TH> <TH> <TH colspan=\"5\">  Best Occs"
                    "<TR align = \"right\"><b><TH> <TH>  "
                    "<TH><span style='color:%s'>A</span>"
                    "<TH><span style='color:%s'>C</span>"
                    "<TH><span style='color:%s'>G</span><TH><span style='color:%s'>T</span>"
                    "<TH> <TH> <TH><span style='color:%s'>A</span>"
                    "<TH><span style='color:%s'>C</span><TH><span style='color:%s'>G</span>"
                    "<TH><span style='color:%s'>T</span></b>",
                    HTML_COLORS[DNA_A], HTML_COLORS[DNA_C],
                    HTML_COLORS[DNA_G], HTML_COLORS[DNA_T],
                    HTML_COLORS[DNA_A], HTML_COLORS[DNA_C],
                    HTML_COLORS[DNA_G], HTML_COLORS[DNA_T]);

            for (pp = 0; pp < strlen(seekpat->next->pat); pp++) {

              sprintf(entry[DNA_A], "%d", seekpat->next->matrix[pp][DNA_A]);
              sprintf(entry[DNA_C], "%d", seekpat->next->matrix[pp][DNA_C]);
              sprintf(entry[DNA_G], "%d", seekpat->next->matrix[pp][DNA_G]);
              sprintf(entry[DNA_T], "%d", seekpat->next->matrix[pp][DNA_T]);

              fprintf(wee_fp, "%d \t%5d %5d %5d %5d\t", pp + 1,
                      seekpat->next->matrix[pp][DNA_A],
                      seekpat->next->matrix[pp][DNA_C],
                      seekpat->next->matrix[pp][DNA_G],
                      seekpat->next->matrix[pp][DNA_T]);
              fprintf(wee_fp, "\t%5d %5d %5d %5d\n",
                      newmatrix[pp][DNA_A],
                      newmatrix[pp][DNA_C],
                      newmatrix[pp][DNA_G],
                      newmatrix[pp][DNA_T]);

              fprintf(html_fp,
                      "<TR ALIGN=\"right\"><TD><b>%d</b> <TD>\t<TD>%d<TD> %d<TD> %d<TD> %d",
                      pp + 1,
                      seekpat->next->matrix[pp][DNA_A],
                      seekpat->next->matrix[pp][DNA_C],
                      seekpat->next->matrix[pp][DNA_G],
                      seekpat->next->matrix[pp][DNA_T]);
              fprintf(html_fp, "<TD>|<TD>\t<TD>%d<TD> %d<TD> %d<TD> %d",
                      newmatrix[pp][DNA_A],
                      newmatrix[pp][DNA_C],
                      newmatrix[pp][DNA_G],
                      newmatrix[pp][DNA_T]);
            }
            fprintf(wee_fp, "\n==========================================\n");
            fprintf(html_fp, "</TABLE><br>");
            fprintf(html_fp, "<br><HR><br>");
          }
	      }
	    }
    }
  }

  if (has_advice) {
    fprintf(wee_fp,  "\nSorry! No advice on this one.....\n");
    fprintf(html_fp, "<br>Sorry! No advice on this one.....<br>");
  }
  fprintf(html_fp, "</body></html>");
  fclose(wee_fp);
  fclose(html_fp);
}

Pattern *mergesrt(Pattern *list, int length)
{
  return (Pattern *) mergesort_linkedlist((BaseNode *) list,
                                          (BaseNode *) &tail_node,
                                          length);
}
