/***************************************/
/*      Weeder 1.4
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

#define DEFAULT_THRESHOLD 85

static char **sequence, **name;
static int *lens, *nuclei, *revnuclei;
static int winlen = 8, totseq = 0, uu, pp, qq, MAXERR = 0;

static BOOL reversed = FALSE;
static FILE *wee_fp, *html_fp;

int main(int argc, char *argv[])
{
  int maah, bestpos, howmany = 0, howmanynew = 0, lastinseq = -1, nchar, len, partiallen = 0,
    threshold = DEFAULT_THRESHOLD;
  int i, j, k, l, jj, nn;
  int **newmatrix, **matrix;
  BOOL flagnope;

  double percentage, max = 0.0, min = 0.0, mintmp, maxtmp, tmpflo, expected;

  char locatepattern[50];
  char buffer[MAXSEQ];
  char entry[4][6];
  char infile[400] = "stdin";
  char weefile[200] = "";
  char htmlfile[200] = "";

  FILE *input_fp;

  for (i = 0; i < 4; i++) bzero(entry[i], 6);
  i = 1;

  lens = (int *) malloc((MAXSEQ * 2) * sizeof(int));

  /* Conteggio lunghezza sequenze */
  if (argc < 4) {
    fprintf(stderr, "\nUsage: %s inputfile oligo substitutions [threshold] [S]\n", argv[0]);
    exit (0);
  }
  sprintf(infile, "%s", argv[1]);
  sprintf(locatepattern, "%s", argv[2]);
  MAXERR = atoi(argv[3]);

  if ((argc > 4) && (strcmp(argv[4], "S") == 0)) reversed = TRUE;
  if ((argc > 5) && (strcmp(argv[5], "S") == 0)) reversed = TRUE;
  if ((argc > 5) && (strcmp(argv[5], "S") != 0)) threshold = atoi(argv[5]);
  if ((argc > 4) && (strcmp(argv[4], "S") != 0)) threshold = atoi(argv[4]);

  if (threshold == 100) threshold = 99;

  if ((threshold < 0) || (threshold > 100)) {
    fprintf(stderr, "\nPercentage threshold should be between 1 and 100\n");
    exit(0);
  }

  input_fp = fopen(infile, "r");
  if (input_fp == NULL) {
    fprintf (stderr, "\nError in opening file %s\n", infile);
    exit (0);
  }

  partiallen = 0;
  l = 0;
  nchar = 0;
  totseq = 0;

  bzero(buffer, MAXSEQ);
  nchar = fscanf(input_fp, "\n%[^\n]", buffer);

  if (buffer[0] != '>') {
    fprintf(stderr, "\nInput file not in FASTA format\nExiting program\n");
    exit(0);
  }

  matrix = (int **) malloc(strlen(locatepattern) * sizeof(int *));

  for (j = 0; j < strlen(locatepattern); j++)
    matrix[j] = (int *) malloc(4 * sizeof(int));

  for (j = 0; j < strlen(locatepattern); j++)
    for (i = 0; i < 4; i++)
      matrix[j][i] = 0;

  while (nchar != EOF) {
    totseq++;
    nchar = fscanf(input_fp, "\n%[^\n]", buffer);

    while ((buffer[0] != '>') && (nchar != EOF)) {
      partiallen += strlen (buffer);
      bzero (buffer, MAXSEQ);
      nchar = fscanf(input_fp, "\n%[^\n]", buffer);
    }
    lens[l++] = partiallen;
    partiallen = 0;
  }

  /* Fine conteggio lunghezza sequenze */
  fclose(input_fp);

  sprintf(weefile,  FORMAT_LOCATOR_WEE_FILE,  infile, locatepattern);
  sprintf(htmlfile, FORMAT_LOCATOR_HTML_FILE, infile, locatepattern);

  fprintf(stderr,
          "\nRunning locator on file %s\nWriting output on files %s and %s\n",
          infile, weefile, htmlfile);
  input_fp = fopen(infile, "r");

  if (!reversed) sequence = (char **) malloc((totseq + 1) * sizeof(char *));
  else           sequence = (char **) malloc((2 * totseq + 1) * sizeof(char *));

  if (!reversed) name = (char **) malloc((totseq + 1) * sizeof(char *));
  else           name = (char **) malloc((2 * totseq + 1) * sizeof(char *));

  for (i = 0; i < totseq + 1; i++) {
    name[i] = (char *) malloc ((80) * (sizeof (char)));

    if (!reversed) sequence[i] = (char *) malloc((lens[i] + 1 + winlen) * (sizeof(char)));
    else           sequence[i] = (char *) malloc(2 * (lens[i] + 1 + winlen) * (sizeof(char)));
  }

  if (reversed) {
    for (i = 0; i < totseq; i++)
      sequence[i + totseq] = (char *) malloc((lens[i] + 1 + winlen) * sizeof(char));
  }

  if (reversed) {
    for (i = 0; i < totseq; i++)
      name[i + totseq] = (char *) malloc(80 * sizeof(char));
  }
  bzero(buffer, MAXSEQ);
  l = 0;
  nchar = fscanf(input_fp, "\n%[^\n]", sequence[l]);

ERE:

  if (sequence[l][0] == '>') {
      strncpy(name[l], sequence[l], 39);

      if (reversed) {
        strcpy(name[l + totseq], "REVCOMP ");
        strcat(name[l + totseq], name[l]);
      }
      bzero(sequence[l], lens[l]);
      bzero(buffer, MAXSEQ);

      nchar = fscanf(input_fp, "\n%[^\n]", buffer);
      while ((buffer[0] != '>') && (nchar != EOF)) {
        strcat(sequence[l], buffer);
        bzero(buffer, MAXSEQ);
        nchar = fscanf(input_fp, "\n%[^\n]", buffer);
      }
  }
  bzero(sequence[l + 1], lens[l + 1]);

  if (buffer[0] == '>') strcat(sequence[l + 1], buffer);

  len = strlen(sequence[l]);

  for (i = 0; i < lens[l]; i++) {
    sequence[l][i] = tolower (sequence[l][i]);
  }
  lens[l] = strlen(sequence[l]);
  l++;

  if (l < totseq) goto ERE; /* can be done with do {} while (...), don't use goto */

  if (reversed) {
    for (l = 0; l < totseq; l++) {
      for (i = lens[l] - 1; i >= 0; i--) {
        sequence[l + totseq][lens[l] - i - 1] = revcomp(sequence[l][i]);
        sequence[l][lens[l]+lens[l]-i-1] =  revcomp(sequence[l][i]);
      }
      lens[l + totseq] = lens[l];
    }
  }
  fclose(input_fp);

  wee_fp  = fopen(weefile, "w");
  html_fp = fopen(htmlfile, "w");

  fprintf(wee_fp, "\n\n**** LOCATING OLIGO :  %s ****\n\n", locatepattern);
  fprintf(html_fp,
          "<html><body bgcolor = silver font = courier>"
          "<span style='font-size:12.0pt;font-family:courier,garamond,arial,sans-serif;color:black;'>"
          "<br><br><big><b>**** LOCATING OLIGO :  %s ****</b></big><br><br>", locatepattern);

  /****************************************************************/
  /*    INIZIO RICERCA ESATTA                                     */
  /****************************************************************/
  winlen = strlen(locatepattern);
  nuclei = (int *) malloc(winlen * (sizeof (int)));
  revnuclei = (int *) malloc(winlen * (sizeof (int)));

  for (uu = 0; uu < totseq; uu++) {
    bestpos = -1;

    if (!reversed) maah = (lens[uu]) - strlen(locatepattern) + 1;
    else maah = (2 * lens[uu]) - strlen(locatepattern) + 1;

    for (j = 0; j < maah; j++) {
      k = 0;
      flagnope = FALSE;

      for (pp = 0; pp < strlen(locatepattern); pp++)
        if (chartoint(tolower(sequence[uu][j + pp])) < 0)
          flagnope = TRUE;

      if (!flagnope) {
			  for (pp = 0; pp < strlen(locatepattern); pp++) {
          nuclei[pp] = chartoint(tolower(sequence[uu][j + pp]));

          if (tolower(sequence[uu][j + pp]) != tolower(locatepattern[pp])) {
            k++;
          }
        }
			} else k = MAXERR + 1;

      if ((k <= MAXERR) &&(!reversed || ((j < lens[uu]- strlen (locatepattern)+1) || (j >= lens[uu])))) {
			  howmany++;
			  nn = strlen(locatepattern) - 1;
			  pp = 0;

			  for (jj = 0; jj <= nn; jj++)
			    pp += nuclei[nn - jj] * pow_4(jj);

			  for (jj = 0; jj <= nn; jj++) {
          matrix[jj][nuclei[jj]]++;
        }
			  if (lastinseq != uu) {
          lastinseq = uu;
          bestpos = j;
        }
			  k = 0;
			}
    }
  }

  for (pp = 0; pp < strlen(locatepattern); pp++) {
    maxtmp = 0;
    mintmp = 1;

    for (qq = 0; qq < 4; qq++) {
      if (((double) matrix[pp][qq] / howmany) > maxtmp)
        maxtmp = (double) matrix[pp][qq] / howmany;

      if (((double) matrix[pp][qq] / howmany) < mintmp)
        mintmp = (double) matrix[pp][qq] / howmany;
    }
    min += log (mintmp + 0.001);
    max += log (maxtmp);
  }
  
  /* QUI SCRIVE */  
  newmatrix = (int **) malloc(strlen(locatepattern) * sizeof(int *));

  for (j = 0; j < strlen(locatepattern); j++)
    newmatrix[j] = (int *) malloc (4 * sizeof (int));

  fprintf(wee_fp, "\n\n%s with %d substitutions and %d percent threshold\n",
          locatepattern, MAXERR, threshold);

  if (reversed) {
    char *revpat = malloc(sizeof(char) * (strlen(locatepattern) + 1));
    fprintf(wee_fp, "[%s]\n ", reverse_pattern(locatepattern, revpat));
    free(revpat);
  }

  fprintf(html_fp,
          "<span style='font-size:14.0pt;font-family:courier,garamond,arial,sans-serif;color:red'>");
  fprintf(html_fp, "<br><br><b>%s %d substitutions and %d percent threshold </b><br>",
          locatepattern, MAXERR, threshold);
  if (reversed) {
    char *revpat = malloc(sizeof(char) * (strlen(locatepattern) + 1));
    fprintf(html_fp, "<b>[%s]</b><br>", reverse_pattern(locatepattern, revpat));
    free(revpat);
  }
  fprintf(html_fp, "</span><br>");
  fprintf(wee_fp, "\n");

  if (howmany == 0) {
    fprintf(wee_fp, "\nPattern not found\n");
    fclose(wee_fp);

    fprintf(html_fp, "<br>Pattern not found<br>");
    fclose(html_fp);
    exit(0);
  }
  fprintf(wee_fp, "\nBest occurrences (match percentage): \n");
  fprintf(html_fp, "<br>Best occurrences (match percentage): <br>");

  for (uu = 0; uu < strlen(locatepattern); uu++)
    for (j = 0; j < 4; j++)
			newmatrix[uu][j] = 0;

  howmanynew = 0;
  lastinseq = -1;
  expected = 0;

  for (uu = 0; uu < totseq; uu++) {
    fprintf(wee_fp, "\n%s\n", name[uu]);
    fprintf(html_fp,"<br>%10s<br>", name[uu]);
    bestpos = -1;

    if (!reversed) maah = (lens[uu]) - strlen(locatepattern) + 1;
    else maah = (2*lens[uu]) - strlen(locatepattern) + 1;
                       
    for (j = 0; j < maah; j++) {
      k = 0;
      flagnope = FALSE;

      for (pp = 0; pp < strlen(locatepattern); pp++)
        if (chartoint(tolower(sequence[uu][j + pp])) < 0) flagnope = TRUE;

      if (!flagnope) {
				for (pp = 0; pp < strlen(locatepattern); pp++) {
          nuclei[pp] = chartoint(tolower(sequence[uu][j + pp]));
        }
      } else k = MAXERR + 1;

      if ((k <= MAXERR) &&
          (!reversed ||
           ((j < lens[uu] - strlen(locatepattern) + 1) || (j >= lens[uu])))) {

				nn = strlen (locatepattern) - 1;
				tmpflo = 0;

				for (jj = 0; jj <= nn; jj++) {
          tmpflo += log((double) matrix[jj][nuclei[jj]]) - log((double) howmany);
        }
				percentage = (100.0 * (tmpflo - min) / (max - min));

				if (percentage >= threshold) {
          howmanynew++;

          if (!reversed ) {
            fprintf(wee_fp, "+ ");
            fprintf(html_fp, "+ ");
          } else {
            if (j < lens[uu]) {
              fprintf(wee_fp, "+ ");
              fprintf(html_fp, "+ ");
            } else {
              fprintf(wee_fp, "- ");
              fprintf(html_fp, "- ");
            }                             
          }

          if (percentage < 90) {
            fprintf(wee_fp, "[");
            fprintf(html_fp, "[");
          } else {
            fprintf(wee_fp, " ");
            fprintf(html_fp, ".");
          }

          for (jj = 0; jj <= nn; jj++) {

            fprintf(wee_fp, "%c", ALPHABET[nuclei[jj]]);
            fprintf(html_fp, "<span style='color:%s'>%c</span>",
                    HTML_COLORS[nuclei[jj]],
                    ALPHABET[nuclei[jj]]);

            if (percentage >= 90) newmatrix[jj][nuclei[jj]]++;
          }

          if (percentage < 90) {
            fprintf(wee_fp, "]");
            fprintf(html_fp, "]");
          } else {
            fprintf(wee_fp, " ");
            fprintf(html_fp, ".");
          }
				   
          if (reversed) {
            fprintf(wee_fp, "\n");
            fprintf(html_fp, "<br>");

            if (percentage < 90) {
					    fprintf(wee_fp, "  [");
					    fprintf(html_fp, "  [");
					  } else {
					    fprintf(wee_fp, "   ");
					    fprintf(html_fp, "  .");
					  }

            for (jj = nn; jj >= 0; jj--) {
					    fprintf(wee_fp, "%c", revcomp(ALPHABET[nuclei[jj]]));
					    fprintf(html_fp, "%c", revcomp(ALPHABET[nuclei[jj]]));
					  }
            if (percentage < 90) {
					    fprintf(wee_fp, "]");
					    fprintf(html_fp, "]");
					  } else {
					    fprintf(wee_fp, " ");
					    fprintf(html_fp, ".");
					  }

            if (j >= lens[uu]) {
              fprintf(wee_fp, " position %d,", (int) (2 * lens[uu] - j - strlen(locatepattern) + 1));
              fprintf(html_fp, " ---</b> position %d,",
                      (int) (2 * lens[uu] - j - strlen (locatepattern) + 1));
            } else {
              fprintf(wee_fp, " position %d,", j + 1);
              fprintf(html_fp, " ---</b> position %d,", j + 1);
            }                                       
          } else {
            fprintf(wee_fp, " position %d,", j + 1);
            fprintf(html_fp, " ---</b> position %d,", j + 1);
          }
          fprintf(wee_fp, " (%.2f)\n", percentage);
          fprintf(html_fp, " (%.2f)<br>", percentage);
        }

				k = 0;
      }
    }
  }

  fprintf(wee_fp,
          "\n\t                   Frequency Matrix               \n"
          "\t      All Occurrences  \t\t      Best Occurrences  \n"
          "\t    A     C     G     T\t\t    A     C     G     T\n");

  fprintf(html_fp,
          "<br><TABLE WIDTH = 250>"
          "<CAPTION><b>Frequency Matrix</b></CAPTION>"
          "<TR><TH><TH colspan=\"5\">  All Occs <TH> <TH> <TH colspan=\"5\">  Best Occs"
          "<TR align = \"right\"><b><TH> <TH>  <TH><span style='color:%s'>A</span>"
          "<TH><span style='color:%s'>C</span><TH><span style='color:%s'>G</span>"
          "<TH><span style='color:%s'>T</span><TH> <TH> "
          "<TH><span style='color:%s'>A</span><TH><span style='color:%s'>C</span>"
          "<TH><span style='color:%s'>G</span><TH><span style='color:%s'>T</span></b>",
          HTML_COLORS[DNA_A], HTML_COLORS[DNA_C],
          HTML_COLORS[DNA_G], HTML_COLORS[DNA_T],
          HTML_COLORS[DNA_A], HTML_COLORS[DNA_C],
          HTML_COLORS[DNA_G], HTML_COLORS[DNA_T]);

  for (pp = 0; pp < strlen (locatepattern); pp++) {
    sprintf(entry[0], "%d", matrix[pp][0]);
    sprintf(entry[1], "%d", matrix[pp][1]);
    sprintf(entry[2], "%d", matrix[pp][2]);
    sprintf(entry[3], "%d", matrix[pp][3]);

    fprintf(wee_fp, "%d \t%5d %5d %5d %5d\t",
            pp + 1, matrix[pp][0], matrix[pp][1], matrix[pp][2], matrix[pp][3]);
    fprintf(wee_fp, "\t%5d %5d %5d %5d\n",
            newmatrix[pp][0], newmatrix[pp][1], newmatrix[pp][2], newmatrix[pp][3]);
    
    fprintf(html_fp, "<TR ALIGN=\"right\"><TD><b>%d</b> <TD>\t<TD>%d<TD> %d<TD> %d<TD> %d",
            pp + 1, matrix[pp][0], matrix[pp][1], matrix[pp][2], matrix[pp][3]);
    fprintf(html_fp, "<TD>|<TD>\t<TD>%d<TD> %d<TD> %d<TD> %d",
            newmatrix[pp][0], newmatrix[pp][1], newmatrix[pp][2], newmatrix[pp][3]);
  }
  fprintf(wee_fp, "\n==========================================\n");
  fprintf(html_fp, "</TABLE><br><br><HR><br></body></html>");

  fclose(wee_fp);
  fclose(html_fp);
  exit(0);
}
