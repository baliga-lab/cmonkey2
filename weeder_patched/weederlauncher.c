/***************************************/
/*      Weeder 1.4
        see LICENSE.txt file for
        terms and conditions of use
*/
#include <sys/types.h>
#include <sys/wait.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <strings.h>
#include "weeder_util.h"

#define DEFAULT_MAX_RESULTS 10
#define WAIT_CHILD_EXIT     wait(NULL)
#define MAX_PARAM_LENGTH    200 /* so we can easily process long paths */
#define INPUT_PATH_LENGTH   200

static BOOL allflag = FALSE, reverseflag = FALSE, multiflag = FALSE;
static char organism[8], analysis[80], inputfile[INPUT_PATH_LENGTH] = "";
static int max_results = DEFAULT_MAX_RESULTS;

void parse_option(const char *option);
void run_small_analysis(void);
void run_medium_analysis(void);
void run_large_analysis(void);
void run_adviser(void);
void run_weeder_tfbs(const char *arg0, const char *arg1, const char *arg2,
                     const char *arg3);

void print_job_info_weefile(int argc, char **argv);
void print_job_info_htmlfile(int argc, char **argv);
void print_job_info(int argc, char **argv);
void print_usage(const char *progname);

void check_6mer_freqfile_exists(void);
void check_8mer_freqfile_exists(void);
void check_inputfile_exists(void);
void create_mix_file(void);

int main(int argc, char *argv[])
{
  /*
    LA LENGTH E' IL TIPO DI ANALISI  
    IL 15 E' IL REVERSED (TRUE / FALSE)
    IL 17 E' IL TITOLO
  */
  if ((argc < 4) || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
    print_usage(argv[0]);
    exit(1);
  }
  strcpy(inputfile, argv[1]);
  strcpy(organism,  argv[2]);
  strcpy(analysis,  argv[3]);

  check_8mer_freqfile_exists();
  check_6mer_freqfile_exists();
  check_inputfile_exists();

  create_mix_file();

  if (argc >= 5) parse_option(argv[4]);
  if (argc >= 6) parse_option(argv[5]);
  if (argc >= 7) parse_option(argv[6]);
  if (argc >= 8) parse_option(argv[7]);

  if (max_results <= 0) {
    fprintf(stderr, "\nI need to report at least one motif per run!\n");
    exit(1);
  }

  /* run specific analysis */
  if (!strcmp(analysis, "small")) {
    print_job_info(argc, argv);
    run_small_analysis();
  } else if (!strcmp(analysis, "medium")) {
    print_job_info(argc, argv);
    run_medium_analysis();
  } else if (!strcmp(analysis, "large") || !strcmp(analysis, "extra")) {
    print_job_info(argc, argv);
    run_large_analysis();
  } else {
    fprintf(stderr,
            "\nUnknown type of analysis %s (please choose small, medium, large or extra)\n",
            analysis);
    exit(1);
  }
  fprintf(stderr, "Weeder exited successfully\n");
  return 0;
}

void print_job_info_weefile(int argc, char **argv)
{
  char filename[800] = "";
  FILE  *fp;
  int i;

  sprintf(filename, FORMAT_WEE_FILE, inputfile);
  fp = fopen(filename, "a");
  fprintf(fp, "\n\nNEW: Starting a %s job\n\n"
          "Organism code: %s\n", analysis, organism);
  if (reverseflag) fprintf(fp, "\nProcessing *both* strands\n");
  for (i = 0; i < argc; i++) fprintf(fp, "%s ", argv[i]);
  fprintf(fp, "\n\n");
  fclose(fp);
}

void print_job_info_htmlfile(int argc, char **argv)
{
  char filename[800] = "";
  FILE *fp;
  int i;

  sprintf(filename, FORMAT_HTML_FILE, inputfile);
  fp = fopen(filename, "a");
  fprintf(fp,
          "<html><body bgcolor = silver font = garamond><head>"
          "<title>Your Weeder Web Results</title></head>"
          "<tt><br><br><b>NEW: Starting a %s job on file %s</b><br><br><b>"
          "<br>Organism code: %s<br>", analysis, inputfile, organism);
  if (reverseflag) fprintf(fp, "</b><br>Processing <b>both</b> strands<br><br><b>");
  for (i = 0; i < argc; i++) fprintf(fp, "%s ", argv[i]);
  fprintf(fp, "</b><br><br>");
  fclose(fp);
}

void print_job_info(int argc, char **argv)
{
  print_job_info_weefile(argc, argv);
  print_job_info_htmlfile(argc, argv);
}

void run_weeder_tfbs(const char *arg0, const char *arg1, const char *arg2,
                     const char *arg3)
{
  int excflag, i;
  char *params[20];

  /* common Weeder TFBS parameters */
  for (i = 0; i < 20; i++) {
    params[i] = (char *) malloc(sizeof(char) * MAX_PARAM_LENGTH);
    bzero(params[i], MAX_PARAM_LENGTH);
  }

  strcpy(params[0], WEEDER_TFBS_CMD);
  strcpy(params[1], "-f");
  sprintf(params[2], "%s", inputfile);

  strcpy(params[3], "-R");
  if (allflag) strcpy(params[4], "100");
  else strcpy(params[4], "50");

  strcpy(params[5], "-O");
  sprintf(params[6], "%s", organism);

  strcpy(params[13], "-T");
  sprintf(params[14], "%d", max_results);

  /* specific parameters */
  strcpy(params[7],  arg0);
  strcpy(params[8],  arg1);
  strcpy(params[9],  arg2);
  strcpy(params[10], arg3);

  if (multiflag) strcpy(params[11], "-M");
  else strcpy(params[11], "-N");

  if (reverseflag) strcpy(params[12], "-S");
  else strcpy(params[12], "-N");
  params[15] = NULL;
  excflag = execvp(WEEDER_TFBS_CMD, params);

  if (excflag == -1)
    fprintf (stderr, "\nCould not start weederTFBS. Make sure it is located in the same directory of weederlauncher. %d\n", errno);
}

void run_small_analysis()
{
  pid_t pid = fork();

  if (pid == 0) run_weeder_tfbs("-W", "6", "-e", "1");
  else {
    WAIT_CHILD_EXIT;
    pid = fork();

    if (pid == 0) run_weeder_tfbs("-W", "8", "-e", "2");
    else run_adviser();
  }
}

void run_medium_analysis()
{
  pid_t pid = fork();

  if (pid == 0) run_weeder_tfbs("-W", "6", "-e", "1");
  else {
    WAIT_CHILD_EXIT;
    pid = fork();

    if (pid == 0) run_weeder_tfbs("-W", "8", "-e", "2");
    else {
      WAIT_CHILD_EXIT;
      pid = fork();

      if (pid == 0) run_weeder_tfbs("-W", "10", "-e", "3");
      else run_adviser();
    }
  }
}

void run_large_analysis()
{
  pid_t pid = fork();

  if (pid == 0) run_weeder_tfbs("-W", "6", "-e", "1");
  else {
    WAIT_CHILD_EXIT;
    pid = fork();

    if (pid == 0) {
      run_weeder_tfbs("-W", "8", "-e", (!strcmp(analysis, "large")) ? "2" : "3");
    } else {
      WAIT_CHILD_EXIT;
      pid = fork();

      if (pid == 0) {
        run_weeder_tfbs("-W", "10", "-e", (!strcmp(analysis, "large")) ? "3" : "4");
      } else {
        WAIT_CHILD_EXIT;
        pid = fork();
        if (pid == 0) run_weeder_tfbs("-W", "12", "-e", "4");
        else run_adviser();
      }
    }
  }
}

void run_adviser()
{
  pid_t pid;
  int excflag;
  WAIT_CHILD_EXIT;
  pid = fork();
  if (pid == 0) {
    char *params[4];
    int i;
    for (i = 0; i < 3; i++) {
      params[i] = (char *) malloc(sizeof(char) * MAX_PARAM_LENGTH);
      bzero(params[i], MAX_PARAM_LENGTH);
    }
    params[3] = NULL;

    strcpy(params[0], ADVISER_CMD);
    strcpy(params[1], inputfile);
    if (!reverseflag) strcpy(params[2], "N");
    else strcpy(params[2], "S");

    fputs("\nRunning adviser....\n", stderr);
    excflag = execvp(ADVISER_CMD, params);

    if (excflag == -1)
      fprintf(stderr, "\nCould not start the adviser. Make sure it is located in the same directory of weederlauncher. %d\n", errno);
  } else {
    WAIT_CHILD_EXIT;
    exit(0);
  }
}

void parse_option(const char *option)
{
  if (!strcmp(option, "A")) allflag     = TRUE;
  if (!strcmp(option, "S")) reverseflag = TRUE;
  if (!strcmp(option, "M")) multiflag   = TRUE;
  if (option[0] == 'T') {
    const char *num_str = &(option[1]);
    max_results = atoi(num_str);
  }
}

void print_usage(const char *progname)
{
  fprintf(stderr,
          "\nWeeder 1.4.2 (ISB Patch)"
          "\nUsage : %s inputfilename speciescode analysistype <options>\n"
          "\nspeciescode: two-letter species code (e.g. HS, MM, RN, SC, and so on)\n"
          "analysistype: small|medium|large|extra\n"
          "options (not required):\n"
          "A: all sequences must contain the motif (default: half)\n"
          "S: process both strands of input sequences (default: single strand)\n"
          "M: multiple motif occurrences in each sequence (default: expect zero or one occurrences per sequence)\n"
          "T<number>: report top <number> motifs of each run\n", progname);
}

void check_8mer_freqfile_exists()
{
  char filename[200];
  sprintf(filename, FORMAT_8MER_FREQ_FILE, organism);
  check_file_exists(filename, "\nMissing frequency file : %s\n");
}

void check_6mer_freqfile_exists()
{
  char filename[200];
  sprintf(filename, FORMAT_6MER_FREQ_FILE, organism);
  check_file_exists(filename, "\nMissing frequency file : %s\n");
}

void check_inputfile_exists()
{
  check_file_exists(inputfile, "\nweederluncher: No such file : %s\n");
}

void create_mix_file()
{
  char mixfile[800] = "";
  FILE *fp;

  sprintf(mixfile, FORMAT_MIX_FILE, inputfile);
  fp = fopen(mixfile, "w");
  fclose(fp);
}
