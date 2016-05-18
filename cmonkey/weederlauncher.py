#!/usr/bin/python
import subprocess as sp
import sys
import argparse

DEFAULT_MAX_RESULTS = 10
WEEDER_TFBS = 'weederTFBS'
ADVISER = 'adviser'

def print_job_info_weefile(inputfile, analysis, orgcode, reverse):
    weefile = inputfile + '.wee'
    with open(weefile, 'a') as outfile:
        outfile.write("\n\nNEW: Starting a %s job\n\nOrganism code: %s\n" % (analysis, orgcode))
        if reverse:
            outfile.write("\nProcessing *both* strands\n")
        outfile.write('weederlauncher %s %s %s' % (inputfile, orgcode, analysis))
        outfile.write("\n\n")

def print_job_info_htmlfile(inputfile, analysis, orgcode, reverse):
    htmlfile = inputfile + '.html'
    with open(htmlfile, 'a') as outfile:
        outfile.write("<html><body bgcolor = silver font = garamond><head><title>Your Weeder Web Results</title></head><tt><br><br><b>NEW: Starting a %s job on file %s</b><br><br><b><br>Organism code: %s<br>" % (analysis, inputfile, orgcode))
        if reverse:
            outfile.write("</b><br>Processing <b>both</b> strands<br><br><b>")
        outfile.write('weederlauncher %s %s %s' % (inputfile, orgcode, analysis))
        outfile.write("</b><br><br>")

def print_job_info(inputfile, analysis, orgcode, reverse):
    print_job_info_weefile(inputfile, analysis, orgcode, reverse)
    print_job_info_htmlfile(inputfile, analysis, orgcode, reverse)


def run_weederTFBS(inputfile, logfile, orgcode, motif_len, num_mutations, max_results, allseqs,
                   reverse, multi, ffdir):
    seq_percentage = '100' if allseqs else '50'
    reverse_param = '-S' if reverse else '-N'
    multi_param = '-M' if multi else '-N'
    command = [WEEDER_TFBS,
               '-f', inputfile,
               '-R', seq_percentage,
               '-O', orgcode,
               '-W', str(motif_len),
               '-e', str(num_mutations),
               multi_param, reverse_param,
               '-T', str(max_results)]
    if ffdir:
        command.extend(['-F', ffdir])
    proc = sp.Popen(command, stdout=logfile, stderr=logfile)
    return proc.wait()


def run_adviser(inputfile, logfile):
    command = [ADVISER, inputfile]
    proc = sp.Popen(command, stdout=logfile, stderr=logfile)
    return proc.wait()


def run_small_analysis(inputfile, orgcode, max_results, reverse, multi, allseqs, ffdir):
    print_job_info(inputfile, 'small', orgcode, reverse)
    with open('weeder.log', 'w') as logfile:
        run_weederTFBS(inputfile, logfile, orgcode, 6, 1, max_results, allseqs,
                       reverse, multi, ffdir)
        run_weederTFBS(inputfile, logfile, orgcode, 8, 2, max_results, allseqs,
                       reverse, multi, ffdir)
        run_adviser(inputfile, logfile)


def run_medium_analysis(inputfile, orgcode, max_results, reverse, multi, allseqs, ffdir):
    print_job_info(inputfile, 'medium', orgcode, reverse)
    with open('weeder.log', 'w') as logfile:
        run_weederTFBS(inputfile, logfile, orgcode, 6, 1, max_results, allseqs,
                       reverse, multi, ffdir)
        run_weederTFBS(inputfile, logfile, orgcode, 8, 2, max_results, allseqs,
                       reverse, multi, ffdir)
        run_weederTFBS(inputfile, logfile, orgcode, 10, 3, max_results, allseqs,
                       reverse, multi, ffdir)
        run_adviser(inputfile, logfile)


def run_large_analysis(inputfile, orgcode, max_results, reverse, multi, allseqs, ffdir):
    print_job_info(inputfile, 'large', orgcode, reverse)
    with open('weeder.log', 'w') as logfile:
        run_weederTFBS(inputfile, logfile, orgcode, 6, 1, max_results, allseqs,
                       reverse, multi, ffdir)
        run_weederTFBS(inputfile, logfile, orgcode, 8, 2, max_results, allseqs,
                       reverse, multi, ffdir)
        run_weederTFBS(inputfile, logfile, orgcode, 10, 3, max_results, allseqs,
                       reverse, multi, ffdir)
        run_weederTFBS(inputfile, logfile, orgcode, 12, 4, max_results, allseqs,
                       reverse, multi, ffdir)
        run_adviser(inputfile, logfile)


def run_extra_analysis(inputfile, orgcode, max_results, reverse, multi, allseqs, ffdir):
    print_job_info(inputfile, 'extra', orgcode, reverse)
    with open('weeder.log', 'w') as logfile:
        run_weederTFBS(inputfile, logfile, orgcode, 6, 1, max_results, allseqs,
                       reverse, multi, ffdir)
        run_weederTFBS(inputfile, logfile, orgcode, 8, 3, max_results, allseqs,
                       reverse, multi, ffdir)
        run_weederTFBS(inputfile, logfile, orgcode, 10, 4, max_results, allseqs,
                       reverse, multi, ffdir)
        run_weederTFBS(inputfile, logfile, orgcode, 12, 4, max_results, allseqs,
                       reverse, multi, ffdir)
        run_adviser(inputfile, logfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="weederlauncher.py")
    parser.add_argument('--input', required=True)
    parser.add_argument('--orgcode', required=True)
    parser.add_argument('--analysis', required=True, help='small|medium|large')
    parser.add_argument('--allseqs', action='store_true',
                        help='all sequences must contain the motif (default: half)')
    parser.add_argument('--reverse', action='store_true',
                        help='process both strands (default: single strand)')
    parser.add_argument('--multi', action='store_true',
                        help='multiple motif occurrences (default: 0 or 1)')
    parser.add_argument('--topresults', type=int, default=DEFAULT_MAX_RESULTS,
                        help='number of results to report')
    parser.add_argument('--ffdir', help='specify alternative FreqFiles directory')
    args = parser.parse_args()

    if args.analysis == 'small':
        run_small_analysis(args.input, args.orgcode, args.topresults,
                           args.reverse, args.multi, args.allseqs, args.ffdir)
    elif args.analysis == 'medium':
        run_medium_analysis(args.input, args.orgcode, args.topresults,
                            args.reverse, args.multi, args.allseqs, args.ffdir)
    elif args.analysis == 'large':
        run_large_analysis(args.input, args.orgcode, args.topresults,
                           args.reverse, args.multi, args.allseqs, args.ffdir)
    elif args.analysis == 'extra':
        run_extra_analysis(args.input, args.orgcode, args.topresults,
                           args.reverse, args.multi, args.allseqs, args.ffdir)
    else:
        print("Analysis type '%s' not supported" % args.analysis)
