"""util.py - extraction helpers

This file is part of cm2meme. Please see README and LICENSE for
more information and licensing details.
"""
import re


def extract_regex(pattern, infoline):
    """generic info line field extraction based on regex"""
    try:
        match = re.search(pattern, infoline)
        return infoline[match.start():match.end()].split('=')[1].strip()
    except:
        logging.error("ERROR in __extract_regex(), pattern: '%s', infoline: '%s'",
                      str(pattern), str(infoline))


def next_regex_index(pat, start_index, lines):
    """finds the line index of the first occurrence of the pattern"""
    line_index = start_index
    pattern = re.compile(pat)
    current_line = lines[line_index]
    while not pattern.match(current_line):
        line_index += 1
        if line_index >= len(lines):
            return -1
        current_line = lines[line_index]
    return line_index
