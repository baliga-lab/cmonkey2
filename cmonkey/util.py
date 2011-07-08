"""cMonkey utility module"""


def next_non_comment_index(lines, comment, line_index):
    """utility method that takes a list of strings and returns the next
    line index of a non-comment line
    """
    if comment:
        line = lines[line_index]
        while line.lstrip().startswith(comment):
            line_index += 1
            line = lines[line_index]
    return line_index


def remove_quotes(astring, quote):
    """removes the quote character from the input string if given"""
    if quote:
        return astring.replace(quote, "")
    else:
        return astring


class DelimitedFile:  # pylint: disable-msg=R0913
    """A file class to read text files that are delimited with certain
    separators. This class offers some flexibility over regular csv
    reading mechanisms by allowing for comments and storing optional
    headers.
    Create a DelimitedFile instance by calling DelimitedFile.read()."""
    def __init__(self, lines, header):
        self.lines = lines
        self.header = header

    @classmethod
    def read(cls, filepath, sep='\t', has_header=False, comment=None,
             quote=None):
        """Creates the reader object"""
        file_header = None
        file_lines = []
        lines = None
        with open(filepath) as inputfile:
            lines = inputfile.readlines()

        line_index = next_non_comment_index(lines, comment, 0)
        if has_header:
            file_header = lines[line_index].rstrip().split(sep)
            file_header = [remove_quotes(elem, quote) for elem in file_header]
            line_index += 1

        num_lines = len(lines)
        while line_index < num_lines:
            line_index = next_non_comment_index(lines, comment, line_index)
            if line_index < num_lines:
                line = lines[line_index].rstrip().split(sep)
                line = [remove_quotes(elem, quote) for elem in line]
                file_lines.append(line)
                line_index += 1

        return DelimitedFile(file_lines, file_header)

    def get_lines(self):
        """returns the lines in the file"""
        return self.lines

    def get_header(self):
        """returns the header in the file"""
        return self.header

__all__ = ['DelimitedFile']
