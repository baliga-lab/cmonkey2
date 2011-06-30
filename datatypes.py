""" Data types for cMonkey """

class DataMatrix:
    """
    A 2 dimensional data matrix class, with optional row and column names
    The matrix is initialized with a fixed dimension size and can not
    be resized after initialization.
    Names and values of a matrix instance can be modified
    DataMatrix does not make any assumptions about its value or name types
    """
    def __init__(self, nrows, ncols, row_names=None, column_names=None):
        """create a DataMatrix instance"""
        self.values = [[0.0 for _ in range(ncols)] for _ in range(nrows)]
        if not row_names:
            self.row_names = ["Row " + str(i) for i in range(nrows)]
        else:
            if len(row_names) != nrows:
                raise ValueError("number of row names should be %d" % nrows)
            self.row_names = row_names

        if not column_names:
            self.column_names = ["Column " + str(i) for i in range(ncols)]
        else:
            if len(column_names) != ncols:
                raise ValueError("number of column names should be %d" % ncols)
            self.column_names = column_names

    def num_rows(self):
        """returns the number of rows"""
        return len(self.values)
    def num_columns(self):
        """returns the number of columns"""
        if self.num_rows() == 0:
            return 0
        else:
            return len(self.values[0])
    def value_at(self, row, column):
        """retrieve the value at the specified position"""
        return self.values[row][column]
    def set_value_at(self, row, column, value):
        """set the value at the specified position"""
        self.values[row][column] = value
    def row_name(self, row):
        """retrieve the name for the specified row"""
        return self.row_names[row]
    def column_name(self, row):
        """retrieve the name for the specified column"""
        return self.column_names[row]
