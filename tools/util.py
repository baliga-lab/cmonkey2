import string

def make_rname(name, flag):
    def is_int(s):
        try:
            int(s)
            return True
        except:
            return False

    """emulate R column replacement that read.delim does by default"""
    if flag:
        if is_int(name):
            return 'X' + name
        return string.replace(name, '-', '.')
    else:
        return name

