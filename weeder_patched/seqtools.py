HAMMING_MAX = 9999

def read_sequences_from_fasta_string(fasta_string):
    """reads the sequences contained in a FASTA string"""
    lines = fasta_string.split('\n')
    sequences = []
    seqbuffer = ""
    seqname = None
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if len(seqbuffer) > 0:
                sequences.append((seqname, seqbuffer))
                seqbuffer = ""
            seqname = line[1:]
        elif line and len(line) > 0:
            seqbuffer += line
    # add the last line
    if len(seqbuffer) > 0:
        sequences.append((seqname, seqbuffer))
    return sequences


def read_sequences_from_fasta_file(filepath):
    """Read the sequences from the specified FASTA file"""
    with open(filepath) as inputfile:
        fasta_string = inputfile.read()
    return read_sequences_from_fasta_string(fasta_string)

def revcomp(sequence):
    """compute the reverse complement of the input string"""
    return "".join([revchar(c) for c in sequence[::-1]])

def overlap(str1, str2, checkreverse):
    result = False
    overlapping = True

    for l in range(1, 3):
        for i in range(len(str1) - l):
            if i >= len(str2) or str1[i + l] != str2[i]:
                overlapping = False
                break
        if overlapping:
            result = True
        overlapping = True

        for i in range(len(str1) - l):
            if (i + l) >= len(str2) or str1[i] != str2[i + l]:
                overlapping = False
                break
        if overlapping:
            result = True
    if checkreverse:
        rev_result = overlap(str1[::-1], str2, False)
        if rev_result:
            result = True
    return result


def hamming_distance(str1, str2, checkreverse):
    dist_forward = 0
    dist_reverse = HAMMING_MAX
    if len(str1) != len(str2) or str1 == str2:
        return HAMMING_MAX
        
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            dist_forward += 1

    if not checkreverse:
        return dist_forward
    else:
        rev = str1[::-1]
        for i in range(len(str1)):
            if rev[i] != str2[i]:
                dist_reverse += 1
    if dist_reverse < dist_forward:
        return dist_reverse
    else:
        return dist_forward

def inside(str1, str2, checkreverse):
    len1 = len(str1)
    len2 = len(str2)
    result = False

    if (len2 - len1) != 2:
        return False

    for i in range(len2 - len1 + 1):
        match = True
        for j in range(i, i + len1):
            if str1[j - i] != str2[j]:
                match = False
                break
        if match:
            result = True

    if checkreverse:
        rev_result = inside(str1[::-1], str2, False)
        if rev_result:
            result = True

    return result

def char_to_int(c):
    c = c.lower()
    if c == 'a':
        return 0;
    elif c == 'c':
        return 1;
    elif c == 'g':
        return 2;
    elif c == 't':
        return 3;
    elif c == '$':
        return 4;
    else:
        return -1;
