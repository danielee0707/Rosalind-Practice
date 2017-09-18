def count_nucleotides(file_name):
    with open(file_name, 'r') as f:
        sequence = f.read().rstrip()
    my_dict = {'adenine': 0, 'cytosine': 0, 'guanine': 0, 'thymine': 0}
    for n in sequence:
        if n == 'A':
            my_dict['adenine'] += 1
        elif n == 'C':
            my_dict['cytosine'] += 1
        elif n == 'G':
            my_dict['guanine'] += 1
        elif n == 'T':
            my_dict['thymine'] += 1
    return my_dict.items()


def transcribe(file_name):
    with open(file_name, 'r') as f:
        sequence = f.read().rstrip()
    sequence = sequence.replace('T', 'U')
    return sequence


def complement(file_name):
    comp_seq = ''
    with open(file_name, 'r') as f:
        sequence = f.read().rstrip()
    sequence = sequence[::-1]
    for n in sequence:
        if n == 'A':
            comp_seq += 'T'
        elif n == 'C':
            comp_seq += 'G'
        elif n == 'G':
            comp_seq += 'C'
        elif n == 'T':
            comp_seq += 'A'
    return comp_seq


# n: n month
# k: k pairs of offspring for each pair
def fib(n, k):
    if n <= 2:
        return 1
    else:
        return fib(n - 1, k) + k * fib(n - 2, k)


def gc_content(file_name):
    sequences = {}
    seq_label = ''
    total_cg = 0
    winner = ''
    record = 0
    with open(file_name, 'r') as f:
        lines = f.read().splitlines()
    for i in range(len(lines)):
        if lines[i][0] == '>':
            seq_label = lines[i][1:]
            sequences[seq_label] = ''
        else:
            sequences[seq_label] += lines[i]
    for item in sequences.items():
        sequence = item[1]
        for n in sequence:
            if n == 'C' or n == 'G':
                total_cg += 1
        if float(total_cg) / len(sequence) >= record:
            record = float(total_cg) / len(sequence)
            winner = item[0]
        total_cg = 0
    return winner + '\n' + str(record * 100)


def count_snp(file_name):
    with open(file_name, 'r') as f:
        sequence1 = f.readline().rstrip()
        sequence2 = f.readline().rstrip()
    count = 0
    for i in range(len(sequence1)):
        if sequence1[i] != sequence2[i]:
            count += 1
    return count


# k: number of homogeneous dominant
# m: number of heterogeneous
# n: number of homogeneous recessive
def mendel(k, m, n):
    k = float(k)
    m = float(m)
    n = float(n)
    total = k + m + n
    a = k * (k - 1) + 2 * k * m + 2 * k * n + m * (
        m - 1) * 0.75 + 2 * m * n * 0.5
    b = total * (total - 1)
    return a / b


def translate(file_name):
    codons = []
    start_index = 0
    peptide = []
    with open(file_name, 'r') as f:
        sequence = f.read()
    for i in range(len(sequence)):
        if sequence[i:(i + 3)] == 'AUG':
            start_index = i
            break
    for i in range(start_index, len(sequence), 3):
        codons.append(sequence[i:(i + 3)])
    if len(codons) == 0:
        return ''
    if len(codons[-1]) != 3:
        codons.pop()
    for codon in codons:
        if codon[0] == 'U':
            if codon in ['UUU', 'UUC']:
                peptide.append('F')
            elif codon in ['UUA', 'UUG']:
                peptide.append('L')
            elif codon[1] == 'C':
                peptide.append('S')
            elif codon in ['UAU', 'UAC']:
                peptide.append('Y')
            elif codon in ['UAA', 'UAG', 'UGA']:
                return peptide.join('')
            elif codon in ['UGU', 'UGC']:
                peptide.append('C')
            elif codon == 'UGG':
                peptide += 'W'
        elif codon[0] == 'C':
            if codon[1] == 'U':
                peptide.append('L')
            elif codon[1] == 'C':
                peptide.append('P')
            elif codon in ['CAU', 'CAC']:
                peptide.append('H')
            elif codon in ['CAA', 'CAG']:
                peptide.append('Q')
            elif codon[1] == 'G':
                peptide.append('R')
        elif codon[0] == 'A':
            if codon in ['AUU', 'AUC', 'AUA']:
                peptide.append('I')
            elif codon == 'AUG':
                peptide.append('M')
            elif codon[1] == 'C':
                peptide.append('T')
            elif codon in ['AAU', 'AAC']:
                peptide.append('N')
            elif codon in ['AAA', 'AAG']:
                peptide.append('K')
            elif codon in ['AGU', 'AGC']:
                peptide.append('S')
            elif codon in ['AGA', 'AGG']:
                peptide.append('R')
        elif codon[0] == 'G':
            if codon[1] == 'U':
                peptide.append('V')
            elif codon[1] == 'C':
                peptide.append('A')
            elif codon in ['GAU', 'GAC']:
                peptide.append('D')
            elif codon in ['GAA', 'GAG']:
                peptide.append('E')
            elif codon[1] == 'G':
                peptide.append('G')
    return peptide.join('')


def find_motif(file_name):
    positions = []
    with open(file_name, 'r') as f:
        sequence = f.readline().rstrip()
        motif = f.readline().rstrip()
    for i in range(len(sequence)):
        if sequence[i:(i + len(motif))] == motif:
            positions.append(str(i + 1))
    return ' '.join(positions)


def profiling(file_name):
    index = -1
    sequences = []
    consensus = []
    result = []
    max_num = 0
    cons_nuc = []
    with open(file_name, 'r') as f:
        lines = f.read().splitlines()
    for i in range(len(lines)):
        if lines[i][0] == '>':
            index += 1
            sequences.append('')
        else:
            sequences[index] += lines[i]
    length = len(sequences[0])
    list_a = [0] * length
    list_t = [0] * length
    list_c = [0] * length
    list_g = [0] * length
    profile = {'A': list_a, 'T': list_t, 'C': list_c, 'G': list_g}
    for i in range(len(sequences[0])):
        for sequence in sequences:
            if sequence[i] == 'A':
                list_a[i] += 1
            elif sequence[i] == 'T':
                list_t[i] += 1
            elif sequence[i] == 'C':
                list_c[i] += 1
            else:
                list_g[i] += 1
        for key, value in profile.items():
            if value[i] > max_num:
                max_num = value[i]
                cons_nuc = key
        consensus.append(cons_nuc)
        max_num = 0
    result.append(''.join(consensus))
    result.append('A: ' + ' '.join(str(x) for x in profile['A']))
    result.append('C: ' + ' '.join(str(x) for x in profile['C']))
    result.append('G: ' + ' '.join(str(x) for x in profile['G']))
    result.append('T: ' + ' '.join(str(x) for x in profile['T']))
    return result


# o: number of overlapping nucleotides
def overlap(file_name, o):
    index = -1
    sequences = []
    seq_labels = []
    result = []
    with open(file_name, 'r') as f:
        lines = f.read().splitlines()
    for i in range(len(lines)):
        if lines[i][0] == '>':
            index += 1
            seq_labels.append(lines[i][1:])
            sequences.append('')
        else:
            sequences[index] += lines[i]
    for j in range(len(sequences)):
        for i in range(len(sequences)):
            if i != j and sequences[i][:o] == sequences[j][-o:]:
                result.append(seq_labels[j] + ' ' + seq_labels[i])
    return result


# parameters: all genotypes
def expected_offspring(a, b, c, d, e, f):
    return 2 * float(a + b + c + 0.75 * d + 0.5 * e + 0 * f)


def shared_motif(file_name):
    index = -1
    sequences = []
    results = []
    motifs = []
    record = 0
    with open(file_name, 'r') as f:
        lines = f.read().splitlines()
    for i in range(len(lines)):
        if lines[i][0] == '>':
            index += 1
            sequences.append('')
        else:
            sequences[index] += lines[i]
    seq_one = sequences[0]
    others = sequences[1:]
    length = len(seq_one)
    for start in range(length - 1):
        for end in range(length - 1, start, -1):
            motifs.append(seq_one[start:end])
    motifs.sort(key=len, reverse=True)
    for motif in motifs:
        if len(motif) < record:
            return results
        try_another = False
        for sequence in others:
            if motif not in sequence:
                try_another = True
        if try_another:
            continue
        results.append(motif)
        record = len(motif)
    return results


# k: kth generation
# n: minimal number for calculation
def probability(k, n):
    from math import factorial
    total = 2 ** k
    result = 0
    for i in range(n, total + 1):
        ncr = factorial(total) / (factorial(i) * factorial(total - i))
        result += ncr * (0.25 ** i) * (0.75 ** (total - i))
    return result


def protein_motif(filename):
    import urllib.request
    address = 'http://www.uniprot.org/uniprot/'
    result = []
    index = -1
    with open(filename, 'r') as f:
        ids = f.read().splitlines()
    for protein_id in ids:
        result.append(protein_id)
        result.append('')
        index += 2
        entry = urllib.request.urlopen(address + protein_id + '.fasta').read()
        raw = entry.decode("utf-8").splitlines()
        protein = ''.join(raw[1: len(raw)])
        for i in range(len(protein) - 3):
            if protein[i] == 'N' and protein[i + 1] != 'P':
                if protein[i + 2] in ('S', 'T') and protein[i + 3] != 'P':
                    if result[index] == '':
                        result[index] = str(i + 1)
                    else:
                        result[index] += ' ' + str(i + 1)
    return result
