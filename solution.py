def count_nucleotides(file_name):
    with open(file_name, 'r') as f:
        sequence = ''.join(f.read().splitlines())
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
        sequence = ''.join(f.read().splitlines())
    sequence = sequence.replace('T', 'U')
    return sequence


def complement(file_name):
    comp_seq = ''
    with open(file_name, 'r') as f:
        sequence = ''.join(f.read().splitlines())
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
        seq1 = f.readline().rstrip()
        seq2 = f.readline().rstrip()
    count = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
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
    start_index = -1
    peptide = []
    with open(file_name, 'r') as f:
        sequence = ''.join(f.read().splitlines())
    for i in range(len(sequence)):
        if sequence[i:(i + 3)] == 'AUG':
            start_index = i
            break
    if start_index == -1:
        return ''
    for i in range(start_index, len(sequence), 3):
        codons.append(sequence[i:(i + 3)])
    if len(codons[-1]) != 3:
        codons.pop()
    for codon in codons:
        aa = codon_table(codon)
        if aa == '':
            break
        peptide.append(aa)
    return ''.join(peptide)


def codon_table(codon):
    if codon[0] == 'U':
        if codon in ['UUU', 'UUC']:
            return 'F'
        elif codon in ['UUA', 'UUG']:
            return'L'
        elif codon[1] == 'C':
            return'S'
        elif codon in ['UAU', 'UAC']:
            return'Y'
        elif codon in ['UAA', 'UAG', 'UGA']:
            return ''
        elif codon in ['UGU', 'UGC']:
            return'C'
        elif codon == 'UGG':
            return 'W'
    elif codon[0] == 'C':
        if codon[1] == 'U':
            return 'L'
        elif codon[1] == 'C':
            return 'P'
        elif codon in ['CAU', 'CAC']:
            return 'H'
        elif codon in ['CAA', 'CAG']:
            return 'Q'
        elif codon[1] == 'G':
            return 'R'
    elif codon[0] == 'A':
        if codon in ['AUU', 'AUC', 'AUA']:
            return 'I'
        elif codon == 'AUG':
            return 'M'
        elif codon[1] == 'C':
            return 'T'
        elif codon in ['AAU', 'AAC']:
            return 'N'
        elif codon in ['AAA', 'AAG']:
            return 'K'
        elif codon in ['AGU', 'AGC']:
            return 'S'
        elif codon in ['AGA', 'AGG']:
            return 'R'
    elif codon[0] == 'G':
        if codon[1] == 'U':
            return 'V'
        elif codon[1] == 'C':
            return 'A'
        elif codon in ['GAU', 'GAC']:
            return 'D'
        elif codon in ['GAA', 'GAG']:
            return 'E'
        elif codon[1] == 'G':
            return 'G'


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


def protein_motif(file_name):
    import urllib.request
    address = 'http://www.uniprot.org/uniprot/'
    result = []
    index = -1
    with open(file_name, 'r') as f:
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


def infer_rna(file_name):
    aa_table = {'A': 4, 'R': 6, 'N': 2, 'D': 2, 'C': 2, 'Q': 2, 'E': 2, 'G': 4,
                'H': 2, 'I': 3, 'L': 6, 'K': 2, 'M': 1, 'F': 2, 'P': 4, 'S': 6,
                'T': 4, 'W': 1, 'Y': 2, 'V': 4}
    result = 1
    with open(file_name, 'r') as f:
        protein = ''.join(f.read().splitlines())
    for aa in protein:
        result *= aa_table[aa]
    return (result * 3) % 1000000


def translate_orf(file_name):
    result = []
    with open(file_name, 'r') as f:
        lines = f.read().splitlines()
    sequence = ''.join(lines[1:])
    comp_seq = []
    rev_seq = sequence[::-1]
    for n in rev_seq:
        if n == 'A':
            comp_seq.append('T')
        elif n == 'C':
            comp_seq.append('G')
        elif n == 'G':
            comp_seq.append('C')
        elif n == 'T':
            comp_seq.append('A')
    comp_seq = ''.join(comp_seq)
    orf = [sequence.replace('T', 'U'),
           sequence[1:].replace('T', 'U'),
           sequence[2:].replace('T', 'U'),
           comp_seq.replace('T', 'U'),
           comp_seq[1:].replace('T', 'U'),
           comp_seq[2:].replace('T', 'U')]
    for frame in orf:
        start_indices = []
        stop_indices = []
        for i in range(0, len(frame), 3):
            if frame[i:(i + 3)] == 'AUG':
                start_indices.append(i)
            if frame[i:(i + 3)] in ['UAA', 'UAG', 'UGA']:
                stop_indices.append(i)
        if not start_indices or not stop_indices:
            continue
        for start in start_indices:
            if start > stop_indices[-1]:
                continue
            codons = []
            peptide = []
            for i in range(start, len(frame), 3):
                codons.append(frame[i:(i + 3)])
            if len(codons[-1]) != 3:
                codons.pop()
            for codon in codons:
                aa = codon_table(codon)
                if aa == '':
                    break
                peptide.append(aa)
            result.append(''.join(peptide))
    return list(set(result))


def permutation(n):
    from itertools import permutations
    perm = permutations(range(1, n + 1))
    p = list(perm)
    for per in p:
        print(' '.join(str(x) for x in per))
    return len(p)


def weight(file_name):
    mass_tab = {'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259,
                'F': 147.06841, 'G': 57.02146, 'H': 137.05891, 'I': 113.08406,
                'K': 128.09496, 'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
                'P': 97.05276, 'Q': 128.05858, 'R': 156.10111, 'S': 87.03203,
                'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333}
    w = 0
    with open(file_name, 'r') as f:
        protein = f.read().splitlines()
        for a in protein[0]:
            w += mass_tab[a]
    return w
