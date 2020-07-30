import networkx as nx
import collections


def sequence_dotbracket_to_graph(seq_info=None, seq_struct=None):
    """Given a sequence and the dotbracket sequence make a graph.
    Parameters
    ----------
    seq_info string
        node labels eg a sequence string
    seq_struct  string
        dotbracket string
    Returns
    -------
        returns a nx.Graph
        secondary struct associated with seq_struct
    """
    graph = nx.Graph()

    lifo = collections.defaultdict(list)
    open_brace_string={")":"(",
                "]":"[",
                ">":"<"}

    for i, (c, b) in enumerate(zip(seq_info, seq_struct)):
        graph.add_node(i, label=c, position=i)
        if i > 0:
            graph.add_edge(i, i - 1, label='-', type='backbone', len=1)
        if b in ['(','[','<']:
            lifo[b].append(i)
        if b in [')',']','>']:
            j = lifo[open_brace_string[b]].pop()
            graph.add_edge(i, j, label='=', type='basepair', len=1)

    return graph


def seq_to_graph(header, sequence):
    """Fold a sequence in a path graph."""
    seq_struct = '.' * len(sequence)
    graph = sequence_dotbracket_to_graph(seq_info=sequence,
                                         seq_struct=seq_struct)
    graph.graph['info'] = 'sequence'
    graph.graph['sequence'] = sequence
    graph.graph['structure'] = seq_struct
    graph.graph['id'] = header
    return graph


def fold(seqs):
    """Fold a list of sequences into path graphs."""
    for header, seq in seqs:
        yield seq_to_graph(header, seq)



def null_modifier(header=None, seq=None, **options):
    """Null modifier."""
    yield header, seq


def load(input, **options):
    """Load sequences."""
    return fasta_to_fasta(input, **options)

def fasta_to_fasta(input, modifier=null_modifier, **options):
    """Take a FASTA file and yield a normalised FASTA file.

    Parameters
    ----------
    input : string
        A pointer to the data source.

    normalize : bool
        If True all characters are uppercased and Ts are replaced by Us
    """
    normalize = options.get('normalize', True)
    iterable = _fasta_to_fasta(input)
    for line in iterable:
        header = line
        seq = iterable.__next__()
        if normalize:
            seq = seq.upper()
            seq = seq.replace('T', 'U')
        seqs = modifier(header=header, seq=seq, **options)
        for seq in seqs:
            yield seq


def _fasta_to_fasta(input):
    seq = ""
    with open(input, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line:
                if line[0] == '>':
                    line = line[1:]
                    if seq:
                        yield seq
                        seq = ""
                    line_str = str(line)
                    yield line_str.strip()
                else:
                    line_str = line.split()
                    if line_str:
                        seq += str(line_str[0]).strip()
    if seq:
        yield seq

