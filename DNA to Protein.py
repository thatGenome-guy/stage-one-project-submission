from collections import Counter

# Define nucleotide bases
NUCLEOTIDE_BASE = {
    "DNA": ["A", "T", "C", "G"],
    "RNA": ["A", "U", "C", "G"]
}

# DNA codon table
DNA_Codons = {
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
}

# RNA codon table
RNA_Codons = {
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UGU": "C", "UGC": "C",
    "GAU": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "UUU": "F", "UUC": "F",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAU": "H", "CAC": "H",
    "AUA": "I", "AUU": "I", "AUC": "I",
    "AAA": "K", "AAG": "K",
    "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUG": "M",
    "AAU": "N", "AAC": "N",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UGG": "W",
    "UAU": "Y", "UAC": "Y",
    "UAA": "_", "UAG": "_", "UGA": "_"
}


class bio_seq:
    def __init__(self, seq, seq_type="DNA"):
        self.seq = seq.upper()
        self.seq_type = seq_type

    def reverse_complement(self):
        """Returns the reverse complement of a DNA or RNA sequence"""
        if self.seq_type == "DNA":
            mapping = str.maketrans("ATCG", "TAGC")
        else:
            mapping = str.maketrans("AUCG", "UAGC")
        return self.seq.translate(mapping)[::-1]

    def translate_seq(self, init_pos=0):
        """Translate a DNA/RNA sequence into amino acids"""
        if self.seq_type == "DNA":
            return [DNA_Codons.get(self.seq[pos:pos + 3], '') for pos in range(init_pos, len(self.seq) - 2, 3)]
        elif self.seq_type == "RNA":
            return [RNA_Codons.get(self.seq[pos:pos + 3], '') for pos in range(init_pos, len(self.seq) - 2, 3)]

    def proteins_from_rf(self, aa_seq):
        """Extract proteins from an amino acid sequence"""
        current_prot = []
        proteins = []
        for aa in aa_seq:
            if aa == "_":
                if current_prot:
                    for p in current_prot:
                        proteins.append(p)
                    current_prot = []
            else:
                if aa == "M":
                    current_prot.append("")
                for i in range(len(current_prot)):
                    current_prot[i] += aa
        return proteins

    def all_proteins_from_orfs(self, ordered=False):
        """Find all proteins from all six reading frames"""
        rfs = []
        rfs.append(self.translate_seq(0))
        rfs.append(self.translate_seq(1))
        rfs.append(self.translate_seq(2))
        tmp_seq = bio_seq(self.reverse_complement(), self.seq_type)
        rfs.append(tmp_seq.translate_seq(0))
        rfs.append(tmp_seq.translate_seq(1))
        rfs.append(tmp_seq.translate_seq(2))

        res = []
        for rf in rfs:
            prots = self.proteins_from_rf(rf)
            for p in prots:
                res.append(p)

        if ordered:
            return sorted(res, key=len, reverse=True)
        return res


# ----------------------------
# Example: test the code here
# ----------------------------
if __name__ == "__main__":
    example_seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
    dna = bio_seq(example_seq, "DNA")

    print("Original DNA Sequence:\n", dna.seq)
    print("\nTranslated Amino Acids (Frame 0):")
    print("".join(dna.translate_seq()))

    print("\nPossible Proteins from All Reading Frames:")
    proteins = dna.all_proteins_from_orfs(ordered=True)
    for p in proteins:
        print(p)
