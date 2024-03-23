#!/usr/bin/env python
# coding: utf-8
# added file helpful variables
aa_mol_weights = {
    "A": 89.09,
    "C": 121.15,
    "D": 133.1,
    "E": 147.13,
    "F": 165.19,
    "G": 75.07,
    "H": 155.16,
    "I": 131.17,
    "K": 146.19,
    "L": 131.17,
    "M": 149.21,
    "N": 132.12,
    "P": 115.13,
    "Q": 146.15,
    "R": 174.2,
    "S": 105.09,
    "T": 119.12,
    "V": 117.15,
    "W": 204.23,
    "X": 0,
    "Y": 181.19,
}

standard_code = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "*",
    "UAG": "*",
    "UGA": "*",
    "UGU": "C",
    "UGC": "C",
    "UGG": "W",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


class seq:
    def __init__(self, name, organism, sequence, type):
        self.name = name
        self.organism = organism
        self.sequence = sequence
        self.type = type

    def info(self):
        print(f"Name: {self.name}")
        print(f"Type: {self.type}")
        print(f"Organism: {self.organism}")
        print(f"Sequence: {self.sequence}")

    def length(self):
        return len(self.sequence)

    def fasta_out(self):
        with open(f"{self.name}.fa", "w") as f:
            f.write(
                ">"
                + self.name
                + "_"
                + self.organism
                + "_"
                + self.type
                + "\n"
                + self.sequence
            )


# class protein
class protein(seq):
    def __init__(self, name, organism, sequence, type, size):
        super().__init__(name, organism, sequence, type)
        self.size = size

    def fasta_out(self):
        with open(f"{self.name}.fa", "w") as f:
            f.write(
                ">"
                + self.name
                + "_"
                + self.organism
                + "_"
                + self.type
                + "_"
                + str(self.size)
                + "kDa\n"
                + self.sequence
            )

    def mol_weight(self):
        mol_weight = 0
        for aa in self.sequence:
            mol_weight += aa_mol_weights.get(aa, 0)
        return mol_weight


# class nucleotide
class nucleotide(seq):
    def __init__(self, name, organism, sequence, type):
        super().__init__(name, organism, sequence, type)

    def gc_content(self):
        g_count = self.sequence.count("G")
        c_count = self.sequence.count("C")
        total_count = len(self.sequence)
        gc_percentage = (g_count + c_count) / total_count * 100
        print(f"GC content: {gc_percentage:.2f}%")


# class DNA
class DNA(nucleotide):
    def __init__(self, name, organism, sequence, type):
        super().__init__(name, organism, sequence, type)

    def transcribe(self):
        rna_sequence = self.sequence.replace("T", "U")
        print(f"Transcribed RNA sequence: {rna_sequence}")
        return rna_sequence

        # added 2 methods

    def reverse_complement(self):
        complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
        rev_comp = "".join(complement.get(base, base) for base in self.sequence[::-1])
        return rev_comp

    def six_frames(self):
        frames = []
        rev_comp_seq = self.reverse_complement()

        for i in range(3):
            frames.append(self.sequence[i:])
            frames.append(rev_comp_seq[i:])

        return frames


# class RNA
class RNA(nucleotide):
    def __init__(self, name, organism, sequence, type):
        super().__init__(name, organism, sequence, type)

    def start(self):
        start_codon_index = self.sequence.find("AUG")
        if start_codon_index != -1:
            print(f"Start codon 'AUG' found at index {start_codon_index}")
        else:
            print("Start codon 'AUG' not found in the sequence.")
            # add a method called translate

    def translate(self):
        start_codon_index = self.sequence.find("AUG")
        if start_codon_index == -1:
            print("Start codon 'AUG' not found in the sequence.")
            return None

        protein_sequence = ""
        codons = [
            self.sequence[start_codon_index + i : start_codon_index + i + 3]
            for i in range(0, len(self.sequence) - start_codon_index, 3)
        ]

        for codon in codons:
            if len(codon) == 3:
                amino_acid = standard_code.get(codon, "_")
                if amino_acid == "_":
                    break
                protein_sequence += amino_acid

        return protein_sequence


# test
def main():
    uidA = DNA(
        name="uidA",
        organism="bacteria",
        sequence="CGCATGTTACGTCCTGTAGAAACCCCAACCCGTGAAATCAAAAAA",
        type="DNA",
    )
    uidA.fasta_out()

    uidA_frames = uidA.six_frames()
    rev_comp_uidA = uidA.reverse_complement()
    print("Six Frames:")
    for frame in uidA_frames:
        print(frame)
    print("Reverse Complement:")
    print(rev_comp_uidA)

    uidA_RNA_sequence = uidA.transcribe()

    uidA_RNA = RNA(
        name="uidA_RNA", organism="bacteria", sequence=uidA_RNA_sequence, type="RNA"
    )
    uidA_RNA.fasta_out()

    uidA_protein_sequence = uidA_RNA.translate()

    uidA_protein = protein(
        name="uidA_protein",
        organism="bacteria",
        sequence=uidA_protein_sequence,
        type="protein",
        size=100,
    )
    uidA_protein.fasta_out()

    mol_weight = uidA_protein.mol_weight()
    print("Molecular Weight of uidA_protein:", mol_weight)


if __name__ == "__main__":
    main()
