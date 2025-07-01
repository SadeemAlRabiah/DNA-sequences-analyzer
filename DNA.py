from Bio.Seq import Seq
import matplotlib.pyplot as plt

def is_vaild_dna(seq):
    return all(base in "ATCG" for base in seq.upper())

def gc_content(seq):
    seq = seq.upper()
    g = seq.count("G")
    c = seq.count("C")
    return round(((g+c)/len(seq))*100, 2)

def base_percentages(seq):
    seq = seq.upper()
    total = len(seq)
    return{base: round(seq.count(base)/total*100, 2) for base in "ATCG"}

def translate_dna(seq):
    return str(Seq(seq).translate(to_stop=True))

def plot_base_composition():
    percentages = base_percentages(seq)
    base = list(percentages.keys())
    value = list(percentages.values())
    colors = ["#FF9999","#66B2FF", "#99FF99", "#FFCC99"]
    plt.figure(figsize=(6, 4))
    plt.bar(base,value,color=colors)
    plt.title("Base Composition in DNA Sequence")
    plt.grid(axis="y", linestyle="--", alpha=0.5)
    plt.tight_layout()
    

seq = "ATGGCGTACGTAGCTAGCGTAGCTAGGCTAG"
if not is_vaild_dna(seq):
    print("invalid DNA sequence. only A, T, C, G, are allowed")
else:
    plot_base_composition()
    print("DNA Sequence Analysis")
    print(f"length:{len(seq)} base")
    print(f"GC content:{gc_content(seq)}%")
    print(base_percentages(seq))
    plt.show()

