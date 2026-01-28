import random
from itertools import combinations
from collections import Counter

INFA2="MALTFALLVALLVLSCKSSCSVGCDLPQTHSLGSRRTLMLLAQMRRISLFSCLKDRHDFGFPQEEFGNQFQKAETIPVLHEMIQQIFNLFSTKDSSAAWDETLLDKFYTELYQQLNDLEACVIQGVGVTETPLMKEDSILAVRKYFQRITLYLKEKKYSPCAWEVVRAEIMRSFSLSTNLQESLRSKE"
CSF3="MAGPATQSPMKLMALQLLLWHSALWTVQEATPLGPASSLPQSFLLKCLEQVRKIQGDGAALQEKLCATYKLCHPEELVLLGHSLGIPWAPLSSCPSQALQLAGCLSQLHSGLFLYQGLLQALEGISPELGPTLDTLQLDVADFATTIWQQMEELGMAPALQPTQGAMPAFASAFQRRAGGVLVASHLQSFLEVSYRVLRHLAQP"
EPO="MGVHECPAWLWLLLSLLSLPLGLPVLGAPPRLICDSRVLERYLLEAKEAENITTGCAEHCSLNENITVPDTKVNFYAWKRMEVGQQAVEVWQGLALLSEAVLRGQALLVNSSQPWEPLQLHVDKAVSGLRSLTTLLRALGAQKEAISPPDAASAAPLRTITADTFRKLFRVYSNFLRGKLKLYTGEACRTGDR"
PLAT="MDAMKRGLCCVLLLCGAVFVSPSQEIHARFRRGARSYQVICRDEKTQMIYQQHQSWLRPVLRSNRVEYCWCNSGRAQCHSVPVKSCSEPRCFNGGTCQQALYFSDFVCQCPEGFAGKCCEIDTRATCYEDQGISYRGTWSTAESGAECTNWNSSALAQKPYSGRRPDAIRLGLGNHNYCRNPDRDSKPWCYVFKAGKYSSEFCSTPACSEGNSDCYFGNGSAYRGTHSLTESGASCLPWNSMILIGKVYTAQNPSAQALGLGKHNYCRNPDGDAKPWCHVLKNRRLTWEYCDVPSCSTCGLRQYSQPQFRIKGGLFADIASHPWQAAIFAKHRRSPGERFLCGGILISSCWILSAAHCFQERFPPHHLTVILGRTYRVVPGEEEQKFEVEKYIVHKEFDDDTYDNDIALLQLKSDSSRCAQESSVVRTVCLPPADLQLPDWTECELSGYGKHEALSPFYSERLKEAHVRLYPSSRCTSQHLLNRTVTDNMLCAGDTRSGGPQANLHDACQGDSGGPLVCLNDGRMTLVGIISWGLGCGQKDVPGVYTKVTNYLDWIRDNMRP"
IGF1="MGKISSLPTQLFKCCFCDFLKVKMHTMSSSHLFYLALCLLTFTSSATAGPETLCGAELVDALQFVCGDRGFYFNKPTGYGSSSRRAPQTGIVDECCFRSCDLRRLEMYCAPLKPAKSARSVRAQRHTDMPKTQKEVHLKNASRGSAGNKNYRM"
CALB="MKLLSLTGVAGVLATCVAATPLVKRLPSGSDPAFSQPKSVLDAGLTCQGASPSSVSKPILLVPGTGTTGPQSFDSNWIPLSTQLGYTPCWISPPPFMLNDTQVNTEYMVNAITALYAGSGNNKLPVLTWSQGGLVAQWGLTFFPSIRSKVDRLMAFAPDYKGTVLAGPLDALAVSAPSVWQQTTGSALTTALRNAGGLTQIVPTTNLYSATDEIVQPQVSNSPLDSSYLFNGKNVQAQAVCGPLFVIDHAGSLTSQFSYVVGRSALRSTTGQARSADYGITDCNPLPANDLTPEQKVAAAALLAPAAAAIVAGPKQNCEPDLMPYARPFAVGKRTCSGIVTP"
mCitrine="MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYGLMCFARYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSYQSKLSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"

# Function to find stretches of homology between two sequences
def find_homologous_stretches(seq1, seq2):
    stretches = []
    start = None  # Starting index of a homology stretch

    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:  # Matching position
            if start is None:
                start = i  # Start a new stretch
        else:
            if start is not None:
                # End the current stretch and record it
                stretches.append((start, i - 1, i - start))
                start = None

    # Check if a stretch ended at the last position
    if start is not None:
        stretches.append((start, len(seq1) - 1, len(seq1) - start))

    # Sort stretches by length in descending order
    stretches.sort(key=lambda x: x[2], reverse=True)

    return stretches

# Function to find all homologous stretches across all pairs and count lengths
def find_all_homologous_stretches_and_count_lengths(sequences):
    # all_stretches = {}
    length_counts = Counter()

    for idx, (seq1, seq2) in enumerate(combinations(sequences, 2)):
        # pair_key = f"Pair {idx + 1} (Seq {sequences.index(seq1) + 1} vs Seq {sequences.index(seq2) + 1})"
        stretches = find_homologous_stretches(seq1, seq2)
        # all_stretches[pair_key] = stretches

        # Update the length count
        for _, _, length in stretches:
            length_counts[length] += 1

    return length_counts

# --- Standard codon table
CODON_TABLE = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'N': ['AAT', 'AAC'],
    'D': ['GAT', 'GAC'],
    'C': ['TGT', 'TGC'],
    'Q': ['CAA', 'CAG'],
    'E': ['GAA', 'GAG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    'H': ['CAT', 'CAC'],
    'I': ['ATT', 'ATC', 'ATA'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'K': ['AAA', 'AAG'],
    'M': ['ATG'],
    'F': ['TTT', 'TTC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'W': ['TGG'],
    'Y': ['TAT', 'TAC'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    '*': ['TAA', 'TAG', 'TGA']  # stop codons
}


def generate_random_dna_sequences(amino_acid_seq: str, n: int):
    """Generate n random DNA sequences encoding the same amino acid sequence."""
    amino_acid_seq = amino_acid_seq.strip().upper()
    sequences = []

    for _ in range(n):
        dna_seq = []
        for aa in amino_acid_seq:
            if aa not in CODON_TABLE:
                raise ValueError(f"Invalid amino acid: '{aa}'")
            dna_seq.append(random.choice(CODON_TABLE[aa]))
        sequences.append(''.join(dna_seq))
    return sequences


protein = [CALB, CSF3, EPO, IGF1, INFA2, mCitrine, PLAT]
names =   ["CALB", "CSF3", "EPO", "IGF1", "INFA2", "mCitrine", "PLAT"]
n =       [10, 10, 10, 10, 10, 10, 10, 6]

for i in range(len(protein)):
    for j in range(1, 4):
        dna_sequences = generate_random_dna_sequences(protein[i], n[i])

        # Find all homologous stretches and count lengths
        length_counts = find_all_homologous_stretches_and_count_lengths(dna_sequences)

        file_name = f"random_{names[i]}_{j}.txt"

        with open(file_name, "w") as f:
            for length, count in sorted(length_counts.items(), reverse=True):
                msg = "Length " + str(length) + ": " + str(count) + " occurrences"
                f.write(msg + "\n")
