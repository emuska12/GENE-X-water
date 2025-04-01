# Assignment 1 Local Sequence Alignment
# Subject: BIOINF_B FIIT STU 2025
# GENE-X scientist: Ema Tuptova 
# Date: 20. March 2025

############################################################################################################

# Assignment
# we found biological samples in the metheorite from the space
# we have extracted DNA sequences from the samples
    # nucleotide sequences
    # amino acid sequences
# after analysis we found out that the sequences are similar to ones on Earth
# we need to find the similarities between the uknown sequences and the Earth protein ones
# new form of life? 
# we need to align the sequences in ordare to compare them and find the similarities
# we need to program the local sequence alignment -- Smith-Waterman algorithm
# we need efficient alignment - dynamic programming

############################################################################################################

import streamlit as st # UI
import numpy as np # quick matrix initialization
from Bio import SeqIO # parse FASTA
from Bio.Align import substitution_matrices # load BLOSUM62, PAM250
from io import StringIO # convert to correct format
import datetime # get current time
import time # compare DP with reccursion memoization
import matplotlib.pyplot as plt

# Smith-Waterman algorithm IMPLEMENTATION
# local sequence alignment - 2 sequences (pairwise)
# similar to global allignment but any negative values are set to 0
# https://github.com/zhukovanadezhda/smith-waterman
# awards and penalties in scoring matrix ... Match (Conservation), Mismatch (Substitution), Gap Extend and Gap Open (Insertion, Deletion)
# local alignment = maximize similarity, align a section of one seq with another, excluding less similar regions

# ALIGNING NUCLEOTIDES (ACTG) SEQUENCES - in DNA
# # ALIGNING AMYNO ACIDS SEQUENCES (20 amino acids) - in proteins
# https://www.youtube.com/watch?v=lu9ScxSejSE
# DP = we take a problem, break it into smaller and simpler subproblems and solve them to derive the solution of the original problem
# DP steps: 1. INITIALIZATION 2. MATRIX FILLING 3. TRACEBACK
# Edit distance = smallest number of operations needed to transform seq 1 to seq 2

def build_alignment_matrices(sequences, M, N, sub_matrix, method):
    # 1. INITIALIZATION
    if method == "Dynamické programovanie":
        scoring_matrix = np.zeros((M + 1, N + 1), dtype=float) # it stores the best score for each position and determines the overall alignment
        # affine gap penalty model - they help distinguish between opening a new gap and extending an existing one
        gap_in_SEQ2 = np.zeros((M + 1, N + 1), dtype=float) 
        gap_in_SEQ1 = np.zeros((M + 1, N + 1), dtype=float) 
    else: # recursive with memoization
        scoring_matrix = np.full((M + 1, N + 1), -1, dtype=float) # -1 indicates uncomputed cells
        gap_in_SEQ2 = np.full((M + 1, N + 1), -1, dtype=float) 
        gap_in_SEQ1 = np.full((M + 1, N + 1), -1, dtype=float)

    traceback_matrix = np.full((M + 1, N + 1), None, dtype=object) # it stores "pointers" for parent cell (direction how we got the value in current cell), later paves the way for alignment reconstruction
    max_value, max_i, max_j = 0, 0, 0 # traceback will start from the cell with the highest score
    
    # 2. MATRIX FILLING - filling up matrix like we did in the lecture
    for i in range(1, M + 1):
        for j in range(1, N + 1):
            element1 = sequences["SEQ_1"][i - 1] # extract nucleotide from first sequence
            element2 = sequences["SEQ_2"][j - 1] # extract nucleotide from second sequence
            if method == "Dynamické programovanie": # DP computes values iteratively
                possible_values = []
                # edna provides scores for match/mismatch between nucleotide pairs
                # EMBOSS water uses EDNA FULL which is 15x15, this one is only 4x4
                # new walue = diagonally before + value in substitional matrix
                match_value = scoring_matrix[i - 1][j - 1] + sub_matrix[element1][element2]

                # calculate delete_value: score for adding a gap in SEQ_2 (deletion in SEQ_1)
                # compare opening a new gap (from scoring_matrix with GAP_OPEN penalty) vs. extending an existing gap in SEQ_2 (from gap_in_SEQ2 with GAP_EXTEND penalty)
                # max chooses the better option; if a gap already exists in SEQ_2, extending is preferred
                gap_in_SEQ2[i][j] = delete_value = max(scoring_matrix[i - 1][j] + modes["GAP_OPEN"], gap_in_SEQ2[i - 1][j] + modes["GAP_EXTEND"])

                # Calculate insert_value: score for adding a gap in SEQ_1 (insertion in SEQ_2)
                # Compare opening a new gap (from scoring_matrix with GAP_OPEN penalty)
                # vs. extending an existing gap in SEQ_1 (from gap_in_SEQ1 with GAP_EXTEND penalty)
                # max chooses the better option; if a gap already exists in SEQ_1, extending is preferred
                gap_in_SEQ1[i][j] = insert_value = max(scoring_matrix[i][j - 1] + modes["GAP_OPEN"], gap_in_SEQ1[i][j - 1] + modes["GAP_EXTEND"])

                # choose maximum value for a cell score, if it is negative --> 0
                possible_values = [match_value, delete_value, insert_value]
                new_value = max(0, max(possible_values))
                scoring_matrix[i][j] = new_value
                # update traceback matrix - store the direction of the parent cell
                if new_value == 0:
                    traceback_matrix[i][j] = None # when we find 0 cell in alignment reconstruction, traceback ends, so no need to save parent cell
                elif new_value == match_value:
                    traceback_matrix[i][j] = (i - 1, j - 1) # we got there diagonally - match/mismatch
                elif new_value == delete_value:
                    traceback_matrix[i][j] = (i - 1, j) # we got there from up cell
                elif new_value == insert_value:
                    traceback_matrix[i][j] = (i, j - 1) # we got there from left cell
            else: # Memoization computes values recursively
                new_value = compute_cell(sequences, i, j, scoring_matrix, gap_in_SEQ2, gap_in_SEQ1, sub_matrix, traceback_matrix, element1, element2)

            # update for max value so we don't need to go through whole matrix in another loop
            if new_value > max_value:
                max_value, max_i, max_j = new_value, i, j

    #print("    " + " ".join(sequences["SEQ_2"]))
    #for i, row in enumerate(scoring_matrix):
    #    row_label = sequences["SEQ_1"][i - 1] if i > 0 else " "
    #    print(row_label, *row)

    if show_matrix_plot:
        plot_matrix(scoring_matrix)

    # 3. TRACEBACK
    traceback_process(max_i, max_j, max_value, scoring_matrix, traceback_matrix, sequences, seq_type)

# helper function for recursive computation
def compute_cell(sequences, i, j, scoring_matrix, gap_in_SEQ2, gap_in_SEQ1, sub_matrix, traceback_matrix, element1, element2):
    # base case is if we reach the edge of the matrix
    if i == 0 or j == 0:
        return 0
    # uncomputed cell is marked with -1, so if already computed --> return the stored value
    if scoring_matrix[i][j] != -1: 
        return scoring_matrix[i][j]
    
    possible_values = []
    # recursiveôy calling for match, delete, and insert value
    # MATCH
    if i > 1:
        prev_element1 = sequences["SEQ_1"][i - 2]
    else:
        prev_element1 = None  # we can't go before the start of the sequence
    
    if j > 1:
        prev_element2 = sequences["SEQ_2"][j - 2]
    else:
        prev_element2 = None  # we can't go before the start of the sequence
    
    # MATCH = diagonally, the score from the diagonal cell and add the substitution score
    diagonal_score = compute_cell(sequences, i - 1, j - 1, scoring_matrix, gap_in_SEQ2, gap_in_SEQ1, sub_matrix, traceback_matrix, prev_element1, prev_element2)
    match_value = diagonal_score + sub_matrix[element1][element2]

    # DELETE = from above cell, add the gap open penalty and compare with the gap extend penalty
    up_score = compute_cell(sequences, i - 1, j, scoring_matrix, gap_in_SEQ2, gap_in_SEQ1, sub_matrix, traceback_matrix, prev_element1, element2)
    gap_in_SEQ2[i][j] = delete_value = max(up_score + modes["GAP_OPEN"], gap_in_SEQ2[i - 1][j] + modes["GAP_EXTEND"])


    # INSERT = from left cell, add the gap open penalty and compare with the gap extend penalty
    left_score = compute_cell(sequences, i, j - 1, scoring_matrix, gap_in_SEQ2, gap_in_SEQ1, sub_matrix, traceback_matrix, element1, prev_element2)
    gap_in_SEQ1[i][j] = insert_value = max(left_score + modes["GAP_OPEN"], gap_in_SEQ1[i][j - 1] + modes["GAP_EXTEND"])
    
    possible_values = [match_value, delete_value, insert_value]
    new_value = max(0, max(possible_values))
    scoring_matrix[i][j] = new_value

    if new_value == 0:
        traceback_matrix[i][j] = None
    elif new_value == match_value:
        traceback_matrix[i][j] = (i - 1, j - 1)
    elif new_value == delete_value:
        traceback_matrix[i][j] = (i - 1, j)
    elif new_value == insert_value:
        traceback_matrix[i][j] = (i, j - 1)
    # traceback matrix is mutable numpy object so no need to return it

    return scoring_matrix[i][j]

def traceback_process(max_i, max_j, score, scoring_matrix, traceback_matrix, sequences, seq_type):
    # initialization of the aligned sequences and traceback path
    ASEQ_1 = "" 
    ASEQ_2 = ""
    traceback_path = "" # ↘, ↓, →

    i, j = max_i, max_j # traceback starts at maximum value in scoring matrix
    alignment_length, matches, gaps = 0, 0, 0 # statistics

    # tracebacking until value 0 is found in cell of scoring matrix
    while scoring_matrix[i][j] != 0:
        
        if traceback_matrix[i][j] is None: # if we are in the end of the scoring matrix, traceback ends
            break
        
        prev_i, prev_j = traceback_matrix[i][j]  # get the coordinates of the parent cell in traceback matrix
        if prev_i is None or prev_j is None:
            break  

        if (prev_i, prev_j) == (i-1, j-1):  # parent matches diagonal cell coordinates - match or mismatch
            ASEQ_1 += sequences["SEQ_1"][i - 1] 
            ASEQ_2 += sequences["SEQ_2"][j - 1]
            if sequences["SEQ_1"][i - 1] == sequences["SEQ_2"][j - 1]:
                matches += 1 
            traceback_path += "↘"
        elif prev_i == i-1: # parent matches upper cell coordinated - gap in SEQ_2 (deletion for SEQ_1)
            ASEQ_1 += sequences["SEQ_1"][i - 1]
            ASEQ_2 += "-"
            gaps += 1
            traceback_path += "↓"
        elif prev_j == j-1: # parent matches left cell coordinated - gap in SEQ_1 (insertion for SEQ_2)
            ASEQ_1 += "-"
            ASEQ_2 += sequences["SEQ_2"][j - 1]
            gaps += 1
            traceback_path += "→"

        alignment_length += 1
        i, j = prev_i, prev_j  # move to that previous cell in next loop

    # reverse the sequences and traceback path because it was built backwards
    ASEQ_1 = ASEQ_1[::-1]
    ASEQ_2 = ASEQ_2[::-1]
    traceback_path = traceback_path[::-1]
    
    # end is basically starting index because it was built backwards
    end_1 = i + 1 
    end_2 = j + 1

    statistics(ASEQ_1, ASEQ_2, traceback_path, seq_type, seq1ID, seq2ID, modes, alignment_length, matches, gaps, score, max_i, max_j, end_1, end_2)


def statistics(ASEQ_1, ASEQ_2, traceback_path, seq_type, seq1ID, seq2ID, modes, alignment_length, matches, gaps, score, max_i, max_j, start_1, start_2):
    similar_matches = 0  # relevant only for proteins
    annotation = ""  

    # visual representation of the alignent based on similarity of the aligned sequences
    for i in range(len(ASEQ_1)):
        x, y = ASEQ_1[i], ASEQ_2[i]

        if x == y:
            annotation += "|"  # MATCH wohou
        elif x == "-" or y == "-":
            annotation += " "  # GAP in one of the sequences
        elif seq_type == "PROTEIN" and sub_matrix[x][y] > 0:
            annotation += ":"  # SIMILAR amino acids (positive value in substitution matrix)
            similar_matches += 1
        else:
            annotation += "."  # MISMATCH

    end_time = time.time()
    execution_time = end_time - start_time

    # print out in streamlit
    st.markdown("########################################")
    st.markdown("**Program:** GENE-X")
    st.markdown(f"**Rundate:** {datetime.datetime.now()}")
    st.markdown(f"**Execution time:** {execution_time:.4f} s")
    st.markdown("**Align_format:** pair")
    st.markdown("**Report_file:** stdout")
    st.markdown("########################################")
    st.markdown("#======================================")
    st.markdown("**Aligned_sequences:** 2")
    st.markdown(f"**Sequence 1:** {seq1ID.split(' ')[0]}")  # identificator of SEQ_1 (firt part so it is not long)
    st.markdown(f"**Sequence 2:** {seq2ID.split(' ')[0]}")  # identificator of SEQ_2
    if seq_type == "DNA":
        st.markdown(f"**Matrix:** EDNA")
    else:
        st.markdown(f"**Matrix:** {matrix_choice}") 
    st.markdown(f"**Gap penalty:** {modes['GAP_OPEN']}")
    st.markdown(f"**Extension penalty:** {modes['GAP_EXTEND']}")
    st.markdown(f"**Alignment Length:** {alignment_length}")
    st.markdown(f"**Identity:** {matches}/{alignment_length} ({(matches / alignment_length)*100:.1f}%)") # exact matches
    st.markdown(f"**Similarity:** {matches + similar_matches}/{alignment_length} ({((matches + similar_matches) / alignment_length)*100:.1f}%)")  # exact + similar matches
    st.markdown(f"**Gaps:** {gaps}/{alignment_length} ({(gaps / alignment_length)*100:.1f}%)") 
    st.markdown(f"**Score (max value):** {score} at ({max_i}, {max_j}) in matrix")

    # Formátovanie a zobrazenie zarovnania
    print_out_alignment(ASEQ_1, ASEQ_2, traceback_path, annotation, start_1, start_2, alignment_length)

def print_out_alignment(aseq1, aseq2, traceback_path, annotation, start_1, start_2, alignment_length):
    output = "#======================================\n\n"
    max_id_length = max(len(seq1ID), len(seq2ID)) # indent based on longest name of sequence so it is not weird to read
    # indexes of sequences also have some width and single digit is not as thick as three digits, so we need to align them according to max digits in number that will occur
    max_num_width = max(len(str(start_1+alignment_length)), len(str(start_2+alignment_length)))
    
    # to watch current alignment position, it is not just +50, +50... because gaps dont count
    current_start_1 = start_1
    current_start_2 = start_2
    
    for i in range(0, alignment_length, max_chars_in_line):
        # get chunks for this block
        line_seq1 = aseq1[i:(i+max_chars_in_line)]
        line_seq1_end = current_start_1 + len([element for element in line_seq1 if element != '-']) - 1
        line_seq2 = aseq2[i:(i+max_chars_in_line)]
        line_seq2_end = current_start_2 + len([element for element in line_seq2 if element != '-']) - 1
        match_line = annotation[i:(i+max_chars_in_line)]
        tr_path_line = traceback_path[i:(i+max_chars_in_line)]
        

        # I used ChatGPT for indenting lines because never heard of this ":<2" syntax and wanted to indent annotation accroding to sequence characters
        output += f"{' ':<{max_id_length + max_num_width + 2}}{tr_path_line}\n"
        output += f"{seq1ID:<{max_id_length}} {current_start_1:>{max_num_width}} {line_seq1} {line_seq1_end:>{max_num_width}}\n"
        output += f"{' ':{max_id_length + max_num_width + 2}}{match_line}\n"
        output += f"{seq2ID:<{max_id_length}} {current_start_2:>{max_num_width}} {line_seq2} {line_seq2_end:>{max_num_width}}\n"
        # to update start position for next block in next loop just +1
        current_start_1 = line_seq1_end + 1
        current_start_2 = line_seq2_end + 1
    output += "\n#======================================\n"
    st.code(output, language="plaintext")

def plot_matrix(scoring_matrix):
    plt.imshow(scoring_matrix, cmap="Grays")
    plt.colorbar(label="Skóre")
    plt.title("Skórovacia matica")
    st.pyplot(plt)
    plt.close()

############################################################################################################

# load BLOSUM62 matrix - because of physical and biological properties of amino acids more complex scoring
b62 = substitution_matrices.load("BLOSUM62")
pam250 = substitution_matrices.load("PAM250")

# EDNA-like full scoring matrix for A, T, G, C https://rosalind.info/glossary/dnafull/
edna= {
    'A': {'A': 5, 'T': -4, 'G': -4, 'C': -4},
    'T': {'A': -4, 'T': 5, 'G': -4, 'C': -4},
    'G': {'A': -4, 'T': -4, 'G': 5, 'C': -4},
    'C': {'A': -4, 'T': -4, 'G': -4, 'C': 5}
}

# User interface 
st.title("🔬 Lokálne zarovnávanie sekvencií")
st.markdown("### Ema Tuptová")
# does the user want to align DNA or PROTEIN sequences, radio buttons
seq_type = st.radio("Typ sekvencií", ["DNA", "PROTEIN"])
# does the user want to align with DP or RECCURSIVE MEMOISATION approach, radio buttons
method = st.radio("Metóda zarovnania", ["Dynamické programovanie", "Rekurzia s memoizáciou"])
# alignemt settings that are pre-set like EMBOSS WATER but user can change them anytime
gap_open_penalty = st.number_input("Penalizácia za otvorenie medzery (Gap Open)", value=-10)
gap_extend_penalty = st.number_input("Penalizácia za predĺženie medzery (Gap Extend)", value=-0.5)
max_chars_in_line = st.number_input("Maximálny počet znakov na riadok vo výstupe", value=50)
modes = {
    #"MATCH": match_score, # not necessary anymore, because we can calculate it from the matrix (even for DNA)
    #"MISS": mismatch_score,
    "GAP_OPEN": gap_open_penalty,
    "GAP_EXTEND": gap_extend_penalty,
}
# if aligning proteins, user can also specify sub_matrix
if seq_type == "PROTEIN":
    matrix_choice = st.radio("Výber matice", ["BLOSUM62", "PAM250"])
    if matrix_choice == "BLOSUM62":
        sub_matrix = b62
    else:
        sub_matrix = pam250
else:
    sub_matrix = edna

show_matrix_plot = st.checkbox("Zobraziť skórovaciu maticu", value=True)

# ok so lets load the sequences from FASTA files
# https://www.youtube.com/watch?v=yj3qld0J5EY --i do not have to write my own parser?
# FASTA file - standard format for storing sequence data
# first identifier line: > sequence_id
# Biopython SeqIO package has parse(file_path, format) which returns seq files as SeqRecord objects containing id (string) and seq
# there is a need to install Bipython (pip install biopython)
# user can download with pip install -r requirements.txt

# ALWAYS listening for user action

# input for seq 1 - textarea or file loader
# if there would be both sequences loaded in one same FASTA -- different approach (no next, but for loop)

seq1 = st.text_area("Sekvencia 1", ">sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens GN=HBA1 PE=1 SV=2\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR" if seq_type == "PROTEIN" else ">test1\nATGAGTCTCTCTGATAAGGACAAGGCTGCTGTGAAAGCCCTATGG")
seq1_file = st.file_uploader("Vybrať súbor pre sekvenciu 1", type=["fasta", "txt"])

if seq1_file: # file
    try:
        # convert to unicode and create text stream
        file_content = seq1_file.getvalue().decode("utf-8")  
        file_handle = StringIO(file_content)  # Vytvor textový stream
        record = next(SeqIO.parse(file_handle, "fasta"))  # BioPython function
        # extract id (first line) and sequence
        seq1 = str(record.seq) 
        seq1ID = str(record.id)
    except Exception as e:
        st.error(f"Chyba pri načítaní súboru: {e}")
        seq1 = ""
else: # textarea
    try:
        file_handle = StringIO(seq1)
        record = next(SeqIO.parse(file_handle, "fasta"))
        seq1 = str(record.seq)
        seq1ID = str(record.id)
    except Exception as e:
        st.error(f"Chyba pri spracovaní sekvencie 1 {e}")


# same for SEQ2
seq2 = st.text_area("Sekvencia 2", ">sp|P01942|HBA_MOUSE Hemoglobin subunit alpha OS=Mus musculus GN=Hba PE=1 SV=2\nMVLSGEDKSNIKAAWGKIGGHGAEYGAEALERMFASFPTTKTYFPHFDVSHGSAQVKGHGKKVADALASAAGHLDDLPGALSALSDLHAHKLRVDPVNFKLLSHCLLVTLASHHPADFTPAVHASLDKFLASVSTVLTSKYR " if seq_type == "PROTEIN" else ">test2\nCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAG")
seq2_file = st.file_uploader("Vybrať súbor pre sekvenciu 2", type=["fasta", "txt"])

if seq2_file:
    try:
        file_content = seq2_file.getvalue().decode("utf-8")
        file_handle = StringIO(file_content)
        record = next(SeqIO.parse(file_handle, "fasta"))
        seq2 = str(record.seq)
        seq2ID = str(record.id)
    except Exception as e:
        st.error(f"Chyba pri načítaní súboru: {e}")
        seq2 = ""
else:
    try:
        file_handle = StringIO(seq2)
        record = next(SeqIO.parse(file_handle, "fasta"))
        seq2 = str(record.seq)
        seq2ID = str(record.id)
    except Exception as e:
        st.error(f"Chyba pri spracovaní sekvencie 2 {e}")


# ALWAYS listen if user clicked on button
if st.button("Zarovnať"):
    if seq1 and seq2: # both sequences need to be provided before continuing
        # now we have the sequences stored in the dictionary
        # we can start implementing the Smith-Waterman algorithm
        # https://www.youtube.com/watch?v=nqqBgWaNoig

        start_time = time.time()

        sequences = {"SEQ_1": seq1, "SEQ_2": seq2} # sequences are saved here
        M = len(sequences["SEQ_1"])
        N = len(sequences["SEQ_2"])

        st.title("📊 Výsledky zarovnania")
        build_alignment_matrices(sequences, M, N, sub_matrix, method)

    else:
        st.error("Prosím zadajte obe sekvencie.")

############################################################################################################


