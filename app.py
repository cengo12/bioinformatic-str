from Bio import Entrez, SeqIO
import streamlit as st
import pandas as pd

Entrez.email = st.secrets["email"]
Entrez.api_key = st.secrets["api_key"]


def get_seq(input_id):
    handle = Entrez.efetch(db="nucleotide", id=input_id, rettype="fasta")
    for record in SeqIO.parse(handle, "fasta"):
        return record.id, str(record.seq)


def get_id(search_term):
    search_term = search_term  # replace with the name of the gene you want to search for

    handle = Entrez.esearch(db="gene", term=search_term)
    record = Entrez.read(handle)

    gene_id = str(record["IdList"][0])  # the NCBI Gene ID for the gene will be the first item in the IdList
    return gene_id


def parse_fasta(fasta):
    first_line = fasta.find('\n')
    return fasta[first_line:]


def find_repeating_substrings(s):
    # Initialize an empty dictionary to store the substrings and their indexes
    substrings = {}

    # Iterate through the string and check for repeating substrings
    for i in range(len(s)):
        for j in range(i + 2, min(i + 7, len(s) + 1)):
            # Get the current substring
            substring = s[i:j]

            # If the substring is already in the dictionary, update its index
            if substring in substrings:
                substrings[substring].append(i)
            # If the substring is not in the dictionary, add it to the dictionary with its index
            else:
                substrings[substring] = [i]

    # Initialize an empty dictionary to store the repeating regions
    repeating_regions = {}

    # Iterate through the substrings and their indexes
    for substring, indexes in substrings.items():
        if len(indexes) >= 4:
            count = 0
            for i in range(len(indexes)-1):
                if indexes[i] + len(substring) == indexes[i+1] and count<(len(indexes)-2):
                    count = count+1
                else:
                    if count >= 3:
                        repeating_regions[(indexes[i-count] +1, indexes[i-count] + (count+1) *len(substring))] = substring
                        count = 0
                    else:
                        count = 0

    # Return the repeating regions
    return repeating_regions


# ----Streamlit page config----
st.set_page_config(page_title="Short Tandem Repeat (STR) Analysis",
                   page_icon=None, layout="wide",
                   initial_sidebar_state="auto",
                   menu_items=None)

# ----Input form----
with st.form(key="id_form"):
    input_type = st.radio("Input Type",
                          ('NCBI ID', 'Name', 'Manual Input'),
                          horizontal=True)
    input_value = st.text_area(label="Enter value below")

    st.form_submit_button("Search")


# ----Output----
if input_type == 'NCBI ID':
    try:
        gene_id, gene_sequence = get_seq(input_value)
        repeats = find_repeating_substrings(gene_sequence)
        st.write(pd.DataFrame(repeats.items(), columns=['Location', 'Sequence']))
    except:
        st.error("Invalid id or server error!")
elif input_type == 'Name':
    try:
        id = get_id(input_value)
        gene_id, gene_sequence = get_seq(id)
        repeats = find_repeating_substrings(gene_sequence)
        st.write(pd.DataFrame(repeats.items(), columns=['Location', 'Sequence']))
    except:
        st.error("Invalid name or server error!")
else:
    try:
        gene_sequence = parse_fasta(input_value)
        repeats = find_repeating_substrings(gene_sequence)
        st.write(pd.DataFrame(repeats.items(), columns=['Location', 'Sequence']))
    except:
        st.error("Invalid input type or server error!")
