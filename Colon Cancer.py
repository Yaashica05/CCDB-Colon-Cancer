import streamlit as st
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
import re  # Import regex module
from collections import Counter
import urllib.parse # Import for URL encoding DOI
import io
import numpy as np # Import NumPy for Dot Plot matrix
from Bio import Restriction # Import Biopython Restriction module
from Bio.Data import CodonTable # Import for ORF finder and Codon Usage

# --- Page Configuration (MUST be the first Streamlit command) ---
st.set_page_config(
    page_title="Colon Cancer Toolbox",
    layout="wide",
    initial_sidebar_state="expanded",
)

# --- Custom CSS Styling ---
# (CSS remains the same as provided)
st.markdown(
    """
    <style>
    /* --- Global Font Settings --- */
    body, .stApp, h1, h2, h3, h4, h5, h6, p, div, span, li, label,
    .stButton>button,
    .stTextInput>div>div>input,
    .stTextArea>div>textarea,
    .stSelectbox>div>div>div,
    table th, table td, table a, /* Include table elements */
    .footer /* Include footer */
    {
        font-family: 'Times New Roman', Times, serif !important; /* Force Times New Roman */
        font-weight: bold !important; /* Force Bold */
    }

    /* --- Specific Element Styling (Adjust colors, backgrounds, etc. as needed) --- */

    /* Body background - a very light, neutral grey */
    body {
        background-color: #f8f8f8;
    }

    /* Main app container background - Linen */
    .stApp {
        background-color: #faf0e6; /* Linen */
        padding: 20px;
        border-radius: 15px;
        box-shadow: 0 4px 8px rgba(0,0,0,0.1);
        padding-bottom: 60px; /* Avoid footer overlap */
        color: #333; /* Default text color */
    }

    /* Headings */
    h1, h2, h3, h4, h5, h6 {
        color: #333; /* Darker grey for headings */
    }

    /* Paragraph text - adjust line height if needed */
    p, div, span, li {
        color: #444; /* Slightly lighter than headings */
        line-height: 1.6;
    }

     /* Labels for input elements */
    .stTextInput>label, .stTextArea>label, .stSelectbox>label, .stNumberInput>label, .stFileUploader>label, .stSlider>label, .stRadio>label, .stCheckbox>label {
        color: #333;
        /* font-family and font-weight already set globally */
    }

    /* Input fields */
    .stTextInput>div>div>input, .stTextArea>div>textarea, .stSelectbox>div>div>div {
        border: 1px solid #ccc;
        border-radius: 4px;
        /* font-family and font-weight already set globally */
    }

    /* Buttons */
    .stButton>button {
        border-radius: 5px;
        background-color: #4682b4; /* SteelBlue */
        color: white;
        border: none;
        padding: 8px 15px;
        transition: background-color 0.3s ease;
        /* font-family and font-weight already set globally */
    }
    .stButton>button:hover {
        background-color: #5a9bd3; /* Lighter SteelBlue on hover */
        color: white;
    }

    /* Table styling */
    table {
        border-collapse: collapse;
        width: 100%;
        margin-bottom: 1rem;
        background-color: #ffffff;
        border: 1px solid #ddd;
    }
    table th {
        text-align: left;
        padding: 8px 12px;
        background-color: #e9ecef; /* Light grey header */
        border-bottom: 2px solid #dee2e6;
        /* font-family and font-weight already set globally */
    }
    table td {
        vertical-align: top;
        text-align: left;
        padding: 8px 12px;
        border-top: 1px solid #eee;
        /* font-family and font-weight already set globally */
    }
    table a { /* Style for links inside tables */
        color: #0056b3; /* Slightly darker blue for links */
        text-decoration: underline;
        /* font-family and font-weight already set globally */
    }
    table a:hover {
        color: #003d80;
        text-decoration: none;
    }
    table tbody tr:nth-of-type(odd) {
        background-color: rgba(0, 0, 0, 0.03);
    }

    /* Footer styling */
    .footer {
        position: fixed;
        left: 0;
        bottom: 0;
        width: 100%;
        background-color: #d2b48c; /* Tan */
        color: #333; /* Dark text for contrast */
        text-align: center;
        padding: 8px;
        font-size: 13px; /* Footer size can be adjusted slightly if needed */
        border-top: 1px solid #a08464; /* Darker tan border */
        z-index: 1000;
        /* font-family and font-weight already set globally */
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# --- Main App Title ---
st.title("🧬 Colon Cancer Database and Bioinformatics Toolbox 🔬")

# --- Sidebar Navigation ---
# (Updated with new tools)
menu = st.sidebar.selectbox(
    "Choose a Feature:",
    [
        "Home",
        "Colon Cancer Database Search",
        "DNA Sequence Analysis",
        "Protein Sequence Analysis",
        "Bioinformatics Tool (Transcription/Translation)",
        "ORF Finder",
        "Primer Design",
        "Restriction Enzyme Analysis",
        "Motif Finder Tool",
        "K-mer Counter",
        "Genome Coverage Plotter",
        "Dot Plot Comparison",
        "Variant Annotation Tool",
        "Codon Usage Analyzer",
    ],
    help="Select a tool or section from the dropdown list."
)

# --- Helper function to create links for the Database Search ---
# (format_link function remains the same)
def format_link(value, column_name):
    """Creates an HTML link based on column name and value.
       Handles specific ID columns and generic URLs.
    """
    if pd.isna(value):
        return "" # Return empty string for missing values

    s_value = str(value).strip() # Convert to string and strip whitespace

    if not s_value: # Handle empty strings after stripping
        return ""

    # --- Specific Column Formatting (Uses standardized Title Case names) ---
    col_name_std = column_name.strip().title() # Ensure check uses standardized name

    if col_name_std == 'Pubmed Id':
        if s_value.replace('.', '', 1).isdigit(): # Allow potential non-integer IDs
            url = f"https://pubmed.ncbi.nlm.nih.gov/{s_value}"
            return f'<a href="{url}" target="_blank">{s_value}</a>'
        elif s_value.startswith("http"):
             return f'<a href="{s_value}" target="_blank">{s_value}</a>'
        else:
            return s_value # Return as text if not a number or URL

    elif col_name_std == 'Doi Id':
        if '/' in s_value and not s_value.startswith("http"):
            try:
                encoded_doi = urllib.parse.quote(s_value, safe='/:')
                url = f"https://doi.org/{encoded_doi}"
                return f'<a href="{url}" target="_blank">{s_value}</a>'
            except Exception:
                 return s_value
        elif s_value.startswith("http"):
             return f'<a href="{s_value}" target="_blank">{s_value}</a>'
        else:
            return s_value

    elif col_name_std == 'Uniprot':
         if re.match(r'^[A-Z0-9_.-]+$', s_value, re.IGNORECASE) and not s_value.startswith("http"):
             url = f"https://www.uniprot.org/uniprotkb/{s_value}/entry"
             return f'<a href="{url}" target="_blank">{s_value}</a>'
         elif s_value.startswith("http"):
             return f'<a href="{s_value}" target="_blank">{s_value}</a>'
         else:
             return s_value

    elif col_name_std in ['Blast', 'Conserved Domain']:
         if s_value.startswith("http"):
              return f'<a href="{s_value}" target="_blank">{s_value}</a>'
         else:
              return s_value

    # --- Default: Handle generic http links ---
    else:
        if s_value.startswith("http"):
            return f'<a href="{s_value}" target="_blank">{s_value}</a>'
        else:
            return s_value

# --- Page Content based on Menu Selection ---

if menu == "Home":
    st.header("Welcome to the Colon Cancer Data Resource & Toolbox")
    st.markdown(
        """
        This application provides access to a curated database of genes implicated in colon cancer
        and a suite of bioinformatics tools for sequence analysis. Our goal is to offer a user-friendly
        platform for researchers, students, and clinicians interested in colon cancer genomics and molecular biology.

        ---

        ### **Understanding Colon Cancer (Colorectal Cancer - CRC)**

        Colon cancer, often grouped with rectal cancer as colorectal cancer (CRC), develops in the large intestine (colon) or the rectum (the final section of the large intestine). It is one of the most common cancers worldwide.

        *   **Development:** Most CRCs begin as noncancerous growths called polyps on the inner lining of the colon or rectum. Over time, some polyps can transform into cancer. Adenomatous polyps are the most common type with malignant potential.
        *   **Spread:** If not detected early, cancer cells can grow into the wall of the colon or rectum, invade nearby lymph nodes or blood vessels, and metastasize (spread) to distant organs like the liver or lungs.

        ### **Key Information:**

        *   **Risk Factors:**
            *   **Age:** Risk increases significantly after age 50.
            *   **Personal/Family History:** Previous CRC, polyps, or inflammatory bowel disease (IBD like Crohn's or Ulcerative Colitis). Family history of CRC or polyps.
            *   **Inherited Syndromes:** Lynch syndrome (HNPCC) and Familial Adenomatous Polyposis (FAP) greatly increase risk.
            *   **Lifestyle:** Diets high in red/processed meats, low fiber, obesity, physical inactivity, smoking, heavy alcohol use.
            *   **Type 2 Diabetes:** Associated with increased risk.
        *   **Symptoms:** *(May be absent in early stages)*
            *   Changes in bowel habits (diarrhea, constipation, narrowing of stool) lasting more than a few days.
            *   Rectal bleeding or blood in the stool.
            *   Persistent abdominal discomfort (cramps, gas, pain).
            *   A feeling that the bowel doesn't empty completely.
            *   Unexplained weight loss, weakness, or fatigue.
        *   **Screening & Early Detection:** Crucial for prevention and better outcomes.
            *   **Colonoscopy:** Visualizes the entire colon, allows polyp removal. (Gold standard)
            *   **Stool-Based Tests:** FIT (Fecal Immunochemical Test), gFOBT (guaiac-based Fecal Occult Blood Test), Stool DNA tests.
            *   **Other Visual Exams:** Sigmoidoscopy, CT Colonography (Virtual Colonoscopy).
            *   *Screening recommendations vary; consult a healthcare provider.*
        *   **Genetic Factors:** Mutations in key genes like *APC*, *KRAS*, *BRAF*, *TP53*, and mismatch repair genes (*MLH1*, *MSH2*, *MSH6*, *PMS2*) are frequently involved in CRC development.
        *   **Treatment Options:** Depend on stage and individual factors. May include surgery, chemotherapy, radiation therapy, targeted drug therapy (e.g., against EGFR or VEGF), and immunotherapy (especially for MSI-High tumors).

        ---

        ### **Available Tools in this Application:**

        Use the sidebar menu to explore these features:

        *   **Colon Cancer Database Search:** Find information on specific genes linked to colon cancer from our curated dataset.
        *   **DNA Sequence Analysis:** Calculate GC content, length, and base composition of a DNA sequence.
        *   **Protein Sequence Analysis:** Determine molecular weight, pI, GRAVY score, amino acid composition, and predicted secondary structure of a protein sequence.
        *   **Bioinformatics Tool (Transcription/Translation):** Transcribe DNA to RNA and translate RNA to protein.
        *   **ORF Finder:** Identify potential Open Reading Frames (protein-coding regions) in a DNA sequence.
        *   **Primer Design:** Generate basic forward and reverse primers from a DNA template.
        *   **Restriction Enzyme Analysis:** Find cut sites for specified restriction enzymes in a DNA sequence.
        *   **Motif Finder Tool:** Locate specific sequence patterns (motifs) within DNA or protein sequences.
        *   **K-mer Counter:** Count occurrences of short subsequences (k-mers) in DNA or protein.
        *   **Genome Coverage Plotter:** Visualize sequencing coverage depth across genomic positions from an uploaded file.
        *   **Dot Plot Comparison:** Generate a dot plot to visualize similarity between two sequences.
        *   **Variant Annotation Tool:** Apply a simple single-base substitution to a reference DNA sequence.
        *   **Codon Usage Analyzer:** Analyze the frequency of codon usage in a coding DNA sequence (CDS).

        ---
        *Disclaimer: This tool is intended for educational and informational purposes only. It is not a substitute for professional medical advice, diagnosis, or treatment. Always consult with a qualified healthcare provider regarding any medical conditions or treatment options.*
        """, unsafe_allow_html=True # Allow HTML like lists
    )

    st.image(
        "https://img.freepik.com/free-vector/colorectal-cancer-crc-infographic-education_1308-47971.jpg",
        caption="Colorectal Cancer Infographic (Source: Freepik - Illustrative purposes)",
        use_container_width=True
    )
    st.markdown("[Learn more about Colorectal Cancer from the American Cancer Society](https://www.cancer.org/cancer/types/colon-rectal-cancer.html)", unsafe_allow_html=True)

    st.subheader("Feedback")
    feedback = st.text_area("Share your feedback about this application:", key="feedback_home")
    if st.button("Submit Feedback", key="submit_feedback_home"):
        if feedback:
            st.success("Thank you for your feedback!")
            # In a real application, you would log this feedback (e.g., write to a file or database)
            print(f"Feedback received: {feedback}")
        else:
            st.warning("Feedback cannot be empty.")

# --- Colon Cancer Database Search ---
elif menu == "Colon Cancer Database Search":
    st.header("Search Colon Cancer Gene Database")
    st.markdown("Enter a gene name (e.g., `APC`, `KRAS`, `TP53`, `MLH1`) to search our curated information.")

    data_path = "Colon Cancer.xlsx" # Ensure this file exists in the same directory or provide the correct path
    query = st.text_input("Enter Gene Name to Search:", key="gene_search_input", placeholder="e.g., KRAS")

    if st.button("Search Database", key="gene_search_button"):
        if not query:
            st.warning("Please enter a gene name to search.")
        else:
            try:
                data = pd.read_excel(data_path)
                # Standardize column names: strip whitespace, convert to Title Case
                data.columns = data.columns.str.strip().str.title()
                # Try to find the gene name column (case-insensitive check after standardization)
                gene_col = None
                possible_gene_cols = ['Gene Name', 'Gene', 'Symbol', 'Gene Symbol'] # Add likely variations
                for col in possible_gene_cols:
                    if col in data.columns:
                        gene_col = col
                        break

                if not gene_col:
                     st.error(f"Error: Could not find a suitable 'Gene Name' column (tried: {possible_gene_cols}) in the Excel file. Columns found: {list(data.columns)}")
                     st.stop() # Stop execution if the key column isn't found

                # Prepare data for case-insensitive search
                data[gene_col] = data[gene_col].astype(str).str.strip().str.lower()
                search_query = query.strip().lower()

                # Perform the search (contains allows partial matches)
                results = data[data[gene_col].str.contains(search_query, na=False, case=False)].copy() # Use .copy() to avoid SettingWithCopyWarning

                if not results.empty:
                    st.success(f"Found {len(results)} result(s) for '{query}'.")
                    st.write("### Search Results:")

                    # Format links in the results DataFrame
                    formatted_results = results.copy() # Work on a copy
                    for col in formatted_results.columns:
                         try:
                             # Apply the link formatting function, passing column name for context
                             formatted_results[col] = formatted_results[col].apply(lambda x: format_link(x, col))
                         except Exception as apply_e:
                             st.warning(f"Could not apply link formatting to column '{col}': {apply_e}. Displaying raw data.")

                    # Convert the formatted DataFrame to HTML
                    html_table = formatted_results.to_html(
                        escape=False, # Allow HTML tags (like <a>)
                        index=False, # Don't show the DataFrame index
                        na_rep='-', # Representation for missing values
                        justify='left', # Align text to the left
                        classes=['dataframe', 'st-table', 'results-table'] # Add CSS classes if needed
                    )
                    st.write(html_table, unsafe_allow_html=True) # Display the HTML table

                else:
                    st.warning(f"No results found matching '{query}'. Please check the gene name or try a different search term.")

            except FileNotFoundError:
                 st.error(f"Error: The database file '{data_path}' was not found. Please ensure the file exists in the correct directory.")
            except ImportError:
                st.error("Error: The 'openpyxl' library is required to read Excel (.xlsx) files. Please install it (`pip install openpyxl`).")
            except Exception as e:
                st.error(f"An unexpected error occurred during the database search:")
                st.error(f"Error Type: {type(e).__name__}")
                st.error(f"Details: {e}")
                st.caption("Check if the Excel file format is correct and accessible.")

# --- DNA Sequence Analysis ---
elif menu == "DNA Sequence Analysis":
    st.header("DNA Sequence Analysis Tool")
    st.markdown("Enter a DNA sequence containing **A, T, C, G, N** characters for basic analysis.")
    sequence = st.text_area("Enter DNA Sequence:", height=150, key="dna_seq_input", placeholder="e.g., ATGCGTAGCTAGCATGCNNNNATGC")

    if st.button("Analyze DNA Sequence", key="dna_analyze_button"):
        if sequence:
            # Clean and validate sequence
            sequence_cleaned = "".join(sequence.split()).upper()
            valid_bases = {'A', 'T', 'C', 'G', 'N'}
            invalid_chars = set(c for c in sequence_cleaned if c not in valid_bases)

            if invalid_chars:
                st.error(f"Invalid characters found in the sequence: {', '.join(sorted(list(invalid_chars)))}. Please use only A, T, C, G, or N.")
            elif not sequence_cleaned:
                 st.error("The entered sequence is empty after removing whitespace.")
            else:
                try:
                    st.subheader("Analysis Results")
                    # --- Basic Stats ---
                    length = len(sequence_cleaned)
                    st.write(f"**Sequence Length:** {length} bp")

                    # GC Content (excluding Ns)
                    seq_no_n = sequence_cleaned.replace('N', '')
                    if len(seq_no_n) > 0:
                        gc = gc_fraction(seq_no_n) * 100
                        st.write(f"**GC Content (excluding Ns):** {gc:.2f}%")
                    else:
                        st.write("**GC Content:** N/A (Sequence contains only 'N' or is empty after removing 'N')")

                    # Base Composition
                    st.write("**Base Composition:**")
                    composition = Counter(sequence_cleaned)
                    comp_data = [{"Base": base, "Count": count, "Percentage": (count/length*100) if length > 0 else 0}
                                 for base, count in sorted(composition.items())]
                    st.dataframe(pd.DataFrame(comp_data).style.format({"Count":"{:,}","Percentage": "{:.2f}%"}), hide_index=True, use_container_width=True)

                    # --- Plotting ---
                    # Plot composition excluding 'N'
                    plot_composition = {k: v for k, v in composition.items() if k != 'N'}
                    if plot_composition:
                        fig, ax = plt.subplots()
                        colors = {'A': '#1f77b4', 'T': '#ff7f0e', 'C': '#2ca02c', 'G': '#d62728'}
                        bar_bases = sorted(plot_composition.keys())
                        bar_counts = [plot_composition[b] for b in bar_bases]
                        bar_colors = [colors.get(base, '#8c564b') for base in bar_bases] # Default color for safety
                        ax.bar(bar_bases, bar_counts, color=bar_colors)
                        ax.set_xlabel("Base")
                        ax.set_ylabel("Count")
                        ax.set_title("Base Composition (excluding N)")
                        plt.tight_layout()
                        st.pyplot(fig)
                    elif 'N' in composition and len(composition) == 1:
                         st.info("Composition plot not shown as the sequence contains only 'N' bases.")
                    elif not composition:
                         st.info("Composition plot not shown as the sequence is empty.")

                except Exception as e:
                    st.error(f"An error occurred during DNA analysis: {e}")
        else:
            st.warning("Please enter a DNA sequence to analyze.")

# --- Protein Sequence Analysis ---
elif menu == "Protein Sequence Analysis":
    st.header("Protein Sequence Analysis Tool")
    st.markdown("Enter a protein sequence using standard single-letter amino acid codes (A-Y). Non-standard characters (e.g., *, X, B, Z) will be ignored in calculations requiring standard AAs.")
    protein_sequence = st.text_area("Enter Protein Sequence:", height=150, key="protein_seq_input", placeholder="e.g., MKTAYIA...")

    if st.button("Analyze Protein Sequence", key="protein_analyze_button"):
        protein_sequence_cleaned = "".join(protein_sequence.split()).upper() # Clean input
        if not protein_sequence_cleaned:
            st.warning("Please enter a protein sequence.")
        else:
            # Define standard AAs and filter the sequence for Biopython analysis
            standard_aa = "ACDEFGHIKLMNPQRSTVWY"
            protein_for_analysis = re.sub(f'[^{standard_aa}]', '', protein_sequence_cleaned) # Remove non-standard

            # Identify non-standard characters present in the input
            all_input_chars = set(protein_sequence_cleaned)
            standard_set = set(standard_aa)
            non_standard_in_input = all_input_chars - standard_set

            if non_standard_in_input:
                st.warning(f"Input sequence contains non-standard characters: `{'`, `'.join(sorted(list(non_standard_in_input)))}`. "
                           f"These are excluded from calculations like Molecular Weight, pI, GRAVY, Instability, and Secondary Structure.")

            if not protein_for_analysis:
                st.error("Sequence contains only non-standard characters after filtering. Cannot perform standard analysis.")
            else:
                try:
                    analyzed_seq = ProteinAnalysis(protein_for_analysis)

                    st.subheader("Analysis Results (Based on Standard AAs)")
                    st.write(f"**Original Input Length:** {len(protein_sequence_cleaned)} aa")
                    st.write(f"**Length Analyzed (Standard AAs):** {len(protein_for_analysis)} aa")
                    st.write(f"**Molecular Weight:** {analyzed_seq.molecular_weight():,.2f} Da")
                    st.write(f"**Isoelectric Point (pI):** {analyzed_seq.isoelectric_point():.2f}")
                    st.write(f"**GRAVY Score (Grand Average of Hydropathicity):** {analyzed_seq.gravy():.3f} "
                             f"*(Positive=Hydrophobic, Negative=Hydrophilic)*")
                    instability = analyzed_seq.instability_index()
                    stability = "Stable" if instability < 40 else "Unstable"
                    st.write(f"**Instability Index:** {instability:.2f} (*{stability}*)")

                    # Amino Acid Composition (based on the original input sequence)
                    st.subheader("Amino Acid Composition (Based on Original Input)")
                    full_composition = Counter(protein_sequence_cleaned)
                    comp_data = [{"Amino Acid": aa, "Count": count, "Percentage": (count/len(protein_sequence_cleaned)*100)}
                                 for aa, count in sorted(full_composition.items())]
                    st.dataframe(pd.DataFrame(comp_data).style.format({'Count': '{:,}', 'Percentage': '{:.2f}%'}), hide_index=True, use_container_width=True)


                    # Secondary Structure Prediction (based on standard AAs)
                    st.subheader("Secondary Structure Fraction (Predicted - Standard AAs only)")
                    try:
                        helix, turn, sheet = analyzed_seq.secondary_structure_fraction()
                        st.write(f"- **Helix:** {helix*100:.2f}%")
                        st.write(f"- **Turn:** {turn*100:.2f}%")
                        st.write(f"- **Sheet:** {sheet*100:.2f}%")
                    except Exception as ss_e:
                        st.warning(f"Could not calculate secondary structure fractions: {ss_e}")

                except ValueError as ve:
                    st.error(f"An error occurred during protein analysis: {ve}. This might happen if the sequence contains characters not recognized even after filtering.")
                except Exception as e:
                    st.error(f"An unexpected error occurred during protein analysis: {e}")

# --- Bioinformatics Tool (Transcription/Translation) ---
elif menu == "Bioinformatics Tool (Transcription/Translation)":
    st.header("DNA Transcription & Translation Tool")
    st.markdown("Enter a DNA sequence (coding strand, **A, T, C, G** only) to transcribe into RNA and translate into protein using the standard genetic code.")
    sequence = st.text_area("Enter DNA Sequence (Coding Strand):", height=150, key="bioinfo_seq_input", placeholder="e.g., ATGCGTAGCTAGCATGC")

    if st.button("Transcribe and Translate", key="bioinfo_submit_button"):
        sequence_cleaned = "".join(sequence.split()).upper()
        if not sequence_cleaned:
            st.warning("Please enter a DNA sequence.")
        else:
            # Validate sequence
            valid_bases = {'A', 'T', 'C', 'G'}
            invalid_chars = set(c for c in sequence_cleaned if c not in valid_bases)
            if invalid_chars:
                st.error(f"Invalid characters found: {', '.join(sorted(list(invalid_chars)))}. Please use only A, T, C, G.")
            else:
                try:
                    seq_obj = Seq(sequence_cleaned)

                    st.subheader("Results")
                    st.write(f"**Input DNA Length:** {len(seq_obj)} bp")
                    st.write(f"**GC Content:** {gc_fraction(seq_obj) * 100:.2f}%")
                    st.markdown("---")

                    # Reverse Complement
                    st.markdown("**DNA Reverse Complement (5' → 3'):**")
                    st.code(str(seq_obj.reverse_complement()), language='text')
                    st.markdown("---")

                    # Transcription
                    st.markdown("**Transcription (DNA → RNA, 5' → 3'):**")
                    rna_transcript = seq_obj.transcribe()
                    st.code(str(rna_transcript), language='text')
                    st.markdown("---")

                    # Translation
                    st.markdown("**Translation (RNA → Protein):**")
                    st.caption("Uses standard genetic code. Translation stops at the first encountered stop codon ('*').")

                    # Check if length is multiple of 3 for translation
                    remainder = len(rna_transcript) % 3
                    if remainder != 0:
                        st.warning(f"RNA sequence length ({len(rna_transcript)} bases) is not a multiple of 3. The last {remainder} base(s) will be ignored during translation.")
                        rna_transcript_trimmed = rna_transcript[:-remainder] # Trim sequence
                    else:
                        rna_transcript_trimmed = rna_transcript

                    if len(rna_transcript_trimmed) >= 3:
                         # Translate, stopping at the first stop codon
                         protein_translation = rna_transcript_trimmed.translate(to_stop=True)
                         st.code(str(protein_translation), language='text')
                         st.write(f"- **Translated Protein Length:** {len(protein_translation)} aa")

                         # Optionally show full translation including internal stops
                         protein_full_translation = rna_transcript_trimmed.translate(to_stop=False)
                         # Only show the full translation if it's different AND contains internal stops
                         # (Don't show if the only '*' is the standard terminal stop)
                         if '*' in protein_full_translation and protein_full_translation != protein_translation and not (len(protein_full_translation) == len(protein_translation)+1 and protein_full_translation.endswith('*')):
                            st.caption("**Full Translation (including internal stop codons):**")
                            st.code(str(protein_full_translation), language='text')

                    elif len(rna_transcript) < 3: # Original length check
                         st.error("Cannot translate: RNA sequence is shorter than 3 bases.")
                    else: # Trimmed length is 0
                         st.error("Cannot translate: RNA sequence length is zero after trimming to a multiple of 3.")

                except Exception as e:
                    st.error(f"An error occurred during processing: {e}")

# --- ORF Finder ---
elif menu == "ORF Finder":
    st.header("Open Reading Frame (ORF) Finder")
    st.markdown("Identifies potential protein-coding regions (ORFs) in a DNA sequence (**A, T, C, G** only) based on standard start (ATG) and stop codons (TAA, TAG, TGA) in all six reading frames.")

    sequence = st.text_area("Enter DNA Sequence:", height=150, key="orf_seq_input", placeholder="e.g., AGATG...")
    min_len_aa = st.number_input("Minimum ORF Length (in Amino Acids):", min_value=10, value=30, step=5, key="orf_min_len", help="Minimum number of amino acids (excluding stop codon) for an ORF to be reported.")

    if st.button("Find ORFs", key="orf_find_button"):
        sequence_cleaned = "".join(sequence.split()).upper()
        if not sequence_cleaned:
            st.warning("Please enter a DNA sequence.")
        else:
            valid_bases = {'A', 'T', 'C', 'G'}
            invalid_chars = set(c for c in sequence_cleaned if c not in valid_bases)
            if invalid_chars:
                st.error(f"Sequence contains invalid characters: {', '.join(sorted(list(invalid_chars)))}. Use only A, T, C, G.")
            else:
                try:
                    st.subheader("ORF Results")
                    seq_obj = Seq(sequence_cleaned)
                    seq_rev_comp = seq_obj.reverse_complement()
                    min_len_bp = min_len_aa * 3 # Minimum length in base pairs (for DNA segment)

                    standard_table = CodonTable.unambiguous_dna_by_id[1]
                    start_codons = standard_table.start_codons
                    stop_codons = standard_table.stop_codons

                    orfs_found = []

                    # Process 3 forward frames and 3 reverse frames
                    for strand, nuc_seq in [("+", seq_obj), ("-", seq_rev_comp)]:
                        for frame in range(3):
                            seq_len = len(nuc_seq)
                            frame_seq = str(nuc_seq[frame:]) # Get sequence for this frame

                            # Simple regex to find potential ORFs (Start ... Stop)
                            # This finds the *longest* ORF between a start and the *first* stop
                            # It doesn't find nested ORFs easily but is simpler
                            for start_match in re.finditer('ATG', frame_seq):
                                start_index_in_frame = start_match.start()
                                # Look for the first stop codon *after* this start codon in the same frame
                                end_index_in_frame = -1
                                min_stop_index = float('inf')

                                for stop_codon in stop_codons:
                                    try:
                                        # Find stop codon index relative to the frame sequence
                                        # Ensure the stop codon search starts *after* the start codon
                                        stop_index = frame_seq.index(stop_codon, start_index_in_frame + 3)
                                        # Ensure the stop codon is in frame
                                        if (stop_index - start_index_in_frame) % 3 == 0:
                                            if stop_index < min_stop_index:
                                                min_stop_index = stop_index
                                    except ValueError:
                                        continue # Stop codon not found after start

                                if min_stop_index != float('inf'):
                                    end_index_in_frame = min_stop_index

                                    # Calculate lengths and original sequence positions
                                    orf_len_bp = end_index_in_frame + 3 - start_index_in_frame # Includes stop codon bp
                                    orf_len_aa = orf_len_bp // 3 - 1 # Excludes stop codon

                                    if orf_len_aa >= min_len_aa:
                                        # Map frame indices back to original sequence indices
                                        if strand == '+':
                                            orf_start_orig = frame + start_index_in_frame + 1 # 1-based
                                            orf_end_orig = frame + end_index_in_frame + 3   # 1-based end position of stop codon
                                        else: # Reverse strand
                                            orf_start_orig = len(seq_obj) - (frame + end_index_in_frame + 3) + 1 # 1-based start
                                            orf_end_orig = len(seq_obj) - (frame + start_index_in_frame) # 1-based end

                                        orf_dna_seq = frame_seq[start_index_in_frame : end_index_in_frame + 3]
                                        orf_prot_seq = str(Seq(orf_dna_seq).translate(to_stop=True)) # Translate ORF DNA

                                        orfs_found.append({
                                            "Strand": strand,
                                            "Frame": frame + 1, # 1, 2, 3
                                            "Start (bp)": orf_start_orig,
                                            "End (bp)": orf_end_orig,
                                            "Length (bp)": orf_len_bp,
                                            "Length (aa)": orf_len_aa,
                                            #"ORF DNA": orf_dna_seq[:60] + ('...' if len(orf_dna_seq)>60 else ''), # Show snippet
                                            "Protein": orf_prot_seq[:60] + ('...' if len(orf_prot_seq)>60 else '') # Show snippet
                                        })
                                        # Optional: break here if we only want the first ORF starting at `start_match.start()`
                                        # break # Move to the next potential start codon


                    if orfs_found:
                        orf_df = pd.DataFrame(orfs_found)
                        orf_df = orf_df.sort_values(by="Length (aa)", ascending=False)
                        st.success(f"Found {len(orfs_found)} potential ORF(s) meeting the minimum length criteria ({min_len_aa} aa).")
                        st.dataframe(orf_df.style.format({"Start (bp)":"{:,}", "End (bp)":"{:,}", "Length (bp)":"{:,}", "Length (aa)":"{:,}"}),
                                     hide_index=True, use_container_width=True)
                    else:
                        st.warning(f"No ORFs found with a minimum length of {min_len_aa} amino acids.")

                except Exception as e:
                    st.error(f"An error occurred during ORF finding: {e}")


# --- Primer Design ---
elif menu == "Primer Design":
    st.header("Basic Primer Design Tool")
    st.markdown("Enter a DNA template sequence (**A, T, C, G** only) to generate simple forward and reverse primers based on sequence ends.")
    st.info("⚠️ This tool provides basic suggestions only. For real-world applications, primers must be validated for melting temperature (Tm), specificity (e.g., using BLAST), potential hairpins, dimers, and other factors using specialized software.")
    sequence = st.text_area("Enter DNA Template Sequence:", height=150, key="primer_seq_input", placeholder="e.g., ATGCGTACGT...")
    primer_length = st.slider("Desired Primer Length:", min_value=15, max_value=30, value=20, step=1, key="primer_len_slider", help="Select the length for both forward and reverse primers.")

    if st.button("Design Basic Primers", key="primer_design_button"):
        sequence_cleaned = "".join(sequence.split()).upper()

        if not sequence_cleaned:
            st.warning("Please enter a DNA template sequence.")
        else:
            # Validate sequence
            valid_bases = {'A', 'T', 'C', 'G'}
            invalid_chars = set(c for c in sequence_cleaned if c not in valid_bases)
            if invalid_chars:
                st.error(f"Invalid characters found: {', '.join(sorted(list(invalid_chars)))}. Please use only A, T, C, G.")
            elif len(sequence_cleaned) < primer_length * 2:
                st.error(f"Sequence is too short ({len(sequence_cleaned)} bp). Needs to be at least twice the primer length ({primer_length * 2} bp).")
            else:
                try:
                    # Forward primer is the start of the sequence
                    forward_primer_seq = Seq(sequence_cleaned[:primer_length])
                    # Reverse primer is the reverse complement of the end of the sequence
                    reverse_primer_template = Seq(sequence_cleaned[-primer_length:])
                    reverse_primer_seq = reverse_primer_template.reverse_complement()

                    st.subheader("Suggested Primers")
                    col1, col2 = st.columns(2)
                    with col1:
                        st.markdown("**Forward Primer (5' → 3'):**")
                        st.code(str(forward_primer_seq), language='text')
                        st.write(f"- Length: {len(forward_primer_seq)} bp")
                        st.write(f"- GC Content: {gc_fraction(forward_primer_seq)*100:.2f}%")
                        # Add basic Tm calculation (optional, very approximate)
                        # from Bio.SeqUtils import MeltingTemp as mt
                        # st.write(f"- Approx Tm (Wallace): {mt.Tm_Wallace(str(forward_primer_seq)):.2f} °C")
                        # st.write(f"- Approx Tm (GC): {mt.Tm_GC(str(forward_primer_seq)):.2f} °C")
                    with col2:
                        st.markdown("**Reverse Primer (5' → 3'):**")
                        st.code(str(reverse_primer_seq), language='text')
                        st.write(f"- Length: {len(reverse_primer_seq)} bp")
                        st.write(f"- GC Content: {gc_fraction(reverse_primer_seq)*100:.2f}%")
                        # st.write(f"- Approx Tm (Wallace): {mt.Tm_Wallace(str(reverse_primer_seq)):.2f} °C")
                        # st.write(f"- Approx Tm (GC): {mt.Tm_GC(str(reverse_primer_seq)):.2f} °C")

                    # Amplicon size is the length of the original template
                    amplicon_length = len(sequence_cleaned)
                    st.write(f"**Expected Amplicon Size:** {amplicon_length} bp (entire template)")

                except Exception as e:
                    st.error(f"An error occurred during primer design: {e}")

# --- Restriction Enzyme Analysis ---
elif menu == "Restriction Enzyme Analysis":
    st.header("Restriction Enzyme Analysis Tool")
    st.markdown("Finds cut sites for specified restriction enzymes in a DNA sequence (**A, T, C, G** only). Uses enzyme definitions from Biopython's REBASE data.")

    dna_sequence = st.text_area("Enter DNA Sequence:", height=150, key="re_dna_input", placeholder="e.g., GAATTC...")

    # Get available enzymes from Biopython
    try:
        all_enzymes = sorted([str(e) for e in Restriction.AllEnzymes])
    except Exception as e:
        st.error(f"Could not load enzyme list from Biopython: {e}")
        all_enzymes = ["EcoRI", "HindIII", "BamHI"] # Fallback examples

    selected_enzymes = st.multiselect(
        "Select Restriction Enzyme(s):",
        options=all_enzymes,
        default=["EcoRI", "HindIII"],
        key="re_enzyme_select",
        help="Choose one or more enzymes from the list (uses REBASE names)."
    )

    is_linear = st.checkbox("Treat DNA as Linear", value=True, key="re_linear_check", help="Uncheck if the DNA sequence represents a circular plasmid.")

    if st.button("Analyze Restriction Sites", key="re_analyze_button"):
        dna_sequence_cleaned = "".join(dna_sequence.split()).upper()
        if not dna_sequence_cleaned:
            st.warning("Please enter a DNA sequence.")
        elif not selected_enzymes:
            st.warning("Please select at least one restriction enzyme.")
        else:
            valid_bases = {'A', 'T', 'C', 'G'}
            invalid_chars = set(c for c in dna_sequence_cleaned if c not in valid_bases)
            if invalid_chars:
                st.error(f"Sequence contains invalid characters: {', '.join(sorted(list(invalid_chars)))}. Use only A, T, C, G.")
            else:
                try:
                    seq_obj = Seq(dna_sequence_cleaned)
                    batch = Restriction.RestrictionBatch(selected_enzymes)
                    analysis_results = batch.search(seq_obj, linear=is_linear)

                    st.subheader("Restriction Analysis Results")
                    st.write(f"**Sequence Length:** {len(seq_obj)} bp ({'Linear' if is_linear else 'Circular'})")

                    results_data = []
                    found_sites = False
                    for enzyme, sites in analysis_results.items():
                        if sites:
                            found_sites = True
                            # Convert 0-based Biopython indices to 1-based for display
                            sites_1based = sorted([s + 1 for s in sites])
                            results_data.append({
                                "Enzyme": str(enzyme),
                                "Recognition Site": enzyme.site,
                                "# Sites": len(sites_1based),
                                "Cut Positions (1-based)": ", ".join(map(str, sites_1based))
                            })
                        else:
                            # Optionally show enzymes with no sites
                            results_data.append({
                                "Enzyme": str(enzyme),
                                "Recognition Site": enzyme.site,
                                "# Sites": 0,
                                "Cut Positions (1-based)": "None found"
                            })

                    if results_data:
                         results_df = pd.DataFrame(results_data)
                         st.dataframe(results_df, hide_index=True, use_container_width=True)
                         if not found_sites:
                              st.info("None of the selected enzymes cut the provided sequence.")
                    else:
                         st.error("Analysis failed or returned no results.") # Should not happen if batch is created

                    # Optional: Fragment calculation (more complex)
                    # try:
                    #    analysis = Restriction.Analysis(batch, seq_obj, linear=is_linear)
                    #    # analysis.print_that() # Prints a lot of info
                    #    fragments = analysis.fragments
                    #    st.subheader("Predicted Fragments (bp)")
                    #    st.write(sorted(fragments, reverse=True))
                    # except Exception as frag_e:
                    #    st.warning(f"Could not calculate fragments: {frag_e}")


                except ValueError as ve:
                    st.error(f"Error during analysis: {ve}. This might be due to an invalid enzyme name if not using the dropdown, or an issue with the sequence.")
                except Exception as e:
                    st.error(f"An unexpected error occurred during restriction analysis: {e}")

# --- Motif Finder Tool ---
elif menu == "Motif Finder Tool":
    st.header("Motif Finder Tool")
    st.markdown("Find occurrences of a motif (exact string or simple regular expression) within a DNA or Protein sequence.")
    st.info("Supports basic regex syntax like `.` (any char), `[AT]` (A or T), `*` (0 or more), `+` (1 or more), `?` (0 or 1), `{n,m}` (n to m repeats). Use `\` to escape special characters.")

    sequence = st.text_area("Enter Sequence (DNA or Protein):", height=150, key="motif_seq_input")
    motif = st.text_input("Enter Motif / Regex Pattern:", key="motif_pattern_input", placeholder="e.g., ATG or Y[DN]R or ^M.*K$")
    # seq_type = st.radio("Sequence Type (for context only):", ("DNA", "Protein"), horizontal=True, key="motif_seq_type") # Context only, doesn't affect search
    case_sensitive = st.checkbox("Case Sensitive Search", value=False, key="motif_case_sensitive")

    if st.button("Find Motif Occurrences", key="motif_find_button"):
        if not sequence or not motif:
            st.warning("Please provide both the sequence and the motif/pattern to search for.")
        else:
            sequence_cleaned = "".join(sequence.split()) # Remove whitespace only
            if not sequence_cleaned:
                st.warning("The sequence is empty after removing whitespace.")
            else:
                 flags = 0 if case_sensitive else re.IGNORECASE
                 try:
                     matches = list(re.finditer(motif, sequence_cleaned, flags=flags))
                     st.subheader("Motif Search Results")

                     if matches:
                         st.success(f"Found {len(matches)} occurrence(s) of the motif.")
                         match_data = [
                             {"Position (1-based Start)": m.start() + 1,
                              "Position (End)": m.end(), # re end is exclusive, so this is correct
                              "Matched String": m.group()}
                             for m in matches
                         ]
                         st.dataframe(match_data, hide_index=True, use_container_width=True)
                     else:
                         st.warning("Motif/pattern not found in the provided sequence with the current settings.")

                 except re.error as e:
                     st.error(f"Invalid regular expression pattern provided: {e}. Please check the syntax.")
                 except Exception as e:
                     st.error(f"An error occurred during the motif search: {e}")

# --- K-mer Counter ---
elif menu == "K-mer Counter":
    st.header("K-mer Counter Tool")
    st.markdown("Counts the occurrences of all overlapping subsequences of length 'k' (k-mers) in a given DNA or Protein sequence.")

    sequence = st.text_area("Enter Sequence (DNA or Protein):", height=150, key="kmer_seq_input")
    k_value = st.number_input("K-mer Length (k):", min_value=1, value=3, step=1, key="kmer_k_value", help="Length of the subsequence to count.")

    if st.button("Count K-mers", key="kmer_count_button"):
        sequence_cleaned = "".join(sequence.split()).upper() # Clean and uppercase
        if not sequence_cleaned:
            st.warning("Please enter a sequence.")
        elif k_value <= 0:
            st.warning("K-mer length (k) must be greater than 0.")
        elif k_value > len(sequence_cleaned):
            st.warning(f"K-mer length ({k_value}) cannot be greater than the sequence length ({len(sequence_cleaned)}).")
        else:
            try:
                st.subheader(f"{k_value}-mer Count Results")
                st.write(f"**Sequence Length:** {len(sequence_cleaned)} bases/residues")
                st.write(f"**K-mer Length (k):** {k_value}")

                # Generate k-mers
                kmers = [sequence_cleaned[i:i + k_value] for i in range(len(sequence_cleaned) - k_value + 1)]
                total_kmers = len(kmers)
                st.write(f"**Total K-mers Found:** {total_kmers:,}")

                if total_kmers > 0:
                    kmer_counts = Counter(kmers)
                    kmer_data = []
                    for kmer, count in kmer_counts.items():
                        frequency = (count / total_kmers) * 100
                        kmer_data.append({"K-mer": kmer, "Count": count, "Frequency (%)": frequency})

                    kmer_df = pd.DataFrame(kmer_data)
                    kmer_df = kmer_df.sort_values(by="Count", ascending=False) # Sort by count

                    st.dataframe(kmer_df.style.format({'Count': '{:,}', 'Frequency (%)': '{:.3f}%'}),
                                 hide_index=True, use_container_width=True)

                    # --- Plotting Option ---
                    st.write("---")
                    show_plot = st.checkbox("Show Plot of Top K-mers", value=True, key="kmer_show_plot")
                    if show_plot:
                        top_n = st.slider("Number of Top K-mers to Plot:", min_value=5, max_value=min(50, len(kmer_df)), value=min(20, len(kmer_df)), step=1, key="kmer_top_n")
                        if top_n > 0 and len(kmer_df) > 0:
                            plot_data = kmer_df.head(top_n)
                            fig, ax = plt.subplots(figsize=(10, 5))
                            ax.bar(plot_data['K-mer'], plot_data['Count'], color='teal')
                            ax.set_xlabel(f"{k_value}-mer")
                            ax.set_ylabel("Count")
                            ax.set_title(f"Top {top_n} Most Frequent {k_value}-mers")
                            plt.xticks(rotation=60, ha='right', fontsize=8)
                            plt.tight_layout()
                            st.pyplot(fig)
                        elif len(kmer_df) == 0:
                            st.info("No k-mers to plot.")

                else:
                    st.info("No k-mers could be generated (sequence might be shorter than k).")

            except Exception as e:
                st.error(f"An error occurred during k-mer counting: {e}")

# --- Genome Coverage Plotter ---
elif menu == "Genome Coverage Plotter":
    st.header("Genome Coverage Plotter")
    st.markdown("Upload a coverage data file (e.g., from `samtools depth`) containing genomic positions and coverage values. The tool expects columns like 'position' and 'coverage'.")
    uploaded_file = st.file_uploader("Upload Coverage Data (CSV, TSV, TXT):", type=['csv', 'tsv', 'txt'], key="coverage_uploader")

    if uploaded_file:
        st.write("---")
        st.subheader("File Settings")
        col1, col2 = st.columns(2)
        with col1:
            # Improved separator detection/selection
            sep_options = {',': 'Comma (,)', '\t': 'Tab (\\t)', ' ': 'Space ( )'}
            # Try to guess separator by reading first few lines
            try:
                 peek_data = uploaded_file.read(1024).decode('utf-8')
                 uploaded_file.seek(0) # Reset file pointer!
                 sniffer = io.StringIO(peek_data)
                 dialect = pd.io.common.sniff_delimiter(sniffer.read(), list(sep_options.keys()))
                 default_sep_index = list(sep_options.keys()).index(dialect) if dialect in sep_options else 1 # Default to tab if guess fails
            except Exception:
                 default_sep_index = 1 # Default to tab on error

            selected_sep_display = st.selectbox("Column Separator:", options=list(sep_options.keys()), index=default_sep_index, format_func=lambda x: sep_options[x], key="cov_sep")
            separator = selected_sep_display if selected_sep_display != ' ' else r'\s+' # Use regex for space to handle multiple spaces
        with col2:
            # Header detection/setting
            use_header = st.checkbox("File has a header row", value=True, key="cov_has_header", help="Check this if the first row contains column names.")
            header_row_num = st.number_input("Header Row Index (if applicable):", min_value=0, value=0, step=1, key="cov_header_num", help="Row number (0-based) containing headers. Only used if 'File has a header row' is checked.")
            header_arg = header_row_num if use_header else None

        st.write("---")
        st.subheader("Data Preview & Column Selection")

        try:
            # Read the CSV/TSV file
            # Add comment='#' to handle potential comment lines (e.g., in some BED files or depth outputs)
            df = pd.read_csv(uploaded_file, sep=separator, header=header_arg, engine='python', comment='#', low_memory=False)

            if df.empty:
                 st.error("File Reading Error: The file appears empty or contains only comments/empty lines after applying settings.")
            else:
                st.write("Preview of uploaded data (first 5 rows):")
                st.dataframe(df.head(), height=200)

                if use_header:
                    available_cols = list(df.columns)
                    # Normalize column names for guessing (lower case, strip spaces)
                    normalized_cols_map = {str(col).strip().lower(): col for col in available_cols}
                    pos_col_guess, cov_col_guess = None, None
                    # Guess position column
                    for p_name in ['position', 'pos', 'coordinate', 'location', 'start', 'chromstart', 'base', 'pos.', 'position(bp)', 'end', 'chromend']: # Added common names
                        if p_name in normalized_cols_map: pos_col_guess = normalized_cols_map[p_name]; break
                    if not pos_col_guess and len(available_cols) > 1: pos_col_guess = available_cols[1] # Fallback guess col 1
                    elif not pos_col_guess: pos_col_guess = available_cols[0] # Fallback guess col 0

                    # Guess coverage column
                    for c_name in ['coverage', 'cov', 'depth', 'reads', 'count', 'value', 'score', 'mean_coverage', 'depth_of_coverage']: # Added common names
                        if c_name in normalized_cols_map: cov_col_guess = normalized_cols_map[c_name]; break
                    if not cov_col_guess and len(available_cols) > 2: cov_col_guess = available_cols[2] # Fallback guess col 2
                    elif not cov_col_guess and len(available_cols) > 1 and available_cols[1] != pos_col_guess: cov_col_guess = available_cols[1] # Fallback guess col 1 if different from pos
                    elif not cov_col_guess: cov_col_guess = available_cols[0] # Fallback guess col 0

                else: # No header
                    available_cols = list(range(df.shape[1]))
                    st.info("No header specified. Please select columns by their index (0, 1, 2...).")
                    # Guess based on common samtools depth output (chr, pos, depth) -> indices 1 and 2
                    pos_col_guess = 1 if len(available_cols) > 1 else 0
                    cov_col_guess = 2 if len(available_cols) > 2 else (1 if len(available_cols) > 1 else 0)

                # Ensure guesses are valid indices/columns
                pos_col_idx = available_cols.index(pos_col_guess) if pos_col_guess in available_cols else 0
                cov_col_idx = available_cols.index(cov_col_guess) if cov_col_guess in available_cols else (1 if len(available_cols)>1 else 0)

                # Column selection widgets
                col_sel1, col_sel2 = st.columns(2)
                with col_sel1: pos_col_selected = st.selectbox("Select Position Column:", available_cols, index=pos_col_idx, key="cov_pos_col")
                with col_sel2: cov_col_selected = st.selectbox("Select Coverage Column:", available_cols, index=cov_col_idx, key="cov_cov_col")

                if pos_col_selected == cov_col_selected and len(available_cols)>1:
                     st.error("Position and Coverage columns cannot be the same.")
                else:
                    st.success(f"Selected: Position='{pos_col_selected}', Coverage='{cov_col_selected}'.")
                    st.write("---")
                    st.subheader("Plotting Options & Results")
                    try:
                        # Select and rename columns for consistency
                        plot_df = df[[pos_col_selected, cov_col_selected]].copy()
                        plot_df.columns = ['Position', 'Coverage']

                        # Convert to numeric, coercing errors to NaN
                        plot_df['Position'] = pd.to_numeric(plot_df['Position'], errors='coerce')
                        plot_df['Coverage'] = pd.to_numeric(plot_df['Coverage'], errors='coerce')

                        # Handle missing/non-numeric values
                        rows_before = len(plot_df)
                        plot_df.dropna(inplace=True)
                        rows_after = len(plot_df)
                        if rows_after < rows_before: st.warning(f"Removed {rows_before - rows_after} rows containing non-numeric or missing data in selected columns.")

                        if plot_df.empty: st.error("No valid numeric data found in the selected columns after cleaning.")
                        else:
                            # Ensure data is sorted by position for plotting
                            plot_df = plot_df.sort_values(by='Position').reset_index(drop=True)

                            # Downsampling for large datasets (optional but recommended)
                            max_points_to_plot = 20000 # Limit points for performance
                            if len(plot_df) > max_points_to_plot:
                                st.info(f"Dataset has {len(plot_df):,} points. Downsampling to {max_points_to_plot:,} for plotting performance.")
                                # Simple uniform sampling
                                indices = np.linspace(0, len(plot_df) - 1, max_points_to_plot, dtype=int)
                                plot_df_display = plot_df.iloc[indices]
                            else:
                                plot_df_display = plot_df

                            # --- Plotting ---
                            fig, ax = plt.subplots(figsize=(12, 5))
                            ax.plot(plot_df_display['Position'], plot_df_display['Coverage'], label='Coverage', lw=0.8, color='steelblue', alpha=0.8) # Thinner line for raw

                            # Fill option
                            fill_plot = st.checkbox("Fill area under curve", value=True, key="cov_fill")
                            if fill_plot: ax.fill_between(plot_df_display['Position'], plot_df_display['Coverage'], alpha=0.3, color='lightblue')

                            # Smoothing option (applied to the *original* sorted data before potential downsampling for display)
                            smoothing_window = st.slider("Smoothing Window Size (0 for none):", 0, min(200, len(plot_df)//10), 0, 5, key="cov_smooth", help="Applies a rolling mean. 0 = no smoothing. Window size is in number of data points.")
                            if smoothing_window > 1:
                                # Calculate smoothed data on the full dataset
                                smoothed_coverage = plot_df['Coverage'].rolling(window=smoothing_window, center=True, min_periods=1).mean()
                                # Select the corresponding smoothed points if downsampling was applied
                                smoothed_coverage_display = smoothed_coverage.iloc[indices] if len(plot_df) > max_points_to_plot else smoothed_coverage
                                # Use the display positions for plotting
                                ax.plot(plot_df_display['Position'], smoothed_coverage_display, label=f'Smoothed ({smoothing_window}pt window)', lw=1.5, color='darkorange', linestyle='-')


                            # Plot formatting
                            ax.set_xlabel(f"Genomic Position ({pos_col_selected})")
                            ax.set_ylabel(f"Coverage Depth ({cov_col_selected})")
                            ax.set_title("Genome Coverage Plot")
                            ax.legend()
                            ax.grid(True, linestyle='--', alpha=0.6)
                            ax.ticklabel_format(style='plain', axis='x', useOffset=False) # Prevent scientific notation on position axis
                            min_pos, max_pos = plot_df_display['Position'].min(), plot_df_display['Position'].max()
                            ax.set_xlim(min_pos, max_pos)
                            ax.set_ylim(bottom=0) # Coverage cannot be negative
                            plt.tight_layout()
                            st.pyplot(fig)

                            # --- Statistics ---
                            st.subheader("Coverage Statistics (Based on Full Dataset)")
                            stats_df = plot_df['Coverage'].describe(percentiles=[.01, .1, .25, .5, .75, .9, .99]).to_frame().T # Add more percentiles
                            stats_df = stats_df.rename(columns={'count': 'Positions Covered', 'mean': 'Mean Coverage', 'std': 'Std Dev', 'min': 'Min Coverage', '1%':'1st Pctl', '10%':'10th Pctl', '25%': '25th Pctl', '50%': 'Median (50th Pctl)', '75%': '75th Pctl', '90%':'90th Pctl', '99%':'99th Pctl', 'max': 'Max Coverage'})
                            st.dataframe(stats_df.style.format("{:,.2f}", na_rep="-"))

                            # Breadth of coverage calculation
                            min_cov_threshold = st.number_input("Calculate % positions with coverage >=", 0, value=10, key="cov_breadth_thresh", help="Minimum coverage depth to be considered 'covered'.")
                            if min_cov_threshold >= 0:
                                 positions_above_thresh = (plot_df['Coverage'] >= min_cov_threshold).sum()
                                 # Calculate total range covered by data, not just number of points, if possible
                                 total_span = plot_df['Position'].max() - plot_df['Position'].min() + 1
                                 # Using number of unique positions as denominator is safer if data is sparse
                                 total_positions_in_data = len(plot_df['Position']) # Or use plot_df['Position'].nunique() if positions can be repeated

                                 if total_positions_in_data > 0:
                                     breadth_pct = (positions_above_thresh / total_positions_in_data) * 100
                                     st.metric(label=f"Breadth of Coverage (>= {min_cov_threshold}x)", value=f"{breadth_pct:.2f}%", delta=f"{positions_above_thresh:,} / {total_positions_in_data:,} positions", help="Percentage of data points in the file meeting the coverage threshold.")
                                 else:
                                     st.metric(label=f"Breadth of Coverage (>= {min_cov_threshold}x)", value="N/A", delta="No positions in data")

                    except Exception as plot_e: st.error(f"An error occurred during data processing or plotting: {plot_e}\nPlease ensure the correct columns are selected and contain numeric data.")

        except pd.errors.EmptyDataError: st.error("File Reading Error: The uploaded file is empty or contains no data.")
        except ValueError as ve: st.error(f"File Reading Error: Possible issue with separator or header row setting. Details: {ve}")
        except Exception as read_e: st.error(f"File Reading Error: Could not parse the file. Please check the format, separator, and header settings. Error: {read_e}")

# --- Dot Plot Comparison ---
elif menu == "Dot Plot Comparison":
    st.header("Dot Plot Sequence Comparison")
    st.markdown("Generates a dot plot to visualize regions of similarity between two sequences (DNA or Protein). A dot is placed where characters match.")
    st.info("For longer sequences, the plot can become dense. Consider using the windowed approach for clearer results on larger inputs.")

    col1, col2 = st.columns(2)
    with col1:
        seq1 = st.text_area("Sequence 1 (Horizontal Axis):", height=150, key="dot_seq1_input", placeholder="Enter first sequence...")
    with col2:
        seq2 = st.text_area("Sequence 2 (Vertical Axis):", height=150, key="dot_seq2_input", placeholder="Enter second sequence...")

    st.write("---")
    st.subheader("Dot Plot Options")
    use_window = st.checkbox("Use Windowed Comparison", value=False, key="dot_use_window", help="Compare windows of characters instead of single characters to reduce noise.")

    window_size, threshold = 1, 1
    if use_window:
        col_w1, col_w2 = st.columns(2)
        with col_w1:
            window_size = st.slider("Window Size:", min_value=2, max_value=20, value=5, step=1, key="dot_window_size", help="Size of the sliding window for comparison.")
        with col_w2:
            threshold = st.slider("Match Threshold:", min_value=1, max_value=window_size, value=max(1, window_size // 2), step=1, key="dot_threshold", help="Minimum number of matches within the window to plot a dot.")

    dot_size = st.slider("Dot Size:", min_value=1, max_value=10, value=2, step=1, key="dot_size", help="Adjust the size of the points on the plot.")

    if st.button("Generate Dot Plot", key="dot_generate_button"):
        seq1_cleaned = "".join(seq1.split()).upper()
        seq2_cleaned = "".join(seq2.split()).upper()

        if not seq1_cleaned or not seq2_cleaned:
            st.warning("Please enter both sequences.")
        elif use_window and (window_size <= 0 or threshold <= 0 or threshold > window_size):
            st.warning("Invalid window/threshold settings. Ensure Window Size > 0, Threshold > 0, and Threshold <= Window Size.")
        else:
            st.subheader("Dot Plot Results")
            len1, len2 = len(seq1_cleaned), len(seq2_cleaned)
            st.write(f"Sequence 1 Length: {len1}")
            st.write(f"Sequence 2 Length: {len2}")

            # Limit sequence length to prevent browser freezing / excessive computation
            max_len_dotplot = 2000
            if len1 > max_len_dotplot or len2 > max_len_dotplot:
                 st.warning(f"One or both sequences exceed the maximum recommended length ({max_len_dotplot}) for direct plotting. Consider using shorter sequences or specialized tools for large comparisons.")
                 # Optionally truncate or disable button here, for now just warn.

            else:
                try:
                    # Generate match coordinates
                    match_coords_x = []
                    match_coords_y = []

                    if not use_window: # Simple direct comparison
                        for i in range(len1):
                            for j in range(len2):
                                if seq1_cleaned[i] == seq2_cleaned[j]:
                                    match_coords_x.append(i)
                                    match_coords_y.append(j)
                    else: # Windowed comparison
                        for i in range(len1 - window_size + 1):
                            for j in range(len2 - window_size + 1):
                                sub_seq1 = seq1_cleaned[i : i + window_size]
                                sub_seq2 = seq2_cleaned[j : j + window_size]
                                matches = sum(c1 == c2 for c1, c2 in zip(sub_seq1, sub_seq2))
                                if matches >= threshold:
                                    # Plot dot at the center of the window
                                    match_coords_x.append(i + window_size // 2)
                                    match_coords_y.append(j + window_size // 2)

                    # --- Plotting ---
                    if not match_coords_x:
                        st.info("No matches found based on the current settings.")
                    else:
                        fig, ax = plt.subplots(figsize=(8, 8)) # Make it square
                        # Use scatter plot for dots
                        ax.scatter(match_coords_x, match_coords_y, s=dot_size, c='black', marker='.')

                        ax.set_xlabel(f"Sequence 1 Position (Length: {len1})")
                        ax.set_ylabel(f"Sequence 2 Position (Length: {len2})")
                        ax.set_title(f"Dot Plot Comparison{' (Windowed)' if use_window else ''}")

                        # Set limits to match sequence lengths (0-based index)
                        ax.set_xlim(0, len1)
                        ax.set_ylim(0, len2)
                        ax.invert_yaxis() # Common convention for dot plots

                        # Improve appearance
                        ax.grid(False) # Grid can be distracting
                        ax.set_aspect('equal', adjustable='box') # Ensure proper aspect ratio if possible

                        plt.tight_layout()
                        st.pyplot(fig)

                except Exception as e:
                    st.error(f"An error occurred while generating the dot plot: {e}")


# --- Variant Annotation Tool ---
elif menu == "Variant Annotation Tool":
    st.header("Simple Variant Substitution Tool")
    st.markdown("Introduce a single nucleotide substitution into a reference DNA sequence (**A, T, C, G** only).")
    ref_seq = st.text_area("Reference DNA Sequence:", height=100, key="var_ref_seq", placeholder="e.g., ATGCGTACGTAGCTAG")

    col1, col2, col3 = st.columns(3)
    with col1:
        variant_position = st.number_input("Variant Position (1-based index):", min_value=1, step=1, key="var_pos", help="The position in the sequence to change (starting from 1).")
    with col2:
        ref_base_display = "?" # Placeholder
        if ref_seq and variant_position > 0 and variant_position <= len("".join(ref_seq.split())):
            ref_base_display = "".join(ref_seq.split()).upper()[variant_position - 1]
        st.text_input("Reference Base (REF):", value=ref_base_display, key="var_ref_display", disabled=True, help="The base currently at the specified position.")
    with col3:
        variant_base = st.selectbox("New Variant Base (ALT):", ["A", "T", "C", "G"], index=0, key="var_base", help="The new base to insert at the specified position.")


    if st.button("Apply Variant Substitution", key="var_apply_button"):
        ref_seq_cleaned = "".join(ref_seq.split()).upper()
        error_flag = False

        # --- Input Validations ---
        if not ref_seq_cleaned:
            st.warning("Please enter a reference DNA sequence.")
            error_flag = True
        if not variant_base: # Should not happen with selectbox, but good practice
            st.warning("Please select a new variant base.")
            error_flag = True
        valid_bases = {'A', 'T', 'C', 'G'}
        invalid_ref_chars = set(c for c in ref_seq_cleaned if c not in valid_bases)
        if invalid_ref_chars:
            st.error(f"Reference sequence contains invalid characters: {', '.join(sorted(list(invalid_ref_chars)))}. Use only A, T, C, G.")
            error_flag = True
        # Check position validity only if sequence is valid
        if not invalid_ref_chars and ref_seq_cleaned and (variant_position <= 0 or variant_position > len(ref_seq_cleaned)):
            st.error(f"Variant position ({variant_position}) is out of range for the sequence length ({len(ref_seq_cleaned)}).")
            error_flag = True
        # --- End Validations ---

        if not error_flag:
            try:
                zero_based_pos = variant_position - 1
                original_base = ref_seq_cleaned[zero_based_pos]

                # Check if the change is redundant
                if original_base == variant_base:
                    st.info(f"The base at position {variant_position} is already '{original_base}'. No change applied.")
                    st.write("**Reference Sequence:**"); st.code(ref_seq_cleaned, language='text')
                else:
                    # Apply the substitution
                    alt_seq_list = list(ref_seq_cleaned)
                    alt_seq_list[zero_based_pos] = variant_base
                    alt_seq = "".join(alt_seq_list)

                    st.subheader("Variant Applied Successfully")
                    st.write(f"**Variant:** Position {variant_position}")
                    st.write(f"**Change:** Reference (`{original_base}`) → Altered (`{variant_base}`)")
                    st.write("**Reference Sequence:**"); st.code(ref_seq_cleaned, language='text')
                    st.write("**Altered Sequence:**"); st.code(alt_seq, language='text')
                    st.markdown("---")

                    # --- Codon Context Analysis ---
                    st.write("**Codon Context & Predicted Effect (Standard Genetic Code):**")
                    # Determine the start of the codon containing the variant (0-based)
                    codon_start_pos_0based = (zero_based_pos // 3) * 3

                    # Check if the full codon is within the sequence boundaries
                    if codon_start_pos_0based + 3 <= len(ref_seq_cleaned):
                        ref_codon = ref_seq_cleaned[codon_start_pos_0based : codon_start_pos_0based + 3]
                        alt_codon = alt_seq[codon_start_pos_0based : codon_start_pos_0based + 3]
                        try:
                            # Use Biopython's CodonTable for translation
                            standard_table = CodonTable.unambiguous_dna_by_id[1]
                            ref_aa = standard_table.forward_table.get(ref_codon) # Returns AA or None
                            if ref_codon in standard_table.stop_codons: ref_aa = '*'
                            if ref_aa is None: ref_aa = '?' # Handle ambiguous codons if they were allowed

                            alt_aa = standard_table.forward_table.get(alt_codon)
                            if alt_codon in standard_table.stop_codons: alt_aa = '*'
                            if alt_aa is None: alt_aa = '?'

                            st.write(f"- **Codon Number:** {(codon_start_pos_0based // 3) + 1}")
                            st.write(f"- **Position within Codon:** {(zero_based_pos % 3) + 1}")
                            st.write(f"- **Reference Codon:** `{ref_codon}` (Starts at base {codon_start_pos_0based + 1}, translates to `{ref_aa}`)")
                            st.write(f"- **Altered Codon:** `{alt_codon}` (Translates to `{alt_aa}`)")

                            # Determine effect
                            if ref_aa == alt_aa:
                                effect = "Silent (Synonymous)"
                                st.success(f"- **Predicted Effect:** **{effect}** (does not change amino acid)")
                            elif alt_aa == '*':
                                effect = "Nonsense"
                                st.error(f"- **Predicted Effect:** **{effect}** (introduces stop codon)")
                            elif ref_aa == '*':
                                effect = "Stop Loss"
                                st.warning(f"- **Predicted Effect:** **{effect}** (removes stop codon)")
                            else: # ref_aa != alt_aa and neither is '*'
                                effect = "Missense (Non-synonymous)"
                                st.warning(f"- **Predicted Effect:** **{effect}** (changes amino acid from {ref_aa} to {alt_aa})")

                        except KeyError as ke:
                            st.warning(f"Could not translate codon '{ke}'. Is it a standard DNA codon?")
                        except Exception as codon_e:
                            st.warning(f"Could not determine codon effect: {codon_e}")
                    else:
                        st.write("- Variant position is too close to the 3' end of the sequence to determine a full codon.")
            except Exception as e:
                st.error(f"An error occurred while applying the variant: {e}")

# --- Codon Usage Analyzer ---
elif menu == "Codon Usage Analyzer":
    st.header("Codon Usage Analyzer")
    st.markdown("Analyze the frequency of codons in a **coding DNA sequence (CDS)**. Sequence must contain only **A, T, C, G** and its length must be a multiple of 3.")
    dna_sequence = st.text_area("Enter Coding DNA Sequence (CDS):", height=150, key="codon_seq_input", placeholder="e.g., ATGCGT...")

    if st.button("Analyze Codon Usage", key="codon_analyze_button"):
        dna_sequence_cleaned = "".join(dna_sequence.split()).upper()
        # --- Input Validation ---
        if not dna_sequence_cleaned:
            st.warning("Please enter a coding DNA sequence.")
        else:
            valid_bases = {'A', 'T', 'C', 'G'}
            invalid_chars = set(c for c in dna_sequence_cleaned if c not in valid_bases)
            if invalid_chars:
                st.error(f"Sequence contains invalid characters: {', '.join(sorted(list(invalid_chars)))}. Use only A, T, C, G.")
            elif len(dna_sequence_cleaned) % 3 != 0:
                st.error(f"Sequence length ({len(dna_sequence_cleaned)} bp) is not a multiple of 3. Please provide a valid CDS (Coding DNA Sequence).")
            # --- End Validation ---
            else:
                try:
                    # Extract codons
                    codons = [dna_sequence_cleaned[i:i + 3] for i in range(0, len(dna_sequence_cleaned), 3)]
                    total_codons = len(codons)
                    codon_counts = Counter(codons)

                    st.subheader("Codon Usage Results")
                    st.write(f"**Total Number of Codons Analyzed:** {total_codons:,}")

                    # Get standard codon table information
                    try:
                        standard_table = CodonTable.unambiguous_dna_by_id[1]
                        all_codons_in_table = list(standard_table.forward_table.keys()) + standard_table.stop_codons
                    except Exception as table_e:
                        st.error(f"Could not load standard codon table from Biopython: {table_e}")
                        st.stop() # Stop if table isn't available

                    # Calculate frequencies and prepare DataFrame
                    codon_data = []
                    for codon in sorted(all_codons_in_table):
                        count = codon_counts.get(codon, 0)
                        freq_percent = (count / total_codons) * 100 if total_codons > 0 else 0
                        # Determine amino acid or stop
                        amino_acid = standard_table.forward_table.get(codon) # AA or None
                        if codon in standard_table.stop_codons: amino_acid = 'Stop'
                        if amino_acid is None: amino_acid = '?' # Should not happen for standard table

                        codon_data.append({"Codon": codon, "Amino Acid": amino_acid, "Count": count, "Frequency (%)": freq_percent})

                    usage_df = pd.DataFrame(codon_data)

                    # Display the main usage table
                    st.dataframe(usage_df.style.format({'Count': '{:,}', 'Frequency (%)': '{:.2f}%'}),
                                 hide_index=True, use_container_width=True)
                    st.write("---")

                    # --- Visualization Options ---
                    st.subheader("Visualize Codon Usage")
                    plot_choice = st.selectbox("Visualize usage for:", ["All Codons", "Specific Amino Acid"], key="codon_plot_choice")

                    if plot_choice == "All Codons":
                         # Plot frequencies for all codons *that are actually used* in the sequence
                         fig, ax = plt.subplots(figsize=(15, 6))
                         plot_df_filtered = usage_df[usage_df['Count'] > 0].sort_values(by='Codon') # Filter and sort for plot
                         if not plot_df_filtered.empty:
                             ax.bar(plot_df_filtered['Codon'], plot_df_filtered['Frequency (%)'], color='skyblue')
                             ax.set_xlabel("Codon")
                             ax.set_ylabel("Frequency (%)")
                             ax.set_title("Frequency of Used Codons in Sequence")
                             plt.xticks(rotation=90, fontsize=8)
                             plt.tight_layout()
                             st.pyplot(fig)
                         else:
                              st.info("No codons were found in the sequence to plot.")

                    else: # Specific Amino Acid
                         amino_acids = sorted(list(set(usage_df['Amino Acid']))) # Get unique AAs from the table
                         selected_aa = st.selectbox("Select Amino Acid (or Stop):", amino_acids, key="codon_aa_select")

                         # Filter data for the selected amino acid
                         aa_df = usage_df[(usage_df['Amino Acid'] == selected_aa) & (usage_df['Count'] > 0)].copy() # Also filter for count>0 here
                         aa_total_count = aa_df['Count'].sum() # Sum counts *only* for codons encoding this AA *present in the sequence*

                         if aa_total_count > 0:
                             # Calculate relative frequency *within that AA*
                             aa_df['Relative Freq (%)'] = (aa_df['Count'] / aa_total_count) * 100

                             # Plot relative frequencies
                             fig, ax = plt.subplots(figsize=(8, 5))
                             ax.bar(aa_df['Codon'], aa_df['Relative Freq (%)'], color='lightcoral')
                             ax.set_xlabel("Codon")
                             ax.set_ylabel("Relative Frequency (%) within AA")
                             ax.set_title(f"Codon Usage for {selected_aa} (Total Count in Seq: {aa_total_count:,})")
                             ax.set_ylim(0, 100) # Relative frequency is 0-100%
                             plt.tight_layout()
                             st.pyplot(fig)

                             # Show table for the selected AA
                             st.dataframe(aa_df[['Codon', 'Count', 'Relative Freq (%)']].style.format({'Count': '{:,}', 'Relative Freq (%)': '{:.2f}%'}),
                                          hide_index=True)
                         else:
                             st.info(f"The amino acid '{selected_aa}' (or its codons) was not found in the provided sequence.")
                except Exception as e:
                    st.error(f"An error occurred during codon usage analysis: {e}")


# --- Footer ---
# (Footer definition remains the same)
st.markdown('<div class="footer">© 2024 Colon Cancer Toolbox | For Educational & Informational Purposes Only | Not for Medical Diagnosis</div>', unsafe_allow_html=True)
