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

# --- Page Configuration (MUST be the first Streamlit command) ---
st.set_page_config(
    page_title="Colon Cancer Toolbox",
    layout="wide",
    initial_sidebar_state="expanded",
)

# --- Custom CSS Styling ---
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
# (Sidebar definition remains the same)
menu = st.sidebar.selectbox(
    "Choose a Feature:",
    [
        "Home",
        "Colon Cancer Database Search",
        "DNA Sequence Analysis",
        "Protein Sequence Analysis",
        "Primer Design",
        "Motif Finder Tool",
        "Bioinformatics Tool (Transcription/Translation)",
        "Genome Coverage Plotter",
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

# (Home page content remains the same)
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
        *   **Primer Design:** Generate basic forward and reverse primers from a DNA template.
        *   **Motif Finder Tool:** Locate specific sequence patterns (motifs) within DNA or protein sequences using exact matches or regular expressions.
        *   **Bioinformatics Tool (Transcription/Translation):** Transcribe DNA to RNA and translate RNA to protein based on the standard genetic code.
        *   **Genome Coverage Plotter:** Visualize sequencing coverage depth across genomic positions from an uploaded file.
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
            print(f"Feedback received: {feedback}")
        else:
            st.warning("Feedback cannot be empty.")

# (Colon Cancer Database Search content remains the same)
elif menu == "Colon Cancer Database Search":
    st.header("Search Colon Cancer Gene Database")
    st.markdown("Enter a gene name (e.g., `APC`, `KRAS`, `TP53`, `MLH1`) to search our curated information.")

    data_path = "Colon Cancer.xlsx"
    query = st.text_input("Enter Gene Name to Search:", key="gene_search_input", placeholder="e.g., KRAS")

    if st.button("Search Database", key="gene_search_button"):
        if not query:
            st.warning("Please enter a gene name to search.")
        else:
            try:
                data = pd.read_excel(data_path)
                data.columns = data.columns.str.strip().str.title()
                gene_col = None
                possible_gene_cols = ['Gene Name', 'Gene', 'Symbol', 'Gene Symbol']
                for col in possible_gene_cols:
                    if col in data.columns:
                        gene_col = col
                        break
                if not gene_col:
                     st.error(f"Error: Could not find a 'Gene Name' column (or similar: {possible_gene_cols}) in the Excel file. Columns found: {list(data.columns)}")
                     st.stop()

                data[gene_col] = data[gene_col].astype(str).str.strip().str.lower()
                search_query = query.strip().lower()
                results = data[data[gene_col].str.contains(search_query, na=False, case=False)].copy()

                if not results.empty:
                    st.success(f"Found {len(results)} result(s) for '{query}'.")
                    st.write("### Search Results:")
                    formatted_results = results.copy()
                    for col in formatted_results.columns:
                         try:
                             formatted_results[col] = formatted_results[col].apply(lambda x: format_link(x, col))
                         except Exception as apply_e:
                             st.warning(f"Could not apply link formatting to column '{col}': {apply_e}. Displaying raw data.")

                    html_table = formatted_results.to_html(
                        escape=False,
                        index=False,
                        na_rep='-',
                        justify='left',
                        classes=['dataframe', 'st-table', 'results-table']
                    )
                    st.write(html_table, unsafe_allow_html=True)

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


# (DNA Sequence Analysis content remains the same)
elif menu == "DNA Sequence Analysis":
    st.header("DNA Sequence Analysis Tool")
    st.markdown("Enter a DNA sequence containing **A, T, C, G, N** characters for basic analysis.")
    sequence = st.text_area("Enter DNA Sequence:", height=150, key="dna_seq_input", placeholder="e.g., ATGCGTAGCTAGCATGCNNNNATGC")

    if st.button("Analyze DNA Sequence", key="dna_analyze_button"):
        if sequence:
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
                    length = len(sequence_cleaned)
                    st.write(f"**Sequence Length:** {length} bp")
                    seq_no_n = sequence_cleaned.replace('N', '')
                    if len(seq_no_n) > 0:
                        gc = gc_fraction(seq_no_n) * 100
                        st.write(f"**GC Content (excluding Ns):** {gc:.2f}%")
                    else:
                        st.write("**GC Content:** N/A (Sequence contains only 'N' or is empty after removing 'N')")

                    st.write("**Base Composition:**")
                    composition = Counter(sequence_cleaned)
                    comp_data = [{"Base": base, "Count": count, "Percentage": (count/length*100) if length > 0 else 0}
                                 for base, count in sorted(composition.items())]
                    st.dataframe(pd.DataFrame(comp_data).style.format({"Count":"{:,}","Percentage": "{:.2f}%"}), hide_index=True, use_container_width=True)

                    plot_composition = {k: v for k, v in composition.items() if k != 'N'}
                    if plot_composition:
                        fig, ax = plt.subplots()
                        colors = {'A': '#1f77b4', 'T': '#ff7f0e', 'C': '#2ca02c', 'G': '#d62728'}
                        bar_bases = sorted(plot_composition.keys())
                        bar_counts = [plot_composition[b] for b in bar_bases]
                        bar_colors = [colors.get(base, '#8c564b') for base in bar_bases]
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


# (Protein Sequence Analysis content remains the same)
elif menu == "Protein Sequence Analysis":
    st.header("Protein Sequence Analysis Tool")
    st.markdown("Enter a protein sequence using standard single-letter amino acid codes (A-Y). Non-standard characters (e.g., *, X, B, Z) will be ignored in calculations requiring standard AAs.")
    protein_sequence = st.text_area("Enter Protein Sequence:", height=150, key="protein_seq_input", placeholder="e.g., MKTAYIA...")

    if st.button("Analyze Protein Sequence", key="protein_analyze_button"):
        protein_sequence_cleaned = "".join(protein_sequence.split()).upper()
        if not protein_sequence_cleaned:
            st.warning("Please enter a protein sequence.")
        else:
            standard_aa = "ACDEFGHIKLMNPQRSTVWY"
            protein_for_analysis = re.sub(f'[^{standard_aa}]', '', protein_sequence_cleaned)
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

                    st.subheader("Amino Acid Composition (Based on Original Input)")
                    full_composition = Counter(protein_sequence_cleaned)
                    comp_data = [{"Amino Acid": aa, "Count": count, "Percentage": (count/len(protein_sequence_cleaned)*100)}
                                 for aa, count in sorted(full_composition.items())]
                    st.dataframe(pd.DataFrame(comp_data).style.format({'Count': '{:,}', 'Percentage': '{:.2f}%'}), hide_index=True, use_container_width=True)

                    st.subheader("Secondary Structure Fraction (Predicted - Standard AAs only)")
                    try:
                        helix, turn, sheet = analyzed_seq.secondary_structure_fraction()
                        st.write(f"- **Helix:** {helix*100:.2f}%")
                        st.write(f"- **Turn:** {turn*100:.2f}%")
                        st.write(f"- **Sheet:** {sheet*100:.2f}%")
                    except Exception as ss_e:
                        st.warning(f"Could not calculate secondary structure fractions: {ss_e}")

                except Exception as e:
                    st.error(f"An error occurred during protein analysis: {e}")

# (Primer Design content remains the same)
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
            valid_bases = {'A', 'T', 'C', 'G'}
            invalid_chars = set(c for c in sequence_cleaned if c not in valid_bases)
            if invalid_chars:
                st.error(f"Invalid characters found: {', '.join(sorted(list(invalid_chars)))}. Please use only A, T, C, G.")
            elif len(sequence_cleaned) < primer_length * 2:
                st.error(f"Sequence is too short ({len(sequence_cleaned)} bp). Needs to be at least twice the primer length ({primer_length * 2} bp).")
            else:
                try:
                    forward_primer_seq = Seq(sequence_cleaned[:primer_length])
                    reverse_primer_template = Seq(sequence_cleaned[-primer_length:])
                    reverse_primer_seq = reverse_primer_template.reverse_complement()

                    st.subheader("Suggested Primers")
                    col1, col2 = st.columns(2)
                    with col1:
                        st.markdown("**Forward Primer (5' → 3'):**")
                        st.code(str(forward_primer_seq), language='text')
                        st.write(f"- Length: {len(forward_primer_seq)} bp")
                        st.write(f"- GC Content: {gc_fraction(forward_primer_seq)*100:.2f}%")
                    with col2:
                        st.markdown("**Reverse Primer (5' → 3'):**")
                        st.code(str(reverse_primer_seq), language='text')
                        st.write(f"- Length: {len(reverse_primer_seq)} bp")
                        st.write(f"- GC Content: {gc_fraction(reverse_primer_seq)*100:.2f}%")

                    amplicon_length = len(sequence_cleaned)
                    st.write(f"**Expected Amplicon Size:** {amplicon_length} bp (entire template)")

                except Exception as e:
                    st.error(f"An error occurred during primer design: {e}")

# (Bioinformatics Tool (Transcription/Translation) content remains the same)
elif menu == "Bioinformatics Tool (Transcription/Translation)":
    st.header("DNA Transcription & Translation Tool")
    st.markdown("Enter a DNA sequence (coding strand, **A, T, C, G** only) to transcribe into RNA and translate into protein using the standard genetic code.")
    sequence = st.text_area("Enter DNA Sequence (Coding Strand):", height=150, key="bioinfo_seq_input", placeholder="e.g., ATGCGTAGCTAGCATGC")

    if st.button("Transcribe and Translate", key="bioinfo_submit_button"):
        sequence_cleaned = "".join(sequence.split()).upper()
        if not sequence_cleaned:
            st.warning("Please enter a DNA sequence.")
        else:
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
                    st.markdown("**DNA Reverse Complement (5' → 3'):**")
                    st.code(str(seq_obj.reverse_complement()), language='text')
                    st.markdown("---")
                    st.markdown("**Transcription (DNA → RNA, 5' → 3'):**")
                    rna_transcript = seq_obj.transcribe()
                    st.code(str(rna_transcript), language='text')
                    st.markdown("---")
                    st.markdown("**Translation (RNA → Protein):**")
                    st.caption("Uses standard genetic code. Translation stops at the first encountered stop codon ('*').")

                    remainder = len(rna_transcript) % 3
                    if remainder != 0:
                        st.warning(f"RNA sequence length ({len(rna_transcript)} bases) is not a multiple of 3. The last {remainder} base(s) will be ignored during translation.")
                        rna_transcript_trimmed = rna_transcript[:-remainder]
                    else:
                        rna_transcript_trimmed = rna_transcript

                    if len(rna_transcript_trimmed) >= 3:
                         protein_translation = rna_transcript_trimmed.translate(to_stop=True)
                         st.code(str(protein_translation), language='text')
                         st.write(f"- **Translated Protein Length:** {len(protein_translation)} aa")
                         protein_full_translation = rna_transcript_trimmed.translate(to_stop=False)
                         if '*' in protein_full_translation and protein_full_translation.endswith('*') and len(protein_full_translation) == len(protein_translation)+1:
                            pass
                         elif '*' in protein_full_translation:
                            st.caption("**Full Translation (including internal stop codons):**")
                            st.code(str(protein_full_translation), language='text')

                    elif len(rna_transcript) < 3:
                         st.error("Cannot translate: RNA sequence is shorter than 3 bases.")
                    else:
                         st.error("Cannot translate: RNA sequence length is zero after trimming to a multiple of 3.")

                except Exception as e:
                    st.error(f"An error occurred during processing: {e}")

# (Motif Finder Tool content remains the same)
elif menu == "Motif Finder Tool":
    st.header("Motif Finder Tool")
    st.markdown("Find occurrences of a motif (exact string or simple regular expression) within a DNA or Protein sequence.")
    st.info("Supports basic regex syntax like `.` (any char), `[AT]` (A or T), `*` (0 or more), `+` (1 or more), `?` (0 or 1), `{n,m}` (n to m repeats). Use `\` to escape special characters.")

    sequence = st.text_area("Enter Sequence (DNA or Protein):", height=150, key="motif_seq_input")
    motif = st.text_input("Enter Motif / Regex Pattern:", key="motif_pattern_input", placeholder="e.g., ATG or Y[DN]R or ^M.*K$")
    seq_type = st.radio("Sequence Type (for context only):", ("DNA", "Protein"), horizontal=True, key="motif_seq_type")
    case_sensitive = st.checkbox("Case Sensitive Search", value=False, key="motif_case_sensitive")

    if st.button("Find Motif Occurrences", key="motif_find_button"):
        if not sequence or not motif:
            st.warning("Please provide both the sequence and the motif/pattern to search for.")
        else:
            sequence_cleaned = "".join(sequence.split())
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
                              "Position (End)": m.end(),
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

# (Genome Coverage Plotter content remains the same)
elif menu == "Genome Coverage Plotter":
    st.header("Genome Coverage Plotter")
    st.markdown("Upload a coverage data file (e.g., from `samtools depth`) containing genomic positions and coverage values. The tool expects columns like 'position' and 'coverage'.")
    uploaded_file = st.file_uploader("Upload Coverage Data (CSV, TSV, TXT):", type=['csv', 'tsv', 'txt'], key="coverage_uploader")

    if uploaded_file:
        st.write("---")
        st.subheader("File Settings")
        col1, col2 = st.columns(2)
        with col1:
            sep_options = {',': 'Comma (,)', '\t': 'Tab (\\t)', ' ': 'Space ( )'}
            selected_sep_display = st.selectbox("Column Separator:", options=list(sep_options.keys()), index=1, format_func=lambda x: sep_options[x], key="cov_sep")
            separator = selected_sep_display if selected_sep_display != ' ' else r'\s+'
        with col2:
            header_row = st.number_input("Header Row Index (0 if first row, None if no header):", min_value=0, value=0, step=1, key="cov_header", help="Enter the row number containing column names (0-based). If no header, plotting will use column indices.")
            use_header = st.checkbox("File has a header row", value=True, key="cov_has_header")
            header_arg = header_row if use_header else None

        st.write("---")
        st.subheader("Data Preview & Column Selection")

        try:
            df = pd.read_csv(uploaded_file, sep=separator, header=header_arg, engine='python', comment='#')
            st.write("Preview of uploaded data (first 5 rows):")
            st.dataframe(df.head(), height=200)

            if use_header:
                available_cols = list(df.columns)
                normalized_cols_map = {col.strip().lower(): col for col in available_cols}
                pos_col_guess, cov_col_guess = None, None
                for p_name in ['position', 'pos', 'coordinate', 'location', 'start', 'chromstart', 'base', 'pos.', 'position(bp)']:
                    if p_name in normalized_cols_map: pos_col_guess = normalized_cols_map[p_name]; break
                for c_name in ['coverage', 'cov', 'depth', 'reads', 'count', 'value']:
                    if c_name in normalized_cols_map: cov_col_guess = normalized_cols_map[c_name]; break
            else:
                available_cols = list(range(df.shape[1]))
                st.info("No header specified. Please select columns by their index (0, 1, 2...).")
                pos_col_guess = 1 if len(available_cols) > 1 else 0
                cov_col_guess = 2 if len(available_cols) > 2 else (1 if len(available_cols) > 1 else 0)

            col_sel1, col_sel2 = st.columns(2)
            with col_sel1: pos_col_selected = st.selectbox("Select Position Column:", available_cols, index=available_cols.index(pos_col_guess) if pos_col_guess in available_cols else 0, key="cov_pos_col")
            with col_sel2: cov_col_selected = st.selectbox("Select Coverage Column:", available_cols, index=available_cols.index(cov_col_guess) if cov_col_guess in available_cols else (1 if len(available_cols)>1 else 0), key="cov_cov_col")

            if pos_col_selected == cov_col_selected:
                 st.error("Position and Coverage columns cannot be the same.")
            else:
                st.success(f"Selected: Position='{pos_col_selected}', Coverage='{cov_col_selected}'.")
                st.write("---")
                st.subheader("Plotting Options & Results")
                try:
                    plot_df = df[[pos_col_selected, cov_col_selected]].copy()
                    plot_df.columns = ['Position', 'Coverage']
                    plot_df['Position'] = pd.to_numeric(plot_df['Position'], errors='coerce')
                    plot_df['Coverage'] = pd.to_numeric(plot_df['Coverage'], errors='coerce')
                    rows_before = len(plot_df)
                    plot_df.dropna(inplace=True)
                    rows_after = len(plot_df)
                    if rows_after < rows_before: st.warning(f"Removed {rows_before - rows_after} rows containing non-numeric data in selected columns.")

                    if plot_df.empty: st.error("No valid numeric data found in the selected columns after cleaning.")
                    else:
                        plot_df = plot_df.sort_values(by='Position')
                        fig, ax = plt.subplots(figsize=(12, 5))
                        ax.plot(plot_df['Position'], plot_df['Coverage'], label='Coverage', lw=1.0, color='steelblue')
                        fill_plot = st.checkbox("Fill area under curve", value=True, key="cov_fill")
                        if fill_plot: ax.fill_between(plot_df['Position'], plot_df['Coverage'], alpha=0.3, color='lightblue')
                        smoothing_window = st.slider("Smoothing Window Size (0 for none):", 0, 100, 0, 5, key="cov_smooth", help="Applies a rolling mean. 0 = no smoothing.")
                        if smoothing_window > 1:
                            plot_df['Smoothed Coverage'] = plot_df['Coverage'].rolling(window=smoothing_window, center=True, min_periods=1).mean()
                            ax.plot(plot_df['Position'], plot_df['Smoothed Coverage'], label=f'Smoothed ({smoothing_window}bp window)', lw=1.5, color='darkorange', linestyle='--')
                        ax.set_xlabel(f"Genomic Position ({pos_col_selected})")
                        ax.set_ylabel(f"Coverage Depth ({cov_col_selected})")
                        ax.set_title("Genome Coverage Plot")
                        ax.legend()
                        ax.grid(True, linestyle='--', alpha=0.6)
                        ax.ticklabel_format(style='plain', axis='x', useOffset=False)
                        plt.tight_layout()
                        st.pyplot(fig)

                        st.subheader("Coverage Statistics")
                        stats_df = plot_df['Coverage'].describe().to_frame().T
                        stats_df = stats_df.rename(columns={'count': 'Positions Covered', 'mean': 'Mean Coverage', 'std': 'Std Dev', 'min': 'Min Coverage', '25%': '25th Pctl', '50%': 'Median (50th Pctl)', '75%': '75th Pctl', 'max': 'Max Coverage'})
                        st.dataframe(stats_df.style.format("{:,.2f}", na_rep="-"))
                        min_cov_threshold = st.number_input("Calculate % positions with coverage >=", 0, value=10, key="cov_breadth_thresh")
                        if min_cov_threshold >= 0:
                             positions_above_thresh = (plot_df['Coverage'] >= min_cov_threshold).sum()
                             total_positions_in_range = plot_df['Position'].nunique()
                             if total_positions_in_range > 0:
                                 breadth_pct = (positions_above_thresh / total_positions_in_range) * 100
                                 st.metric(label=f"Breadth of Coverage (>= {min_cov_threshold}x)", value=f"{breadth_pct:.2f}%", delta=f"{positions_above_thresh:,} / {total_positions_in_range:,} positions")
                             else:
                                 st.metric(label=f"Breadth of Coverage (>= {min_cov_threshold}x)", value="N/A", delta="No positions in range")
                except Exception as plot_e: st.error(f"An error occurred during data processing or plotting: {plot_e}\nPlease ensure the correct columns are selected and contain numeric data.")

        except pd.errors.EmptyDataError: st.error("File Reading Error: The uploaded file is empty or contains no data.")
        except ValueError as ve: st.error(f"File Reading Error: Possible issue with separator or header row setting. Details: {ve}")
        except Exception as read_e: st.error(f"File Reading Error: Could not parse the file. Please check the format, separator, and header settings. Error: {read_e}")

# (Variant Annotation Tool content remains the same)
elif menu == "Variant Annotation Tool":
    st.header("Simple Variant Substitution Tool")
    st.markdown("Introduce a single nucleotide substitution into a reference DNA sequence (**A, T, C, G** only).")
    ref_seq = st.text_area("Reference DNA Sequence:", height=100, key="var_ref_seq", placeholder="e.g., ATGCGTACGTAGCTAG")
    col1, col2, col3 = st.columns(3)
    with col1: variant_position = st.number_input("Variant Position (1-based index):", min_value=1, step=1, key="var_pos", help="The position in the sequence to change (starting from 1).")
    with col2: variant_base = st.selectbox("New Variant Base (ALT):", ["A", "T", "C", "G"], index=0, key="var_base", help="The new base to insert at the specified position.")

    if st.button("Apply Variant Substitution", key="var_apply_button"):
        ref_seq_cleaned = "".join(ref_seq.split()).upper()
        error_flag = False
        if not ref_seq_cleaned: st.warning("Please enter a reference DNA sequence."); error_flag = True
        if not variant_base: st.warning("Please select a new variant base."); error_flag = True
        valid_bases = {'A', 'T', 'C', 'G'}
        invalid_ref_chars = set(c for c in ref_seq_cleaned if c not in valid_bases)
        if invalid_ref_chars: st.error(f"Reference sequence contains invalid characters: {', '.join(sorted(list(invalid_ref_chars)))}. Use only A, T, C, G."); error_flag = True
        if not error_flag and (variant_position <= 0 or variant_position > len(ref_seq_cleaned)): st.error(f"Variant position ({variant_position}) is out of range for the sequence length ({len(ref_seq_cleaned)})."); error_flag = True

        if not error_flag:
            try:
                zero_based_pos = variant_position - 1
                original_base = ref_seq_cleaned[zero_based_pos]
                if original_base == variant_base:
                    st.info(f"The base at position {variant_position} is already '{original_base}'. No change applied.")
                    st.write("**Reference Sequence:**"); st.code(ref_seq_cleaned, language='text')
                else:
                    alt_seq_list = list(ref_seq_cleaned); alt_seq_list[zero_based_pos] = variant_base; alt_seq = "".join(alt_seq_list)
                    st.subheader("Variant Applied Successfully")
                    st.write(f"**Variant:** Position {variant_position}")
                    st.write(f"**Change:** Reference (`{original_base}`) → Altered (`{variant_base}`)")
                    st.write("**Reference Sequence:**"); st.code(ref_seq_cleaned, language='text')
                    st.write("**Altered Sequence:**"); st.code(alt_seq, language='text')
                    st.markdown("---")
                    st.write("**Codon Context (if applicable):**")
                    codon_start_pos_0based = (zero_based_pos // 3) * 3
                    if codon_start_pos_0based + 3 <= len(ref_seq_cleaned):
                        ref_codon = ref_seq_cleaned[codon_start_pos_0based : codon_start_pos_0based + 3]
                        alt_codon = alt_seq[codon_start_pos_0based : codon_start_pos_0based + 3]
                        try:
                            from Bio.Seq import Seq; from Bio.Data import CodonTable
                            standard_table = CodonTable.unambiguous_dna_by_id[1]
                            ref_aa = standard_table.forward_table.get(ref_codon, '*' if ref_codon in standard_table.stop_codons else '?')
                            alt_aa = standard_table.forward_table.get(alt_codon, '*' if alt_codon in standard_table.stop_codons else '?')
                            st.write(f"- **Reference Codon:** `{ref_codon}` (Starts at base {codon_start_pos_0based + 1}, translates to `{ref_aa}`)")
                            st.write(f"- **Altered Codon:** `{alt_codon}` (Translates to `{alt_aa}`)")
                            if ref_aa != alt_aa:
                                if alt_aa == '*': st.write("- **Effect:** **Nonsense** mutation (introduces stop codon)")
                                else: st.write("- **Effect:** **Missense** mutation (changes amino acid)")
                            else: st.write("- **Effect:** **Silent (Synonymous)** mutation (does not change amino acid)")
                        except ImportError: st.warning("Could not determine codon effect: Biopython library not fully available.")
                        except Exception as codon_e: st.warning(f"Could not determine codon effect: {codon_e}")
                    else: st.write("- Variant position is too close to the end of the sequence to determine a full codon.")
            except Exception as e: st.error(f"An error occurred while applying the variant: {e}")

# (Codon Usage Analyzer content remains the same)
elif menu == "Codon Usage Analyzer":
    st.header("Codon Usage Analyzer")
    st.markdown("Analyze the frequency of codons in a **coding DNA sequence (CDS)**. Sequence must contain only **A, T, C, G** and its length must be a multiple of 3.")
    dna_sequence = st.text_area("Enter Coding DNA Sequence (CDS):", height=150, key="codon_seq_input", placeholder="e.g., ATGCGT...")

    if st.button("Analyze Codon Usage", key="codon_analyze_button"):
        dna_sequence_cleaned = "".join(dna_sequence.split()).upper()
        if not dna_sequence_cleaned: st.warning("Please enter a coding DNA sequence.")
        else:
            valid_bases = {'A', 'T', 'C', 'G'}
            invalid_chars = set(c for c in dna_sequence_cleaned if c not in valid_bases)
            if invalid_chars: st.error(f"Sequence contains invalid characters: {', '.join(sorted(list(invalid_chars)))}. Use only A, T, C, G.")
            elif len(dna_sequence_cleaned) % 3 != 0: st.error(f"Sequence length ({len(dna_sequence_cleaned)} bp) is not a multiple of 3. Please provide a valid CDS.")
            else:
                try:
                    codons = [dna_sequence_cleaned[i:i + 3] for i in range(0, len(dna_sequence_cleaned), 3)]
                    total_codons = len(codons); codon_counts = Counter(codons)
                    st.subheader("Codon Usage Results"); st.write(f"**Total Number of Codons Analyzed:** {total_codons:,}")
                    from Bio.Data import CodonTable
                    try: standard_table = CodonTable.unambiguous_dna_by_id[1]; all_codons_in_table = list(standard_table.forward_table.keys()) + standard_table.stop_codons
                    except ImportError: st.error("Biopython library needed for codon table. Please install it (`pip install biopython`)."); st.stop()
                    except Exception as table_e: st.error(f"Could not load standard codon table: {table_e}"); st.stop()

                    codon_data = []
                    for codon in sorted(all_codons_in_table):
                        count = codon_counts.get(codon, 0); freq_percent = (count / total_codons) * 100 if total_codons > 0 else 0
                        amino_acid = standard_table.forward_table.get(codon, 'Stop' if codon in standard_table.stop_codons else '?')
                        codon_data.append({"Codon": codon, "Amino Acid": amino_acid, "Count": count, "Frequency (%)": freq_percent})
                    usage_df = pd.DataFrame(codon_data)
                    st.dataframe(usage_df.style.format({'Count': '{:,}', 'Frequency (%)': '{:.2f}%'}), hide_index=True, use_container_width=True)

                    plot_choice = st.selectbox("Visualize usage for:", ["All Codons", "Specific Amino Acid"], key="codon_plot_choice")
                    if plot_choice == "All Codons":
                         fig, ax = plt.subplots(figsize=(15, 6))
                         plot_df_filtered = usage_df[usage_df['Count'] > 0]
                         ax.bar(plot_df_filtered['Codon'], plot_df_filtered['Frequency (%)'], color='skyblue')
                         ax.set_xlabel("Codon"); ax.set_ylabel("Frequency (%)"); ax.set_title("Frequency of Used Codons")
                         plt.xticks(rotation=90, fontsize=8); plt.tight_layout(); st.pyplot(fig)
                    else:
                         amino_acids = sorted(list(set(usage_df['Amino Acid'])))
                         selected_aa = st.selectbox("Select Amino Acid (or Stop):", amino_acids, key="codon_aa_select")
                         aa_df = usage_df[usage_df['Amino Acid'] == selected_aa].copy()
                         aa_total_count = aa_df['Count'].sum()
                         if aa_total_count > 0:
                             aa_df['Relative Freq (%)'] = (aa_df['Count'] / aa_total_count) * 100
                             fig, ax = plt.subplots(figsize=(8, 5))
                             ax.bar(aa_df['Codon'], aa_df['Relative Freq (%)'], color='lightcoral')
                             ax.set_xlabel("Codon"); ax.set_ylabel("Relative Frequency (%) within AA")
                             ax.set_title(f"Codon Usage for {selected_aa} (Total Count: {aa_total_count:,})")
                             ax.set_ylim(0, 100); plt.tight_layout(); st.pyplot(fig)
                             st.dataframe(aa_df[['Codon', 'Count', 'Relative Freq (%)']].style.format({'Count': '{:,}', 'Relative Freq (%)': '{:.2f}%'}), hide_index=True)
                         else: st.info(f"The amino acid '{selected_aa}' was not found in the translated sequence.")
                except Exception as e: st.error(f"An error occurred during codon usage analysis: {e}")


# --- Footer ---
# (Footer definition remains the same)
st.markdown('<div class="footer">© 2024 Colon Cancer Toolbox | For Educational & Informational Purposes Only | Not for Medical Diagnosis</div>', unsafe_allow_html=True)
