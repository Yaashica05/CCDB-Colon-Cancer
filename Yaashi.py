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
# (Includes styles for table links, fonts, footer, etc. - kept from previous version)
st.markdown(
    """
    <style>
    body {
        background-color: #f0f8ff; /* AliceBlue */
        font-family: 'Times New Roman', Times, serif;
    }
    .stApp {
        background-color: #f0f8ff;
        padding: 15px;
        border-radius: 10px;
    }
    h1, h2, h3, h4, h5, h6, p, div, span, label, .stButton>button, .stTextInput>div>div>input, .stTextArea>div>textarea, .stSelectbox>div>div>div {
        font-family: 'Times New Roman', Times, serif;
        font-weight: bold;
    }
    .stTextInput>label, .stTextArea>label, .stSelectbox>label, .stNumberInput>label, .stFileUploader>label, .stSlider>label, .stRadio>label {
        font-weight: bold;
        font-family: 'Times New Roman', Times, serif;
    }
    table a { /* Style for links inside the generated table */
        color: #0000FF; /* Blue */
        text-decoration: underline;
        font-weight: normal;
    }
    table a:hover {
        color: #00008B; /* DarkBlue on hover */
    }
    table td {
        font-weight: normal;
        font-family: 'Times New Roman', Times, serif;
        vertical-align: top;
        text-align: left; /* Ensure text aligns left */
        padding: 4px 8px; /* Add some padding */
    }
    table th {
        font-weight: bold;
        font-family: 'Times New Roman', Times, serif;
        text-align: left;
        padding: 4px 8px;
    }
    .footer {
        position: fixed;
        left: 0;
        bottom: 0;
        width: 100%;
        background-color: #e6e6fa; /* Lavender */
        color: black;
        text-align: center;
        padding: 5px;
        font-size: 12px;
        font-family: 'Times New Roman', Times, serif;
        font-weight: normal;
        border-top: 1px solid #d3d3d3;
        z-index: 1000;
    }
     .stApp > div:first-child {
        padding-bottom: 50px; /* Avoid overlap with footer */
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# --- Main App Title ---
st.title("Colon Cancer Database and Bioinformatics Toolbox")

# --- Sidebar Navigation ---
menu = st.sidebar.selectbox(
    "Choose a Feature:",
    [
        "Home",
        "Colon Cancer Database Search",
        "DNA Sequence Analysis",
        "Primer Design",
        "Bioinformatics Tool",
        "Protein Sequence Analysis",
        "Motif Finder Tool",
        "Genome Coverage Plotter",
        "Variant Annotation Tool",
        "Codon Usage Analyzer",
    ],
    help="Select a tool or section from the dropdown list."
)

# --- Helper function to create links for the Database Search ---
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
    # These names MUST match the column names *after* .str.strip().str.title()
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
         # Basic check for typical Uniprot ID format (can be complex)
         # or if it's already a URL
         if re.match(r'^[A-Z0-9_.-]+$', s_value, re.IGNORECASE) and not s_value.startswith("http"):
             # Check length constraints if needed, e.g., usually 6 or 10 chars for accession
             # if len(s_value) == 6 or len(s_value) == 10: # Example check
             url = f"https://www.uniprot.org/uniprotkb/{s_value}/entry"
             return f'<a href="{url}" target="_blank">{s_value}</a>'
         elif s_value.startswith("http"):
             return f'<a href="{s_value}" target="_blank">{s_value}</a>'
         else:
             return s_value

    # --- Formatting for BLAST, Conserved Domain (Make links ONLY if they are URLs) ---
    elif col_name_std in ['Blast', 'Conserved Domain']:
         if s_value.startswith("http"):
              return f'<a href="{s_value}" target="_blank">{s_value}</a>'
         else:
              # Return as plain text if not a URL (e.g., 'CD', 'Pfam:XYZ', 'Yes')
              return s_value

    # --- Default: Handle generic http links in ANY other column ---
    else:
        # This catches URLs in columns not explicitly handled above
        if s_value.startswith("http"):
            return f'<a href="{s_value}" target="_blank">{s_value}</a>'
        else:
            # Otherwise, just return the value as plain text
            return s_value

# --- Page Content based on Menu Selection ---

# Home page
if menu == "Home":
    st.header("Welcome to the Colon Cancer Data Resource")
    st.markdown(
        """
        This application serves as a central hub for exploring colon cancer-related gene information and utilizing various bioinformatics tools for sequence analysis.

        **Colon Cancer Overview:**
        Colon cancer, also known as colorectal cancer, originates in the large intestine (colon) or the rectum. It often begins as small, noncancerous clumps of cells called polyps. Early detection through screening is crucial.

        ### **Key Information:**
        - **Risk Factors**: Age, history, diet, lifestyle, genetics.
        - **Symptoms**: Changes in bowel habits, bleeding, abdominal pain, weight loss, fatigue. *Early stages may have no symptoms.*
        - **Prevention & Screening**: Healthy lifestyle, regular screening (e.g., colonoscopy, FIT).
        - **Treatment Options**: Surgery, chemotherapy, radiation, targeted therapy, immunotherapy.

        Use the sidebar menu to navigate the tools.
        ---
        *Disclaimer: For informational purposes only. Consult healthcare professionals for medical advice.*
        """
    )
    st.image(
        "https://img.freepik.com/free-vector/colorectal-cancer-crc-infographic-education_1308-47971.jpg",
        caption="Colorectal Cancer Infographic (Source: Freepik)",
        use_container_width=True
    )
    st.markdown("[Learn more about Colorectal Cancer on Wikipedia](https://en.wikipedia.org/wiki/Colorectal_cancer)", unsafe_allow_html=True)

    st.subheader("Feedback")
    feedback = st.text_area("Share your feedback about this application:")
    if st.button("Submit Feedback"):
        if feedback:
            st.success("Thank you for your feedback!")
            print(f"Feedback received: {feedback}")
        else:
            st.warning("Feedback cannot be empty.")

# Colon Cancer Database Search
elif menu == "Colon Cancer Database Search":
    st.header("Search Colon Cancer Gene Database")
    st.markdown("Enter a gene name (e.g., APC, KRAS, TP53) to search the database.")

    # !!! IMPORTANT: Update this path to your actual Excel file location !!!
    data_path = "Colon Cancer.xlsx"  # Assumes file is in the same folder as the script

    query = st.text_input("Enter Gene Name to Search:", key="gene_search_input")

    if st.button("Search", key="gene_search_button"):
        if not query:
            st.warning("Please enter a gene name to search.")
        else:
            try:
                # Load data
                data = pd.read_excel(data_path)

                # --- Data Cleaning and Preparation ---
                # 1. Standardize column names
                data.columns = data.columns.str.strip().str.title()

                # 2. Identify the Gene Name column
                gene_col = None
                possible_gene_cols = ['Gene Name', 'Gene', 'Symbol']
                for col in possible_gene_cols:
                    if col in data.columns:
                        gene_col = col
                        break
                if not gene_col:
                     st.error(f"Error: Could not find a 'Gene Name' column (or similar) in the Excel file. Standardized columns found: {list(data.columns)}")
                     st.stop()

                # 3. Clean Gene Name column for searching
                data[gene_col] = data[gene_col].astype(str).str.strip().str.lower()

                # --- Perform Search ---
                search_query = query.strip().lower()
                # Make a copy to avoid SettingWithCopyWarning when applying format_link
                results = data[data[gene_col].str.contains(search_query, na=False)].copy()

                if not results.empty:
                    st.success(f"Found {len(results)} result(s) for '{query}'.")
                    st.write("### Search Results:")

                    # --- Apply formatting using the helper function ---
                    # Iterate through ALL columns in the filtered results
                    for col in results.columns:
                         # Apply the format_link function to each cell in the column
                         # The function itself handles which columns/values get linked
                         try:
                             results[col] = results[col].apply(lambda x: format_link(x, col))
                         except Exception as apply_e:
                             st.warning(f"Could not apply formatting to column '{col}': {apply_e}")
                             # Continue to next column if one fails

                    # --- Convert the MODIFIED DataFrame (with HTML strings) to HTML ---
                    html_table = results.to_html(
                        escape=False,        # Crucial: Allow embedded HTML (<a> tags)
                        index=False,         # Hide DataFrame index
                        na_rep='',           # Show missing values as empty
                        justify='left',      # Align text left
                        classes=['dataframe', 'st-table'] # Apply CSS classes if needed
                    )

                    # --- Display the generated HTML table ---
                    st.write(html_table, unsafe_allow_html=True) # Crucial: Render the HTML

                else:
                    st.warning(f"No results found matching '{query}'. Please check the gene name or try a broader term.")

            except FileNotFoundError:
                 st.error(f"Error: The database file was not found at: '{data_path}'. Ensure the path is correct.")
            except ImportError:
                st.error("Error: 'openpyxl' needed for .xlsx files. Install: pip install openpyxl")
            except Exception as e:
                st.error(f"An unexpected error occurred:")
                st.error(f"Type: {type(e).__name__}, Details: {e}")
                import traceback
                st.code(traceback.format_exc())


# DNA Sequence Analysis
elif menu == "DNA Sequence Analysis":
    st.header("DNA Sequence Analysis Tool")
    st.markdown("Enter a DNA sequence (A, T, C, G, N) for basic analysis.")
    sequence = st.text_area("Enter DNA Sequence:", height=150, key="dna_seq_input", placeholder="e.g., ATGCGTAGCTAGCATGCNNNNATGC")

    if st.button("Analyze Sequence", key="dna_analyze_button"):
        if sequence:
            sequence_cleaned = "".join(sequence.split()).upper()
            valid_bases = {'A', 'T', 'C', 'G', 'N'}
            invalid_chars = set(c for c in sequence_cleaned if c not in valid_bases)

            if invalid_chars:
                st.error(f"Invalid characters found: {', '.join(sorted(list(invalid_chars)))}. Use only A, T, C, G, N.")
            elif not sequence_cleaned:
                 st.error("Sequence is empty after cleaning.")
            else:
                try:
                    st.subheader("Analysis Results")
                    length = len(sequence_cleaned)
                    st.write(f"**Sequence Length:** {length} bp")

                    seq_no_n = sequence_cleaned.replace('N', '')
                    gc = gc_fraction(seq_no_n) * 100 if len(seq_no_n) > 0 else 0
                    st.write(f"**GC Content (excluding Ns):** {gc:.2f}%" if len(seq_no_n) > 0 else "**GC Content:** N/A (Only Ns?)")

                    st.write("**Base Composition:**")
                    composition = Counter(sequence_cleaned)
                    comp_data = [{"Base": base, "Count": count, "Percentage": (count/length*100) if length > 0 else 0}
                                 for base, count in sorted(composition.items())]
                    st.dataframe(pd.DataFrame(comp_data).style.format({"Percentage": "{:.2f}%"}), hide_index=True)

                    plot_composition = {k: v for k, v in composition.items() if k != 'N'}
                    if plot_composition:
                        fig, ax = plt.subplots()
                        colors = {'A': 'blue', 'T': 'red', 'C': 'orange', 'G': 'green'}
                        bar_colors = [colors.get(base, 'gray') for base in plot_composition.keys()]
                        ax.bar(plot_composition.keys(), plot_composition.values(), color=bar_colors)
                        ax.set_xlabel("Base"); ax.set_ylabel("Count"); ax.set_title("Base Composition (excluding N)")
                        st.pyplot(fig)
                    elif 'N' in composition and len(composition) == 1:
                         st.info("Plot not shown (sequence contains only 'N' bases).")

                except Exception as e:
                    st.error(f"Analysis error: {e}")
        else:
            st.warning("Please enter a DNA sequence.")

# --- Other Menu Sections (Primer Design, Bioinformatics Tool, etc.) ---
# Keep the code for the other sections as it was in the previous complete version.
# They are omitted here for brevity but should be included in your final file.

# Primer Design Tool
elif menu == "Primer Design":
    st.header("Basic Primer Design Tool")
    st.markdown("Enter a DNA template sequence (A, T, C, G only) to generate simple forward and reverse primers.")
    sequence = st.text_area("Enter DNA Template Sequence:", height=150, key="primer_seq_input", placeholder="e.g., ATGCGTACGT...")
    primer_length = st.slider("Desired Primer Length:", min_value=15, max_value=30, value=20, key="primer_len_slider")

    if st.button("Design Primers", key="primer_design_button"):
        sequence_cleaned = "".join(sequence.split()).upper()

        if not sequence_cleaned:
            st.warning("Please enter a DNA template sequence.")
        else:
            valid_bases = {'A', 'T', 'C', 'G'}
            invalid_chars = set(c for c in sequence_cleaned if c not in valid_bases)
            if invalid_chars:
                st.error(f"Invalid characters found: {', '.join(sorted(list(invalid_chars)))}. Use only A, T, C, G.")
            elif len(sequence_cleaned) < primer_length * 2:
                st.error(f"Sequence too short ({len(sequence_cleaned)} bp). Needs >= {primer_length * 2} bp.")
            else:
                try:
                    forward_primer = sequence_cleaned[:primer_length]
                    reverse_primer_template = sequence_cleaned[-primer_length:]
                    reverse_primer = str(Seq(reverse_primer_template).reverse_complement())

                    st.subheader("Suggested Primers")
                    col1, col2 = st.columns(2)
                    with col1:
                        st.markdown("**Forward Primer (5' -> 3'):**")
                        st.code(forward_primer, language='text')
                        st.write(f"- Length: {len(forward_primer)} bp")
                        st.write(f"- GC Content: {gc_fraction(forward_primer)*100:.2f}%")
                    with col2:
                        st.markdown("**Reverse Primer (5' -> 3'):**")
                        st.code(reverse_primer, language='text')
                        st.write(f"- Length: {len(reverse_primer)} bp")
                        st.write(f"- GC Content: {gc_fraction(reverse_primer)*100:.2f}%")
                    st.info("Note: Basic suggestions only. Check Tm, specificity, secondary structures, etc. for real applications.")
                except Exception as e:
                    st.error(f"Primer design error: {e}")

# Bioinformatics Tool (Transcription/Translation)
elif menu == "Bioinformatics Tool":
    st.header("DNA Transcription & Translation Tool")
    st.markdown("Enter a DNA sequence (coding strand, A, T, C, G only) to transcribe and translate.")
    sequence = st.text_area("Enter DNA Sequence (Coding Strand):", height=150, key="bioinfo_seq_input", placeholder="e.g., ATGCGTAGCTAGCATGC")

    if st.button("Transcribe and Translate", key="bioinfo_submit_button"):
        sequence_cleaned = "".join(sequence.split()).upper()
        if not sequence_cleaned:
            st.warning("Please enter a DNA sequence.")
        else:
            valid_bases = {'A', 'T', 'C', 'G'}
            invalid_chars = set(c for c in sequence_cleaned if c not in valid_bases)
            if invalid_chars:
                st.error(f"Invalid characters found: {', '.join(sorted(list(invalid_chars)))}. Use only A, T, C, G.")
            else:
                try:
                    seq_obj = Seq(sequence_cleaned)
                    st.subheader("Results")
                    st.write(f"**DNA Length:** {len(seq_obj)} bp | **GC Content:** {gc_fraction(seq_obj) * 100:.2f}%")
                    st.markdown("**Reverse Complement (DNA):**"); st.code(str(seq_obj.reverse_complement()), language='text')
                    st.markdown("**Transcription (DNA -> RNA):**"); rna_transcript = seq_obj.transcribe(); st.code(str(rna_transcript), language='text')
                    st.markdown("**Translation (RNA -> Protein):**"); st.caption("Standard genetic code, stops at first stop codon ('*').")

                    remainder = len(rna_transcript) % 3
                    if remainder != 0:
                        st.warning(f"RNA length ({len(rna_transcript)}) not multiple of 3. Last {remainder} base(s) ignored.")
                        rna_transcript_trimmed = rna_transcript[:-remainder]
                    else:
                        rna_transcript_trimmed = rna_transcript

                    if len(rna_transcript_trimmed) > 0:
                         protein_translation = rna_transcript_trimmed.translate(to_stop=True)
                         st.code(str(protein_translation), language='text')
                         st.write(f"- Protein Length: {len(protein_translation)} aa")
                    else:
                         st.error("Cannot translate: Sequence length 0 after trimming.")
                except Exception as e:
                    st.error(f"Processing error: {e}")


# Protein Sequence Analysis Tool
elif menu == "Protein Sequence Analysis":
    st.header("Protein Sequence Analysis Tool")
    st.markdown("Enter a protein sequence (single-letter codes) for property analysis.")
    protein_sequence = st.text_area("Enter Protein Sequence:", height=150, key="protein_seq_input", placeholder="e.g., MKTAYIA...")

    if st.button("Analyze Protein", key="protein_analyze_button"):
        protein_sequence_cleaned = "".join(protein_sequence.split()).upper()
        if not protein_sequence_cleaned:
            st.warning("Please enter a protein sequence.")
        else:
            valid_aa = set("ACDEFGHIKLMNPQRSTVWY") # Standard 20 AA for ProteinAnalysis
            protein_for_analysis = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', protein_sequence_cleaned)
            invalid_chars_in_input = set(protein_sequence_cleaned) - set("ACDEFGHIKLMNPQRSTVWYX*")

            if invalid_chars_in_input:
                st.warning(f"Input sequence contains non-standard characters: {', '.join(sorted(list(invalid_chars_in_input)))}. Analysis performed only on standard AAs (A-Y).")
            if not protein_for_analysis:
                st.error("Sequence contains only non-standard characters after filtering. Cannot perform analysis.")
            else:
                try:
                    analyzed_seq = ProteinAnalysis(protein_for_analysis)
                    st.subheader("Analysis Results (Standard AAs only)")
                    st.write(f"**Analyzed Length:** {len(protein_for_analysis)} aa")
                    st.write(f"**Molecular Weight:** {analyzed_seq.molecular_weight():.2f} Da")
                    st.write(f"**Isoelectric Point (pI):** {analyzed_seq.isoelectric_point():.2f}")
                    st.write(f"**GRAVY Score:** {analyzed_seq.gravy():.3f}")
                    instability = analyzed_seq.instability_index(); stability = "Stable" if instability < 40 else "Unstable"
                    st.write(f"**Instability Index:** {instability:.2f} ({stability})")

                    st.subheader("Amino Acid Composition (%)")
                    aa_percent = analyzed_seq.get_amino_acids_percent()
                    aa_df = pd.DataFrame([{"Amino Acid": aa, "Percentage": p * 100} for aa, p in aa_percent.items()])
                    st.dataframe(aa_df.sort_values(by='Percentage', ascending=False).style.format({'Percentage': '{:.2f}%'}), hide_index=True)

                    st.subheader("Secondary Structure Fraction (Predicted)")
                    try: helix, turn, sheet = analyzed_seq.secondary_structure_fraction(); st.write(f"- Helix: {helix*100:.2f}%, Turn: {turn*100:.2f}%, Sheet: {sheet*100:.2f}%")
                    except Exception as ss_e: st.warning(f"Could not calculate secondary structure: {ss_e}")
                except Exception as e: st.error(f"Analysis error: {e}")


# Motif Finder Tool
elif menu == "Motif Finder Tool":
    st.header("Motif Finder Tool")
    st.markdown("Find occurrences of a motif (exact string or simple regex) in a DNA/Protein sequence.")
    sequence = st.text_area("Enter DNA/Protein Sequence:", height=150, key="motif_seq_input")
    motif = st.text_input("Enter Motif to Find:", key="motif_pattern_input", placeholder="e.g., ATG or Y[DN]R or peptide")
    seq_type = st.radio("Sequence Type:", ("DNA", "Protein"), horizontal=True, key="motif_seq_type")
    case_sensitive = st.checkbox("Case Sensitive Search", key="motif_case_sensitive")

    if st.button("Find Motif", key="motif_find_button"):
        if not sequence or not motif:
            st.warning("Please enter both sequence and motif.")
        else:
            sequence_cleaned = "".join(sequence.split())
            if not sequence_cleaned: st.warning("Sequence is empty.")
            else:
                 flags = 0 if case_sensitive else re.IGNORECASE
                 try:
                     matches = list(re.finditer(motif, sequence_cleaned, flags=flags))
                     st.subheader("Motif Search Results")
                     if matches:
                         st.success(f"Found {len(matches)} occurrence(s).")
                         match_data = [{"Position (1-based Start)": m.start() + 1, "Match": m.group()} for m in matches]
                         st.dataframe(match_data, hide_index=True)
                         # Highlighting logic (optional, kept from previous)
                     else: st.warning("Motif not found.")
                 except re.error as e: st.error(f"Invalid motif pattern (Regex Error): {e}.")
                 except Exception as e: st.error(f"Motif search error: {e}")


# Genome Coverage Plotter
elif menu == "Genome Coverage Plotter":
    st.header("Genome Coverage Plotter")
    st.markdown("Upload a coverage file (CSV/TSV) with 'position' and 'coverage' columns.")
    uploaded_file = st.file_uploader("Upload Coverage Data:", type=['csv', 'tsv', 'txt'], key="coverage_uploader")

    if uploaded_file:
        col1, col2 = st.columns(2)
        with col1: sep = st.selectbox("Separator:", (",", "\\t", " "), index=0, key="cov_sep", format_func=lambda x: {' ': 'Space', ',': 'Comma', '\\t': 'Tab'}.get(x))
        with col2: header_row = st.number_input("Header Row (0-index):", 0, value=0, key="cov_header")
        separator = sep if sep != ' ' else r'\s+'

        try:
            df = pd.read_csv(uploaded_file, sep=separator, header=header_row, engine='python')
            st.dataframe(df.head(), height=200)

            pos_col, cov_col = None, None
            normalized_cols = {col.strip().lower(): col for col in df.columns}
            for p in ['position', 'pos', 'coordinate', 'location', 'start']:
                if p in normalized_cols: pos_col = normalized_cols[p]; break
            for c in ['coverage', 'cov', 'depth', 'reads']:
                if c in normalized_cols: cov_col = normalized_cols[c]; break

            if not pos_col or not cov_col: st.error(f"Could not detect 'position'/'coverage' columns. Found: {list(df.columns)}. Check settings.")
            else:
                st.success(f"Using: Pos='{pos_col}', Cov='{cov_col}'.")
                try:
                    df[pos_col] = pd.to_numeric(df[pos_col], errors='coerce')
                    df[cov_col] = pd.to_numeric(df[cov_col], errors='coerce')
                    rows_before = len(df)
                    df.dropna(subset=[pos_col, cov_col], inplace=True)
                    if len(df) < rows_before: st.warning(f"Removed {rows_before - len(df)} rows with non-numeric data.")

                    if df.empty: st.error("No valid numeric data found.")
                    else:
                        st.subheader("Genome Coverage Plot")
                        df = df.sort_values(by=pos_col)
                        fig, ax = plt.subplots(figsize=(12, 6))
                        ax.plot(df[pos_col], df[cov_col], label='Coverage', lw=1.5)
                        # Smoothing option (kept from previous)
                        ax.set_xlabel(f"Position ({pos_col})"); ax.set_ylabel(f"Coverage ({cov_col})"); ax.set_title("Coverage Plot"); ax.legend(); ax.grid(True, alpha=0.5); ax.ticklabel_format(style='plain', axis='x')
                        st.pyplot(fig)
                        st.subheader("Coverage Statistics"); st.dataframe(df[cov_col].describe().to_frame().T.style.format("{:,.2f}"))
                except Exception as plot_e: st.error(f"Plotting/Processing Error: {plot_e}")
        except Exception as read_e: st.error(f"File Reading Error: {read_e}. Check format/settings.")


# Variant Annotation Tool (Simple Substitution)
elif menu == "Variant Annotation Tool":
    st.header("Simple Variant Substitution Tool")
    st.markdown("Enter reference DNA (A,T,C,G) and specify a single base substitution.")
    ref_seq = st.text_area("Reference DNA Sequence:", height=100, key="var_ref_seq", placeholder="e.g., ATGCGTACGTAGCTAG")
    col1, col2 = st.columns(2)
    with col1: variant_position = st.number_input("Variant Position (1-based):", 1, step=1, key="var_pos")
    with col2: variant_base = st.text_input("New Variant Base (A/T/C/G):", max_chars=1, key="var_base").upper()

    if st.button("Apply Variant Substitution", key="var_apply_button"):
        ref_seq_cleaned = "".join(ref_seq.split()).upper()
        # Validation... (kept from previous, ensures inputs are valid)
        error_flag = False
        if not ref_seq_cleaned or not variant_base: st.warning("Need sequence and variant base."); error_flag = True
        valid_bases = {'A', 'T', 'C', 'G'}
        if any(b not in valid_bases for b in ref_seq_cleaned): st.error("Invalid ref base."); error_flag = True
        if variant_base not in valid_bases and variant_base: st.error("Invalid variant base."); error_flag = True
        if not error_flag and (variant_position > len(ref_seq_cleaned) or variant_position <= 0) : st.error("Position out of range."); error_flag = True

        if not error_flag:
            try:
                zero_based_pos = variant_position - 1
                original_base = ref_seq_cleaned[zero_based_pos]
                if original_base == variant_base: st.info(f"Base at {variant_position} is already '{original_base}'.")
                else:
                    alt_seq = list(ref_seq_cleaned); alt_seq[zero_based_pos] = variant_base; alt_seq = "".join(alt_seq)
                    st.subheader("Variant Applied")
                    st.write(f"**Variant:** Pos {variant_position} (Ref: {original_base} -> Alt: {variant_base})")
                    st.write("**Reference:**"); st.code(ref_seq_cleaned, language='text')
                    st.write("**Altered:**"); st.code(alt_seq, language='text')
                    # Context and Codon Impact (optional, kept from previous)
            except Exception as e: st.error(f"Error applying variant: {e}")


# Codon Usage Analyzer
elif menu == "Codon Usage Analyzer":
    st.header("Codon Usage Analyzer")
    st.markdown("Enter a coding DNA sequence (CDS - multiple of 3, A,T,C,G only).")
    dna_sequence = st.text_area("Coding DNA Sequence (CDS):", height=150, key="codon_seq_input", placeholder="e.g., ATGCGT...")

    if st.button("Analyze Codon Usage", key="codon_analyze_button"):
        dna_sequence_cleaned = "".join(dna_sequence.split()).upper()
        if not dna_sequence_cleaned: st.warning("Please enter sequence.")
        else:
            valid_bases = {'A', 'T', 'C', 'G'}
            if any(b not in valid_bases for b in dna_sequence_cleaned): st.error("Invalid bases found.")
            elif len(dna_sequence_cleaned) % 3 != 0: st.error(f"Length ({len(dna_sequence_cleaned)}) not multiple of 3.")
            else:
                try:
                    codons = [dna_sequence_cleaned[i:i + 3] for i in range(0, len(dna_sequence_cleaned), 3)]
                    total_codons = len(codons)
                    codon_counts = Counter(codons)
                    st.subheader("Codon Usage Results")
                    st.write(f"**Total Codons:** {total_codons}")

                    from Bio.Data import CodonTable
                    try: standard_table = CodonTable.unambiguous_dna_by_id[1] # Standard code
                    except Exception as table_e: st.error(f"Codon table error: {table_e}"); st.stop()

                    codon_data = []
                    for codon in sorted(list(standard_table.forward_table.keys()) + standard_table.stop_codons):
                        count = codon_counts.get(codon, 0)
                        freq = (count / total_codons) * 100 if total_codons > 0 else 0
                        aa = standard_table.forward_table.get(codon, '*' if codon in standard_table.stop_codons else '?')
                        codon_data.append({"Codon": codon, "AA": aa, "Count": count, "Freq (%)": freq})

                    usage_df = pd.DataFrame(codon_data)
                    st.dataframe(usage_df.style.format({'Freq (%)': '{:.2f}%'}), hide_index=True)
                    # Plotting (optional, kept from previous)
                except Exception as e: st.error(f"Codon analysis error: {e}")


# --- Footer ---
st.markdown('<div class="footer">© 2024 Colon Cancer Toolbox | For Educational & Informational Purposes Only</div>', unsafe_allow_html=True)
