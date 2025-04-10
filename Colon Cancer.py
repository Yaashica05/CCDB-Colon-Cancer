import streamlit as st
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, MeltingTemp
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO, pairwise2
from Bio.Data import CodonTable, IUPACData
from Bio.Align import substitution_matrices
from Bio.Restriction import RestrictionBatch, Analysis
import matplotlib.pyplot as plt
import re
from collections import Counter
import urllib.parse
import io
import os # For checking file existence

# --- Page Configuration (MUST be the first Streamlit command) ---
st.set_page_config(
    page_title="CCDB: Colon Cancer Database and Bioinformatics Tools",
    layout="wide",
    initial_sidebar_state="expanded",
)

# --- Custom CSS Styling (Refined Dark Theme + Times New Roman + Bold) ---
st.markdown(
    """
    <style>
    /* --- Global Font Settings --- */
    body, .stApp, h1, h2, h3, h4, h5, h6, p, div, span, li, label,
    .stButton>button,
    .stTextInput>div>div>input,
    .stTextArea>div>textarea,
    .stSelectbox>div>div>div, /* Covers selected item */
    .stSelectbox [data-baseweb="select"] > div:first-child, /* Dropdown display */
    .stMultiSelect [data-baseweb="tag"] span, /* Multiselect tags */
    .stSlider label, .stNumberInput label, .stRadio label, .stCheckbox label,
    table th, table td, table a,
    .stAlert,
    .footer,
    /* Ensure plot text elements inherit if possible (often needs Python-side control) */
    .plot-container text
    {
        font-family: 'Times New Roman', Times, serif !important;
        font-weight: bold !important;
    }

    /* Remove underlines globally, except on hover for links */
    a, a:link, a:visited {
        text-decoration: none !important;
    }
    a:hover {
        text-decoration: underline !important;
    }

    /* --- Specific Element Styling (REFINED DARK THEME) --- */

    /* Body background - Subtle Dark Gradient *** UPDATED *** */
    body {
        background-color: #0b0c10; /* Fallback */
        background-image: linear-gradient(to bottom right, #0b0c10, #1a1a2e); /* Dark blue/purple gradient */
        background-attachment: fixed; /* Keep gradient fixed during scroll */
        min-height: 100vh; /* Ensure gradient covers full height */
    }

    /* Main app container - Dark Blue-Grey */
    .stApp {
        background-color: #1f2833;
        padding: 30px; border-radius: 10px;
        border: 1px solid #45a29e; /* Teal border */
        box-shadow: 0 5px 15px rgba(0, 0, 0, 0.4);
        padding-bottom: 80px; /* Space for footer */
        color: #d0d0d0; /* Default text color - light grey */
    }

    /* Sidebar Styling - Slightly different shade */
    [data-testid="stSidebar"] {
        background-color: #1a1a2e; /* Dark Blue/Purple */
        border-right: 2px solid #45a29e; /* Teal separator */
        padding-top: 2rem; /* Add some top padding */
    }
    /* Sidebar Header & Widget Labels */
    [data-testid="stSidebar"] h1,
    [data-testid="stSidebar"] h2,
    [data-testid="stSidebar"] h3,
    [data-testid="stSidebar"] label /* General label targeting */
    {
         color: #66fcf1 !important; /* Bright Cyan */
    }
    /* Sidebar Radio Button Options Text */
    [data-testid="stSidebar"] .stRadio div[role="radiogroup"] > label {
         color: #e0e0e0 !important; /* Light grey for options */
         font-weight: bold !important; /* Ensure boldness */
         font-family: 'Times New Roman', Times, serif !important;
    }
     /* Sidebar expander */
    [data-testid="stSidebar"] [data-testid="stExpander"] summary {
        color: #97e8e1; /* Lighter Cyan */
    }

    /* Headings - Bright Cyan */
    h1, h2, h3, h4, h5, h6 { color: #66fcf1; }
    h1 { border-bottom: 2px solid #45a29e; padding-bottom: 10px; margin-bottom: 20px; }
    h2 { padding-bottom: 5px; margin-top: 30px; margin-bottom: 15px;}
    h3 { margin-top: 25px; margin-bottom: 10px; color: #97e8e1; } /* Lighter Cyan for H3 */


    /* Paragraph text - Light Grey */
    p, div, span, li { color: #e0e0e0; line-height: 1.7; font-size: 16px; }

     /* Labels for Widgets - Lighter Cyan */
    .stTextInput>label, .stTextArea>label, .stSelectbox>label, .stNumberInput>label,
    .stFileUploader>label, .stSlider>label, .stRadio>label, .stCheckbox>label,
    .stMultiSelect>label
    {
        color: #97e8e1; font-size: 1.1em; margin-bottom: 5px;
    }

    /* Input fields, Select boxes, Text Areas */
    .stTextInput>div>div>input,
    .stTextArea>div>textarea,
    .stSelectbox>div>div>div, /* Container for selected value */
    .stMultiSelect>div>div /* Container for multiselect */
    {
        border: 1px solid #4f5b66; border-radius: 5px; background-color: #2c3e50; /* Darker input bg */
        color: #e1e1e1; padding: 9px 12px;
    }
    .stTextArea>div>textarea { min-height: 150px; }
    .stSelectbox [data-baseweb="select"] > div { background-color: #2c3e50; } /* Ensure dropdown bg matches */
    .stMultiSelect [data-baseweb="tag"] { background-color: #45a29e; color: #0b0c10; } /* Teal tags */

    /* Input Placeholders */
    .stTextInput>div>div>input::placeholder, .stTextArea>div>textarea::placeholder { color: #a0a0a0; opacity: 0.8; }

    /* Buttons - Cyan */
    .stButton>button {
        border-radius: 5px; background-color: #66fcf1; color: #0b0c10; /* Cyan bg, very dark text */
        border: 1px solid #45a29e; padding: 10px 22px;
        transition: background-color 0.3s ease, transform 0.1s ease, box-shadow 0.2s ease;
        font-size: 1.05em; box-shadow: 0 2px 4px rgba(0,0,0,0.2);
    }
    .stButton>button:hover {
        background-color: #48dbcd; color: #0b0c10; transform: translateY(-2px);
        box-shadow: 0 4px 8px rgba(102, 252, 241, 0.2); border-color: #66fcf1;
    }
    .stButton>button:active {
        background-color: #45a29e; transform: translateY(0px);
        box-shadow: inset 0 2px 4px rgba(0,0,0,0.3); border-color: #31726e;
    }
     /* Disabled button styling */
    .stButton>button:disabled {
        background-color: #4f5b66; color: #999; border-color: #4f5b66;
        cursor: not-allowed;
    }

    /* Tables - Consistent Dark Theme */
    table {
        border-collapse: separate; border-spacing: 0; width: 100%; margin-bottom: 1.5rem;
        background-color: #25303a; /* Slightly lighter than main bg */
        border: 1px solid #4f5b66; border-radius: 6px; overflow: hidden;
    }
    table th {
        text-align: left; padding: 12px 15px; background-color: #2c3e50; /* Header bg */
        color: #97e8e1; /* Header text - Lighter Cyan */
        border-bottom: 2px solid #45a29e; /* Teal bottom border */
    }
    table td {
        vertical-align: middle; text-align: left; padding: 12px 15px;
        border-top: 1px solid #4f5b66; /* Internal border */
        color: #d0d0d0; /* Cell text */
    }
    /* Table links */
    table a { color: #fca311; /* Bright Orange links */ }
    table a:hover { color: #ffbf00; /* Lighter Orange on hover */ }
    /* Zebra striping */
    table tbody tr:nth-of-type(odd) { background-color: #25303a; } /* Base row color */
    table tbody tr:nth-of-type(even) { background-color: #2c3e50; } /* Slightly different even row */
    table tbody tr:hover { background-color: #3a4b5a; } /* Hover row color */

    /* Alerts - Themed */
    .stAlert { border-radius: 5px; padding: 15px; font-size: 1.05em; border: 1px solid; }
    .stAlert > div[data-testid="stMarkdownContainer"] > p { color: inherit !important; font-weight: bold !important; font-family: 'Times New Roman', Times, serif !important; }
    div[data-baseweb="alert"][role="alert"] { border-left: 5px solid; }
    /* Info */
    div[data-testid="stInfo"] { border-color: #66fcf1; background-color: #2a3b4c; color: #d1ffff; }
    /* Success */
    div[data-testid="stSuccess"] { border-color: #50fa7b; background-color: #283a2f; color: #c3e88d; }
    /* Warning */
    div[data-testid="stWarning"] { border-color: #ffbf00; background-color: #4d3c1a; color: #ffd54f; }
    /* Error */
    div[data-testid="stError"] { border-color: #ff5555; background-color: #4d2a2a; color: #ff8a80; }

    /* Footer - Fixed at bottom */
    .footer {
        position: fixed; left: 0; bottom: 0; width: 100%;
        background-color: #0b0c10; /* Match body gradient start */
        background-image: linear-gradient(to right, #0b0c10, #1a1a2e); /* Match body gradient */
        color: #97e8e1; /* Lighter Cyan text */
        text-align: center; padding: 12px; font-size: 14px;
        border-top: 3px solid #45a29e; /* Teal top border */
        z-index: 1000;
        font-family: 'Times New Roman', Times, serif !important; font-weight: bold !important;
    }

    /* Ensure Plotly/Matplotlib backgrounds are transparent or match theme */
    .plot-container.plotly, .stPlotlyChart, .stVegaLiteChart {
        background-color: transparent !important;
    }
    /* Style Streamlit expander */
    [data-testid="stExpander"] {
        border: 1px solid #4f5b66;
        border-radius: 5px;
        background-color: #25303a; /* Slightly lighter than main bg */
    }
    [data-testid="stExpander"] summary {
        background-color: #2c3e50; /* Header part */
        color: #97e8e1; /* Lighter cyan text */
        border-radius: 5px 5px 0 0;
        padding: 10px 15px;
    }
    [data-testid="stExpander"] summary:hover {
        background-color: #3a4b5a;
    }

    /* Style Streamlit Divider */
    hr {
        border-top: 1px solid #4f5b66; /* Subtle divider */
        margin-top: 1rem;
        margin-bottom: 1rem;
    }

    </style>
    """,
    unsafe_allow_html=True,
)


# --- Global Variables & Constants ---
DATA_PATH = "Colon Cancer.xlsx"
STANDARD_DNA_BASES = {'A', 'T', 'C', 'G'}
EXTENDED_DNA_BASES = {'A', 'T', 'C', 'G', 'N'}
STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"
AMBIGUOUS_AA = "BZXJ*"

# --- Helper function to parse FASTA files ---
def parse_fasta(uploaded_file):
    """Parses the first sequence from an uploaded FASTA file."""
    sequence = None
    seq_id = None
    if uploaded_file is not None:
        try:
            # Reset file pointer just in case it was read before
            uploaded_file.seek(0)
            content_bytes = uploaded_file.getvalue()
            try: content_str = content_bytes.decode("utf-8")
            except UnicodeDecodeError: content_str = content_bytes.decode("latin-1")

            stringio = io.StringIO(content_str)
            records = list(SeqIO.parse(stringio, "fasta"))
            if records:
                first_record = records[0]
                sequence = str(first_record.seq).upper()
                seq_id = first_record.id
                st.success(f"Read FASTA ID: **{seq_id}** ({len(sequence):,} bases/AAs).")
            else: st.error("FASTA file empty or invalid.")
        except Exception as e:
            st.error(f"FASTA Read Error: {e}")
    return sequence, seq_id

# --- Helper function to create links for the Database Search ---
def format_link(value, column_name):
    """Creates an HTML link based on column name and value with styling."""
    if pd.isna(value): return ""
    s_value = str(value).strip()
    if not s_value: return ""

    col_name_std = column_name.strip().title()
    # Using the orange link color defined in CSS (fca311)
    link_style = 'style="color:#fca311;"' # Bold is global, underline on hover handled by CSS

    try:
        if col_name_std == 'Pubmed Id':
            if s_value.replace('.', '', 1).isdigit(): url = f"https://pubmed.ncbi.nlm.nih.gov/{s_value}"
            elif s_value.startswith("http"): url = s_value
            else: return s_value # Return as text if not a number or URL
            return f'<a href="{url}" target="_blank" {link_style}>{s_value}</a>'
        elif col_name_std == 'Doi Id':
            if '/' in s_value and not s_value.startswith("http"):
                encoded_doi = urllib.parse.quote(s_value, safe='/:')
                url = f"https://doi.org/{encoded_doi}"
            elif s_value.startswith("http"): url = s_value
            else: return s_value
            return f'<a href="{url}" target="_blank" {link_style}>{s_value}</a>'
        elif col_name_std == 'Uniprot':
             # Simple check for common UniProtKB ID format or URL
             if re.match(r'^[A-Z0-9_.\-]+$', s_value, re.IGNORECASE) and not s_value.startswith("http"):
                 url = f"https://www.uniprot.org/uniprotkb/{s_value}/entry"
             elif s_value.startswith("http"): url = s_value
             else: return s_value
             return f'<a href="{url}" target="_blank" {link_style}>{s_value}</a>'
        # General URL detection for specified columns or columns containing 'link', 'url', etc.
        elif col_name_std in ['Blast', 'Conserved Domain', 'Link', 'Url', 'Reference'] or 'link' in col_name_std.lower() or 'url' in col_name_std.lower():
             if s_value.startswith("http"):
                  display_url = s_value if len(s_value) < 50 else s_value[:47] + "..."
                  return f'<a href="{s_value}" target="_blank" {link_style} title="{s_value}">{display_url}</a>'
             else: return s_value # Return non-HTTP values as text
        else:
            # Fallback: If any other column contains a value starting with http, treat it as a link
            if s_value.startswith("http"):
                display_url = s_value if len(s_value) < 50 else s_value[:47] + "..."
                return f'<a href="{s_value}" target="_blank" {link_style} title="{s_value}">{display_url}</a>'
            else: return s_value # Default: return value as text
    except Exception: # Catch any unexpected formatting error
        return s_value # Return original value if linking fails

# --- Helper function for Matplotlib Plot Styling ---
def style_plot(fig, ax, title="", xlabel="", ylabel=""):
    """Applies consistent dark theme styling to a Matplotlib plot."""
    fig.patch.set_facecolor('#1f2833') # Match app background
    fig.patch.set_alpha(1.0) # Ensure figure background is opaque
    ax.set_facecolor('#2c3e50') # Slightly different background for plot area

    # Titles and Labels
    ax.set_title(title, color="#66fcf1", fontweight='bold', fontsize=14, fontfamily='Times New Roman')
    ax.set_xlabel(xlabel, color="#e0e0e0", fontweight='bold', fontsize=12, fontfamily='Times New Roman')
    ax.set_ylabel(ylabel, color="#e0e0e0", fontweight='bold', fontsize=12, fontfamily='Times New Roman')

    # Ticks
    ax.tick_params(axis='x', colors='#d0d0d0', labelsize=10)
    ax.tick_params(axis='y', colors='#d0d0d0', labelsize=10)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontfamily('Times New Roman')
        label.set_fontweight('bold')

    # Grid
    ax.grid(True, linestyle=':', alpha=0.4, color='#888888')

    # Spines (borders of the plot area)
    for spine in ax.spines.values():
        spine.set_edgecolor('#4f5b66')

    plt.tight_layout()

# --- Main App Title ---
st.title("🧬 CCDB: Colon Cancer Database and Bioinformatics Tools 🔬")
# (Removed the markdown separator here)

# --- Sidebar Navigation ---
with st.sidebar:
    st.header("🛠️ Toolbox Menu")
    menu = st.radio(
        "Select Feature:",
        [
            "Home", "Colon Cancer Database Search", "DNA Sequence Analysis",
            "Protein Sequence Analysis", "Primer Design", "Restriction Enzyme Analysis",
            "Pairwise Sequence Alignment", "Motif Finder Tool",
            "Bioinformatics Tool (Transcription/Translation)", "Genome Coverage Plotter",
            "Variant Annotation Tool", "Codon Usage Analyzer",
        ],
        key="main_menu", help="Choose a tool from the list."
    )
    st.divider()
    st.info("Ensure input files (FASTA, coverage data, etc.) are readily available for upload.")
    # Optionally add more sidebar content, links, or info here


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
        *   **Primer Design:** Generate basic forward and reverse primers from a DNA template.
        *   **Restriction Enzyme Analysis:** Find restriction enzyme cut sites in a DNA sequence.
        *   **Pairwise Sequence Alignment:** Align two DNA or protein sequences using global or local algorithms.
        *   **Motif Finder Tool:** Locate specific sequence patterns (motifs) within DNA or protein sequences using exact matches or regular expressions.
        *   **Bioinformatics Tool (Transcription/Translation):** Transcribe DNA to RNA and translate RNA to protein based on the standard genetic code.
        *   **Genome Coverage Plotter:** Visualize sequencing coverage depth across genomic positions from an uploaded file.
        *   **Variant Annotation Tool:** Apply a simple single-base substitution to a reference DNA sequence and predict its codon effect.
        *   **Codon Usage Analyzer:** Analyze the frequency of codon usage in a coding DNA sequence (CDS).

        ---
        *Disclaimer: This tool is intended for educational and informational purposes only. It is not a substitute for professional medical advice, diagnosis, or treatment. Always consult with a qualified healthcare provider regarding any medical conditions or treatment options.*
        """, unsafe_allow_html=True
    )

    # Using the Freepik Infographic image
    st.image(
        "https://img.freepik.com/free-vector/colorectal-cancer-crc-infographic-education_1308-47971.jpg",
        caption="Colorectal Cancer Infographic (Source: Freepik - Illustrative purposes)",
        use_container_width=True
    )
    st.markdown("""<a href="https://www.cancer.org/cancer/types/colon-rectal-cancer.html" target="_blank" style="color:#fca311;">Learn more about Colorectal Cancer from the American Cancer Society</a>""", unsafe_allow_html=True)

    st.divider()
    st.subheader("Feedback")
    feedback = st.text_area("Share your feedback about this application:", key="feedback_home_page", placeholder="Your thoughts, suggestions, or issues...")
    if st.button("Submit Feedback", key="submit_feedback_home_page"):
        if feedback:
            st.success("Thank you for your feedback!")
            # In a real app, you'd save this feedback (e.g., to a file or database)
            print(f"Feedback received (Home Page): {feedback}")
        else:
            st.warning("Feedback cannot be empty.")


elif menu == "Colon Cancer Database Search":
    st.header("Colon Cancer Gene Database Search")
    st.markdown(f"Search for information on genes associated with Colorectal Cancer using a local database file (`{DATA_PATH}`). Please enter a gene symbol.")

    # Check if the database file exists
    if not os.path.exists(DATA_PATH):
        st.error(f"**Database File Not Found:** The required file `{DATA_PATH}` could not be located in the application directory.")
        st.info("Please ensure the Excel database file is present in the same folder as the Streamlit script.")
        st.stop() # Stop execution if the file is missing

    query = st.text_input("Enter Gene Symbol:", key="gene_search_input", placeholder="e.g., APC, KRAS, MLH1")

    if st.button("Search Database", key="gene_search_button"):
        if not query:
            st.warning("Please enter a gene symbol to search for.")
        else:
            try:
                # Load the Excel file
                data = pd.read_excel(DATA_PATH)
                # Standardize column names (strip whitespace, title case)
                data.columns = data.columns.str.strip().str.title()

                # Identify the gene symbol column (flexible matching)
                gene_col = None
                possible_gene_cols = ['Gene Symbol', 'Gene Name', 'Symbol', 'Gene']
                for col in possible_gene_cols:
                    if col in data.columns:
                        gene_col = col
                        break

                if not gene_col:
                     st.error(f"Could not find a suitable gene symbol column in the database. Expected one of: {possible_gene_cols}. Found columns: {list(data.columns)}")
                     st.stop()

                # Prepare search query and perform search (case-insensitive exact match first)
                data[gene_col] = data[gene_col].astype(str).str.strip() # Ensure string type and remove whitespace
                search_query = query.strip()
                # Use str.fullmatch for exact matching (case-insensitive)
                results = data[data[gene_col].str.fullmatch(search_query, case=False, na=False)].copy()

                # Fallback to contains search if exact match yields no results
                if results.empty:
                     st.info(f"No exact match found for '{search_query}'. Trying a broader search...")
                     # Using regex=False for simple substring containment might be safer if gene names have special chars
                     results = data[data[gene_col].str.contains(search_query, case=False, na=False, regex=False)].copy()

                if not results.empty:
                    st.success(f"Found {len(results)} result(s) containing '{query}'.")
                    st.write("### Search Results:")
                    # Apply link formatting to relevant columns
                    formatted_results = results.copy()
                    for col in formatted_results.columns:
                         if col in ['Pubmed Id', 'Doi Id', 'Uniprot', 'Blast', 'Conserved Domain', 'Link', 'Url', 'Reference'] or 'link' in col.lower() or 'url' in col.lower():
                             try:
                                 formatted_results[col] = formatted_results[col].apply(lambda x: format_link(x, col))
                             except Exception as apply_e:
                                 st.warning(f"Link formatting failed for column '{col}': {apply_e}. Displaying raw data.")
                    # Display results as an HTML table to render links correctly
                    html_table = formatted_results.to_html(escape=False, index=False, na_rep='-', justify='left', classes=['dataframe', 'st-table']) # Add st-table class for potential Streamlit overrides
                    st.write(html_table, unsafe_allow_html=True)
                else:
                    st.warning(f"No results found matching or containing '{query}' in the '{gene_col}' column.")

            except ImportError:
                st.error("The `openpyxl` library is required to read Excel files. Please install it (`pip install openpyxl`).")
            except FileNotFoundError:
                 st.error(f"**Database File Not Found:** The file `{DATA_PATH}` was not found during the read attempt.")
            except Exception as e:
                st.error(f"An error occurred during the database search:")
                st.exception(e) # Display the full error traceback for debugging


elif menu == "DNA Sequence Analysis":
    st.header("DNA Sequence Analysis")
    st.markdown("Analyze basic properties of a DNA sequence (accepts A, T, C, G, N). Input sequence by pasting or uploading a FASTA file.")

    input_method = st.radio("Input Method:", ("Paste", "FASTA"), key="dna_input", horizontal=True)
    sequence, seq_id = "", "Pasted Sequence"
    if input_method == "Paste":
        sequence = st.text_area("Enter DNA Sequence:", height=150, key="dna_seq_in", placeholder="Paste your DNA sequence here (e.g., ATGC...). Non-ATCGN characters will be ignored.")
    else:
        uploaded_file = st.file_uploader("Upload FASTA File:", type=['fasta','fa','fna'], key="dna_fasta_up")
        if uploaded_file:
            sequence, seq_id = parse_fasta(uploaded_file)
            if sequence:
                st.text_area("Sequence Preview (from FASTA):", value=f">{seq_id}\n{sequence[:100]}{'...' if len(sequence) > 100 else ''}", height=100, key="dna_disp", disabled=True)

    st.divider()
    if st.button("Analyze DNA", key="dna_analyze_btn"):
        if not sequence:
            st.warning("Please provide a DNA sequence via paste or FASTA upload.")
        else:
            # Clean the sequence: keep only A, T, C, G, N (case-insensitive)
            seq_clean = "".join(re.findall(r'[ATCGN]', sequence.upper()))
            if not seq_clean:
                st.error("The provided input does not contain any valid DNA bases (A, T, C, G, N).")
            else:
                try:
                    st.subheader(f"Analysis Results for: '{seq_id}'")
                    seq_obj = Seq(seq_clean)
                    length = len(seq_obj)
                    st.metric("Total Length", f"{length:,} bp")

                    # Calculate GC/AT content excluding 'N' bases
                    seq_no_n = seq_clean.replace('N', '')
                    len_no_n = len(seq_no_n)
                    if len_no_n > 0:
                        gc_val = gc_fraction(seq_no_n) * 100
                        at_val = 100.0 - gc_val
                        st.metric("GC Content (excluding N)", f"{gc_val:.1f}%")
                        st.metric("AT Content (excluding N)", f"{at_val:.1f}%")
                    else:
                        st.info("GC/AT content could not be calculated (sequence contains only N or is empty after cleaning).")

                    # Base Composition Calculation
                    st.divider()
                    st.write("#### Base Composition (including N):")
                    composition = Counter(seq_clean)
                    comp_data = [{"Base": b, "Count": c, "Percentage (%)": (c/length*100) if length > 0 else 0}
                                 for b, c in sorted(composition.items())]
                    comp_df = pd.DataFrame(comp_data)
                    st.dataframe(comp_df.style.format({"Count":"{:,}", "Percentage (%)":"{:.1f}%"}),
                                 hide_index=True, use_container_width=True)

                    # Plot Base Composition (excluding N)
                    plot_comp = {k: v for k, v in composition.items() if k in STANDARD_DNA_BASES}
                    if plot_comp:
                        fig, ax = plt.subplots(figsize=(6, 4))
                        bases = sorted(plot_comp.keys())
                        counts = [plot_comp[b] for b in bases]
                        # Define colors consistent with theme accents
                        base_colors = {'A': '#66fcf1', 'T': '#fca311', 'C': '#ffbf00', 'G': '#ff5555'}
                        bar_colors = [base_colors.get(b, '#888888') for b in bases]

                        ax.bar(bases, counts, color=bar_colors)
                        style_plot(fig, ax, title="Base Counts (excluding N)", xlabel="Base", ylabel="Count")
                        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ','))) # Format y-axis ticks
                        st.pyplot(fig)
                    elif 'N' in composition and len(composition) == 1:
                        st.info("Composition plot omitted (sequence contains only N).")
                    elif not plot_comp:
                         st.info("Composition plot omitted (no standard bases A,T,C,G found).")

                    # Basic ORF Finder
                    st.divider()
                    st.write("#### Basic Open Reading Frame (ORF) Finder")
                    min_orf_aa = st.number_input("Minimum ORF Length (Amino Acids):", min_value=10, max_value=500, value=30, step=5, key="orf_len")
                    orfs = []
                    st.info("Searching for ORFs (Met 'M' to Stop '*') on the forward strand (3 frames)...")
                    with st.spinner("Scanning for ORFs... This might take a moment for long sequences."):
                         for frame in range(3):
                            # Translate from current frame, stopping at the first stop codon encountered
                            try:
                                translation = str(seq_obj[frame:].translate(to_stop=True))
                            except Exception as trans_err:
                                st.warning(f"Translation error in frame {frame+1}: {trans_err}. Skipping frame.")
                                continue

                            current_pos_in_translation = 0
                            while 'M' in translation[current_pos_in_translation:]:
                                # Find the next potential start codon 'M'
                                start_codon_index = translation.find('M', current_pos_in_translation)
                                if start_codon_index == -1:
                                    break # No more 'M's in the remaining translation

                                # The ORF is the sequence from 'M' onwards (as translate already stopped at '*')
                                potential_orf = translation[start_codon_index:]

                                if len(potential_orf) >= min_orf_aa:
                                    # Calculate DNA coordinates (1-based)
                                    start_dna_pos = frame + (start_codon_index * 3) + 1
                                    end_dna_pos = start_dna_pos + (len(potential_orf) * 3) - 1
                                    orfs.append({
                                        "Frame": frame + 1,
                                        "Start (DNA)": start_dna_pos,
                                        "End (DNA)": end_dna_pos,
                                        "Length (AA)": len(potential_orf),
                                        "Protein Sequence (Preview)": potential_orf[:40] + "..." if len(potential_orf) > 40 else potential_orf
                                    })

                                # Move search position past the current 'M' to find subsequent ORFs
                                current_pos_in_translation = start_codon_index + 1

                    if orfs:
                         st.success(f"Found {len(orfs)} potential ORF(s) ≥ {min_orf_aa} AA.")
                         orf_df = pd.DataFrame(orfs).sort_values(["Frame", "Start (DNA)"])
                         st.dataframe(orf_df.style.format({"Start (DNA)":"{:,}", "End (DNA)":"{:,}", "Length (AA)":"{:,}"}),
                                      hide_index=True, use_container_width=True)
                    else:
                         st.info(f"No potential ORFs found with a minimum length of {min_orf_aa} amino acids.")
                except Exception as e:
                    st.error("An error occurred during DNA sequence analysis:")
                    st.exception(e)


elif menu == "Protein Sequence Analysis":
    st.header("Protein Sequence Analysis")
    st.markdown("Analyze biochemical properties of a protein sequence. Input sequence by pasting or uploading a FASTA file.")

    input_method = st.radio("Input Method:", ("Paste", "FASTA"), key="prot_in", horizontal=True)
    sequence, seq_id = "", "Pasted Sequence"
    if input_method == "Paste":
        sequence = st.text_area("Enter Protein Sequence:", height=150, key="prot_seq", placeholder="Paste your protein sequence here (e.g., MKT...). Non-standard/ambiguous characters will be handled.")
    else:
        uploaded_file = st.file_uploader("Upload FASTA File:", type=['fasta','fa','faa'], key="prot_fasta")
        if uploaded_file:
            sequence, seq_id = parse_fasta(uploaded_file)
            if sequence:
                 st.text_area("Sequence Preview (from FASTA):", value=f">{seq_id}\n{sequence[:100]}{'...' if len(sequence) > 100 else ''}", height=100, key="prot_disp", disabled=True)

    st.divider()
    if st.button("Analyze Protein", key="prot_analyze_btn"):
        if not sequence:
            st.warning("Please provide a protein sequence via paste or FASTA upload.")
        else:
            # Clean the sequence: remove whitespace and convert to uppercase
            seq_clean = "".join(sequence.split()).upper()
            all_chars = set(seq_clean)
            standard_set = set(STANDARD_AA)
            ambiguous_set = set(AMBIGUOUS_AA)

            # Identify non-standard and ambiguous characters
            non_standard_chars = all_chars - standard_set - ambiguous_set
            ambiguous_chars_present = all_chars.intersection(ambiguous_set)

            # Sequence used for calculations (only standard AAs)
            seq_for_calc = "".join(c for c in seq_clean if c in standard_set)

            # Display warnings if non-standard or ambiguous characters are found
            if non_standard_chars:
                st.warning(f"Ignoring unknown characters: `{'`, `'.join(sorted(list(non_standard_chars)))}`")
            if ambiguous_chars_present:
                st.warning(f"Ignoring ambiguous amino acid codes for calculations: `{'`, `'.join(sorted(list(ambiguous_chars_present)))}`")

            if not seq_for_calc:
                st.error("No standard amino acids found in the sequence. Cannot perform analysis.")
            else:
                try:
                    st.subheader(f"Analysis Results for: '{seq_id}'")
                    pa = ProteinAnalysis(seq_for_calc)

                    # Display key metrics
                    col1, col2 = st.columns(2)
                    with col1:
                        st.metric("Input Length (incl. non-std)", f"{len(seq_clean):,} aa")
                        st.metric("Analyzed Length (standard AA)", f"{len(seq_for_calc):,} aa")
                        st.metric("Molecular Weight", f"{pa.molecular_weight():,.1f} Da")
                        st.metric("Isoelectric Point (pI)", f"{pa.isoelectric_point():.2f}")
                    with col2:
                        # GRAVY Score
                        gravy_score = pa.gravy()
                        hydrophobicity = "Hydrophobic" if gravy_score > 0 else ("Hydrophilic" if gravy_score < 0 else "Neutral")
                        st.metric("GRAVY Index", f"{gravy_score:.3f}", delta=hydrophobicity, delta_color="off") # delta_color="off" prevents default up/down arrows

                        # Instability Index
                        instability_idx = pa.instability_index()
                        stability = "Likely Stable" if instability_idx < 40 else "Likely Unstable"
                        # Use delta_color to indicate stability (normal=good=stable)
                        st.metric("Instability Index", f"{instability_idx:.2f}", delta=stability,
                                  delta_color="normal" if stability == "Likely Stable" else "inverse")

                    # Amino Acid Composition (based on the full input sequence)
                    st.divider()
                    st.write("#### Amino Acid Composition (Full Input Sequence):")
                    composition = Counter(seq_clean)
                    comp_data = [{"Amino Acid": aa, "Count": c, "Percentage (%)": (c/len(seq_clean)*100) if seq_clean else 0}
                                 for aa, c in sorted(composition.items())]
                    comp_df = pd.DataFrame(comp_data)
                    st.dataframe(comp_df.style.format({'Count':'{:,}', 'Percentage (%)':'{:.1f}%'}),
                                 hide_index=True, use_container_width=True)

                    # Secondary Structure Prediction
                    st.divider()
                    st.write("#### Predicted Secondary Structure Fractions:")
                    try:
                        helix, turn, sheet = pa.secondary_structure_fraction()
                        labels = 'Helix', 'Turn', 'Sheet'
                        sizes = [helix * 100, turn * 100, sheet * 100] # Convert to percentages
                        # Define colors
                        colors = ['#66fcf1', '#fca311', '#ffbf00'] # Cyan, Orange, Amber

                        fig, ax = plt.subplots()
                        # Use wedgeprops for styling, autopct for percentages
                        wedges, texts, autotexts = ax.pie(sizes, labels=labels, startangle=90, colors=colors, autopct='%1.1f%%',
                                              wedgeprops={'edgecolor': '#1f2833', 'linewidth': 1},
                                              textprops={'color': '#e0e0e0', 'fontweight': 'bold', 'fontfamily': 'Times New Roman'})

                        # Improve autopct text visibility
                        plt.setp(autotexts, size=10, weight="bold", color='#0b0c10', family='Times New Roman') # Dark text on light wedges

                        ax.axis('equal') # Equal aspect ratio ensures a circular pie chart

                        # Style the plot using the helper function (sets background, etc.)
                        style_plot(fig, ax, title="Predicted Secondary Structure")
                        # Remove axis labels for pie chart
                        ax.set_xlabel("")
                        ax.set_ylabel("")

                        # Add a legend (optional, as labels are on wedges)
                        # ax.legend(wedges, labels, title="Structure", loc="center left", bbox_to_anchor=(1.05, 0.5),
                        #           labelcolor='#e0e0e0', facecolor='#2c3e50', edgecolor='#4f5b66',
                        #           prop={'family':'Times New Roman', 'weight':'bold', 'size': 10})

                        st.pyplot(fig)

                    except ZeroDivisionError:
                         st.warning("Could not calculate secondary structure fractions (possibly due to sequence length or composition).")
                    except Exception as ss_e:
                        st.warning(f"Secondary structure prediction failed: {ss_e}")

                except Exception as e:
                    st.error("An error occurred during protein sequence analysis:")
                    st.exception(e)


elif menu == "Primer Design":
    st.header("Basic Primer Design Tool")
    st.markdown("Generates simple forward and reverse primers from the ends of a DNA template (accepts A, T, C, G only). Estimates Melting Temperature (Tm).")
    st.warning("⚠️ **Educational Tool:** This performs very basic primer selection from sequence ends. It **does not** check for specificity (e.g., BLAST), secondary structures (hairpins, dimers), primer-dimers, or compatibility. Use dedicated primer design software for experimental work.", icon="🔬")

    input_method = st.radio("Input Method:", ("Paste", "FASTA"), key="p_in", horizontal=True)
    sequence, seq_id = "", "Pasted Sequence"
    if input_method == "Paste":
        sequence = st.text_area("Enter DNA Template Sequence:", height=150, key="p_seq", placeholder="Paste your DNA sequence here (e.g., ATGC...). Only ATCG characters will be used.")
    else:
        uploaded_file = st.file_uploader("Upload FASTA File:", type=['fasta','fa','fna'], key="p_fasta")
        if uploaded_file:
            sequence, seq_id = parse_fasta(uploaded_file)
            if sequence:
                st.text_area("Sequence Preview (from FASTA):", value=f">{seq_id}\n{sequence[:100]}{'...' if len(sequence) > 100 else ''}", height=100, key="p_disp", disabled=True)

    # Primer Design Parameters
    st.divider()
    st.subheader("Primer Parameters")
    col1, col2 = st.columns(2)
    with col1:
        primer_len = st.slider("Desired Primer Length (bp):", min_value=15, max_value=35, value=20, step=1, key="p_len")
    with col2:
        # Tm Calculation Method Selection
        tm_method = st.selectbox("Melting Temperature (Tm) Method:",
                                 ["Tm_NN", "Tm_GC", "Tm_Wallace"],
                                 index=0, key="p_tm_m",
                                 help="Tm_NN (Nearest Neighbor) is generally more accurate but requires concentration inputs. Tm_GC and Tm_Wallace are simpler approximations.")
        # Conditional inputs for Tm_NN
        dna_nM, salt_mM = 50.0, 50.0 # Default values
        if tm_method == "Tm_NN":
            with st.expander("Tm_NN Parameters (Optional)"):
                dna_nM = st.number_input("Primer Concentration (nM):", min_value=1.0, max_value=1000.0, value=50.0, step=1.0, key="p_dna_c", format="%.1f")
                salt_mM = st.number_input("Monovalent Cation Conc. (mM, e.g., Na+):", min_value=1.0, max_value=200.0, value=50.0, step=1.0, key="p_salt_c", format="%.1f")
        else:
            # Keep defaults if not Tm_NN, but they won't be used by Tm_GC/Tm_Wallace
            dna_nM, salt_mM = 50.0, 50.0

    st.divider()
    if st.button("Design Primers", key="p_design_btn"):
        if not sequence:
            st.warning("Please provide a DNA template sequence.")
        else:
            # Clean sequence: keep only A, T, C, G (case-insensitive)
            seq_clean = "".join(re.findall(r'[ATCG]', sequence.upper()))
            if not seq_clean:
                st.error("The provided input does not contain any valid DNA bases (A, T, C, G) required for primer design.")
            elif len(seq_clean) < primer_len * 2:
                st.error(f"Template sequence is too short ({len(seq_clean):,} bp). It needs to be at least twice the primer length ({primer_len * 2} bp).")
            else:
                try:
                    # Design primers from the ends
                    fw_p = Seq(seq_clean[:primer_len])
                    rv_template_segment = Seq(seq_clean[-primer_len:])
                    rv_p = rv_template_segment.reverse_complement()

                    # Calculate GC content
                    fw_gc = gc_fraction(fw_p) * 100
                    rv_gc = gc_fraction(rv_p) * 100

                    # Calculate Tm
                    fw_tm, rv_tm = "N/A", "N/A" # Initialize as Not Available
                    try:
                        tm_args = {'strict': False} # Allow approximate calculations if needed
                        if tm_method == "Tm_NN":
                            # Convert nM to M for dnac1 if needed (check Bio.SeqUtils docs - typically expects M)
                            # Assuming MeltingTemp expects nM directly based on prior code structure. If errors occur, conversion might be needed.
                            # tm_args.update({'dnac1': dna_nM / 1e9, 'Na': salt_mM}) # Example if M is needed
                            tm_args.update({'dnac1': dna_nM, 'Na': salt_mM}) # Pass nM and mM as is (verify expected units)
                        else:
                             # Tm_GC and Tm_Wallace don't typically use these args, but pass defaults just in case
                             tm_args.update({'dnac1': 50.0, 'Na': 50.0})

                        # Get the specific Tm function from MeltingTemp module
                        tm_func = getattr(MeltingTemp, tm_method)
                        fw_tm = tm_func(fw_p, **tm_args)
                        rv_tm = tm_func(rv_p, **tm_args)
                    except ValueError as tm_ve:
                        st.warning(f"Tm calculation issue: {tm_ve}. Tm values might be inaccurate or unavailable.")
                    except AttributeError:
                         st.error(f"Could not find the Tm calculation method '{tm_method}' in Bio.SeqUtils.MeltingTemp.")
                    except Exception as tm_e:
                         st.error(f"An unexpected error occurred during Tm calculation: {tm_e}")

                    # Display Results
                    st.subheader("Suggested Primers (from sequence ends)")
                    p_col1, p_col2 = st.columns(2)
                    with p_col1:
                        st.markdown("#### Forward Primer (5' → 3')")
                        st.code(str(fw_p), language='text')
                        st.metric("Length", f"{len(fw_p)} bp")
                        st.metric("GC Content", f"{fw_gc:.1f}%")
                        st.metric(f"Estimated Tm ({tm_method})", f"{fw_tm:.1f} °C" if isinstance(fw_tm, (float, int)) else "N/A")
                    with p_col2:
                        st.markdown("#### Reverse Primer (5' → 3')")
                        st.code(str(rv_p), language='text')
                        st.metric("Length", f"{len(rv_p)} bp")
                        st.metric("GC Content", f"{rv_gc:.1f}%")
                        st.metric(f"Estimated Tm ({tm_method})", f"{rv_tm:.1f} °C" if isinstance(rv_tm, (float, int)) else "N/A")

                    st.divider()
                    st.metric(f"Resulting Amplicon Size (based on full template '{seq_id}')", f"{len(seq_clean):,} bp")

                    # Add caption about Tm parameters if Tm_NN was used and successful
                    if tm_method == "Tm_NN" and isinstance(fw_tm, (float, int)):
                        st.caption(f"Tm ({tm_method}) calculated using: Primer Conc. = {dna_nM:.1f} nM, Salt Conc. = {salt_mM:.1f} mM.")
                    elif isinstance(fw_tm, (float, int)):
                        st.caption(f"Tm calculated using the '{tm_method}' formula/method.")

                except Exception as e:
                    st.error("An error occurred during primer design:")
                    st.exception(e)


elif menu == "Restriction Enzyme Analysis":
    st.header("Restriction Enzyme Analysis")
    st.markdown("Identify restriction enzyme cut sites within a DNA sequence (accepts A, T, C, G only). Input sequence by pasting or uploading a FASTA file.")
    st.info("This tool uses enzyme data from Biopython's REBASE interface. Analysis is performed on the provided DNA strand.")

    input_method = st.radio("Input Method:", ("Paste", "FASTA"), key="re_in", horizontal=True)
    sequence, seq_id = "", "Pasted Sequence"
    if input_method == "Paste":
        sequence = st.text_area("Enter DNA Sequence:", height=150, key="re_seq", placeholder="Paste your DNA sequence here (e.g., GAATTC...). Only ATCG characters will be used.")
    else:
        uploaded_file = st.file_uploader("Upload FASTA File:", type=['fasta','fa','fna'], key="re_fasta")
        if uploaded_file:
            sequence, seq_id = parse_fasta(uploaded_file)
            if sequence:
                st.text_area("Sequence Preview (from FASTA):", value=f">{seq_id}\n{sequence[:100]}{'...' if len(sequence) > 100 else ''}", height=100, key="re_disp", disabled=True)

    # Enzyme Selection
    st.divider()
    st.subheader("Enzyme Selection")
    # List of common enzymes for convenience
    common_enzymes = sorted(['EcoRI','BamHI','HindIII','NotI','XhoI','PstI','SacI','SalI','SpeI','SphI','KpnI','SmaI','EcoRV','HaeIII','HhaI','AluI','BglII','NcoI','NdeI','XbaI','ClaI','DpnI','MboI','Sau3AI','TaqI'])
    selected_common = st.multiselect("Select common enzymes:", common_enzymes, default=['EcoRI','BamHI','HindIII'], key="re_sel")

    custom_enzymes_str = st.text_input("Enter custom enzymes (comma-separated, case-sensitive):", key="re_cust", placeholder="e.g., BsaI, BsmBI, PacI")
    # Combine selected and custom enzymes
    all_selected_enzymes = set(selected_common)
    if custom_enzymes_str:
        custom_list = {e.strip() for e in custom_enzymes_str.split(',') if e.strip()}
        all_selected_enzymes.update(custom_list)

    final_enzyme_list = sorted(list(all_selected_enzymes))

    # Topology option
    is_linear = st.checkbox("Assume Linear DNA", value=True, key="re_lin", help="Uncheck if the sequence represents circular DNA (e.g., a plasmid).")

    st.divider()
    if st.button("Analyze Restriction Sites", key="re_analyze_btn"):
        if not sequence:
            st.warning("Please provide a DNA sequence.")
        elif not final_enzyme_list:
            st.warning("Please select at least one restriction enzyme.")
        else:
            # Clean sequence: keep only A, T, C, G
            seq_clean = "".join(re.findall(r'[ATCG]', sequence.upper()))
            if not seq_clean:
                st.error("The provided input does not contain any valid DNA bases (A, T, C, G).")
            else:
                try:
                    seq_obj = Seq(seq_clean)
                    st.info(f"Analyzing sequence '{seq_id}' ({len(seq_obj):,} bp, {'Linear' if is_linear else 'Circular'}).")

                    # Validate selected enzymes against Biopython's known enzymes
                    valid_enzymes_obj = []
                    invalid_enzyme_names = []
                    from Bio.Restriction import AllEnzymes # Load the enzyme list
                    with st.spinner("Validating selected enzymes..."):
                        for enz_name in final_enzyme_list:
                            if enz_name in AllEnzymes:
                                valid_enzymes_obj.append(AllEnzymes.get(enz_name)) # Get the enzyme object
                            else:
                                invalid_enzyme_names.append(enz_name)

                    if invalid_enzyme_names:
                        st.warning(f"The following enzyme names were not recognized and will be ignored: `{'`, `'.join(invalid_enzyme_names)}`")

                    if not valid_enzymes_obj:
                        st.error("No valid restriction enzymes were selected or found. Please check the enzyme names.")
                        st.stop()

                    # Perform the restriction analysis
                    st.write(f"**Analyzing with:** {', '.join(map(str, valid_enzymes_obj))}")
                    with st.spinner("Searching for cut sites..."):
                        rb = RestrictionBatch(valid_enzymes_obj) # Use the validated enzyme objects
                        analysis = Analysis(rb, seq_obj, linear=is_linear)

                    st.subheader("Restriction Analysis Results")
                    st.write(f"**Sequence:** '{seq_id}' | **Length:** {len(seq_obj):,} bp | **Topology:** {'Linear' if is_linear else 'Circular'}")
                    st.divider()

                    # Get the results and summarize
                    results_dict = analysis.full() # Dictionary: {EnzymeObject: [site_positions]}
                    summary_data = []
                    total_cuts = 0
                    sites_found = False

                    for enz_obj, sites in results_dict.items():
                        cut_count = len(sites)
                        total_cuts += cut_count
                        if cut_count > 0:
                            sites_found = True
                        summary_data.append({
                            "Enzyme": str(enz_obj),
                            "Recognition Site": str(enz_obj.site),
                            "Number of Cuts": cut_count,
                            "Cut Positions (1-based)": ", ".join(map(str, sites)) if sites else "None"
                        })

                    if not sites_found:
                        st.success("✅ No cut sites were found for the selected valid enzymes in this sequence.")
                    else:
                        st.write("#### Cut Site Summary:")
                        summary_df = pd.DataFrame(summary_data)
                        # Sort by number of cuts (descending), then by enzyme name (ascending)
                        summary_df = summary_df.sort_values(by=["Number of Cuts", "Enzyme"], ascending=[False, True])
                        st.dataframe(summary_df.style.format({"Number of Cuts":"{:,}"}),
                                     hide_index=True, use_container_width=True)

                        # Fragment Analysis
                        st.divider()
                        st.write("#### Predicted DNA Fragments:")
                        try:
                            fragments = analysis.fragments # Get fragment sequences (as strings)
                            fragment_lengths = sorted([len(f) for f in fragments], reverse=True) # Sort lengths descending

                            f_col1, f_col2 = st.columns(2)
                            with f_col1:
                                st.metric("Total Cuts Detected", f"{total_cuts:,}")
                            with f_col2:
                                st.metric("Number of Fragments", f"{len(fragments):,}")

                            st.write("**Fragment Lengths (bp, sorted descending):**")
                            if fragment_lengths:
                                if len(fragment_lengths) > 30: # Show first few and count rest if too many
                                    st.text(", ".join(f"{l:,}" for l in fragment_lengths[:15]) + f", ... ({len(fragment_lengths)-15} more fragments)")
                                else:
                                    st.text(", ".join(f"{l:,}" for l in fragment_lengths))
                            else:
                                # Should only happen if 0 cuts on linear, or 1 cut on circular (results in one linear piece)
                                st.text(f"Single fragment of {len(seq_obj):,} bp (or no fragments if analysis failed).")

                        except Exception as frag_e:
                            st.warning(f"Could not calculate or display fragment lengths: {frag_e}")

                        # Option to show text map
                        st.divider()
                        show_map = st.checkbox("Show Text-Based Restriction Map", value=False, key="re_map")
                        if show_map:
                            st.write("#### Restriction Map (Text Representation):")
                            try:
                                # Use print_that to get map string, then display in code block
                                map_output = io.StringIO()
                                analysis.print_that(out=map_output, title=f"Map for {seq_id}", top=True) # Render map to string buffer
                                map_string = map_output.getvalue()
                                map_output.close()
                                st.code(map_string, language='text')
                            except Exception as map_e:
                                st.error(f"Failed to generate text map: {map_e}")

                except ValueError as ve:
                    st.error(f"Restriction analysis failed. Check sequence and enzyme names. Error: {ve}")
                except Exception as e:
                    st.error("An unexpected error occurred during restriction analysis:")
                    st.exception(e)


elif menu == "Pairwise Sequence Alignment":
    st.header("Pairwise Sequence Alignment")
    st.markdown("Align two DNA or protein sequences using global (Needleman-Wunsch) or local (Smith-Waterman) algorithms via Biopython's `pairwise2`.")
    st.info("Scores depend heavily on the chosen parameters (match/mismatch scores or substitution matrix, gap penalties). Ensure parameters are appropriate for your sequence type (DNA/Protein).")

    # Input areas for two sequences
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("#### Sequence 1")
        input_method1 = st.radio("Input Source:",("Paste","FASTA"), key="al_in1", horizontal=True, label_visibility="collapsed")
        seq1, seq_id1 = "", "Sequence 1"
        if input_method1 == "Paste":
            seq1 = st.text_area("Enter Sequence 1:", height=120, key="al_s1", placeholder="Paste sequence here...")
        else:
            uploaded_file1 = st.file_uploader("Upload FASTA File 1:", type=['fasta','fa','fna','faa'], key="al_f1")
            if uploaded_file1:
                seq1, seq_id1 = parse_fasta(uploaded_file1)
        # Display preview if sequence is loaded
        if seq1:
            st.text_area("Preview Seq 1:", value=f">{seq_id1}\n{seq1[:80]}{'...' if len(seq1)>80 else ''}", height=75, key="al_d1", disabled=True)

    with col2:
        st.markdown("#### Sequence 2")
        input_method2 = st.radio("Input Source:",("Paste","FASTA"), key="al_in2", horizontal=True, label_visibility="collapsed")
        seq2, seq_id2 = "", "Sequence 2"
        if input_method2 == "Paste":
            seq2 = st.text_area("Enter Sequence 2:", height=120, key="al_s2", placeholder="Paste sequence here...")
        else:
            uploaded_file2 = st.file_uploader("Upload FASTA File 2:", type=['fasta','fa','fna','faa'], key="al_f2")
            if uploaded_file2:
                seq2, seq_id2 = parse_fasta(uploaded_file2)
        # Display preview if sequence is loaded
        if seq2:
             st.text_area("Preview Seq 2:", value=f">{seq_id2}\n{seq2[:80]}{'...' if len(seq2)>80 else ''}", height=75, key="al_d2", disabled=True)

    # Alignment Parameters
    st.divider()
    st.subheader("Alignment Parameters")
    pcol1, pcol2, pcol3 = st.columns(3)
    with pcol1:
        alignment_mode = st.selectbox("Alignment Mode:", ["Global (Needleman-Wunsch)", "Local (Smith-Waterman)"], key="al_m")
        sequence_type = st.radio("Sequence Type:", ("DNA", "Protein"), key="al_t", horizontal=True)

    with pcol2:
        st.write("**Scoring Parameters:**")
        substitution_matrix = None
        match_score, mismatch_penalty = None, None
        if sequence_type == "DNA":
            match_score = st.number_input("Match Score:", value=2.0, step=0.5, key="al_mt", help="Score for identical characters.")
            mismatch_penalty = st.number_input("Mismatch Penalty:", value=-1.0, step=-0.5, key="al_ms", help="Penalty for non-identical characters (should be negative or zero).")
            mismatch_penalty = min(0.0, mismatch_penalty) # Ensure penalty is not positive
        else: # Protein
            available_matrices = sorted(substitution_matrices.list_matrices())
            default_matrix_index = available_matrices.index('BLOSUM62') if 'BLOSUM62' in available_matrices else 0
            selected_matrix_name = st.selectbox("Substitution Matrix:", available_matrices, index=default_matrix_index, key="al_mat")
            try:
                substitution_matrix = substitution_matrices.load(selected_matrix_name)
                st.caption(f"Using **{selected_matrix_name}** matrix.")
            except Exception as e:
                st.error(f"Error loading substitution matrix '{selected_matrix_name}': {e}")
                substitution_matrix = None # Ensure it's None if loading failed

    with pcol3:
        st.write("**Gap Penalties:**")
        gap_open_penalty = st.number_input("Gap Open Penalty:", value=-10.0, step=-0.5, key="al_go", help="Penalty for starting a gap (should be negative or zero).")
        gap_extend_penalty = st.number_input("Gap Extend Penalty:", value=-0.5, step=-0.1, key="al_ge", help="Penalty for extending a gap (should be negative or zero).")
        gap_open_penalty = min(0.0, gap_open_penalty) # Ensure penalties are not positive
        gap_extend_penalty = min(0.0, gap_extend_penalty)

    st.divider()
    if st.button("Align Sequences", key="al_btn"):
        # --- Input Validation ---
        valid_input = True
        if not seq1:
            st.warning("Sequence 1 is missing.")
            valid_input = False
        if not seq2:
            st.warning("Sequence 2 is missing.")
            valid_input = False
        if sequence_type == "Protein" and not substitution_matrix:
            st.error("A valid protein substitution matrix must be selected and loaded.")
            valid_input = False
        if not valid_input:
            st.stop() # Halt execution if inputs are invalid

        # --- Sequence Cleaning (Remove non-alphanumeric, keep standard chars + '-') ---
        # Define valid characters based on type
        valid_chars_pattern = re.compile(f"[^{re.escape(IUPACData.protein_letters if sequence_type == 'Protein' else IUPACData.unambiguous_dna_letters)}]")
        seq1_clean = valid_chars_pattern.sub('', seq1.upper())
        seq2_clean = valid_chars_pattern.sub('', seq2.upper())

        if not seq1_clean:
            st.error("Sequence 1 contains no valid characters after cleaning.")
            valid_input = False
        if not seq2_clean:
            st.error("Sequence 2 contains no valid characters after cleaning.")
            valid_input = False
        if not valid_input:
            st.stop()

        # --- Perform Alignment ---
        try:
            st.subheader("Alignment Results")
            # Determine alignment function and parameters based on mode and type
            mode_prefix = "global" if alignment_mode.startswith("Global") else "local"
            params = {}
            align_function_name = ""

            if substitution_matrix: # Protein alignment
                align_function_name = f"{mode_prefix}ds" # Use matrix (d) and default scores (s) - Biopython convention
                params = {'matrix': substitution_matrix, 'open': gap_open_penalty, 'extend': gap_extend_penalty}
            else: # DNA alignment
                align_function_name = f"{mode_prefix}ms" # Use match/mismatch scores (m) and default scores (s)
                params = {'match': match_score, 'mismatch': mismatch_penalty, 'open': gap_open_penalty, 'extend': gap_extend_penalty}

            # Check if the alignment function exists in pairwise2.align
            if hasattr( pairwise2.align, align_function_name):
                 align_func = getattr( pairwise2.align, align_function_name)
            else:
                st.error(f"Internal error: Could not find alignment function '{align_function_name}' in `pairwise2.align`.")
                st.stop()

            # Execute alignment
            st.info(f"Performing {alignment_mode} {sequence_type} alignment...")
            with st.spinner("Aligning sequences... This may take time for long sequences."):
                # Request only the best alignment for simplicity
                alignments = align_func(seq1_clean, seq2_clean, **params, one_alignment_only=True)

            if not alignments:
                st.warning("No alignment could be generated with the given sequences and parameters.")
            else:
                # Process the best alignment
                best_alignment = alignments[0]
                # Unpack alignment tuple
                aligned1, aligned2, score, begin, end = best_alignment

                st.metric("Alignment Score", f"{score:.2f}")
                if mode_prefix == "local":
                    st.write(f"**Aligned Region (0-based index):** Start = {begin}, End = {end}")

                # Format and display the alignment
                st.divider()
                st.write("#### Best Alignment:")

                # Call format_alignment WITHOUT name1/name2 arguments
                formatted_alignment_string = pairwise2.format_alignment(
                    aligned1, aligned2, score, begin, end,
                    full_sequences=(mode_prefix == 'global') # Show full length for global
                )

                # Manually construct the output string to include names
                display_text = f"Sequence 1: {seq_id1}\nSequence 2: {seq_id2}\n\n{formatted_alignment_string}"

                # Display the combined text
                st.code(display_text, language='text')

                st.caption("Alignment key: '|' = Match, '.' = Mismatch (or matrix score > 0), ' ' = Negative score mismatch, '-' = Gap.")

                # Calculate Identity and Gaps
                identity_count = 0
                gap_count = 0
                alignment_length = len(aligned1) # Length of the aligned sequences (including gaps)
                aligned_pairs_count = 0 # Number of positions where neither sequence has a gap

                for i in range(alignment_length):
                    res1 = aligned1[i]
                    res2 = aligned2[i]
                    if res1 == '-' or res2 == '-':
                        gap_count += 1
                    else:
                        aligned_pairs_count += 1
                        if res1 == res2:
                            identity_count += 1

                # Calculate percentages
                identity_percent = (identity_count / aligned_pairs_count * 100) if aligned_pairs_count > 0 else 0
                gap_percent = (gap_count / alignment_length * 100) if alignment_length > 0 else 0

                # Display Statistics
                st.divider()
                st.write("#### Alignment Statistics:")
                scol1, scol2, scol3 = st.columns(3)
                with scol1:
                    st.metric("Alignment Length", f"{alignment_length:,}")
                with scol2:
                    # Identity is calculated over non-gapped positions
                    st.metric("Percent Identity", f"{identity_percent:.1f}%", f"{identity_count:,} / {aligned_pairs_count:,} pairs")
                with scol3:
                    st.metric("Gap Percentage", f"{gap_percent:.1f}%", f"{gap_count:,} / {alignment_length:,} total")
                st.caption("Percent Identity is calculated as (Identical Matches) / (Aligned Pairs excluding gaps).")

        except Exception as e:
            st.error("An error occurred during sequence alignment:")
            st.exception(e) # Display the full traceback for debugging


elif menu == "Motif Finder Tool":
    st.header("Sequence Motif Finder")
    st.markdown("Search for specific patterns (motifs) within a DNA or protein sequence using exact matching or regular expressions (regex).")
    st.info("Regex examples: `GAATTC` (exact), `G[AT]ATTC` (A or T at pos 2), `M[A-Z]{5}K` (M, 5 any AA, K), `^ATG` (starts with ATG), `TAG$` (ends with TAG). Use `\\` to escape special characters like `.` or `*`.")

    input_method = st.radio("Input Method:",("Paste","FASTA"), key="mo_in", horizontal=True)
    sequence, seq_id="","Pasted Sequence"
    if input_method=="Paste":
        sequence=st.text_area("Enter Sequence:", height=150, key="mo_seq", placeholder="Paste DNA or protein sequence here...")
    else:
        uploaded_file=st.file_uploader("Upload FASTA File:", type=['fasta','fa','fna','faa'], key="mo_f")
        if uploaded_file:
            sequence, seq_id = parse_fasta(uploaded_file)
    # Display preview if sequence is loaded
    if sequence:
        st.text_area("Sequence Preview:", value=f">{seq_id}\n{sequence[:80]}{'...' if len(sequence)>80 else ''}", height=75, key="mo_d", disabled=True)

    # Motif Input and Options
    st.divider()
    motif_pattern = st.text_input("Enter Motif or Regular Expression (Regex):", key="mo_pat", placeholder="e.g., GAATTC or ^M[A-Z]+")

    mcol1, mcol2 = st.columns(2)
    with mcol1:
        case_sensitive = st.checkbox("Case Sensitive Search", value=False, key="mo_case", help="If checked, 'A' will not match 'a'.")
    with mcol2:
        allow_overlap = st.checkbox("Allow Overlapping Matches", value=True, key="mo_ov", help="If unchecked, the search continues after the end of a found match.")

    st.divider()
    if st.button("Find Motifs", key="mo_find"):
        if not sequence:
            st.warning("Please provide a sequence to search within.")
        elif not motif_pattern:
            st.warning("Please enter a motif or regular expression to search for.")
        else:
            # Clean sequence (remove newlines/whitespace)
            seq_cleaned = "".join(sequence.split())
            if not seq_cleaned:
                st.warning("The provided sequence is empty after cleaning.")
            else:
                try:
                    # Set regex flags based on options
                    regex_flags = 0 if case_sensitive else re.IGNORECASE

                    found_matches = []
                    search_position = 0
                    st.info(f"Searching for motif '{motif_pattern}' in '{seq_id}'...")
                    with st.spinner("Searching..."):
                        # Use finditer to get match objects with positions
                        for match in re.finditer(motif_pattern, seq_cleaned, flags=regex_flags):
                            # Check for overlap condition
                            if allow_overlap or match.start() >= search_position:
                                found_matches.append(match)
                                # Update search position for non-overlapping search
                                if not allow_overlap:
                                    search_position = match.end()

                    st.subheader(f"Search Results for '{motif_pattern}' in '{seq_id}'")
                    if found_matches:
                        st.success(f"Found {len(found_matches):,} match(es).")
                        # Prepare data for table display
                        match_data = [{
                            "#": i + 1,
                            "Start (1-based)": m.start() + 1,
                            "End (1-based)": m.end(),
                            "Length": m.end() - m.start(),
                            "Matched Sequence": m.group()
                        } for i, m in enumerate(found_matches)]

                        match_df = pd.DataFrame(match_data)
                        st.dataframe(match_df.style.format({
                            "Start (1-based)":"{:,}",
                            "End (1-based)":"{:,}",
                            "Length":"{:,}"
                        }), hide_index=True, use_container_width=True)

                        # Option to highlight matches in the sequence (limited length)
                        st.divider()
                        show_highlight = st.checkbox("Highlight Matches in Sequence (first 2000 chars)", value=False, key="mo_hi")
                        if show_highlight:
                             highlight_limit = 2000
                             sequence_to_highlight = seq_cleaned[:highlight_limit]
                             highlighted_html = ""
                             last_end = 0
                             # Define highlighting style
                             mark_style = "background-color:#ffbf00; padding: 1px 3px; border-radius: 3px; color:#1a1a2e; font-weight:bold;" # Amber highlight, dark text
                             mark_open = f"<mark style='{mark_style}'>"
                             mark_close = "</mark>"

                             for match in found_matches:
                                 if match.start() >= highlight_limit:
                                     break # Stop if match starts beyond limit
                                 start, end = match.start(), min(match.end(), highlight_limit) # Cap end at limit
                                 if start >= last_end: # Ensure we don't process overlapping parts multiple times in highlight
                                     highlighted_html += sequence_to_highlight[last_end:start] # Add text before match
                                     highlighted_html += mark_open + sequence_to_highlight[start:end] + mark_close # Add highlighted match
                                     last_end = end

                             highlighted_html += sequence_to_highlight[last_end:] # Add remaining text after last match

                             st.markdown(f"**Highlighted Sequence (first {highlight_limit:,} characters):**")
                             # Use a div with monospace font for better alignment
                             st.markdown(f"<div style='font-family:monospace; word-wrap:break-word; line-height:1.6; border: 1px solid #4f5b66; padding: 10px; border-radius: 5px; background-color: #2c3e50;'>{highlighted_html}{'...' if len(seq_cleaned) > highlight_limit else ''}</div>", unsafe_allow_html=True)

                    else:
                        st.info("The specified motif/pattern was not found in the sequence.")

                except re.error as regex_e:
                    st.error(f"Invalid Regular Expression: '{motif_pattern}'. Please check the syntax.")
                    st.exception(regex_e)
                except Exception as e:
                    st.error("An error occurred during the motif search:")
                    st.exception(e)


elif menu == "Bioinformatics Tool (Transcription/Translation)":
    st.header("DNA Transcription and Translation Tool")
    st.markdown("Transcribe a DNA sequence (coding strand, A, T, C, G only) into RNA, and then translate the RNA into a protein sequence using the standard genetic code.")

    input_method = st.radio("Input Method:",("Paste","FASTA"), key="tr_in", horizontal=True)
    dna_sequence, seq_id="","Pasted Sequence"
    if input_method=="Paste":
        dna_sequence=st.text_area("Enter DNA Sequence (Coding Strand, 5' to 3'):", height=120, key="tr_seq", placeholder="Paste your DNA sequence here (e.g., ATGCGT...). Only ATCG characters will be used.")
    else:
        uploaded_file=st.file_uploader("Upload FASTA File (containing DNA sequence):", type=['fasta','fa','fna'], key="tr_f")
        if uploaded_file:
            dna_sequence, seq_id = parse_fasta(uploaded_file)
            if dna_sequence:
                 st.text_area("Sequence Preview (from FASTA):", value=f">{seq_id}\n{dna_sequence[:80]}{'...' if len(dna_sequence)>80 else ''}", height=75, key="tr_d", disabled=True)

    st.divider()
    if st.button("Transcribe and Translate", key="tr_btn"):
        if not dna_sequence:
            st.warning("Please provide a DNA sequence.")
        else:
            # Clean DNA sequence: keep only A, T, C, G
            dna_cleaned = "".join(re.findall(r'[ATCG]', dna_sequence.upper()))
            if not dna_cleaned:
                st.error("The provided input does not contain any valid DNA bases (A, T, C, G).")
            else:
                try:
                    st.subheader(f"Transcription & Translation Results for: '{seq_id}'")
                    dna_seq_obj = Seq(dna_cleaned)
                    dna_length = len(dna_seq_obj)

                    # Display DNA Info
                    col1, col2 = st.columns(2)
                    with col1:
                        st.metric("Input DNA Length", f"{dna_length:,} bp")
                        st.write("**Input DNA (Coding Strand, 5'→3'):**")
                        st.code(str(dna_seq_obj), language='text')
                    with col2:
                        st.metric("GC Content", f"{gc_fraction(dna_seq_obj)*100:.1f}%")
                        st.write("**Reverse Complement DNA (5'→3'):**")
                        st.code(str(dna_seq_obj.reverse_complement()), language='text')

                    # Transcription
                    st.divider()
                    st.write("#### Transcription (DNA → RNA)")
                    rna_seq_obj = dna_seq_obj.transcribe()
                    st.write("**Resulting RNA Sequence (5'→3'):**")
                    st.code(str(rna_seq_obj), language='text')
                    st.metric("RNA Length", f"{len(rna_seq_obj):,} bases")

                    # Translation
                    st.divider()
                    st.write("#### Translation (RNA → Protein)")
                    st.caption("Translation uses the standard genetic code and starts from the beginning of the RNA sequence. '*' indicates a STOP codon.")

                    # Check if RNA length is a multiple of 3
                    remainder = len(rna_seq_obj) % 3
                    rna_for_translation = rna_seq_obj
                    if remainder != 0:
                        st.warning(f"RNA sequence length ({len(rna_seq_obj):,}) is not a multiple of 3. The last {remainder} base(s) will be ignored for translation.")
                        rna_for_translation = rna_seq_obj[:-remainder]

                    if len(rna_for_translation) >= 3:
                        # Translate, stopping at the first stop codon
                        protein_seq_to_stop = rna_for_translation.translate(to_stop=True)
                        st.write("**Protein Sequence (Translated until first STOP codon):**")
                        st.code(str(protein_seq_to_stop), language='text')
                        st.metric("Protein Length (until STOP)", f"{len(protein_seq_to_stop):,} aa")

                        # Option to show full translation including internal stop codons
                        st.divider()
                        show_full_translation = st.checkbox("Show Full Translation (including internal STOP codons)", value=False, key="tr_full")
                        if show_full_translation:
                            protein_seq_full = rna_for_translation.translate(to_stop=False) # Use default table=1 (standard code)
                            st.write("**Full Protein Sequence (includes internal '*' = STOP):**")
                            st.code(str(protein_seq_full), language='text')
                            st.metric("Full Translation Length", f"{len(protein_seq_full):,} aa")

                    elif len(rna_seq_obj) < 3:
                        st.error("Cannot translate: The RNA sequence is shorter than 3 bases.")
                    else: # Length was multiple of 3, but 0 after trimming (shouldn't happen if initial length >= 3)
                         st.error("Cannot translate: RNA length became 0 after trimming non-codon bases.")

                except Exception as e:
                    st.error("An error occurred during processing:")
                    st.exception(e)


elif menu == "Genome Coverage Plotter":
    st.header("Genome Coverage Plotter")
    st.markdown("Visualize sequencing coverage depth across genomic positions. Upload a file containing position and coverage data (e.g., output from `samtools depth`, BEDGRAPH, or a simple CSV/TSV).")
    st.info("The file should contain at least two columns: one for genomic position and one for coverage depth at that position.")

    uploaded_file = st.file_uploader("Upload Coverage Data File:", type=['csv','tsv','txt','bed', 'bedgraph'], key="cov_up")

    if uploaded_file:
        st.divider()
        st.subheader("File Parsing Settings")
        # Guess separator based on extension
        file_extension = os.path.splitext(uploaded_file.name)[-1].lower()
        default_separator = '\t' # Default to tab
        if file_extension == '.csv':
            default_separator = ','
        elif file_extension in ['.bed', '.bedgraph']:
            default_separator = '\t'

        col1, col2 = st.columns(2)
        with col1:
            # Separator selection with user-friendly names
            separator_options = {',': 'Comma (,)', '\t': 'Tab (\\t)', ' ': 'Whitespace ( )'}
            selected_separator_key = st.selectbox("Column Separator:", list(separator_options.keys()),
                                                  index=list(separator_options.keys()).index(default_separator) if default_separator in separator_options else 1,
                                                  format_func=lambda x: separator_options[x], key="cov_sep")
            # Use regex for whitespace to handle multiple spaces/tabs robustly
            separator_regex = r'\s+' if selected_separator_key == ' ' else selected_separator_key
        with col2:
            # Header detection and comment character
            default_has_header = file_extension not in ['.bed', '.bedgraph'] # BED/BEDGRAPH typically don't have headers
            has_header = st.checkbox("File has header row", value=default_has_header, key="cov_head")
            header_arg = 'infer' if has_header else None
            comment_char = st.text_input("Comment Character (lines starting with this ignored):", '#', max_chars=1, key="cov_comm", help="Leave blank if no comment lines.")

        st.divider()
        st.subheader("Data Loading and Column Selection")
        try:
            # Read the file content
            file_content = uploaded_file.getvalue().decode('utf-8')
            # Use pandas to read the data
            # *** CORRECTED LINE: Removed low_memory=False ***
            df_coverage = pd.read_csv(io.StringIO(file_content),
                                     sep=separator_regex,
                                     header=header_arg,
                                     engine='python', # 'python' engine handles regex separators better
                                     comment=comment_char if comment_char else None)


            if df_coverage.empty:
                st.error("The uploaded file appears to be empty or contains only comments.")
                st.stop()

            # Assign generic column names if no header
            if not has_header:
                num_cols = df_coverage.shape[1]
                if num_cols >= 2:
                    df_coverage.columns = [f'Column_{i+1}' for i in range(num_cols)]
                    st.info(f"No header detected; generic column names assigned: {list(df_coverage.columns)}")
                else:
                    st.error("Could not parse columns. Ensure the correct separator is selected and the file has at least two columns.")
                    st.stop()

            st.write("Data Preview (first 5 rows):")
            st.dataframe(df_coverage.head(), height=150, use_container_width=True)

            # Column selection for position and coverage
            available_columns = list(df_coverage.columns)
            # Attempt to guess columns (case-insensitive matching)
            normalized_columns = {str(c).strip().lower(): c for c in available_columns}
            pos_guess, cov_guess = None, None
            # Common column names to check
            pos_keywords = ['pos', 'position', 'start', 'coordinate', 'loc', 'location']
            cov_keywords = ['depth', 'cov', 'coverage', 'score', 'value', 'count']

            # Find best guess based on keywords (first match wins)
            for key in pos_keywords:
                 if key in normalized_columns: pos_guess = normalized_columns[key]; break
            for key in cov_keywords:
                 if key in normalized_columns: cov_guess = normalized_columns[key]; break

            # Set default indices based on guesses or simple column order
            pos_default_index = available_columns.index(pos_guess) if pos_guess in available_columns else (1 if len(available_columns) > 1 else 0) # Often 2nd col in BED/depth (0-based index 1)
            cov_default_index = available_columns.index(cov_guess) if cov_guess in available_columns else (2 if len(available_columns) > 2 else (1 if len(available_columns) > 1 else 0)) # Often 3rd col in depth (index 2)

            col_s1, col_s2 = st.columns(2)
            with col_s1:
                position_column = st.selectbox("Select Position Column:", available_columns, index=pos_default_index, key="cov_pos")
            with col_s2:
                coverage_column = st.selectbox("Select Coverage Column:", available_columns, index=cov_default_index, key="cov_cov")

            if position_column == coverage_column:
                st.error("Position and Coverage columns cannot be the same. Please select different columns.")
            else:
                st.success(f"Using column **'{position_column}'** for Position and **'{coverage_column}'** for Coverage.")
                st.divider()
                st.subheader("Plotting Options")
                try:
                    # Prepare data for plotting
                    plot_data = df_coverage[[position_column, coverage_column]].copy()
                    plot_data.columns = ['Position', 'Coverage'] # Standardize names

                    # Convert columns to numeric, coercing errors to NaN
                    plot_data['Position'] = pd.to_numeric(plot_data['Position'], errors='coerce')
                    plot_data['Coverage'] = pd.to_numeric(plot_data['Coverage'], errors='coerce')

                    # Handle missing/invalid values
                    rows_before_drop = len(plot_data)
                    plot_data.dropna(subset=['Position', 'Coverage'], inplace=True)
                    rows_after_drop = len(plot_data)
                    if rows_after_drop < rows_before_drop:
                        st.warning(f"Removed {rows_before_drop - rows_after_drop:,} rows with invalid non-numeric data in selected columns.")

                    if plot_data.empty:
                        st.error("No valid numeric data available for plotting after cleaning. Check column selections and file content.")
                    else:
                        # Sort data by position for correct plotting
                        plot_data = plot_data.sort_values(by='Position').reset_index(drop=True)
                        min_pos, max_pos = plot_data['Position'].min(), plot_data['Position'].max()
                        st.write(f"Plotting coverage data from position **{min_pos:,.0f}** to **{max_pos:,.0f}**.")

                        # Plotting customization
                        col_p1, col_p2 = st.columns(2)
                        with col_p1:
                            fill_area = st.checkbox("Fill Area Under Curve", value=True, key="cov_fill")
                            log_scale_y = st.checkbox("Use Logarithmic Y-axis", value=False, key="cov_log", help="Useful if coverage varies greatly.")
                        with col_p2:
                            # Smoothing (Rolling Mean)
                            smoothing_window = st.slider("Smoothing Window Size (0 = None):", min_value=0, max_value=201, value=0, step=2, key="cov_smooth", help="Applies a rolling mean. Use odd numbers for symmetrical window (e.g., 5, 11, 51). 0 disables smoothing.")
                            # Ensure window size is odd if > 0
                            if smoothing_window > 0 and smoothing_window % 2 == 0:
                                smoothing_window += 1
                                st.caption(f"Adjusted smoothing window to odd number: {smoothing_window}")

                        # Generate Plot
                        fig, ax = plt.subplots(figsize=(14, 5)) # Wider figure

                        # Plot raw coverage
                        ax.plot(plot_data['Position'], plot_data['Coverage'], label='Raw Coverage',
                                linewidth=0.7, color='#66fcf1', alpha=0.8) # Cyan, slightly thinner line
                        if fill_area:
                            ax.fill_between(plot_data['Position'], plot_data['Coverage'], alpha=0.25, color='#97e8e1') # Lighter cyan fill

                        # Plot smoothed coverage if window size > 0
                        if smoothing_window > 1: # Need at least 3 for rolling mean to make sense
                            plot_data['Smoothed_Coverage'] = plot_data['Coverage'].rolling(window=smoothing_window, center=True, min_periods=1).mean()
                            ax.plot(plot_data['Position'], plot_data['Smoothed_Coverage'], label=f'Smoothed ({smoothing_window}bp window)',
                                    linewidth=1.8, color='#fca311', linestyle='-') # Orange, thicker line

                        # Apply styling
                        style_plot(fig, ax, title=f"Genome Coverage Plot ({uploaded_file.name})",
                                   xlabel=f"Genomic Position ({position_column})",
                                   ylabel=f"Coverage Depth ({coverage_column}){ ' (Log Scale)' if log_scale_y else ''}")

                        # Set y-axis scale
                        if log_scale_y:
                            ax.set_yscale('log')
                            # Adjust y-min for log scale to avoid issues with zero coverage
                            min_positive_coverage = plot_data[plot_data['Coverage'] > 0]['Coverage'].min()
                            ax.set_ylim(bottom=max(0.1, min_positive_coverage * 0.5) if pd.notna(min_positive_coverage) else 0.1)

                        # Format axes tick labels
                        ax.ticklabel_format(style='plain', axis='x', useOffset=False) # Prevent scientific notation on x-axis
                        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ','))) # Comma separators for y-axis

                        # Add Legend
                        legend = ax.legend(facecolor='#2c3e50', labelcolor='#e0e0e0', framealpha=0.8, fontsize=10)
                        for text in legend.get_texts():
                            text.set_fontfamily('Times New Roman')
                            text.set_fontweight('bold')

                        # Display the plot
                        st.pyplot(fig)

                        # Basic Statistics
                        st.divider()
                        st.subheader("Coverage Statistics (Raw Data)")
                        coverage_stats = plot_data['Coverage'].describe()
                        scol1, scol2, scol3, scol4 = st.columns(4)
                        with scol1: st.metric("Mean Coverage", f"{coverage_stats.get('mean', 0):,.1f}x")
                        with scol2: st.metric("Median Coverage", f"{coverage_stats.get('50%', 0):,.1f}x")
                        with scol3: st.metric("Minimum Coverage", f"{coverage_stats.get('min', 0):,.0f}x")
                        with scol4: st.metric("Maximum Coverage", f"{coverage_stats.get('max', 0):,.0f}x")

                        # Calculate Breadth of Coverage
                        max_cov_int = int(coverage_stats.get('max', 1000))
                        coverage_threshold = st.number_input("Calculate Breadth of Coverage at Depth ≥:",
                                                             min_value=0, max_value=max_cov_int, value=10, step=1, key="cov_b_thresh")
                        if coverage_threshold >= 0:
                            positions_above_threshold = (plot_data['Coverage'] >= coverage_threshold).sum()
                            total_positions = len(plot_data)
                            if total_positions > 0:
                                breadth_percent = (positions_above_threshold / total_positions * 100)
                                st.metric(f"Breadth of Coverage (≥ {coverage_threshold}x)",
                                          f"{breadth_percent:.1f}%",
                                          f"{positions_above_threshold:,} / {total_positions:,} positions")
                            else:
                                st.metric(f"Breadth of Coverage (≥ {coverage_threshold}x)", "N/A", "No data points")

                except Exception as plot_e:
                    st.error("An error occurred during plotting:")
                    st.exception(plot_e)

        except pd.errors.EmptyDataError:
            st.error("File Read Error: The file appears to be empty or could not be parsed with the selected settings.")
        except ValueError as ve:
            # Make the error message more specific to the potential low_memory/engine conflict
            st.error(f"File Read Error: Could not parse the file. Check separator, header, and comment settings. Details: {ve}")
        except Exception as read_e:
            st.error("An unexpected error occurred while reading the file:")
            st.exception(read_e)


elif menu == "Variant Annotation Tool":
    st.header("Simple Variant Substitution Tool")
    st.markdown("Introduce a single nucleotide polymorphism (SNP) into a reference DNA sequence (A, T, C, G only) and observe the change, including the effect on the corresponding codon.")

    input_method = st.radio("Input Reference DNA:",("Paste","FASTA"), key="var_in", horizontal=True)
    ref_sequence, seq_id="","Reference Sequence"
    if input_method=="Paste":
        ref_sequence=st.text_area("Enter Reference DNA Sequence:", height=100, key="var_ref", placeholder="Paste reference DNA (e.g., ATGC...). Only ATCG characters used.")
    else:
        uploaded_file=st.file_uploader("Upload Reference FASTA File:", type=['fasta','fa','fna'], key="var_f")
        if uploaded_file:
            ref_sequence, seq_id = parse_fasta(uploaded_file)
            if ref_sequence:
                st.text_area("Reference Preview:", value=f">{seq_id}\n{ref_sequence[:80]}{'...' if len(ref_sequence)>80 else ''}", height=75, key="var_d", disabled=True)

    # Variant Information Input
    st.divider()
    st.subheader("Variant Details")
    vcol1, vcol2 = st.columns(2)
    with vcol1:
        # Using 1-based indexing for user input, as is common in genomics
        variant_pos = st.number_input("Variant Position (1-based index):", min_value=1, step=1, value=1, key="var_p", help="The position in the sequence where the change occurs.")
    with vcol2:
        variant_base = st.selectbox("New Base (Alternate Allele):", ["A", "T", "C", "G"], key="var_b", help="The nucleotide to introduce at the specified position.")

    st.divider()
    if st.button("Apply Variant and Annotate", key="var_btn"):
        if not ref_sequence:
            st.warning("Please provide a reference DNA sequence.")
        else:
            # Clean reference sequence
            ref_cleaned = "".join(re.findall(r'[ATCG]', ref_sequence.upper()))
            error_occurred = False
            if not ref_cleaned:
                st.error("The reference sequence contains no valid DNA bases (A, T, C, G).")
                error_occurred = True

            ref_length = len(ref_cleaned)
            # Validate variant position (adjust to 0-based index internally)
            zero_based_pos = variant_pos - 1
            if not error_occurred and (zero_based_pos < 0 or zero_based_pos >= ref_length):
                st.error(f"Variant position {variant_pos:,} is outside the valid range of the sequence (1 to {ref_length:,}).")
                error_occurred = True

            if not error_occurred:
                try:
                    # Get the original base at the variant position
                    original_base = ref_cleaned[zero_based_pos]

                    # Check if the variant base is the same as the original
                    if original_base == variant_base:
                        st.info(f"No change needed: The base at position {variant_pos:,} in the reference ('{seq_id}') is already '{original_base}'.")
                        st.write(f"**Reference Sequence ('{seq_id}'):**")
                        st.code(ref_cleaned, language='text')
                    else:
                        # Create the altered sequence
                        alt_list = list(ref_cleaned)
                        alt_list[zero_based_pos] = variant_base
                        alt_sequence = "".join(alt_list)

                        st.subheader("Variant Applied Successfully")
                        st.write(f"**Sequence:** '{seq_id}' | **Variant:** Position {variant_pos:,}")
                        st.success(f"**Change:** Reference Base (`{original_base}`) → Alternate Base (`{variant_base}`)")

                        # Display Reference vs. Altered sequences side-by-side
                        s_col1, s_col2 = st.columns(2)
                        with s_col1:
                            st.write("**Reference Sequence:**")
                            st.code(ref_cleaned, language='text')
                        with s_col2:
                            st.write("**Altered Sequence:**")
                            st.code(alt_sequence, language='text')

                        # Codon Context and Annotation
                        st.divider()
                        st.write("#### Codon Context and Predicted Effect:")
                        # Determine the codon boundaries (0-based)
                        codon_start_index = (zero_based_pos // 3) * 3
                        position_in_codon = (zero_based_pos % 3) + 1 # 1, 2, or 3

                        # Check if a full codon can be extracted
                        if codon_start_index + 3 <= ref_length:
                            ref_codon_str = ref_cleaned[codon_start_index : codon_start_index + 3]
                            alt_codon_str = alt_sequence[codon_start_index : codon_start_index + 3]

                            st.write(f"The variant at position {variant_pos:,} occurs at position **{position_in_codon}** within the codon starting at base **{codon_start_index + 1:,}**.")

                            try:
                                # Translate codons using standard table
                                ref_codon_obj = Seq(ref_codon_str)
                                alt_codon_obj = Seq(alt_codon_str)
                                ref_aa = str(ref_codon_obj.translate(table=1)) # Standard code
                                alt_aa = str(alt_codon_obj.translate(table=1))

                                # Display codon change
                                c_col1, c_col2 = st.columns(2)
                                with c_col1:
                                    st.write("**Reference Codon:**")
                                    st.code(f"Codon: {ref_codon_str}\nAA:    {ref_aa}", language='text')
                                with c_col2:
                                    st.write("**Altered Codon:**")
                                    st.code(f"Codon: {alt_codon_str}\nAA:    {alt_aa}", language='text')

                                # Predict the effect
                                effect = "Unknown"
                                if ref_aa == alt_aa:
                                    effect = "**Silent (Synonymous)**"
                                elif alt_aa == '*':
                                    effect = "**Nonsense (Stop Gained)**"
                                elif ref_aa == '*':
                                    effect = "**Stop-Lost**"
                                else: # ref_aa != alt_aa and neither is stop (or only ref was stop)
                                    effect = "**Missense (Non-synonymous)**"

                                st.markdown(f"**Predicted Effect:** {effect}")

                            except Exception as codon_translate_e:
                                st.warning(f"Could not translate codons or determine effect: {codon_translate_e}")

                        else:
                            st.info("The variant occurs too close to the end of the sequence to determine the full codon context.")

                except IndexError:
                     st.error(f"Internal error: Position {variant_pos:,} seems invalid despite checks.")
                except Exception as e:
                    st.error("An error occurred while applying the variant:")
                    st.exception(e)


elif menu == "Codon Usage Analyzer":
    st.header("Codon Usage Analyzer")
    st.markdown("Analyze the frequency of codon usage within a provided Coding DNA Sequence (CDS). The sequence must contain only A, T, C, G and its length must be a multiple of 3.")

    input_method = st.radio("Input CDS:",("Paste","FASTA"), key="cod_in", horizontal=True)
    cds_sequence, seq_id="","CDS Sequence"
    if input_method=="Paste":
        cds_sequence=st.text_area("Enter Coding DNA Sequence (CDS):", height=120, key="cod_seq", placeholder="Paste CDS here (e.g., ATGCGT...). Must be multiple of 3, ATCG only.")
    else:
        uploaded_file=st.file_uploader("Upload FASTA File (containing CDS):", type=['fasta','fa','fna'], key="cod_f")
        if uploaded_file:
            cds_sequence, seq_id = parse_fasta(uploaded_file)
            if cds_sequence:
                st.text_area("CDS Preview:", value=f">{seq_id}\n{cds_sequence[:80]}{'...' if len(cds_sequence)>80 else ''}", height=75, key="cod_d", disabled=True)

    st.divider()
    if st.button("Analyze Codon Usage", key="cod_btn"):
        if not cds_sequence:
            st.warning("Please provide a Coding DNA Sequence (CDS).")
        else:
            # Clean and validate sequence
            cds_cleaned = "".join(re.findall(r'[ATCG]', cds_sequence.upper()))
            error_occurred = False
            if not cds_cleaned:
                st.error("The provided sequence contains no valid DNA bases (A, T, C, G).")
                error_occurred = True
            elif len(cds_cleaned) % 3 != 0:
                st.error(f"Sequence length ({len(cds_cleaned):,}) is not a multiple of 3. Please provide a valid CDS.")
                error_occurred = True

            if not error_occurred:
                try:
                    st.subheader(f"Codon Usage Analysis for: '{seq_id}'")
                    # Split sequence into codons
                    codons_list = [cds_cleaned[i:i+3] for i in range(0, len(cds_cleaned), 3)]
                    total_codons = len(codons_list)
                    codon_counts = Counter(codons_list)
                    st.metric("Total Number of Codons Analyzed", f"{total_codons:,}")

                    # Get standard codon table information
                    try:
                        standard_table = CodonTable.unambiguous_dna_by_id[1]
                        all_possible_codons = list(standard_table.forward_table.keys()) + standard_table.stop_codons
                        aa_map = standard_table.forward_table # Codon -> AA
                        stop_codons = set(standard_table.stop_codons)
                    except Exception as table_e:
                        st.error(f"Failed to load standard codon table: {table_e}")
                        st.stop()

                    # Calculate frequencies
                    usage_data = []
                    for codon in sorted(all_possible_codons):
                        count = codon_counts.get(codon, 0)
                        frequency = (count / total_codons * 100) if total_codons > 0 else 0
                        amino_acid = aa_map.get(codon, 'Stop' if codon in stop_codons else '?') # Assign AA or 'Stop'
                        usage_data.append({
                            "Codon": codon,
                            "Amino Acid (AA)": amino_acid,
                            "Count": count,
                            "Frequency (%)": frequency
                        })

                    usage_df = pd.DataFrame(usage_data)
                    st.write("#### Codon Usage Table (Standard Genetic Code):")
                    st.dataframe(usage_df.style.format({'Count':'{:,}', 'Frequency (%)':'{:.1f}%'}),
                                 hide_index=True, use_container_width=True)

                    # Visualization Options
                    st.divider()
                    st.subheader("Codon Usage Visualization")
                    plot_choice = st.selectbox("Select Visualization Type:",
                                               ["Frequency of All Codons", "Relative Usage for Specific Amino Acid"],
                                               key="cod_plot")

                    # Filter out codons that are not present in the sequence for plotting
                    present_codons_df = usage_df[usage_df['Count'] > 0].copy()

                    if present_codons_df.empty:
                        st.info("No codons found in the sequence to visualize.")
                    else:
                        if plot_choice == "Frequency of All Codons":
                            plot_df_all = present_codons_df.sort_values("Codon")
                            fig, ax = plt.subplots(figsize=(16, 6)) # Wide figure for all codons

                            # Color bars by amino acid type (simple grouping)
                            unique_aas = sorted(list(plot_df_all['Amino Acid (AA)'].unique()))
                            # Use a perceptually uniform colormap if many AAs, or define manually
                            cmap = plt.get_cmap('tab20', len(unique_aas)) # 'tab20' has 20 distinct colors
                            aa_colors = {aa: cmap(i) for i, aa in enumerate(unique_aas)}
                            bar_colors = [aa_colors.get(aa, '#888888') for aa in plot_df_all['Amino Acid (AA)']]

                            ax.bar(plot_df_all['Codon'], plot_df_all['Frequency (%)'], color=bar_colors)

                            style_plot(fig, ax, title=f"Codon Frequency Distribution for '{seq_id}'",
                                       xlabel="Codon", ylabel="Frequency (%)")
                            # Rotate x-axis labels for better readability
                            ax.tick_params(axis='x', rotation=90, labelsize=9)
                            st.pyplot(fig)

                        elif plot_choice == "Relative Usage for Specific Amino Acid":
                             # Get list of AAs/Stop present in the sequence
                             available_aas = sorted(list(present_codons_df['Amino Acid (AA)'].unique()))
                             if not available_aas:
                                 st.info("No amino acids encoded by the provided sequence.")
                             else:
                                 selected_aa = st.selectbox("Select Amino Acid or 'Stop':", available_aas, key="cod_aa")

                                 # Filter data for the selected AA/Stop
                                 aa_specific_df = present_codons_df[present_codons_df['Amino Acid (AA)'] == selected_aa].copy()
                                 total_count_for_aa = aa_specific_df['Count'].sum()

                                 if total_count_for_aa > 0 and not aa_specific_df.empty:
                                     # Calculate relative frequency within this AA's codons
                                     aa_specific_df['Relative Freq (%)'] = (aa_specific_df['Count'] / total_count_for_aa * 100)
                                     aa_specific_df = aa_specific_df.sort_values("Codon")

                                     # Plot relative usage
                                     fig, ax = plt.subplots(figsize=(max(6, len(aa_specific_df)*1.5), 5)) # Adjust width based on # codons
                                     ax.bar(aa_specific_df['Codon'], aa_specific_df['Relative Freq (%)'], color='#fca311') # Use theme's orange

                                     style_plot(fig, ax, title=f"Relative Codon Usage for '{selected_aa}' in '{seq_id}' (Total Count: {total_count_for_aa:,})",
                                                xlabel="Codon", ylabel="Relative Frequency (%)")
                                     ax.set_ylim(0, 105) # Set y-limit slightly above 100%
                                     st.pyplot(fig)

                                     # Show table for selected AA
                                     st.write(f"**Usage Details for '{selected_aa}':**")
                                     st.dataframe(aa_specific_df[['Codon', 'Count', 'Relative Freq (%)']].style.format({'Count':'{:,}', 'Relative Freq (%)':'{:.1f}%'}),
                                                  hide_index=True, use_container_width=True)
                                 else:
                                     st.info(f"No codons encoding '{selected_aa}' were found in this sequence.")

                except Exception as e:
                    st.error("An error occurred during codon usage analysis:")
                    st.exception(e)


# --- Footer ---
# Use st.markdown with unsafe_allow_html=True to apply the CSS class
st.markdown("---") # Visual separator line
st.markdown('<div class="footer">Copyright © 2024 CCDB Tools | For Educational and Informational Purposes Only</div>', unsafe_allow_html=True)
