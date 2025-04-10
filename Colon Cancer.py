# -*- coding: utf-8 -*-
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

# --- Custom CSS Styling (LIGHT BLUE THEME + Times New Roman + Bold) ---
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
    /* Sidebar specific text elements */
    [data-testid="stSidebar"] h1,
    [data-testid="stSidebar"] h2,
    [data-testid="stSidebar"] h3,
    [data-testid="stSidebar"] label,
    [data-testid="stSidebar"] .stRadio div[role="radiogroup"] > label, /* Sidebar Radio OPTIONS text */
    [data-testid="stSidebar"] [data-testid="stExpander"] summary,
    /* Ensure plot text elements inherit if possible (often needs Python-side control) */
    .plot-container text, /* General plot text target */
    /* Target common streamlit elements containing text */
    [data-testid="stMarkdownContainer"], [data-testid="stText"],
    [data-testid="stMetricLabel"], [data-testid="stMetricValue"], [data-testid="stMetricDelta"]
    {
        font-family: 'Times New Roman', Times, serif !important;
        font-weight: bold !important;
        color: #1E2A3A; /* Default DARK text color */
    }

    /* Remove underlines globally, except on hover for links */
    a, a:link, a:visited {
        text-decoration: none !important;
        color: #0056b3 !important; /* Dark Blue link color */
    }
    a:hover {
        text-decoration: underline !important;
        color: #003d80 !important; /* Darker blue on hover */
    }

    /* --- Specific Element Styling (LIGHT BLUE THEME) --- */

    /* Body background - Subtle Light Blue Gradient */
    body {
        background-color: #e6f7ff; /* Fallback light blue */
        background-image: linear-gradient(to bottom right, #e6f7ff, #cceeff); /* Light blue gradient */
        background-attachment: fixed;
        min-height: 100vh;
    }

    /* Main app container - Off-white/Very Light Grey */
    .stApp {
        background-color: #f8f9fa; /* Very light grey/off-white */
        padding: 30px; border-radius: 10px;
        border: 1px solid #b3d9ff; /* Lighter blue border */
        box-shadow: 0 5px 15px rgba(0, 0, 0, 0.1); /* Lighter shadow */
        padding-bottom: 80px; /* Space for footer */
        color: #1E2A3A; /* Default dark text color for app container */
    }

    /* Sidebar Styling - Slightly different light shade */
    [data-testid="stSidebar"] {
        background-color: #e0f2f7; /* Slightly different light blue/cyan */
        border-right: 2px solid #b3d9ff; /* Light blue separator */
        padding-top: 2rem;
    }
    /* Sidebar Header & Widget Labels */
    [data-testid="stSidebar"] h1,
    [data-testid="stSidebar"] h2,
    [data-testid="stSidebar"] h3,
    [data-testid="stSidebar"] label
    {
         color: #004080 !important; /* Darker Blue */
    }
    /* Sidebar Radio Button Options Text */
    [data-testid="stSidebar"] .stRadio div[role="radiogroup"] > label {
         color: #333333 !important; /* Dark Grey for options */
    }
     /* Sidebar expander */
    [data-testid="stSidebar"] [data-testid="stExpander"] summary {
        color: #0056b3; /* Dark Blue */
    }

    /* Headings - Dark Blue */
    h1, h2, h3, h4, h5, h6 { color: #004080; }
    h1 { border-bottom: 2px solid #80ccff; padding-bottom: 10px; margin-bottom: 20px; } /* Lighter blue border */
    h2 { padding-bottom: 5px; margin-top: 30px; margin-bottom: 15px;}
    h3 { margin-top: 25px; margin-bottom: 10px; color: #0056b3; } /* Slightly lighter dark blue for H3 */


    /* Paragraph text - Dark Grey/Black */
    p, div, span, li, [data-testid="stMarkdownContainer"] p {
        color: #333333; line-height: 1.7; font-size: 16px;
    }

     /* Labels for Widgets - Dark Blue */
    .stTextInput>label, .stTextArea>label, .stSelectbox>label, .stNumberInput>label,
    .stFileUploader>label, .stSlider>label, .stRadio>label, .stCheckbox>label,
    .stMultiSelect>label, [data-testid="stMetricLabel"]
    {
        color: #0056b3 !important; font-size: 1.1em; margin-bottom: 5px;
    }

    /* Input fields, Select boxes, Text Areas */
    .stTextInput>div>div>input,
    .stTextArea>div>textarea,
    .stSelectbox>div>div>div,
    .stMultiSelect>div>div
    {
        border: 1px solid #ced4da; /* Standard grey border */
        border-radius: 5px; background-color: #ffffff; /* White background */
        color: #212529; /* Dark text */
        padding: 9px 12px;
    }
    .stTextArea>div>textarea { min-height: 150px; }
    .stSelectbox [data-baseweb="select"] > div { background-color: #ffffff; } /* Ensure dropdown bg matches */
    .stMultiSelect [data-baseweb="tag"] { background-color: #007bff; color: #ffffff; } /* Primary blue tags, white text */

    /* Input Placeholders */
    .stTextInput>div>div>input::placeholder, .stTextArea>div>textarea::placeholder {
        color: #6c757d; opacity: 0.8; /* Grey placeholder */
        font-weight: normal !important;
    }

    /* Buttons - Good Contrast Blue */
    .stButton>button {
        border-radius: 5px; background-color: #007bff; color: #ffffff; /* Primary Blue bg, WHITE text */
        border: 1px solid #0069d9; padding: 10px 22px;
        transition: background-color 0.3s ease, transform 0.1s ease, box-shadow 0.2s ease;
        font-size: 1.05em; box-shadow: 0 2px 4px rgba(0,0,0,0.15);
    }
    .stButton>button:hover {
        background-color: #0056b3; color: #ffffff; transform: translateY(-2px); /* Darker blue on hover */
        box-shadow: 0 4px 8px rgba(0, 123, 255, 0.3); border-color: #0056b3;
    }
    .stButton>button:active {
        background-color: #004085; transform: translateY(0px); /* Even darker blue */
        box-shadow: inset 0 2px 4px rgba(0,0,0,0.2); border-color: #00376b;
    }
     /* Disabled button styling */
    .stButton>button:disabled {
        background-color: #adb5bd; color: #6c757d; border-color: #adb5bd; /* Greyed out */
        cursor: not-allowed;
    }

    /* Tables - Consistent Light Theme */
    table {
        border-collapse: separate; border-spacing: 0; width: 100%; margin-bottom: 1.5rem;
        background-color: #ffffff; /* White background */
        border: 1px solid #dee2e6; /* Lighter grey border */
        border-radius: 6px; overflow: hidden;
    }
    table th {
        text-align: left; padding: 12px 15px; background-color: #e9ecef; /* Light grey header bg */
        color: #212529; /* Dark header text */
        border-bottom: 2px solid #dee2e6;
    }
    table td {
        vertical-align: middle; text-align: left; padding: 12px 15px;
        border-top: 1px solid #dee2e6; /* Internal border */
        color: #333333; /* Dark cell text */
    }
    /* Table links - Use global link color */
    table a, table a:link, table a:visited {
        color: #0056b3 !important;
        font-weight: bold !important;
        text-decoration: none !important;
    }
    table a:hover {
        color: #003d80 !important;
        text-decoration: underline !important;
    }
    /* Zebra striping */
    table tbody tr:nth-of-type(odd) { background-color: #ffffff; } /* White row */
    table tbody tr:nth-of-type(even) { background-color: #f8f9fa; } /* Slightly off-white even row */
    table tbody tr:hover { background-color: #e9ecef; } /* Light grey hover row */

    /* Alerts - Themed for Light Background */
    .stAlert { border-radius: 5px; padding: 15px; font-size: 1.05em; border: 1px solid; }
    .stAlert > div[data-testid="stMarkdownContainer"] > p {
        color: inherit !important; /* Inherit text color from parent alert div */
        font-weight: bold !important;
        font-family: 'Times New Roman', Times, serif !important;
    }
    div[data-baseweb="alert"][role="alert"] { border-left: 5px solid; }
    /* Info */
    div[data-testid="stInfo"] { border-color: #17a2b8; background-color: #d1ecf1; color: #0c5460; } /* Cyan */
    /* Success */
    div[data-testid="stSuccess"] { border-color: #28a745; background-color: #d4edda; color: #155724; } /* Green */
    /* Warning */
    div[data-testid="stWarning"] { border-color: #ffc107; background-color: #fff3cd; color: #856404; } /* Yellow */
    /* Error */
    div[data-testid="stError"] { border-color: #dc3545; background-color: #f8d7da; color: #721c24; } /* Red */

    /* Footer - Fixed at bottom */
    .footer {
        position: fixed; left: 0; bottom: 0; width: 100%;
        background-color: #e0f2f7; /* Match sidebar background */
        /* background-image: linear-gradient(to right, #e0f2f7, #cceeff); */ /* Optional gradient */
        color: #004080; /* Dark Blue text */
        text-align: center; padding: 12px; font-size: 14px;
        border-top: 3px solid #b3d9ff; /* Light blue top border */
        z-index: 1000;
        font-family: 'Times New Roman', Times, serif !important; font-weight: bold !important;
    }

    /* Ensure Plotly/Matplotlib backgrounds are transparent or match theme */
    .plot-container.plotly, .stPlotlyChart, .stVegaLiteChart {
        background-color: transparent !important;
    }
    /* Style Streamlit expander */
    [data-testid="stExpander"] {
        border: 1px solid #ced4da;
        border-radius: 5px;
        background-color: #ffffff; /* White background */
        margin-bottom: 1rem;
    }
    [data-testid="stExpander"] summary {
        background-color: #e9ecef; /* Light grey header */
        color: #0056b3; /* Dark blue text */
        border-radius: 5px 5px 0 0;
        padding: 10px 15px;
    }
    [data-testid="stExpander"] summary:hover {
        background-color: #d6dade; /* Slightly darker grey */
    }
    /* Expander content */
    [data-testid="stExpander"] [data-testid="stVerticalBlock"] {
        padding: 15px;
        background-color: #f8f9fa; /* Off-white content background */
    }


    /* Style Streamlit Divider */
    hr {
        border-top: 1px solid #dee2e6; /* Light grey divider */
        margin-top: 1rem;
        margin-bottom: 1rem;
    }

    /* Style Streamlit Metric */
    [data-testid="stMetric"] {
        background-color: #e9ecef; /* Light grey background */
        border: 1px solid #dee2e6;
        border-radius: 5px;
        padding: 15px;
        margin-bottom: 10px;
    }
    [data-testid="stMetricValue"] {
        color: #1E2A3A !important; /* Dark metric value */
        font-size: 1.6em !important;
    }
    [data-testid="stMetricDelta"] {
         color: #333333 !important; /* Dark grey delta */
    }

     /* Ensure code blocks match theme */
    .stCodeBlock code {
        background-color: #e9ecef; /* Light grey background */
        color: #212529; /* Dark text */
        border: 1px solid #dee2e6;
        border-radius: 5px;
        padding: 10px;
        font-family: monospace !important; /* Keep monospace for code */
        font-weight: normal !important; /* Code usually not bold */
        font-size: 0.95em;
        white-space: pre-wrap; /* Allow wrapping */
        word-wrap: break-word;
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
    # Use the dark blue link color defined in CSS (#0056b3)
    link_style = 'style="color:#0056b3 !important;"'

    try:
        # 1. Check if value is already a full URL
        if s_value.startswith("http"):
             display_url = s_value if len(s_value) < 50 else s_value[:47] + "..."
             return f'<a href="{s_value}" target="_blank" {link_style} title="{s_value}">{display_url}</a>'

        # 2. Handle specific ID types that need URL construction
        elif col_name_std == 'Pubmed Id':
            if s_value.replace('.', '', 1).isdigit():
                url = f"https://pubmed.ncbi.nlm.nih.gov/{s_value}"
                return f'<a href="{url}" target="_blank" {link_style}>{s_value}</a>'

        elif col_name_std == 'Doi Id':
            if '/' in s_value: # Simple check for potential DOI structure
                encoded_doi = urllib.parse.quote(s_value, safe='/:')
                url = f"https://doi.org/{encoded_doi}"
                return f'<a href="{url}" target="_blank" {link_style}>{s_value}</a>'

        elif col_name_std == 'Uniprot':
             if re.match(r'^[A-Z0-9_.\-]+$', s_value, re.IGNORECASE): # Check for UniProt ID format
                 url = f"https://www.uniprot.org/uniprotkb/{s_value}/entry"
                 return f'<a href="{url}" target="_blank" {link_style}>{s_value}</a>'

        # 3. For other columns, display as plain text since we already checked for 'http'
        else:
            return s_value

    except Exception: # Catch any unexpected formatting error
        return s_value

# --- Helper function for Matplotlib Plot Styling ---
def style_plot(fig, ax, title="", xlabel="", ylabel=""):
    """Applies consistent LIGHT theme styling (Times New Roman Bold) to a Matplotlib plot."""
    fig.patch.set_facecolor('#f8f9fa') # Match app background (off-white)
    fig.patch.set_alpha(1.0)
    ax.set_facecolor('#ffffff') # White plot area background

    font_props = {'family': 'Times New Roman', 'weight': 'bold'}
    text_color = '#1E2A3A'; accent_color = '#004080'
    ax.set_title(title, color=accent_color, fontsize=14, **font_props)
    ax.set_xlabel(xlabel, color=text_color, fontsize=12, **font_props)
    ax.set_ylabel(ylabel, color=text_color, fontsize=12, **font_props)

    tick_color = '#333333'
    ax.tick_params(axis='x', colors=tick_color, labelsize=10)
    ax.tick_params(axis='y', colors=tick_color, labelsize=10)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontfamily('Times New Roman'); label.set_fontweight('bold'); label.set_color(tick_color)

    ax.grid(True, linestyle=':', alpha=0.6, color='#cccccc')
    for spine in ax.spines.values(): spine.set_edgecolor('#dee2e6')
    plt.tight_layout()

# --- Main App Title ---
st.title("🧬 CCDB: Colon Cancer Database and Bioinformatics Tools 🔬")

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


# --- Page Content based on Menu Selection ---

if menu == "Home":
    st.header("Welcome to the Colon Cancer Data Resource & Toolbox")
    st.markdown( # Content unchanged
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
    st.image( # Content unchanged
        "https://img.freepik.com/free-vector/colorectal-cancer-crc-infographic-education_1308-47971.jpg",
        caption="Colorectal Cancer Infographic (Source: Freepik - Illustrative purposes)",
        use_container_width=True
    )
    st.markdown("""<a href="https://www.cancer.org/cancer/types/colon-rectal-cancer.html" target="_blank" style="color:#0056b3;">Learn more about Colorectal Cancer from the American Cancer Society</a>""", unsafe_allow_html=True)
    st.divider()
    st.subheader("Feedback")
    feedback = st.text_area("Share your feedback about this application:", key="feedback_home_page", placeholder="Your thoughts, suggestions, or issues...")
    if st.button("Submit Feedback", key="submit_feedback_home_page"):
        if feedback: st.success("Thank you for your feedback!")
        else: st.warning("Feedback cannot be empty.")


elif menu == "Colon Cancer Database Search":
    st.header("Colon Cancer Gene Database Search")
    st.markdown(f"Search for information on genes associated with Colorectal Cancer using a local database file (`{DATA_PATH}`). Please enter a gene symbol.")
    if not os.path.exists(DATA_PATH):
        st.error(f"**Database File Not Found:** `{DATA_PATH}` missing."); st.stop()
    query = st.text_input("Enter Gene Symbol:", key="gene_search_input", placeholder="e.g., APC, KRAS, MLH1")
    if st.button("Search Database", key="gene_search_button"):
        if not query: st.warning("Please enter a gene symbol.")
        else:
            try:
                data = pd.read_excel(DATA_PATH); data.columns = data.columns.str.strip().str.title()
                gene_col = None; possible_gene_cols = ['Gene Symbol', 'Gene Name', 'Symbol', 'Gene']
                for col in possible_gene_cols:
                    if col in data.columns: gene_col = col; break
                if not gene_col: st.error(f"Gene symbol column not found. Found: {list(data.columns)}"); st.stop()
                data[gene_col] = data[gene_col].astype(str).str.strip(); search_query = query.strip()
                results = data[data[gene_col].str.fullmatch(search_query, case=False, na=False)].copy()
                if results.empty:
                     st.info(f"No exact match for '{search_query}'. Trying broader search..."); results = data[data[gene_col].str.contains(search_query, case=False, na=False, regex=False)].copy()
                if not results.empty:
                    st.success(f"Found {len(results)} result(s) containing '{query}'."); st.write("### Search Results:")
                    formatted_results = results.copy()
                    # UPDATED list of columns where link formatting should be attempted if value looks like a URL
                    link_potential_columns = [
                        'Pubmed Id', 'Doi Id', 'Uniprot', 'Blast', 'Conserved Domain',
                        'Link', 'Url', 'Reference', 'Entry', 'History', 'Variant Viewer',
                        'Feature Viewer', 'Genomic Coordinates', 'Publications',
                        'External Links', 'Sequence'
                        ]
                    link_potential_columns.extend([col for col in formatted_results.columns if ('link' in col.lower() or 'url' in col.lower()) and col not in link_potential_columns])
                    link_potential_columns = list(set(link_potential_columns))
                    for col in formatted_results.columns:
                        if col in link_potential_columns: # Check if column name suggests a link
                            try: formatted_results[col] = formatted_results[col].apply(lambda x: format_link(x, col)) # format_link handles URL check
                            except Exception as apply_e: st.warning(f"Link formatting failed for '{col}': {apply_e}.")
                    html_table = formatted_results.to_html(escape=False, index=False, na_rep='-', justify='left', classes=['st-table'])
                    st.write(html_table, unsafe_allow_html=True)
                else: st.warning(f"No results found matching/containing '{query}' in '{gene_col}'.")
            except ImportError: st.error("`openpyxl` needed for Excel. `pip install openpyxl`.")
            except FileNotFoundError: st.error(f"Database file `{DATA_PATH}` not found.")
            except Exception as e: st.error(f"Database search error: {e}"); st.exception(e)

# --- DNA Sequence Analysis (Unchanged logic, styling updates via CSS/style_plot) ---
elif menu == "DNA Sequence Analysis":
    st.header("DNA Sequence Analysis")
    st.markdown("Analyze basic properties of a DNA sequence (A/T/C/G/N).")
    input_method = st.radio("Input Method:", ("Paste", "FASTA"), key="dna_input", horizontal=True)
    sequence, seq_id = "", "Pasted Sequence"
    if input_method == "Paste": sequence = st.text_area("Enter DNA Sequence:", height=150, key="dna_seq_in", placeholder="Paste DNA sequence...")
    else: uploaded_file = st.file_uploader("Upload FASTA File:", type=['fasta','fa','fna'], key="dna_fasta_up"); sequence, seq_id = parse_fasta(uploaded_file) if uploaded_file else (None, None)
    if sequence: st.text_area("Sequence Preview:", value=f">**{seq_id}**\n{sequence[:100]}{'...' if len(sequence) > 100 else ''}", height=100, key="dna_disp", disabled=True)
    st.divider()
    if st.button("Analyze DNA", key="dna_analyze_btn"):
        if not sequence: st.warning("Please provide DNA sequence.")
        else:
            seq_clean = "".join(re.findall(r'[ATCGN]', sequence.upper()))
            if not seq_clean: st.error("Input contains no valid DNA bases (A/T/C/G/N).")
            else:
                try:
                    st.subheader(f"Analysis Results for: '**{seq_id}**'")
                    seq_obj = Seq(seq_clean); length = len(seq_obj); st.metric("Total Length", f"{length:,} bp")
                    seq_no_n = seq_clean.replace('N', ''); len_no_n = len(seq_no_n)
                    if len_no_n > 0: gc_val = gc_fraction(seq_no_n)*100; at_val = 100.0 - gc_val; st.metric("GC% (no N)", f"{gc_val:.1f}%"); st.metric("AT% (no N)", f"{at_val:.1f}%")
                    else: st.info("GC/AT% not calculated (only N?).")
                    st.divider(); st.write("#### Base Composition (incl. N):")
                    composition = Counter(seq_clean); comp_data = [{"Base": b, "Count": c, "Percentage (%)": (c/length*100) if length > 0 else 0} for b, c in sorted(composition.items())]
                    comp_df = pd.DataFrame(comp_data); st.dataframe(comp_df.style.format({"Count":"{:,}", "Percentage (%)":"{:.1f}%"}), hide_index=True, use_container_width=True)
                    plot_comp = {k: v for k, v in composition.items() if k in STANDARD_DNA_BASES}
                    if plot_comp:
                        fig, ax = plt.subplots(figsize=(6, 4)); bases = sorted(plot_comp.keys()); counts = [plot_comp[b] for b in bases]
                        base_colors = {'A': '#007bff', 'T': '#ff7f0e', 'C': '#2ca02c', 'G': '#d62728'} # Light theme colors
                        bar_colors = [base_colors.get(b, '#888888') for b in bases]; ax.bar(bases, counts, color=bar_colors)
                        style_plot(fig, ax, title="Base Counts (no N)", xlabel="Base", ylabel="Count"); ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ','))); st.pyplot(fig)
                    elif 'N' in composition and len(composition) == 1: st.info("Plot omitted (only N).")
                    elif not plot_comp: st.info("Plot omitted (no standard bases).")
                    st.divider(); st.write("#### Basic ORF Finder"); min_orf_aa = st.number_input("Min ORF Length (AA):", 10, 500, 30, 5, key="orf_len"); orfs = []
                    st.info("Searching ORFs (M to *) on forward strand...");
                    with st.spinner("Scanning..."):
                         for frame in range(3):
                            try: translation = str(seq_obj[frame:].translate(table=1, to_stop=True, cds=False))
                            except: translation = "" # Ignore translation errors for simple ORF scan
                            if not translation: continue
                            current_pos = 0
                            while 'M' in translation[current_pos:]:
                                start_idx = translation.find('M', current_pos)
                                if start_idx == -1: break
                                potential_orf = translation[start_idx:]
                                if len(potential_orf) >= min_orf_aa:
                                    start_dna = frame + (start_idx * 3) + 1; end_dna = start_dna + (len(potential_orf) * 3) - 1
                                    if end_dna <= length and not any(o['Start (DNA)'] == start_dna and o['Frame'] == frame+1 for o in orfs):
                                        orfs.append({"Frame": frame+1, "Start (DNA)": start_dna, "End (DNA)": end_dna, "Length (AA)": len(potential_orf), "Protein (Preview)": potential_orf[:40] + ("..." if len(potential_orf)>40 else "")})
                                current_pos = start_idx + 1
                    if orfs: st.success(f"Found {len(orfs)} ORF(s) ≥ {min_orf_aa} AA."); orf_df = pd.DataFrame(orfs).sort_values(["Frame", "Start (DNA)"]); st.dataframe(orf_df.style.format({"Start (DNA)":"{:,}", "End (DNA)":"{:,}", "Length (AA)":"{:,}"}), hide_index=True, use_container_width=True)
                    else: st.info(f"No ORFs ≥ {min_orf_aa} AA found.")
                except Exception as e: st.error(f"DNA analysis error: {e}"); st.exception(e)

# --- Protein Sequence Analysis (Unchanged logic, styling updates via CSS/style_plot) ---
elif menu == "Protein Sequence Analysis":
    st.header("Protein Sequence Analysis")
    st.markdown("Analyze biochemical properties of a protein sequence.")
    input_method = st.radio("Input Method:", ("Paste", "FASTA"), key="prot_in", horizontal=True)
    sequence, seq_id = "", "Pasted Sequence"
    if input_method == "Paste": sequence = st.text_area("Enter Protein Sequence:", height=150, key="prot_seq", placeholder="Paste protein sequence...")
    else: uploaded_file = st.file_uploader("Upload FASTA File:", type=['fasta','fa','faa'], key="prot_fasta"); sequence, seq_id = parse_fasta(uploaded_file) if uploaded_file else (None, None)
    if sequence: st.text_area("Sequence Preview:", value=f">**{seq_id}**\n{sequence[:100]}{'...' if len(sequence) > 100 else ''}", height=100, key="prot_disp", disabled=True)
    st.divider()
    if st.button("Analyze Protein", key="prot_analyze_btn"):
        if not sequence: st.warning("Please provide protein sequence.")
        else:
            seq_clean = "".join(sequence.split()).upper(); all_chars = set(seq_clean); standard_set = set(STANDARD_AA); ambiguous_set = set(AMBIGUOUS_AA)
            non_std = all_chars - standard_set - ambiguous_set; ambiguous_present = all_chars.intersection(ambiguous_set)
            seq_for_calc = "".join(c for c in seq_clean if c in standard_set)
            if non_std: st.warning(f"Ignoring unknown: `{'`, `'.join(sorted(list(non_std)))}`")
            if ambiguous_present: st.warning(f"Ignoring ambiguous: `{'`, `'.join(sorted(list(ambiguous_present)))}`")
            if not seq_for_calc: st.error("No standard amino acids found.")
            else:
                try:
                    st.subheader(f"Analysis Results for: '**{seq_id}**'")
                    pa = ProteinAnalysis(seq_for_calc)
                    col1, col2 = st.columns(2)
                    with col1: st.metric("Input Len", f"{len(seq_clean):,} aa"); st.metric("Analyzed Len", f"{len(seq_for_calc):,} aa"); st.metric("Mol. Weight", f"{pa.molecular_weight():,.1f} Da"); st.metric("pI", f"{pa.isoelectric_point():.2f}")
                    with col2:
                        gravy = pa.gravy(); hydro = "Hydrophobic" if gravy > 0 else ("Hydrophilic" if gravy < 0 else "Neutral"); st.metric("GRAVY", f"{gravy:.3f}", delta=hydro, delta_color="off")
                        instab = pa.instability_index(); stability = "Stable" if instab < 40 else "Unstable"; st.metric("Instability", f"{instab:.2f}", delta=stability, delta_color="normal" if stability == "Stable" else "inverse")
                    st.divider(); st.write("#### AA Composition (Full Input):")
                    composition = Counter(seq_clean); comp_data = [{"AA": aa, "Count": c, "%": (c/len(seq_clean)*100) if seq_clean else 0} for aa, c in sorted(composition.items())]
                    comp_df = pd.DataFrame(comp_data); st.dataframe(comp_df.style.format({'Count':'{:,}', '%':'{:.1f}%'}), hide_index=True, use_container_width=True)
                    st.divider(); st.write("#### Predicted Secondary Structure:")
                    try:
                        helix, turn, sheet = pa.secondary_structure_fraction(); labels = 'Helix', 'Turn', 'Sheet'; sizes = [helix*100, turn*100, sheet*100]
                        if sum(sizes) > 0.1:
                            colors = ['#007bff', '#ff7f0e', '#2ca02c'] # Light theme colors
                            fig, ax = plt.subplots(figsize=(6, 4)); wedges, texts, autotexts = ax.pie(sizes, labels=labels, startangle=90, colors=colors, autopct='%1.1f%%', wedgeprops={'edgecolor': '#ffffff', 'linewidth': 1.5}, textprops={'color': '#333333'})
                            for autotext in autotexts: autotext.set_color('#ffffff'); autotext.set_fontfamily('Times New Roman'); autotext.set_fontweight('bold'); autotext.set_fontsize(10)
                            for text in texts: text.set_fontfamily('Times New Roman'); text.set_fontweight('bold')
                            ax.axis('equal'); style_plot(fig, ax, title="Predicted Secondary Structure"); ax.set_xlabel(""); ax.set_ylabel(""); st.pyplot(fig)
                        else: st.info("Structure fractions near zero.")
                    except Exception as ss_e: st.warning(f"Structure prediction failed: {ss_e}")
                except Exception as e: st.error(f"Protein analysis error: {e}"); st.exception(e)

# --- (Remaining tool sections - unchanged logic, styling updates via CSS/style_plot) ---

elif menu == "Primer Design":
    st.header("Basic Primer Design Tool")
    st.markdown("Generates simple forward/reverse primers from DNA template ends (A/T/C/G only) & estimates Tm.")
    st.warning("⚠️ **Educational Tool:** Basic selection, no checks for specificity, secondary structures, etc. Use dedicated tools for experiments.", icon="🔬")
    input_method = st.radio("Input:", ("Paste", "FASTA"), key="p_in", horizontal=True)
    sequence, seq_id = "", "Pasted Sequence"
    if input_method == "Paste": sequence = st.text_area("DNA Template:", height=150, key="p_seq", placeholder="Paste DNA sequence...")
    else: uploaded_file = st.file_uploader("FASTA:", type=['fasta','fa','fna'], key="p_fasta"); sequence, seq_id = parse_fasta(uploaded_file) if uploaded_file else (None, None)
    if sequence: st.text_area("Preview:", value=f">**{seq_id}**\n{sequence[:100]}{'...' if len(sequence) > 100 else ''}", height=100, key="p_disp", disabled=True)
    st.divider(); st.subheader("Parameters")
    col1, col2 = st.columns(2)
    with col1: primer_len = st.slider("Primer Length (bp):", 15, 35, 20, 1, key="p_len")
    with col2: tm_method = st.selectbox("Tm Method:", ["Tm_NN", "Tm_GC", "Tm_Wallace"], 0, key="p_tm_m"); dna_nM, salt_mM, mg_mM, dNTPs_mM = 50.0, 50.0, 0.0, 0.0
        if tm_method == "Tm_NN":
            with st.expander("Tm_NN Params"):
                dna_nM = st.number_input("Primer (nM):", 1.0, 1000.0, 50.0, 1.0, key="p_dna_c", format="%.1f"); salt_mM = st.number_input("Na+ (mM):", 1.0, 200.0, 50.0, 1.0, key="p_salt_c", format="%.1f")
                mg_mM = st.number_input("Mg++ (mM):", 0.0, 50.0, 0.0, 0.1, key="p_mg_c", format="%.1f"); dNTPs_mM = st.number_input("dNTPs (mM):", 0.0, 10.0, 0.0, 0.1, key="p_dntp_c", format="%.1f")
    st.divider()
    if st.button("Design Primers", key="p_design_btn"):
        if not sequence: st.warning("Provide DNA template.")
        else:
            seq_clean = "".join(re.findall(r'[ATCG]', sequence.upper()))
            if not seq_clean: st.error("No valid DNA bases found.")
            elif len(seq_clean) < primer_len * 2: st.error(f"Sequence too short ({len(seq_clean):,} bp) for primers of {primer_len} bp.")
            else:
                try:
                    fw_p=Seq(seq_clean[:primer_len]); rv_p=Seq(seq_clean[-primer_len:]).reverse_complement(); fw_gc=gc_fraction(fw_p)*100; rv_gc=gc_fraction(rv_p)*100; fw_tm, rv_tm="N/A", "N/A"; tm_params_display=f"Method: {tm_method}"
                    try:
                        tm_args={'strict': False};
                        if tm_method=="Tm_NN": tm_args.update({'Na': salt_mM,'Mg': mg_mM,'dNTPs': dNTPs_mM,'dnac1': dna_nM,'nn_table': MeltingTemp.DNA_NN4}); tm_params_display+=f", P={dna_nM:.1f}nM, Na={salt_mM:.1f}mM, Mg={mg_mM:.1f}mM, dNTP={dNTPs_mM:.1f}mM"
                        if hasattr(MeltingTemp, tm_method): tm_func=getattr(MeltingTemp, tm_method); fw_tm=tm_func(fw_p, **tm_args); rv_tm=tm_func(rv_p, **tm_args)
                        else: st.error(f"Tm method '{tm_method}' not found.")
                    except Exception as tm_e: st.warning(f"Tm calc failed: {tm_e}")
                    st.subheader(f"Suggested Primers ('**{seq_id}**')")
                    p_col1, p_col2=st.columns(2)
                    with p_col1: st.markdown("#### Forward (5'→3')"); st.code(str(fw_p)); st.metric("Len", f"{len(fw_p)} bp"); st.metric("GC%", f"{fw_gc:.1f}%"); st.metric("Tm", f"{fw_tm:.1f}°C" if isinstance(fw_tm,(float,int)) else "N/A")
                    with p_col2: st.markdown("#### Reverse (5'→3')"); st.code(str(rv_p)); st.metric("Len", f"{len(rv_p)} bp"); st.metric("GC%", f"{rv_gc:.1f}%"); st.metric("Tm", f"{rv_tm:.1f}°C" if isinstance(rv_tm,(float,int)) else "N/A")
                    st.divider(); st.metric(f"Amplicon Size ('**{seq_id}**')", f"{len(seq_clean):,} bp")
                    if isinstance(fw_tm,(float,int)): st.caption(f"Tm calculated using: {tm_params_display}")
                except Exception as e: st.error(f"Primer design error: {e}")

elif menu == "Restriction Enzyme Analysis":
    st.header("Restriction Enzyme Analysis")
    st.markdown("Identify restriction enzyme cut sites in DNA (A/T/C/G only).")
    st.info("Uses Biopython REBASE data.")
    input_method=st.radio("Input:",("Paste","FASTA"),key="re_in",horizontal=True)
    sequence,seq_id="","Pasted Sequence"
    if input_method=="Paste": sequence=st.text_area("DNA Sequence:",height=150,key="re_seq",placeholder="Paste DNA...")
    else: uploaded_file=st.file_uploader("FASTA:",type=['fasta','fa','fna'],key="re_fasta"); sequence,seq_id=parse_fasta(uploaded_file) if uploaded_file else (None,None)
    if sequence: st.text_area("Preview:",value=f">**{seq_id}**\n{sequence[:100]}{'...' if len(sequence)>100 else ''}",height=100,key="re_disp",disabled=True)
    st.divider(); st.subheader("Enzyme Selection")
    common_enzymes=sorted(['EcoRI','BamHI','HindIII','NotI','XhoI','PstI','SacI','SalI','SpeI','KpnI','SmaI','EcoRV','HaeIII','AluI','BglII','NcoI','NdeI','XbaI'])
    selected_common=st.multiselect("Common enzymes:",common_enzymes,default=['EcoRI','BamHI','HindIII'],key="re_sel")
    custom_enzymes_str=st.text_input("Custom enzymes (comma-sep):",key="re_cust",placeholder="e.g., BsaI")
    all_selected=set(selected_common); if custom_enzymes_str: all_selected.update({e.strip() for e in custom_enzymes_str.split(',') if e.strip()})
    final_list=sorted(list(all_selected)); is_linear=st.checkbox("Linear DNA",value=True,key="re_lin",help="Uncheck for circular.")
    st.divider()
    if st.button("Analyze Sites",key="re_analyze_btn"):
        if not sequence: st.warning("Provide DNA.")
        elif not final_list: st.warning("Select enzyme(s).")
        else:
            seq_clean="".join(re.findall(r'[ATCG]',sequence.upper()))
            if not seq_clean: st.error("No valid DNA bases found.")
            else:
                try:
                    seq_obj=Seq(seq_clean); st.info(f"Analyzing '**{seq_id}**' ({len(seq_obj):,} bp, {'Linear' if is_linear else 'Circular'}).")
                    valid_obj,invalid_names=[],[]
                    try: from Bio.Restriction import AllEnzymes
                    except ImportError: st.error("Cannot import RE data."); st.stop()
                    with st.spinner("Validating..."):
                        for name in final_list:
                            try: enz=AllEnzymes.get(name); valid_obj.append(enz) if enz else invalid_names.append(name)
                            except ValueError: invalid_names.append(name)
                    if invalid_names: st.warning(f"Ignored unknown: `{'`, `'.join(invalid_names)}`")
                    if not valid_obj: st.error("No valid enzymes."); st.stop()
                    st.write(f"**Using:** {', '.join(map(str, valid_obj))}")
                    with st.spinner("Searching..."): rb=RestrictionBatch(valid_obj); analysis=Analysis(rb,seq_obj,linear=is_linear)
                    st.subheader("Results"); st.write(f"**Seq:** '**{seq_id}**' | **Len:** {len(seq_obj):,} bp | **Type:** {'Linear' if is_linear else 'Circular'}"); st.divider()
                    results_dict=analysis.full(); summary,total_cuts,sites_found=[],0,False
                    for enz,sites in results_dict.items():
                        count=len(sites); total_cuts+=count; sites_found = sites_found or (count>0)
                        summary.append({"Enzyme": str(enz),"Site": str(enz.site),"Cuts": count,"Positions (1-based)": ", ".join(f"{s:,}" for s in sites) if sites else "None"})
                    if not sites_found: st.success("✅ No cut sites found.")
                    else:
                        st.write("#### Cut Summary:"); summary_df=pd.DataFrame(summary).sort_values(["Cuts","Enzyme"],ascending=[False,True]); st.dataframe(summary_df.style.format({"Cuts":"{:,}"}),hide_index=True,use_container_width=True)
                        st.divider(); st.write("#### Fragments:")
                        try:
                            frag_lens=sorted(analysis.fragments(),reverse=True) # Assumes returns lengths
                            f_col1,f_col2=st.columns(2); f_col1.metric("Cuts",f"{total_cuts:,}"); f_col2.metric("Fragments",f"{len(frag_lens):,}")
                            st.write("**Lengths (bp):**")
                            if frag_lens: st.text(", ".join(f"{l:,}" for l in frag_lens[:15]) + (f", ... ({len(frag_lens)-15} more)" if len(frag_lens)>30 else ""))
                            else: st.text("Fragment info unavailable.")
                        except Exception as frag_e: st.warning(f"Cannot get fragments: {frag_e}")
                        st.divider();
                        if st.checkbox("Show Text Map",False,key="re_map"):
                            st.write("#### Map:")
                            try: map_out=io.StringIO(); map_w=min(100,max(60,len(seq_obj)//10)); analysis.print_that(out=map_out,title=f"Map: {seq_id}",top=True,rc=True,nc=map_w); st.code(map_out.getvalue()); map_out.close()
                            except Exception as map_e: st.error(f"Map failed: {map_e}")
                except Exception as e: st.error(f"Analysis error: {e}")

elif menu == "Pairwise Sequence Alignment": # Logic mostly unchanged, FIX applied
    st.header("Pairwise Sequence Alignment")
    st.markdown("Align two DNA or protein sequences using global (Needleman-Wunsch) or local (Smith-Waterman) algorithms.")
    st.info("Scores depend heavily on parameters (match/mismatch/matrix, gap penalties).")

    col1, col2 = st.columns(2)
    with col1:
        st.markdown("#### Sequence 1"); input_method1 = st.radio("Input:",("Paste","FASTA"), key="al_in1", horizontal=True, label_visibility="collapsed")
        seq1, seq_id1 = "", "Sequence 1"
        if input_method1 == "Paste": seq1 = st.text_area("Seq 1:", height=120, key="al_s1", placeholder="Paste sequence...")
        else: uploaded_file1 = st.file_uploader("FASTA 1:", type=['fasta','fa','fna','faa'], key="al_f1"); seq1, seq_id1 = parse_fasta(uploaded_file1) if uploaded_file1 else (None, None)
        if seq1: st.text_area("Preview 1:", value=f">**{seq_id1}**\n{seq1[:80]}{'...' if len(seq1)>80 else ''}", height=75, key="al_d1", disabled=True)
    with col2:
        st.markdown("#### Sequence 2"); input_method2 = st.radio("Input:",("Paste","FASTA"), key="al_in2", horizontal=True, label_visibility="collapsed")
        seq2, seq_id2 = "", "Sequence 2"
        if input_method2 == "Paste": seq2 = st.text_area("Seq 2:", height=120, key="al_s2", placeholder="Paste sequence...")
        else: uploaded_file2 = st.file_uploader("FASTA 2:", type=['fasta','fa','fna','faa'], key="al_f2"); seq2, seq_id2 = parse_fasta(uploaded_file2) if uploaded_file2 else (None, None)
        if seq2: st.text_area("Preview 2:", value=f">**{seq_id2}**\n{seq2[:80]}{'...' if len(seq2)>80 else ''}", height=75, key="al_d2", disabled=True)

    st.divider(); st.subheader("Alignment Parameters")
    pcol1, pcol2, pcol3 = st.columns(3)
    with pcol1: alignment_mode = st.selectbox("Mode:", ["Global (Needleman-Wunsch)", "Local (Smith-Waterman)"], key="al_m"); sequence_type = st.radio("Type:", ("DNA", "Protein"), key="al_t", horizontal=True)
    with pcol2:
        st.write("**Scoring:**"); substitution_matrix, match_score, mismatch_penalty = None, None, None
        if sequence_type == "DNA":
            match_score = st.number_input("Match:", value=2.0, step=0.5, key="al_mt"); mismatch_penalty = st.number_input("Mismatch:", value=-1.0, step=-0.5, key="al_ms"); mismatch_penalty = min(0.0, mismatch_penalty)
            st.caption(f"Match={match_score}, Mismatch={mismatch_penalty}")
        else:
            available_matrices = sorted(substitution_matrices.list_matrices()); default_matrix_index = available_matrices.index('BLOSUM62') if 'BLOSUM62' in available_matrices else 0
            selected_matrix_name = st.selectbox("Matrix:", available_matrices, index=default_matrix_index, key="al_mat")
            try: substitution_matrix = substitution_matrices.load(selected_matrix_name); st.caption(f"Using **{selected_matrix_name}**")
            except Exception as e: st.error(f"Error loading matrix '{selected_matrix_name}': {e}")
    with pcol3:
        st.write("**Gap Penalties:**"); gap_open_penalty = st.number_input("Open:", value=-10.0, step=-0.5, key="al_go"); gap_extend_penalty = st.number_input("Extend:", value=-0.5, step=-0.1, key="al_ge")
        gap_open_penalty = min(0.0, gap_open_penalty); gap_extend_penalty = min(0.0, gap_extend_penalty)
        st.caption(f"Open={gap_open_penalty}, Extend={gap_extend_penalty}")

    st.divider()
    if st.button("Align Sequences", key="al_btn"):
        valid_input = True
        if not seq1: st.warning("Seq 1 missing."); valid_input = False
        if not seq2: st.warning("Seq 2 missing."); valid_input = False
        if sequence_type == "Protein" and not substitution_matrix: st.error("Valid protein matrix needed."); valid_input = False
        if not valid_input: st.stop()
        seq1_clean, seq2_clean = "".join(seq1.split()).upper(), "".join(seq2.split()).upper()
        if not seq1_clean: st.error("Seq 1 empty."); valid_input = False
        if not seq2_clean: st.error("Seq 2 empty."); valid_input = False
        if not valid_input: st.stop()

        try:
            st.subheader("Alignment Results")
            mode_prefix = "global" if alignment_mode.startswith("Global") else "local"
            gap_args = {'open': gap_open_penalty, 'extend': gap_extend_penalty}; align_func = None; params = {}
            if substitution_matrix: align_function_name = f"{mode_prefix}dx"; params = {'matrix': substitution_matrix, **gap_args}
            else: align_function_name = f"{mode_prefix}mx"; params = {'match': match_score, 'mismatch': mismatch_penalty, **gap_args}
            if hasattr(pairwise2.align, align_function_name): align_func = getattr(pairwise2.align, align_function_name)
            else: st.error(f"Could not find alignment function '{align_function_name}'. Biopython version issue?"); st.stop()

            st.info(f"Performing {alignment_mode} {sequence_type} alignment for '**{seq_id1}**' vs '**{seq_id2}**'...")
            with st.spinner("Aligning..."): alignments = align_func(seq1_clean, seq2_clean, **params, one_alignment_only=True)

            if not alignments: st.warning("No alignment generated.")
            else:
                aligned1, aligned2, score, begin, end = alignments[0]
                st.metric("Score", f"{score:.2f}")
                if mode_prefix == "local": st.write(f"**Aligned Region (0-based):** Start={begin}, End={end}")
                st.divider(); st.write("#### Best Alignment:")
                formatted_alignment_string = pairwise2.format_alignment(aligned1, aligned2, score, begin, end, full_sequences=(mode_prefix == 'global'))
                display_text = f"Sequence 1: {seq_id1}\nSequence 2: {seq_id2}\n\n{formatted_alignment_string}"
                st.code(display_text, language='text')
                st.caption("Key: '|'=Match, '.'=Mismatch(+), ' '=Mismatch(-), '-'=Gap.")

                identity_count, gap_count1, gap_count2, aligned_pairs_count = 0, 0, 0, 0
                alignment_length = len(aligned1)
                for i in range(alignment_length):
                    res1, res2 = aligned1[i], aligned2[i]; is_gap1, is_gap2 = (res1 == '-'), (res2 == '-')
                    # *** SYNTAX FIX APPLIED HERE ***
                    if is_gap1: gap_count1 += 1
                    if is_gap2: gap_count2 += 1
                    # *** END FIX ***
                    if not is_gap1 and not is_gap2: aligned_pairs_count += 1; if res1 == res2: identity_count += 1
                total_gaps = gap_count1 + gap_count2
                identity_percent = (identity_count / aligned_pairs_count * 100) if aligned_pairs_count > 0 else 0
                gap_percent = (total_gaps / (alignment_length * 2) * 100) if alignment_length > 0 else 0

                st.divider(); st.write("#### Alignment Stats:")
                scol1, scol2, scol3 = st.columns(3)
                with scol1: st.metric("Length", f"{alignment_length:,}")
                with scol2: st.metric("Identity", f"{identity_percent:.1f}%", f"{identity_count:,}/{aligned_pairs_count:,}")
                with scol3: st.metric("Gaps", f"{gap_percent:.1f}%", f"{total_gaps:,}/{alignment_length*2:,}")
                st.caption("Identity = Matches / Aligned Columns (no gaps). Gaps = Total Gap Chars / (2 * Length).")
        except Exception as e: st.error(f"Alignment error: {e}"); st.exception(e)


elif menu == "Motif Finder Tool": # Logic unchanged, styling updates via CSS
    st.header("Sequence Motif Finder")
    st.markdown("Search for patterns (motifs) using exact matching or regular expressions (regex).")
    st.info("Regex examples: `GAATTC`, `G[AT]ATTC`, `M[A-Z]{5}K`, `^ATG`, `TAG$`. Escape special chars with `\\`.")
    input_method = st.radio("Input:",("Paste","FASTA"), key="mo_in", horizontal=True)
    sequence, seq_id="","Pasted Sequence"
    if input_method=="Paste": sequence=st.text_area("Sequence:", height=150, key="mo_seq", placeholder="Paste DNA or protein sequence...")
    else: uploaded_file=st.file_uploader("FASTA:", type=['fasta','fa','fna','faa'], key="mo_f"); sequence, seq_id = parse_fasta(uploaded_file) if uploaded_file else (None, None)
    if sequence: st.text_area("Preview:", value=f">**{seq_id}**\n{sequence[:80]}{'...' if len(sequence)>80 else ''}", height=75, key="mo_d", disabled=True)
    st.divider(); motif_pattern = st.text_input("Motif / Regex:", key="mo_pat", placeholder="e.g., GAATTC or ^M[A-Z]+")
    mcol1, mcol2 = st.columns(2); with mcol1: case_sensitive = st.checkbox("Case Sensitive", value=False, key="mo_case"); with mcol2: allow_overlap = st.checkbox("Allow Overlaps", value=True, key="mo_ov")
    st.divider()
    if st.button("Find Motifs", key="mo_find"):
        if not sequence: st.warning("Provide sequence."); elif not motif_pattern: st.warning("Enter motif/regex.")
        else:
            seq_cleaned = "".join(sequence.split())
            if not seq_cleaned: st.warning("Sequence empty after cleaning.")
            else:
                try:
                    regex_flags = 0 if case_sensitive else re.IGNORECASE; found_matches = []
                    st.info(f"Searching for '**{motif_pattern}**' in '**{seq_id}**'...")
                    with st.spinner("Searching..."):
                        compiled_pattern = re.compile(motif_pattern, flags=regex_flags)
                        if allow_overlap: found_matches = list(compiled_pattern.finditer(seq_cleaned))
                        else: current_pos = 0; while current_pos < len(seq_cleaned): match = compiled_pattern.search(seq_cleaned, current_pos); if match: found_matches.append(match); current_pos = match.end(); else: break
                    st.subheader(f"Results for '**{motif_pattern}**' in '**{seq_id}**'")
                    if found_matches:
                        st.success(f"Found **{len(found_matches):,}** match(es).")
                        match_data = [{"#": i + 1, "Start (1-based)": m.start() + 1, "End (1-based)": m.end(), "Length": m.end() - m.start(), "Match": m.group()} for i, m in enumerate(found_matches)]
                        match_df = pd.DataFrame(match_data); st.dataframe(match_df.style.format({"Start (1-based)":"{:,}", "End (1-based)":"{:,}", "Length":"{:,}"}), hide_index=True, use_container_width=True)
                        st.divider()
                        if st.checkbox("Highlight Matches (first 2000 chars)", value=False, key="mo_hi"):
                             limit = 2000; seq_hl = seq_cleaned[:limit]; html = ""; last = 0
                             mark_style = "background-color:#fff3cd; padding: 1px 3px; border-radius: 3px; color:#856404; font-weight:bold;" # Light theme highlight
                             mark_open = f"<mark style='{mark_style}'>"; mark_close = "</mark>"; sorted_m = sorted(found_matches, key=lambda m: m.start())
                             for m in sorted_m:
                                 if m.start() >= limit: break; start, end = m.start(), min(m.end(), limit)
                                 if start >= last: txt = seq_hl[last:start].replace("&","&").replace("<","<").replace(">",">"); html += txt; match_txt = seq_hl[start:end].replace("&","&").replace("<","<").replace(">",">"); html += mark_open + match_txt + mark_close; last = end
                                 elif end > last: overlap = seq_hl[last:end].replace("&","&").replace("<","<").replace(">",">"); html += mark_open + overlap + mark_close; last = end
                             rem = seq_hl[last:].replace("&","&").replace("<","<").replace(">",">"); html += rem
                             st.markdown(f"**Highlighted ({limit:,} chars):**"); st.markdown(f"<div style='font-family:monospace; word-wrap:break-word; border:1px solid #dee2e6; padding:10px; border-radius:5px; background-color:#f8f9fa; color:#333;'>{html}{'...' if len(seq_cleaned)>limit else ''}</div>", unsafe_allow_html=True)
                    else: st.info("Motif/pattern not found.")
                except re.error as regex_e: st.error(f"Invalid Regex: {regex_e}")
                except Exception as e: st.error(f"Motif search error: {e}"); st.exception(e)

elif menu == "Bioinformatics Tool (Transcription/Translation)": # Logic unchanged, styling updates via CSS
    st.header("DNA Transcription and Translation Tool")
    st.markdown("Transcribe DNA (coding strand, A/T/C/G) to RNA, translate RNA to protein.")
    input_method = st.radio("Input:",("Paste","FASTA"), key="tr_in", horizontal=True)
    dna_sequence, seq_id="","Pasted Sequence"
    if input_method=="Paste": dna_sequence=st.text_area("DNA Sequence (Coding Strand):", height=120, key="tr_seq", placeholder="Paste DNA sequence...")
    else: uploaded_file=st.file_uploader("FASTA:", type=['fasta','fa','fna'], key="tr_f"); dna_sequence, seq_id = parse_fasta(uploaded_file) if uploaded_file else (None, None)
    if dna_sequence: st.text_area("Preview:", value=f">**{seq_id}**\n{dna_sequence[:80]}{'...' if len(dna_sequence)>80 else ''}", height=75, key="tr_d", disabled=True)
    st.divider()
    if st.button("Transcribe & Translate", key="tr_btn"):
        if not dna_sequence: st.warning("Provide DNA.")
        else:
            dna_cleaned = "".join(re.findall(r'[ATCG]', dna_sequence.upper()))
            if not dna_cleaned: st.error("No valid DNA bases.")
            else:
                try:
                    st.subheader(f"Results for: '**{seq_id}**'")
                    dna_obj=Seq(dna_cleaned); dna_len=len(dna_obj)
                    col1,col2=st.columns(2)
                    with col1: st.metric("DNA Len",f"{dna_len:,} bp"); st.write("**Input DNA:**"); st.code(str(dna_obj))
                    with col2: gc=gc_fraction(dna_obj)*100 if dna_len>0 else 0; st.metric("GC%",f"{gc:.1f}%"); st.write("**Rev Comp:**"); st.code(str(dna_obj.reverse_complement()))
                    st.divider(); st.write("#### Transcription"); rna_obj=dna_obj.transcribe(); st.write("**RNA:**"); st.code(str(rna_obj)); st.metric("RNA Len",f"{len(rna_obj):,} bases")
                    st.divider(); st.write("#### Translation"); st.caption("Standard code. '*' = STOP.")
                    rem=len(rna_obj)%3; rna_trans=rna_obj; warn=False
                    if len(rna_obj)<3: st.error("RNA < 3 bases."); rna_trans=None
                    elif rem!=0: warn=True; rna_trans=rna_obj[:-rem]; if not rna_trans: st.error("Effective RNA len zero."); rna_trans=None
                    if rna_trans:
                        prot_stop=rna_trans.translate(table=1,to_stop=True,cds=False); st.write("**Protein (to STOP):**"); st.code(str(prot_stop)); st.metric("Len (to STOP)",f"{len(prot_stop):,} aa")
                        if warn: st.warning(f"RNA len not multiple of 3; last {rem} base(s) ignored.")
                        st.divider()
                        if st.checkbox("Show Full Translation",False,key="tr_full"): prot_full=rna_trans.translate(table=1,to_stop=False,cds=False); st.write("**Full Protein:**"); st.code(str(prot_full)); st.metric("Full Len",f"{len(prot_full):,} aa")
                except Exception as e: st.error(f"Processing error: {e}")

elif menu == "Genome Coverage Plotter": # Logic unchanged, styling updates via CSS/style_plot
    st.header("Genome Coverage Plotter")
    st.markdown("Visualize sequencing coverage depth. Upload position/coverage data.")
    st.info("Requires columns for position and coverage depth.")
    uploaded_file=st.file_uploader("Upload Coverage Data:",type=['csv','tsv','txt','bed','bedgraph','depth'],key="cov_up")
    if uploaded_file:
        st.divider(); st.subheader("Parsing Settings")
        ext=os.path.splitext(uploaded_file.name)[-1].lower(); def_sep='\t' if ext!='.csv' else ','
        c1,c2=st.columns(2)
        with c1: sep_opts={',':'Comma (,)','\t':'Tab (\\t)',' ':'Whitespace ( )'}; sep_keys=list(sep_opts.keys()); def_idx=sep_keys.index(def_sep) if def_sep in sep_keys else 1; sel_key=st.selectbox("Separator:",sep_keys,index=def_idx,format_func=lambda x:sep_opts[x],key="cov_sep"); sep_re=r'\s+' if sel_key==' ' else sel_key
        with c2: def_head=ext in ['.csv','.tsv','.txt']; head=st.checkbox("Has header",value=def_head,key="cov_head"); head_arg=0 if head else None; comm=st.text_input("Comment Char:",'#',max_chars=1,key="cov_comm")
        st.divider(); st.subheader("Data Loading & Columns")
        try:
            df_cov=pd.read_csv(io.BytesIO(uploaded_file.getvalue()),sep=sep_re,header=head_arg,engine='python',comment=comm if comm else None,low_memory=False)
            if df_cov.empty: st.error("File empty/only comments."); st.stop()
            if head_arg is None and all(isinstance(c,int) for c in df_cov.columns): n_cols=df_cov.shape[1]; df_cov.columns=[f'Col_{i+1}' for i in range(n_cols)] if n_cols>=2 else df_cov.columns; st.info(f"Generic names assigned: {list(df_cov.columns)}")
            if df_cov.shape[1]<2: st.error("Need >= 2 columns."); st.stop()
            st.write("Preview:"); st.dataframe(df_cov.head(),height=200,use_container_width=True)
            cols=list(df_cov.columns); norm_cols={str(c).strip().lower():c for c in cols}; p_g,c_g,ch_g=None,None,None; ch_k=['chr','chrom','#chr']; p_k=['pos','position','start']; c_k=['depth','cov','coverage','score']
            for k in ch_k:   if k in norm_cols: ch_g=norm_cols[k]; break
            for k in p_k:   if k in norm_cols: p_g=norm_cols[k]; break
            for k in c_k:   if k in norm_cols: c_g=norm_cols[k]; break
            ch_idx=cols.index(ch_g) if ch_g else 0; p_idx=cols.index(p_g) if p_g else 1; c_idx=cols.index(c_g) if c_g else (3 if len(cols)>3 and ext in ['.bedgraph','.bed'] else 2); ch_idx,p_idx,c_idx=min(ch_idx,len(cols)-1),min(p_idx,len(cols)-1),min(c_idx,len(cols)-1)
            cs1,cs2,cs3=st.columns(3);
            with cs1: chr_col=st.selectbox("Chrom Col (Opt):",[None]+cols,index=(cols.index(ch_g)+1 if ch_g else 0),key="cov_chr")
            with cs2: pos_col=st.selectbox("Pos Col:",cols,index=p_idx,key="cov_pos")
            with cs3: cov_col=st.selectbox("Cov Col:",cols,index=c_idx,key="cov_cov")
            if pos_col==cov_col or (chr_col and (chr_col==pos_col or chr_col==cov_col)): st.error("Columns must differ.")
            else:
                st.success(f"Using '{pos_col}' (Pos), '{cov_col}' (Cov)" + (f", '{chr_col}' (Chrom)" if chr_col else ".")); st.divider(); st.subheader("Filtering & Plotting")
                try:
                    use_c=[pos_col,cov_col]; if chr_col: use_c.insert(0,chr_col); plot_df_full=df_cov[use_c].copy(); new_n=['Position','Coverage']; if chr_col: new_n.insert(0,'Chromosome'); plot_df_full.columns=new_n; sel_chr=None
                    if chr_col: u_chrs=sorted(plot_df_full['Chromosome'].astype(str).unique()); opts=['All']+u_chrs; sel_chr=st.selectbox("Filter Chrom:",opts,0,key="cov_filter_chr"); plot_data=plot_df_full[plot_df_full['Chromosome'].astype(str)==str(sel_chr)].copy() if sel_chr!='All' else plot_df_full.copy(); st.info(f"Filtered for **{sel_chr}**.") if sel_chr!='All' else None
                    else: plot_data=plot_df_full.copy()
                    plot_data['Position']=pd.to_numeric(plot_data['Position'],errors='coerce'); plot_data['Coverage']=pd.to_numeric(plot_data['Coverage'],errors='coerce'); r_bef=len(plot_data); plot_data.dropna(subset=['Position','Coverage'],inplace=True); r_aft=len(plot_data);
                    if r_aft<r_bef: st.warning(f"Removed **{r_bef-r_aft:,}** rows with invalid data.")
                    if plot_data.empty: st.error("No valid data left."); st.stop()
                    plot_data=plot_data.sort_values('Position').reset_index(drop=True); mn_p,mx_p=plot_data['Position'].min(),plot_data['Position'].max(); t_sfx=f" ({sel_chr})" if sel_chr and sel_chr!='All' else f" (All)" if chr_col else ""; st.write(f"Plotting **{mn_p:,.0f}** to **{mx_p:,.0f}**{t_sfx}.")
                    cp1,cp2=st.columns(2); with cp1: fill=st.checkbox("Fill",True,key="cov_fill"); log_y=st.checkbox("Log Y",False,key="cov_log")
                    with cp2: mx_sm=min(201,max(1,len(plot_data)//2)); sm_w=st.slider("Smoothing (0=None):",0,mx_sm,0,2,key="cov_smooth"); if sm_w>0 and sm_w%2==0: sm_w+=1; st.caption(f"Smooth win: {sm_w}")
                    fig,ax=plt.subplots(figsize=(14,5)); ax.plot(plot_data['Position'],plot_data['Coverage'],label='Raw',lw=0.7,color='#007bff',alpha=0.8); if fill: ax.fill_between(plot_data['Position'],plot_data['Coverage'],alpha=0.2,color='#80ccff')
                    if sm_w>1: plot_data['Smoothed']=plot_data['Coverage'].rolling(sm_w,center=True,min_periods=1).mean(); ax.plot(plot_data['Position'],plot_data['Smoothed'],label=f'Smooth ({sm_w}bp)',lw=1.8,color='#ff7f0e',ls='-')
                    p_t=f"Coverage: {uploaded_file.name}{t_sfx}"; p_x=f"Pos ({pos_col})"; p_y=f"Depth ({cov_col}){ ' (Log)' if log_y else ''}"; style_plot(fig,ax,title=p_t,xlabel=p_x,ylabel=p_y)
                    if log_y: ax.set_yscale('log'); mn_pos_c=plot_data[plot_data['Coverage']>0]['Coverage'].min(); ax.set_ylim(bottom=max(0.1,mn_pos_c*0.5) if pd.notna(mn_pos_c) and mn_pos_c>0 else 0.1)
                    else: ax.set_ylim(bottom=0)
                    ax.ticklabel_format(style='plain',axis='x',useOffset=False); ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x,p:format(int(x),','))); leg=ax.legend(facecolor='#ffffff',labelcolor='#333333'); [t.set_fontfamily('Times New Roman') or t.set_fontweight('bold') for t in leg.get_texts()]; st.pyplot(fig)
                    st.divider(); st.subheader(f"Stats (Raw{t_sfx})"); stats=plot_data['Coverage'].describe(); sc1,sc2,sc3,sc4=st.columns(4); sc1.metric("Mean",f"{stats.get('mean',0):,.1f}x"); sc2.metric("Median",f"{stats.get('50%',0):,.1f}x"); sc3.metric("Min",f"{stats.get('min',0):,.0f}x"); sc4.metric("Max",f"{stats.get('max',0):,.0f}x")
                    mx_c=int(stats.get('max',1000)); df_th=min(10,max(1,int(stats.get('mean',10)))); mx_th=max(1,mx_c); cov_th=st.number_input("Breadth Depth ≥:",0,mx_th,df_th,1,key="cov_b_thresh")
                    if cov_th>=0: above=(plot_data['Coverage']>=cov_th).sum(); total=len(plot_data); if total>0: brd=(above/total*100); st.metric(f"Breadth (≥{cov_th}x)",f"{brd:.1f}%",f"{above:,}/{total:,} pos"); else: st.metric(f"Breadth (≥{cov_th}x)","N/A")
                except Exception as plot_e: st.error(f"Plotting error: {plot_e}")
        except Exception as read_e: st.error(f"Read error: {read_e}"); st.exception(read_e)

elif menu == "Variant Annotation Tool": # Logic unchanged, styling updates via CSS
    st.header("Simple Variant Substitution Tool"); st.markdown("Introduce SNP into DNA (A/T/C/G) & see codon effect.")
    in_meth=st.radio("Input Ref:",("Paste","FASTA"),key="var_in",horizontal=True); ref_seq,seq_id="","Ref Seq"
    if in_meth=="Paste": ref_seq=st.text_area("Reference DNA:",height=100,key="var_ref",placeholder="Paste ref DNA...")
    else: up_f=st.file_uploader("Ref FASTA:",type=['fasta','fa','fna'],key="var_f"); ref_seq,seq_id=parse_fasta(up_f) if up_f else (None,None)
    if ref_seq: st.text_area("Ref Preview:",value=f">**{seq_id}**\n{ref_seq[:80]}{'...' if len(ref_seq)>80 else ''}",height=75,key="var_d",disabled=True)
    st.divider(); st.subheader("Variant"); v1,v2=st.columns(2); with v1: var_pos=st.number_input("Position (1-based):",1,step=1,value=1,key="var_p"); with v2: var_base=st.selectbox("New Base (Alt):",["A","T","C","G"],key="var_b")
    st.divider()
    if st.button("Apply & Annotate",key="var_btn"):
        if not ref_seq: st.warning("Provide ref DNA.")
        else:
            ref_cl="".join(re.findall(r'[ATCG]',ref_seq.upper())); err=False
            if not ref_cl: st.error("Ref has no valid DNA."); err=True
            else: ref_len=len(ref_cl); z_pos=var_pos-1; if z_pos<0 or z_pos>=ref_len: st.error(f"Pos {var_pos:,} outside range (1-{ref_len:,})."); err=True
            if not err:
                try:
                    orig=ref_cl[z_pos]
                    if orig==var_base: st.info(f"No change: Base at **{var_pos:,}** ('**{seq_id}**') already '{orig}'."); st.write("**Reference:**"); st.code(ref_cl)
                    else:
                        alt_l=list(ref_cl); alt_l[z_pos]=var_base; alt_seq="".join(alt_l)
                        st.subheader("Variant Applied"); st.write(f"**Seq:** '**{seq_id}**' | **Var:** Pos **{var_pos:,}**"); st.markdown(f"**Change:** Ref (`{orig}`) → Alt (`{var_base}`)")
                        s1,s2=st.columns(2); with s1: st.write("**Ref:**"); st.code(ref_cl); with s2: st.write("**Alt:**"); st.code(alt_seq)
                        st.divider(); st.write("#### Codon Context & Effect:")
                        c_start=(z_pos//3)*3; pos_in_c=(z_pos%3)+1
                        if c_start+3<=ref_len:
                            ref_c,alt_c=ref_cl[c_start:c_start+3],alt_seq[c_start:c_start+3]; st.write(f"Var at pos **{var_pos:,}** is pos **{pos_in_c}** in codon at **{c_start+1:,}**.")
                            try:
                                ref_aa,alt_aa=str(Seq(ref_c).translate(1,cds=False)),str(Seq(alt_c).translate(1,cds=False));
                                cc1,cc2=st.columns(2); with cc1: st.write("**Ref Codon:**"); st.code(f"Codon:{ref_c}\nAA:   {ref_aa}"); with cc2: st.write("**Alt Codon:**"); st.code(f"Codon:{alt_c}\nAA:   {alt_aa}")
                                eff="?"; if ref_aa==alt_aa: eff="✅ Silent"; elif alt_aa=='*': eff="🛑 Nonsense" if ref_aa!='*' else "❓ Stop Retained"; elif ref_aa=='*': eff="➡️ Stop-Lost"; else: eff="🔄 Missense"
                                st.markdown(f"**Effect:** {eff}")
                            except Exception as ce: st.warning(f"Cannot translate: {ce}")
                        else: st.info("Variant near end; full codon unavailable.")
                except Exception as e: st.error(f"Error applying variant: {e}")

elif menu == "Codon Usage Analyzer": # Logic unchanged, styling updates via CSS/style_plot
    st.header("Codon Usage Analyzer"); st.markdown("Analyze codon frequency in CDS (A/T/C/G, len multiple of 3).")
    in_m=st.radio("Input CDS:",("Paste","FASTA"),key="cod_in",horizontal=True); cds_seq,seq_id="","CDS Seq"
    if in_m=="Paste": cds_seq=st.text_area("CDS:",height=120,key="cod_seq",placeholder="Paste CDS...")
    else: up_f=st.file_uploader("FASTA:",type=['fasta','fa','fna'],key="cod_f"); cds_seq,seq_id=parse_fasta(up_f) if up_f else (None,None)
    if cds_seq: st.text_area("Preview:",value=f">**{seq_id}**\n{cds_seq[:80]}{'...' if len(cds_seq)>80 else ''}",height=75,key="cod_d",disabled=True)
    st.divider()
    if st.button("Analyze Usage",key="cod_btn"):
        if not cds_seq: st.warning("Provide CDS.")
        else:
            cds_cl="".join(re.findall(r'[ATCG]',cds_seq.upper())); err=False
            if not cds_cl: st.error("No valid DNA."); err=True; elif len(cds_cl)==0: st.error("Empty."); err=True; elif len(cds_cl)%3!=0: st.error(f"Length ({len(cds_cl):,}) not multiple of 3."); err=True
            if not err:
                try:
                    st.subheader(f"Codon Usage for: '**{seq_id}**'")
                    codons=[cds_cl[i:i+3] for i in range(0,len(cds_cl),3)]; tot_c=len(codons); counts=Counter(codons); st.metric("Total Codons",f"{tot_c:,}")
                    try: std_t=CodonTable.unambiguous_dna_by_id[1]; all_p=list(std_t.forward_table.keys())+std_t.stop_codons; aa_map=std_t.forward_table; stops=set(std_t.stop_codons)
                    except Exception as te: st.error(f"Cannot load table: {te}"); st.stop()
                    usage=[]
                    for cdn in sorted(all_p): cnt=counts.get(cdn,0); freq=(cnt/tot_c*100) if tot_c>0 else 0; aa=aa_map.get(cdn,'Stop' if cdn in stops else '?'); usage.append({"Codon":cdn,"AA":aa,"Count":cnt,"Freq (%)":freq})
                    usage_df=pd.DataFrame(usage); st.write("#### Usage Table:"); st.dataframe(usage_df.style.format({'Count':'{:,}','Freq (%)':'{:.1f}%'}),hide_index=True,use_container_width=True)
                    st.divider(); st.subheader("Visualization"); plot_c=st.selectbox("Plot:",["Frequency (All)","Relative Usage (AA)"],key="cod_plot")
                    present_df=usage_df[usage_df['Count']>0].copy()
                    if present_df.empty: st.info("No codons to visualize.")
                    else:
                        if plot_c=="Frequency (All)":
                            df_all=present_df.sort_values("Codon"); fig,ax=plt.subplots(figsize=(16,6)); u_aas=sorted(list(df_all['AA'].unique())); cmap=plt.get_cmap('tab20',len(u_aas)); aa_cols={aa:cmap(i) for i,aa in enumerate(u_aas)}; bar_cols=[aa_cols.get(aa,'#888') for aa in df_all['AA']]; ax.bar(df_all['Codon'],df_all['Freq (%)'],color=bar_cols)
                            style_plot(fig,ax,title=f"Codon Freq for '{seq_id}'",xlabel="Codon",ylabel="Freq (%)"); ax.tick_params(axis='x',rotation=90,labelsize=9); [l.set_fontfamily('Times New Roman') or l.set_fontweight('bold') for l in ax.get_xticklabels()]; st.pyplot(fig)
                        elif plot_c=="Relative Usage (AA)":
                            avail=[aa for aa in present_df['AA'].unique() if aa!='?'];
                            if not avail: st.info("No standard AAs/Stop.")
                            else:
                                sel_aa=st.selectbox("Select AA/Stop:",sorted(avail),key="cod_aa"); aa_df=present_df[present_df['AA']==sel_aa].copy(); tot_aa=aa_df['Count'].sum()
                                if tot_aa>0 and not aa_df.empty:
                                    aa_df['Rel Freq (%)']=(aa_df['Count']/tot_aa*100); aa_df=aa_df.sort_values("Codon"); fig,ax=plt.subplots(figsize=(max(6,len(aa_df)*1.5),5)); ax.bar(aa_df['Codon'],aa_df['Rel Freq (%)'],color='#ff7f0e') # Orange bars
                                    style_plot(fig,ax,title=f"Relative Usage for '{sel_aa}' ({tot_aa:,} total)",xlabel="Codon",ylabel="Rel Freq (%)"); ax.set_ylim(0,105); st.pyplot(fig)
                                    st.write(f"**Details for '{sel_aa}':**"); st.dataframe(aa_df[['Codon','Count','Rel Freq (%)']].style.format({'Count':'{:,}','Rel Freq (%)':'{:.1f}%'}),hide_index=True,use_container_width=True)
                                else: st.info(f"No codons for '{sel_aa}'.")
                except Exception as e: st.error(f"Codon usage error: {e}"); st.exception(e)

# --- Footer ---
st.markdown("---")
st.markdown('<div class="footer">Copyright © 2024 CCDB Tools | For Educational and Informational Purposes Only</div>', unsafe_allow_html=True)
