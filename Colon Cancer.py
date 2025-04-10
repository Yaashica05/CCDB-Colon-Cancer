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
            # else return s_value (if not number and not http -> plain text)

        elif col_name_std == 'Doi Id':
            if '/' in s_value: # Simple check for potential DOI structure
                encoded_doi = urllib.parse.quote(s_value, safe='/:')
                url = f"https://doi.org/{encoded_doi}"
                return f'<a href="{url}" target="_blank" {link_style}>{s_value}</a>'
            # else return s_value

        elif col_name_std == 'Uniprot':
             if re.match(r'^[A-Z0-9_.\-]+$', s_value, re.IGNORECASE): # Check for UniProt ID format
                 url = f"https://www.uniprot.org/uniprotkb/{s_value}/entry"
                 return f'<a href="{url}" target="_blank" {link_style}>{s_value}</a>'
             # else return s_value

        # 3. For other columns (including the newly requested ones):
        # Since we already checked for 'http' start, if we reach here,
        # the value is likely plain text or an ID we don't have a specific
        # URL structure for. Display as plain text.
        else:
            return s_value

    except Exception: # Catch any unexpected formatting error
        return s_value # Return original value if linking fails


# --- Helper function for Matplotlib Plot Styling ---
def style_plot(fig, ax, title="", xlabel="", ylabel=""):
    """Applies consistent LIGHT theme styling (Times New Roman Bold) to a Matplotlib plot."""
    fig.patch.set_facecolor('#f8f9fa') # Match app background (off-white)
    fig.patch.set_alpha(1.0)
    ax.set_facecolor('#ffffff') # White plot area background

    # Titles and Labels (Times New Roman Bold, Dark Color)
    font_props = {'family': 'Times New Roman', 'weight': 'bold'}
    text_color = '#1E2A3A' # Dark text
    accent_color = '#004080' # Dark Blue for title
    ax.set_title(title, color=accent_color, fontsize=14, **font_props)
    ax.set_xlabel(xlabel, color=text_color, fontsize=12, **font_props)
    ax.set_ylabel(ylabel, color=text_color, fontsize=12, **font_props)

    # Ticks (Times New Roman Bold, Dark Color)
    tick_color = '#333333' # Dark grey ticks
    ax.tick_params(axis='x', colors=tick_color, labelsize=10)
    ax.tick_params(axis='y', colors=tick_color, labelsize=10)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontfamily('Times New Roman')
        label.set_fontweight('bold')
        label.set_color(tick_color) # Ensure label color matches ticks

    # Grid (Lighter grey)
    ax.grid(True, linestyle=':', alpha=0.6, color='#cccccc')

    # Spines (borders of the plot area) - Lighter grey
    for spine in ax.spines.values():
        spine.set_edgecolor('#dee2e6')

    plt.tight_layout() # Adjust layout

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

    st.image(
        "https://img.freepik.com/free-vector/colorectal-cancer-crc-infographic-education_1308-47971.jpg",
        caption="Colorectal Cancer Infographic (Source: Freepik - Illustrative purposes)",
        use_container_width=True
    )
    # Apply link styling consistent with the theme
    st.markdown("""<a href="https://www.cancer.org/cancer/types/colon-rectal-cancer.html" target="_blank" style="color:#0056b3;">Learn more about Colorectal Cancer from the American Cancer Society</a>""", unsafe_allow_html=True)

    st.divider()
    st.subheader("Feedback")
    feedback = st.text_area("Share your feedback about this application:", key="feedback_home_page", placeholder="Your thoughts, suggestions, or issues...")
    if st.button("Submit Feedback", key="submit_feedback_home_page"):
        if feedback:
            st.success("Thank you for your feedback!")
            # print(f"Feedback received (Home Page): {feedback}")
        else:
            st.warning("Feedback cannot be empty.")


elif menu == "Colon Cancer Database Search":
    st.header("Colon Cancer Gene Database Search")
    st.markdown(f"Search for information on genes associated with Colorectal Cancer using a local database file (`{DATA_PATH}`). Please enter a gene symbol.")

    if not os.path.exists(DATA_PATH):
        st.error(f"**Database File Not Found:** The required file `{DATA_PATH}` could not be located in the application directory.")
        st.info("Please ensure the Excel database file is present in the same folder as the Streamlit script.")
        st.stop()

    query = st.text_input("Enter Gene Symbol:", key="gene_search_input", placeholder="e.g., APC, KRAS, MLH1")

    if st.button("Search Database", key="gene_search_button"):
        if not query:
            st.warning("Please enter a gene symbol to search for.")
        else:
            try:
                data = pd.read_excel(DATA_PATH)
                data.columns = data.columns.str.strip().str.title() # Standardize column names

                gene_col = None
                possible_gene_cols = ['Gene Symbol', 'Gene Name', 'Symbol', 'Gene']
                for col in possible_gene_cols:
                    if col in data.columns:
                        gene_col = col
                        break

                if not gene_col:
                     st.error(f"Could not find a suitable gene symbol column in the database. Expected one of: {possible_gene_cols}. Found columns: {list(data.columns)}")
                     st.stop()

                data[gene_col] = data[gene_col].astype(str).str.strip()
                search_query = query.strip()
                results = data[data[gene_col].str.fullmatch(search_query, case=False, na=False)].copy()

                if results.empty:
                     st.info(f"No exact match found for '{search_query}'. Trying a broader search...")
                     results = data[data[gene_col].str.contains(search_query, case=False, na=False, regex=False)].copy()

                if not results.empty:
                    st.success(f"Found {len(results)} result(s) matching or containing '{query}'.")
                    st.write("### Search Results:")

                    formatted_results = results.copy()

                    # *** UPDATED: Add user-requested columns to the list for link formatting ***
                    link_potential_columns = [
                        'Pubmed Id', 'Doi Id', 'Uniprot',
                        'Blast', 'Conserved Domain', 'Link', 'Url', 'Reference', # Original potential links
                        'Entry', 'History', 'Variant Viewer', 'Feature Viewer', # New requests
                        'Genomic Coordinates', 'Publications', 'External Links', 'Sequence' # New requests
                        ]
                    # Add any other columns that contain 'link' or 'url' in their name (case-insensitive)
                    link_potential_columns.extend([col for col in formatted_results.columns
                                                   if ('link' in col.lower() or 'url' in col.lower())
                                                   and col not in link_potential_columns])
                    # Ensure uniqueness and handle potential case differences from Excel file
                    link_potential_columns = list(set(link_potential_columns))

                    # Iterate through all columns in the results
                    for col in formatted_results.columns:
                        # Check if the *standardized* column name is in our list
                        # (Standardization happens when reading the Excel file now)
                        if col in link_potential_columns:
                            try:
                                # Pass the current column name ('col') to format_link
                                formatted_results[col] = formatted_results[col].apply(lambda x: format_link(x, col))
                            except Exception as apply_e:
                                st.warning(f"Link formatting attempt failed for column '{col}': {apply_e}. Displaying raw data.")

                    # Display results as an HTML table
                    html_table = formatted_results.to_html(escape=False, index=False, na_rep='-', justify='left', classes=['st-table'])
                    st.write(html_table, unsafe_allow_html=True)
                else:
                    st.warning(f"No results found matching or containing '{query}' in the '{gene_col}' column.")

            except ImportError:
                st.error("The `openpyxl` library is required to read Excel files. Please install it (`pip install openpyxl`).")
            except FileNotFoundError:
                 st.error(f"**Database File Not Found:** The file `{DATA_PATH}` was not found.")
            except Exception as e:
                st.error(f"An error occurred during the database search:")
                st.exception(e)

# --- Other Tool Sections (Unchanged from previous state, but reflecting light theme) ---

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
                st.text_area("Sequence Preview (from FASTA):", value=f">**{seq_id}**\n{sequence[:100]}{'...' if len(sequence) > 100 else ''}", height=100, key="dna_disp", disabled=True)

    st.divider()
    if st.button("Analyze DNA", key="dna_analyze_btn"):
        if not sequence:
            st.warning("Please provide a DNA sequence via paste or FASTA upload.")
        else:
            seq_clean = "".join(re.findall(r'[ATCGN]', sequence.upper()))
            if not seq_clean:
                st.error("The provided input does not contain any valid DNA bases (A, T, C, G, N).")
            else:
                try:
                    st.subheader(f"Analysis Results for: '**{seq_id}**'")
                    seq_obj = Seq(seq_clean)
                    length = len(seq_obj)
                    st.metric("Total Length", f"{length:,} bp")

                    seq_no_n = seq_clean.replace('N', '')
                    len_no_n = len(seq_no_n)
                    if len_no_n > 0:
                        gc_val = gc_fraction(seq_no_n) * 100
                        at_val = 100.0 - gc_val
                        st.metric("GC Content (excluding N)", f"{gc_val:.1f}%")
                        st.metric("AT Content (excluding N)", f"{at_val:.1f}%")
                    else:
                        st.info("GC/AT content could not be calculated.")

                    st.divider()
                    st.write("#### Base Composition (including N):")
                    composition = Counter(seq_clean)
                    comp_data = [{"Base": b, "Count": c, "Percentage (%)": (c/length*100) if length > 0 else 0}
                                 for b, c in sorted(composition.items())]
                    comp_df = pd.DataFrame(comp_data)
                    st.dataframe(comp_df.style.format({"Count":"{:,}", "Percentage (%)":"{:.1f}%"}),
                                 hide_index=True, use_container_width=True)

                    plot_comp = {k: v for k, v in composition.items() if k in STANDARD_DNA_BASES}
                    if plot_comp:
                        fig, ax = plt.subplots(figsize=(6, 4))
                        bases = sorted(plot_comp.keys())
                        counts = [plot_comp[b] for b in bases]
                        # Define colors suitable for light theme (adjust as needed)
                        base_colors = {'A': '#007bff', 'T': '#ff7f0e', 'C': '#2ca02c', 'G': '#d62728'} # Blue, Orange, Green, Red
                        bar_colors = [base_colors.get(b, '#888888') for b in bases]

                        ax.bar(bases, counts, color=bar_colors)
                        style_plot(fig, ax, title="Base Counts (excluding N)", xlabel="Base", ylabel="Count") # Apply light theme style
                        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ',')))
                        st.pyplot(fig)
                    elif 'N' in composition and len(composition) == 1: st.info("Composition plot omitted (sequence contains only N).")
                    elif not plot_comp: st.info("Composition plot omitted (no standard bases A,T,C,G found).")

                    st.divider()
                    st.write("#### Basic Open Reading Frame (ORF) Finder")
                    min_orf_aa = st.number_input("Minimum ORF Length (Amino Acids):", min_value=10, max_value=500, value=30, step=5, key="orf_len")
                    orfs = []
                    st.info("Searching for ORFs (Met 'M' to Stop '*') on the forward strand (3 frames)...")
                    with st.spinner("Scanning for ORFs..."):
                         for frame in range(3):
                            try:
                                translation = str(seq_obj[frame:].translate(table=1, to_stop=True, cds=False))
                            except CodonTable.TranslationError as trans_err:
                                st.warning(f"Translation warning in frame {frame+1}: {trans_err}.")
                                translation = ""
                            except Exception as trans_err:
                                st.warning(f"General translation error in frame {frame+1}: {trans_err}.")
                                continue

                            if not translation: continue

                            current_pos_in_translation = 0
                            while 'M' in translation[current_pos_in_translation:]:
                                start_codon_index = translation.find('M', current_pos_in_translation)
                                if start_codon_index == -1: break

                                potential_orf = translation[start_codon_index:]
                                if len(potential_orf) >= min_orf_aa:
                                    start_dna_pos = frame + (start_codon_index * 3) + 1
                                    end_dna_pos = start_dna_pos + (len(potential_orf) * 3) - 1
                                    if end_dna_pos <= length:
                                        if not any(o['Start (DNA)'] == start_dna_pos and o['Frame'] == frame+1 for o in orfs):
                                            orfs.append({
                                                "Frame": frame + 1, "Start (DNA)": start_dna_pos, "End (DNA)": end_dna_pos,
                                                "Length (AA)": len(potential_orf),
                                                "Protein Sequence (Preview)": potential_orf[:40] + "..." if len(potential_orf) > 40 else potential_orf
                                            })
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
                 st.text_area("Sequence Preview (from FASTA):", value=f">**{seq_id}**\n{sequence[:100]}{'...' if len(sequence) > 100 else ''}", height=100, key="prot_disp", disabled=True)

    st.divider()
    if st.button("Analyze Protein", key="prot_analyze_btn"):
        if not sequence:
            st.warning("Please provide a protein sequence via paste or FASTA upload.")
        else:
            seq_clean = "".join(sequence.split()).upper()
            all_chars = set(seq_clean)
            standard_set = set(STANDARD_AA)
            ambiguous_set = set(AMBIGUOUS_AA)

            non_standard_chars = all_chars - standard_set - ambiguous_set
            ambiguous_chars_present = all_chars.intersection(ambiguous_set)
            seq_for_calc = "".join(c for c in seq_clean if c in standard_set)

            if non_standard_chars: st.warning(f"Ignoring unknown characters: `{'`, `'.join(sorted(list(non_standard_chars)))}`")
            if ambiguous_chars_present: st.warning(f"Ignoring ambiguous amino acid codes for calculations: `{'`, `'.join(sorted(list(ambiguous_chars_present)))}`")

            if not seq_for_calc:
                st.error("No standard amino acids found in the sequence. Cannot perform analysis.")
            else:
                try:
                    st.subheader(f"Analysis Results for: '**{seq_id}**'")
                    pa = ProteinAnalysis(seq_for_calc)

                    col1, col2 = st.columns(2)
                    with col1:
                        st.metric("Input Length (incl. non-std)", f"{len(seq_clean):,} aa")
                        st.metric("Analyzed Length (standard AA)", f"{len(seq_for_calc):,} aa")
                        st.metric("Molecular Weight", f"{pa.molecular_weight():,.1f} Da")
                        st.metric("Isoelectric Point (pI)", f"{pa.isoelectric_point():.2f}")
                    with col2:
                        gravy_score = pa.gravy()
                        hydrophobicity = "Hydrophobic" if gravy_score > 0 else ("Hydrophilic" if gravy_score < 0 else "Neutral")
                        st.metric("GRAVY Index", f"{gravy_score:.3f}", delta=hydrophobicity, delta_color="off")

                        instability_idx = pa.instability_index()
                        stability = "Likely Stable" if instability_idx < 40 else "Likely Unstable"
                        st.metric("Instability Index", f"{instability_idx:.2f}", delta=stability, delta_color="normal" if stability == "Likely Stable" else "inverse")

                    st.divider()
                    st.write("#### Amino Acid Composition (Full Input Sequence):")
                    composition = Counter(seq_clean)
                    comp_data = [{"Amino Acid": aa, "Count": c, "Percentage (%)": (c/len(seq_clean)*100) if seq_clean else 0}
                                 for aa, c in sorted(composition.items())]
                    comp_df = pd.DataFrame(comp_data)
                    st.dataframe(comp_df.style.format({'Count':'{:,}', 'Percentage (%)':'{:.1f}%'}), hide_index=True, use_container_width=True)

                    st.divider()
                    st.write("#### Predicted Secondary Structure Fractions:")
                    try:
                        helix, turn, sheet = pa.secondary_structure_fraction()
                        labels = 'Helix', 'Turn', 'Sheet'
                        sizes = [helix * 100, turn * 100, sheet * 100]

                        if sum(sizes) > 0.1:
                            # Colors suitable for light theme
                            colors = ['#007bff', '#ff7f0e', '#2ca02c'] # Blue, Orange, Green

                            fig, ax = plt.subplots(figsize=(6, 4))
                            wedges, texts, autotexts = ax.pie(sizes, labels=labels, startangle=90, colors=colors, autopct='%1.1f%%',
                                                wedgeprops={'edgecolor': '#ffffff', 'linewidth': 1.5}, # White edge
                                                textprops={'color': '#333333'}) # Dark text for labels

                            # Improve autopct text visibility (bold, potentially white on dark wedges)
                            for autotext in autotexts:
                                autotext.set_color('#ffffff') # White text on wedges
                                autotext.set_fontfamily('Times New Roman')
                                autotext.set_fontweight('bold')
                                autotext.set_fontsize(10)

                            for text in texts:
                                text.set_fontfamily('Times New Roman')
                                text.set_fontweight('bold')

                            ax.axis('equal')
                            style_plot(fig, ax, title="Predicted Secondary Structure") # Apply light theme style
                            ax.set_xlabel("")
                            ax.set_ylabel("")
                            st.pyplot(fig)
                        else:
                            st.info("Secondary structure fractions are essentially zero.")

                    except ZeroDivisionError: st.warning("Could not calculate secondary structure fractions.")
                    except Exception as ss_e: st.warning(f"Secondary structure prediction failed: {ss_e}")

                except Exception as e:
                    st.error("An error occurred during protein sequence analysis:")
                    st.exception(e)

# --- (Continue with other menus, ensuring plots use style_plot for light theme) ---

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
                st.text_area("Sequence Preview (from FASTA):", value=f">**{seq_id}**\n{sequence[:100]}{'...' if len(sequence) > 100 else ''}", height=100, key="p_disp", disabled=True)

    st.divider()
    st.subheader("Primer Parameters")
    col1, col2 = st.columns(2)
    with col1:
        primer_len = st.slider("Desired Primer Length (bp):", min_value=15, max_value=35, value=20, step=1, key="p_len")
    with col2:
        tm_method = st.selectbox("Melting Temperature (Tm) Method:", ["Tm_NN", "Tm_GC", "Tm_Wallace"], index=0, key="p_tm_m", help="Tm_NN is generally more accurate.")
        dna_nM, salt_mM, mg_mM, dNTPs_mM = 50.0, 50.0, 0.0, 0.0
        if tm_method == "Tm_NN":
            with st.expander("Tm_NN Parameters (Advanced)"):
                dna_nM = st.number_input("Primer Concentration (nM - 'dnac1'):", min_value=1.0, max_value=1000.0, value=50.0, step=1.0, key="p_dna_c", format="%.1f")
                salt_mM = st.number_input("Monovalent Cation Conc. (mM - 'Na'):", min_value=1.0, max_value=200.0, value=50.0, step=1.0, key="p_salt_c", format="%.1f")
                mg_mM = st.number_input("Mg2+ Concentration (mM - 'Mg'):", min_value=0.0, max_value=50.0, value=0.0, step=0.1, key="p_mg_c", format="%.1f")
                dNTPs_mM = st.number_input("dNTP Concentration (mM - 'dNTPs'):", min_value=0.0, max_value=10.0, value=0.0, step=0.1, key="p_dntp_c", format="%.1f")

    st.divider()
    if st.button("Design Primers", key="p_design_btn"):
        if not sequence: st.warning("Please provide a DNA template sequence.")
        else:
            seq_clean = "".join(re.findall(r'[ATCG]', sequence.upper()))
            if not seq_clean: st.error("Input contains no valid DNA bases (A, T, C, G).")
            elif len(seq_clean) < primer_len * 2: st.error(f"Template sequence too short ({len(seq_clean):,} bp) for primer length {primer_len} bp.")
            else:
                try:
                    fw_p = Seq(seq_clean[:primer_len])
                    rv_p = Seq(seq_clean[-primer_len:]).reverse_complement()
                    fw_gc = gc_fraction(fw_p) * 100
                    rv_gc = gc_fraction(rv_p) * 100
                    fw_tm, rv_tm = "N/A", "N/A"
                    tm_params_display = f"Method: {tm_method}"
                    try:
                        tm_args = {'strict': False}
                        if tm_method == "Tm_NN":
                             tm_args.update({'Na': salt_mM, 'Mg': mg_mM, 'dNTPs': dNTPs_mM, 'dnac1': dna_nM, 'nn_table': MeltingTemp.DNA_NN4})
                             tm_params_display += f", Primer={dna_nM:.1f}nM, Na={salt_mM:.1f}mM, Mg={mg_mM:.1f}mM, dNTPs={dNTPs_mM:.1f}mM"
                        if hasattr(MeltingTemp, tm_method):
                            tm_func = getattr(MeltingTemp, tm_method)
                            fw_tm = tm_func(fw_p, **tm_args)
                            rv_tm = tm_func(rv_p, **tm_args)
                        else: st.error(f"Could not find Tm method '{tm_method}'.")
                    except Exception as tm_e: st.warning(f"Tm calculation failed: {tm_e}")

                    st.subheader(f"Suggested Primers (from ends of '**{seq_id}**')")
                    p_col1, p_col2 = st.columns(2)
                    with p_col1:
                        st.markdown("#### Forward Primer (5' → 3')")
                        st.code(str(fw_p), language='text')
                        st.metric("Length", f"{len(fw_p)} bp"); st.metric("GC Content", f"{fw_gc:.1f}%"); st.metric(f"Est. Tm", f"{fw_tm:.1f} °C" if isinstance(fw_tm, (float, int)) else "N/A")
                    with p_col2:
                        st.markdown("#### Reverse Primer (5' → 3')")
                        st.code(str(rv_p), language='text')
                        st.metric("Length", f"{len(rv_p)} bp"); st.metric("GC Content", f"{rv_gc:.1f}%"); st.metric(f"Est. Tm", f"{rv_tm:.1f} °C" if isinstance(rv_tm, (float, int)) else "N/A")
                    st.divider()
                    st.metric(f"Resulting Amplicon Size (full template '**{seq_id}**')", f"{len(seq_clean):,} bp")
                    if isinstance(fw_tm, (float, int)): st.caption(f"Tm calculated using: {tm_params_display}")
                except Exception as e: st.error(f"Error during primer design: {e}")


elif menu == "Restriction Enzyme Analysis":
    st.header("Restriction Enzyme Analysis")
    st.markdown("Identify restriction enzyme cut sites within a DNA sequence (A, T, C, G only).")
    st.info("Uses enzyme data from Biopython's REBASE interface.")

    input_method = st.radio("Input Method:", ("Paste", "FASTA"), key="re_in", horizontal=True)
    sequence, seq_id = "", "Pasted Sequence"
    if input_method == "Paste":
        sequence = st.text_area("Enter DNA Sequence:", height=150, key="re_seq", placeholder="Paste DNA sequence here...")
    else:
        uploaded_file = st.file_uploader("Upload FASTA File:", type=['fasta','fa','fna'], key="re_fasta")
        if uploaded_file:
            sequence, seq_id = parse_fasta(uploaded_file)
            if sequence:
                st.text_area("Sequence Preview (from FASTA):", value=f">**{seq_id}**\n{sequence[:100]}{'...' if len(sequence) > 100 else ''}", height=100, key="re_disp", disabled=True)

    st.divider()
    st.subheader("Enzyme Selection")
    common_enzymes = sorted(['EcoRI','BamHI','HindIII','NotI','XhoI','PstI','SacI','SalI','SpeI','SphI','KpnI','SmaI','EcoRV','HaeIII','HhaI','AluI','BglII','NcoI','NdeI','XbaI','ClaI'])
    selected_common = st.multiselect("Select common enzymes:", common_enzymes, default=['EcoRI','BamHI','HindIII'], key="re_sel")
    custom_enzymes_str = st.text_input("Enter custom enzymes (comma-separated, case-sensitive):", key="re_cust", placeholder="e.g., BsaI, PacI")
    all_selected_enzymes = set(selected_common)
    if custom_enzymes_str: all_selected_enzymes.update({e.strip() for e in custom_enzymes_str.split(',') if e.strip()})
    final_enzyme_list = sorted(list(all_selected_enzymes))
    is_linear = st.checkbox("Assume Linear DNA", value=True, key="re_lin", help="Uncheck for circular DNA (e.g., plasmid).")

    st.divider()
    if st.button("Analyze Restriction Sites", key="re_analyze_btn"):
        if not sequence: st.warning("Please provide a DNA sequence.")
        elif not final_enzyme_list: st.warning("Please select at least one enzyme.")
        else:
            seq_clean = "".join(re.findall(r'[ATCG]', sequence.upper()))
            if not seq_clean: st.error("Input contains no valid DNA bases (A, T, C, G).")
            else:
                try:
                    seq_obj = Seq(seq_clean)
                    st.info(f"Analyzing sequence '**{seq_id}**' ({len(seq_obj):,} bp, {'Linear' if is_linear else 'Circular'}).")
                    valid_enzymes_obj, invalid_enzyme_names = [], []
                    try: from Bio.Restriction import AllEnzymes
                    except ImportError: st.error("Could not import Restriction enzyme data."); st.stop()
                    with st.spinner("Validating enzymes..."):
                        for enz_name in final_enzyme_list:
                            try:
                                enz_obj = AllEnzymes.get(enz_name)
                                if enz_obj: valid_enzymes_obj.append(enz_obj)
                                else: invalid_enzyme_names.append(enz_name)
                            except ValueError: invalid_enzyme_names.append(enz_name)
                    if invalid_enzyme_names: st.warning(f"Unrecognized enzymes ignored: `{'`, `'.join(invalid_enzyme_names)}`")
                    if not valid_enzymes_obj: st.error("No valid enzymes selected."); st.stop()

                    st.write(f"**Analyzing with:** {', '.join(map(str, valid_enzymes_obj))}")
                    with st.spinner("Searching for cut sites..."):
                        rb = RestrictionBatch(valid_enzymes_obj)
                        analysis = Analysis(rb, seq_obj, linear=is_linear)

                    st.subheader("Restriction Analysis Results")
                    st.write(f"**Sequence:** '**{seq_id}**' | **Length:** {len(seq_obj):,} bp | **Topology:** {'Linear' if is_linear else 'Circular'}")
                    st.divider()
                    results_dict = analysis.full()
                    summary_data, total_cuts, sites_found = [], 0, False
                    for enz_obj, sites in results_dict.items():
                        cut_count = len(sites); total_cuts += cut_count
                        if cut_count > 0: sites_found = True
                        summary_data.append({"Enzyme": str(enz_obj), "Recognition Site": str(enz_obj.site), "Number of Cuts": cut_count, "Cut Positions (1-based)": ", ".join(f"{s:,}" for s in sites) if sites else "None"})

                    if not sites_found: st.success("✅ No cut sites found for selected enzymes.")
                    else:
                        st.write("#### Cut Site Summary:")
                        summary_df = pd.DataFrame(summary_data).sort_values(by=["Number of Cuts", "Enzyme"], ascending=[False, True])
                        st.dataframe(summary_df.style.format({"Number of Cuts":"{:,}"}), hide_index=True, use_container_width=True)
                        st.divider(); st.write("#### Predicted DNA Fragments:")
                        try:
                            fragment_lengths = sorted(analysis.fragments(), reverse=True) # Assumes fragments() returns lengths
                            f_col1, f_col2 = st.columns(2)
                            with f_col1: st.metric("Total Cuts", f"{total_cuts:,}")
                            with f_col2: st.metric("Number of Fragments", f"{len(fragment_lengths):,}")
                            st.write("**Fragment Lengths (bp, sorted):**")
                            if fragment_lengths:
                                if len(fragment_lengths) > 30: st.text(", ".join(f"{l:,}" for l in fragment_lengths[:15]) + f", ... ({len(fragment_lengths)-15} more)")
                                else: st.text(", ".join(f"{l:,}" for l in fragment_lengths))
                            elif total_cuts == 0 and is_linear: st.text(f"Single fragment of {len(seq_obj):,} bp.")
                            elif total_cuts == 1 and not is_linear: st.text(f"Single fragment of {len(seq_obj):,} bp.")
                            else: st.text("Fragment info unavailable.")
                        except Exception as frag_e: st.warning(f"Could not calculate fragment lengths: {frag_e}")

                        st.divider()
                        if st.checkbox("Show Text-Based Map", value=False, key="re_map"):
                            st.write("#### Restriction Map:")
                            try:
                                map_output = io.StringIO()
                                map_width = min(100, max(60, len(seq_obj) // 10))
                                analysis.print_that(out=map_output, title=f"Map for {seq_id}", top=True, rc=True, nc=map_width)
                                st.code(map_output.getvalue(), language='text')
                                map_output.close()
                            except Exception as map_e: st.error(f"Failed to generate map: {map_e}")
                except Exception as e: st.error(f"Error during restriction analysis: {e}")


elif menu == "Pairwise Sequence Alignment":
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
                    if is_gap1: gap_count1 += 1; if is_gap2: gap_count2 += 1
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
        except Exception as e: st.error(f"Alignment error: {e}")


elif menu == "Motif Finder Tool":
    st.header("Sequence Motif Finder")
    st.markdown("Search for patterns (motifs) using exact matching or regular expressions (regex).")
    st.info("Regex examples: `GAATTC`, `G[AT]ATTC`, `M[A-Z]{5}K`, `^ATG`, `TAG$`. Escape special chars with `\\`.")

    input_method = st.radio("Input:",("Paste","FASTA"), key="mo_in", horizontal=True)
    sequence, seq_id="","Pasted Sequence"
    if input_method=="Paste": sequence=st.text_area("Sequence:", height=150, key="mo_seq", placeholder="Paste DNA or protein sequence...")
    else: uploaded_file=st.file_uploader("FASTA:", type=['fasta','fa','fna','faa'], key="mo_f"); sequence, seq_id = parse_fasta(uploaded_file) if uploaded_file else (None, None)
    if sequence: st.text_area("Preview:", value=f">**{seq_id}**\n{sequence[:80]}{'...' if len(sequence)>80 else ''}", height=75, key="mo_d", disabled=True)

    st.divider(); motif_pattern = st.text_input("Motif / Regex:", key="mo_pat", placeholder="e.g., GAATTC or ^M[A-Z]+")
    mcol1, mcol2 = st.columns(2)
    with mcol1: case_sensitive = st.checkbox("Case Sensitive", value=False, key="mo_case")
    with mcol2: allow_overlap = st.checkbox("Allow Overlaps", value=True, key="mo_ov")

    st.divider()
    if st.button("Find Motifs", key="mo_find"):
        if not sequence: st.warning("Please provide a sequence.")
        elif not motif_pattern: st.warning("Please enter a motif/regex.")
        else:
            seq_cleaned = "".join(sequence.split())
            if not seq_cleaned: st.warning("Sequence is empty after cleaning.")
            else:
                try:
                    regex_flags = 0 if case_sensitive else re.IGNORECASE; found_matches = []
                    st.info(f"Searching for '**{motif_pattern}**' in '**{seq_id}**'...")
                    with st.spinner("Searching..."):
                        compiled_pattern = re.compile(motif_pattern, flags=regex_flags)
                        if allow_overlap: found_matches = list(compiled_pattern.finditer(seq_cleaned))
                        else:
                            current_pos = 0
                            while current_pos < len(seq_cleaned):
                                match = compiled_pattern.search(seq_cleaned, current_pos)
                                if match: found_matches.append(match); current_pos = match.end()
                                else: break
                    st.subheader(f"Search Results for '**{motif_pattern}**' in '**{seq_id}**'")
                    if found_matches:
                        st.success(f"Found **{len(found_matches):,}** match(es).")
                        match_data = [{"#": i + 1, "Start (1-based)": m.start() + 1, "End (1-based)": m.end(), "Length": m.end() - m.start(), "Matched Sequence": m.group()} for i, m in enumerate(found_matches)]
                        match_df = pd.DataFrame(match_data)
                        st.dataframe(match_df.style.format({"Start (1-based)":"{:,}", "End (1-based)":"{:,}", "Length":"{:,}"}), hide_index=True, use_container_width=True)

                        st.divider()
                        if st.checkbox("Highlight Matches (first 2000 chars)", value=False, key="mo_hi"):
                             highlight_limit = 2000; sequence_to_highlight = seq_cleaned[:highlight_limit]; highlighted_html = ""; last_end = 0
                             mark_style = "background-color:#fff3cd; padding: 1px 3px; border-radius: 3px; color:#856404; font-weight:bold;" # Match warning alert style
                             mark_open = f"<mark style='{mark_style}'>" ; mark_close = "</mark>"
                             sorted_matches = sorted(found_matches, key=lambda m: m.start())
                             for match in sorted_matches:
                                 if match.start() >= highlight_limit: break
                                 start, end = match.start(), min(match.end(), highlight_limit)
                                 if start >= last_end:
                                     text_segment = sequence_to_highlight[last_end:start].replace("&","&").replace("<","<").replace(">",">")
                                     highlighted_html += text_segment
                                     match_content = sequence_to_highlight[start:end].replace("&","&").replace("<","<").replace(">",">")
                                     highlighted_html += mark_open + match_content + mark_close
                                     last_end = end
                                 elif end > last_end:
                                     overlap_content = sequence_to_highlight[last_end:end].replace("&","&").replace("<","<").replace(">",">")
                                     highlighted_html += mark_open + overlap_content + mark_close
                                     last_end = end
                             remaining_text = sequence_to_highlight[last_end:].replace("&","&").replace("<","<").replace(">",">")
                             highlighted_html += remaining_text
                             st.markdown(f"**Highlighted Sequence (first {highlight_limit:,} chars):**")
                             st.markdown(f"<div style='font-family: monospace, Courier, mono; word-wrap:break-word; line-height:1.6; border: 1px solid #dee2e6; padding: 10px; border-radius: 5px; background-color: #f8f9fa; color: #333333;'>{highlighted_html}{'...' if len(seq_cleaned) > highlight_limit else ''}</div>", unsafe_allow_html=True)
                    else: st.info("Motif/pattern not found.")
                except re.error as regex_e: st.error(f"Invalid Regex: {regex_e}")
                except Exception as e: st.error(f"Error during motif search: {e}")


elif menu == "Bioinformatics Tool (Transcription/Translation)":
    st.header("DNA Transcription and Translation Tool")
    st.markdown("Transcribe DNA (coding strand, A/T/C/G only) to RNA, and translate RNA to protein (standard code).")

    input_method = st.radio("Input:",("Paste","FASTA"), key="tr_in", horizontal=True)
    dna_sequence, seq_id="","Pasted Sequence"
    if input_method=="Paste": dna_sequence=st.text_area("DNA Sequence (Coding Strand, 5'→3'):", height=120, key="tr_seq", placeholder="Paste DNA sequence...")
    else: uploaded_file=st.file_uploader("FASTA:", type=['fasta','fa','fna'], key="tr_f"); dna_sequence, seq_id = parse_fasta(uploaded_file) if uploaded_file else (None, None)
    if dna_sequence: st.text_area("Preview:", value=f">**{seq_id}**\n{dna_sequence[:80]}{'...' if len(dna_sequence)>80 else ''}", height=75, key="tr_d", disabled=True)

    st.divider()
    if st.button("Transcribe and Translate", key="tr_btn"):
        if not dna_sequence: st.warning("Please provide DNA sequence.")
        else:
            dna_cleaned = "".join(re.findall(r'[ATCG]', dna_sequence.upper()))
            if not dna_cleaned: st.error("Input contains no valid DNA bases (A, T, C, G).")
            else:
                try:
                    st.subheader(f"Results for: '**{seq_id}**'")
                    dna_seq_obj = Seq(dna_cleaned); dna_length = len(dna_seq_obj)
                    col1, col2 = st.columns(2)
                    with col1: st.metric("DNA Length", f"{dna_length:,} bp"); st.write("**Input DNA:**"); st.code(str(dna_seq_obj), language='text')
                    with col2:
                        try: gc_cont = gc_fraction(dna_seq_obj)*100 if dna_length > 0 else 0; st.metric("GC Content", f"{gc_cont:.1f}%")
                        except: st.metric("GC Content", "N/A")
                        st.write("**Reverse Comp.:**"); st.code(str(dna_seq_obj.reverse_complement()), language='text')

                    st.divider(); st.write("#### Transcription (DNA → RNA)")
                    rna_seq_obj = dna_seq_obj.transcribe()
                    st.write("**RNA Sequence:**"); st.code(str(rna_seq_obj), language='text'); st.metric("RNA Length", f"{len(rna_seq_obj):,} bases")

                    st.divider(); st.write("#### Translation (RNA → Protein)")
                    st.caption("Standard genetic code. '*' = STOP.")
                    remainder = len(rna_seq_obj) % 3; rna_for_translation = rna_seq_obj; warn_incomplete = False
                    if len(rna_seq_obj) < 3: st.error("Cannot translate: RNA < 3 bases."); rna_for_translation = None
                    elif remainder != 0: warn_incomplete = True; rna_for_translation = rna_seq_obj[:-remainder]; if not rna_for_translation: st.error("Cannot translate: Effective RNA length zero."); rna_for_translation = None

                    if rna_for_translation:
                        protein_seq_to_stop = rna_for_translation.translate(table=1, to_stop=True, cds=False)
                        st.write("**Protein (until STOP):**"); st.code(str(protein_seq_to_stop), language='text'); st.metric("Length (until STOP)", f"{len(protein_seq_to_stop):,} aa")
                        if warn_incomplete: st.warning(f"RNA length not multiple of 3; last {remainder} base(s) ignored.")
                        st.divider()
                        if st.checkbox("Show Full Translation (incl. internal STOPs)", value=False, key="tr_full"):
                            protein_seq_full = rna_for_translation.translate(table=1, to_stop=False, cds=False)
                            st.write("**Full Protein:**"); st.code(str(protein_seq_full), language='text'); st.metric("Full Length", f"{len(protein_seq_full):,} aa")
                except Exception as e: st.error(f"Processing error: {e}")


elif menu == "Genome Coverage Plotter":
    st.header("Genome Coverage Plotter")
    st.markdown("Visualize sequencing coverage depth. Upload position/coverage data (e.g., `samtools depth`, BEDGRAPH, CSV/TSV).")
    st.info("Requires columns for position and coverage depth.")

    uploaded_file = st.file_uploader("Upload Coverage Data:", type=['csv','tsv','txt','bed', 'bedgraph', 'depth'], key="cov_up")

    if uploaded_file:
        st.divider(); st.subheader("Parsing Settings")
        file_extension = os.path.splitext(uploaded_file.name)[-1].lower()
        default_separator = '\t' if file_extension != '.csv' else ','
        col1, col2 = st.columns(2)
        with col1:
            separator_options = {',': 'Comma (,)', '\t': 'Tab (\\t)', ' ': 'Whitespace ( )'}; sep_keys = list(separator_options.keys())
            default_sep_index = sep_keys.index(default_separator) if default_separator in sep_keys else 1
            selected_separator_key = st.selectbox("Separator:", sep_keys, index=default_sep_index, format_func=lambda x: separator_options[x], key="cov_sep")
            separator_regex = r'\s+' if selected_separator_key == ' ' else selected_separator_key
        with col2:
            default_has_header = file_extension in ['.csv', '.tsv', '.txt']; has_header = st.checkbox("Has header row", value=default_has_header, key="cov_head")
            header_arg = 0 if has_header else None; comment_char = st.text_input("Comment Char:", '#', max_chars=1, key="cov_comm")

        st.divider(); st.subheader("Data Loading & Columns")
        try:
            bytes_data = uploaded_file.getvalue(); s_buf = io.BytesIO(bytes_data)
            df_coverage = pd.read_csv(s_buf, sep=separator_regex, header=header_arg, engine='python', comment=comment_char if comment_char else None, low_memory=False)
            if df_coverage.empty: st.error("File empty or only comments/header."); st.stop()
            if header_arg is None and all(isinstance(c, int) for c in df_coverage.columns):
                num_cols = df_coverage.shape[1]
                if num_cols >= 2: df_coverage.columns = [f'Column_{i+1}' for i in range(num_cols)]; st.info(f"Generic column names assigned: {list(df_coverage.columns)}")
                else: st.error("Need >= 2 columns."); st.stop()
            elif df_coverage.shape[1] < 2: st.error("Need >= 2 columns."); st.stop()

            st.write("Data Preview:"); st.dataframe(df_coverage.head(), height=200, use_container_width=True)
            available_columns = list(df_coverage.columns)
            normalized_columns = {str(c).strip().lower(): c for c in available_columns}
            pos_guess, cov_guess, chr_guess = None, None, None
            chr_keywords = ['chr', 'chrom', '#chr', '#chrom', 'contig']; pos_keywords = ['pos', 'position', 'start', 'base']; cov_keywords = ['depth', 'cov', 'coverage', 'score', 'reads']
            for key in chr_keywords:   if key in normalized_columns: chr_guess = normalized_columns[key]; break
            for key in pos_keywords:   if key in normalized_columns: pos_guess = normalized_columns[key]; break
            for key in cov_keywords:   if key in normalized_columns: cov_guess = normalized_columns[key]; break
            chr_default_index = available_columns.index(chr_guess) if chr_guess in available_columns else 0
            pos_default_index = available_columns.index(pos_guess) if pos_guess in available_columns else 1
            cov_default_index = available_columns.index(cov_guess) if cov_guess in available_columns else (3 if len(available_columns) > 3 and file_extension in ['.bedgraph', '.bed'] else 2)
            chr_default_index = min(chr_default_index, len(available_columns)-1); pos_default_index = min(pos_default_index, len(available_columns)-1); cov_default_index = min(cov_default_index, len(available_columns)-1)

            col_s1, col_s2, col_s3 = st.columns(3)
            with col_s1: chromosome_column = st.selectbox("Chrom/Contig Col (Optional):", [None] + available_columns, index=(available_columns.index(chr_guess)+1 if chr_guess else 0), key="cov_chr")
            with col_s2: position_column = st.selectbox("Position Col:", available_columns, index=pos_default_index, key="cov_pos")
            with col_s3: coverage_column = st.selectbox("Coverage Col:", available_columns, index=cov_default_index, key="cov_cov")

            if position_column == coverage_column: st.error("Position/Coverage cols must differ.")
            elif chromosome_column and (chromosome_column == position_column or chromosome_column == coverage_column): st.error("Chrom/Pos/Coverage cols must differ.")
            else:
                st.success(f"Using '{position_column}' (Pos), '{coverage_column}' (Cov)" + (f", '{chromosome_column}' (Chrom)" if chromosome_column else "."))
                st.divider(); st.subheader("Filtering & Plotting")
                try:
                    cols_to_use = [position_column, coverage_column]; if chromosome_column: cols_to_use.insert(0, chromosome_column)
                    plot_data_full = df_coverage[cols_to_use].copy()
                    new_names = ['Position', 'Coverage']; if chromosome_column: new_names.insert(0, 'Chromosome')
                    plot_data_full.columns = new_names; selected_chr = None
                    if chromosome_column:
                        unique_chrs = sorted(plot_data_full['Chromosome'].astype(str).unique()); filter_options = ['All'] + unique_chrs
                        selected_chr = st.selectbox("Filter by Chromosome:", filter_options, index=0, key="cov_filter_chr")
                        if selected_chr != 'All': plot_data = plot_data_full[plot_data_full['Chromosome'].astype(str) == str(selected_chr)].copy(); st.info(f"Filtered for **{selected_chr}**.")
                        else: plot_data = plot_data_full.copy()
                    else: plot_data = plot_data_full.copy()

                    plot_data['Position'] = pd.to_numeric(plot_data['Position'], errors='coerce'); plot_data['Coverage'] = pd.to_numeric(plot_data['Coverage'], errors='coerce')
                    rows_before = len(plot_data); plot_data.dropna(subset=['Position', 'Coverage'], inplace=True); rows_after = len(plot_data)
                    if rows_after < rows_before: st.warning(f"Removed **{rows_before - rows_after:,}** rows with invalid Pos/Cov data.")

                    if plot_data.empty: st.error("No valid data left for plotting after cleaning/filtering."); st.stop()

                    plot_data = plot_data.sort_values(by='Position').reset_index(drop=True)
                    min_pos, max_pos = plot_data['Position'].min(), plot_data['Position'].max()
                    title_suffix = f" ({selected_chr})" if selected_chr and selected_chr != 'All' else f" (All Chrom)" if chromosome_column else ""
                    st.write(f"Plotting positions **{min_pos:,.0f}** to **{max_pos:,.0f}**{title_suffix}.")

                    col_p1, col_p2 = st.columns(2)
                    with col_p1: fill_area = st.checkbox("Fill Area", value=True, key="cov_fill"); log_scale_y = st.checkbox("Log Y-axis", value=False, key="cov_log")
                    with col_p2:
                        max_smooth = min(201, max(1, len(plot_data)//2)); smoothing_window = st.slider("Smoothing Window (0=None):", 0, max_smooth, 0, 2, key="cov_smooth", help="Rolling mean window.")
                        if smoothing_window > 0 and smoothing_window % 2 == 0: smoothing_window += 1; st.caption(f"Adjusted smoothing to odd: {smoothing_window}")

                    fig, ax = plt.subplots(figsize=(14, 5))
                    ax.plot(plot_data['Position'], plot_data['Coverage'], label='Raw Coverage', linewidth=0.7, color='#007bff', alpha=0.8) # Blue line
                    if fill_area: ax.fill_between(plot_data['Position'], plot_data['Coverage'], alpha=0.2, color='#80ccff') # Light blue fill
                    if smoothing_window > 1:
                        plot_data['Smoothed_Coverage'] = plot_data['Coverage'].rolling(smoothing_window, center=True, min_periods=1).mean()
                        ax.plot(plot_data['Position'], plot_data['Smoothed_Coverage'], label=f'Smoothed ({smoothing_window}bp)', linewidth=1.8, color='#ff7f0e', linestyle='-') # Orange smooth line

                    plot_title = f"Genome Coverage: {uploaded_file.name}{title_suffix}"; plot_xlabel = f"Position ({position_column})"; plot_ylabel = f"Depth ({coverage_column}){ ' (Log)' if log_scale_y else ''}"
                    style_plot(fig, ax, title=plot_title, xlabel=plot_xlabel, ylabel=plot_ylabel) # Apply light theme style

                    if log_scale_y: ax.set_yscale('log'); min_pos_cov = plot_data[plot_data['Coverage'] > 0]['Coverage'].min(); ax.set_ylim(bottom=max(0.1, min_pos_cov * 0.5) if pd.notna(min_pos_cov) and min_pos_cov > 0 else 0.1)
                    else: ax.set_ylim(bottom=0)
                    ax.ticklabel_format(style='plain', axis='x', useOffset=False); ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ',')))
                    legend = ax.legend(facecolor='#ffffff', labelcolor='#333333', framealpha=0.8, fontsize=10)
                    for text in legend.get_texts(): text.set_fontfamily('Times New Roman'); text.set_fontweight('bold')
                    st.pyplot(fig)

                    st.divider(); st.subheader(f"Coverage Stats (Raw{title_suffix})")
                    coverage_stats = plot_data['Coverage'].describe()
                    scol1, scol2, scol3, scol4 = st.columns(4)
                    with scol1: st.metric("Mean", f"{coverage_stats.get('mean', 0):,.1f}x")
                    with scol2: st.metric("Median", f"{coverage_stats.get('50%', 0):,.1f}x")
                    with scol3: st.metric("Min", f"{coverage_stats.get('min', 0):,.0f}x")
                    with scol4: st.metric("Max", f"{coverage_stats.get('max', 0):,.0f}x")

                    max_cov_int = int(coverage_stats.get('max', 1000)); default_thresh = min(10, max(1, int(coverage_stats.get('mean', 10)))); max_thresh = max(1, max_cov_int)
                    coverage_threshold = st.number_input("Breadth of Coverage Depth ≥:", 0, max_thresh, default_thresh, 1, key="cov_b_thresh")
                    if coverage_threshold >= 0:
                        positions_above = (plot_data['Coverage'] >= coverage_threshold).sum(); total_positions = len(plot_data)
                        if total_positions > 0: breadth_percent = (positions_above / total_positions * 100); st.metric(f"Breadth (≥ {coverage_threshold}x)", f"{breadth_percent:.1f}%", f"{positions_above:,} / {total_positions:,} pos")
                        else: st.metric(f"Breadth (≥ {coverage_threshold}x)", "N/A", "No data")
                except Exception as plot_e: st.error(f"Filtering/Plotting error: {plot_e}")
        except Exception as read_e: st.error(f"File reading error: {read_e}")


elif menu == "Variant Annotation Tool":
    st.header("Simple Variant Substitution Tool")
    st.markdown("Introduce a single nucleotide polymorphism (SNP) into reference DNA (A/T/C/G only) and see codon effect.")

    input_method = st.radio("Input Ref DNA:",("Paste","FASTA"), key="var_in", horizontal=True)
    ref_sequence, seq_id="","Reference Sequence"
    if input_method=="Paste": ref_sequence=st.text_area("Reference DNA:", height=100, key="var_ref", placeholder="Paste reference DNA...")
    else: uploaded_file=st.file_uploader("Reference FASTA:", type=['fasta','fa','fna'], key="var_f"); ref_sequence, seq_id = parse_fasta(uploaded_file) if uploaded_file else (None, None)
    if ref_sequence: st.text_area("Ref Preview:", value=f">**{seq_id}**\n{ref_sequence[:80]}{'...' if len(ref_sequence)>80 else ''}", height=75, key="var_d", disabled=True)

    st.divider(); st.subheader("Variant Details")
    vcol1, vcol2 = st.columns(2)
    with vcol1: variant_pos = st.number_input("Variant Position (1-based):", 1, step=1, value=1, key="var_p")
    with vcol2: variant_base = st.selectbox("New Base (Alt):", ["A", "T", "C", "G"], key="var_b")

    st.divider()
    if st.button("Apply Variant and Annotate", key="var_btn"):
        if not ref_sequence: st.warning("Please provide reference DNA.")
        else:
            ref_cleaned = "".join(re.findall(r'[ATCG]', ref_sequence.upper())); error_occurred = False
            if not ref_cleaned: st.error("Ref sequence has no valid DNA bases."); error_occurred = True
            else:
                ref_length = len(ref_cleaned); zero_based_pos = variant_pos - 1
                if zero_based_pos < 0 or zero_based_pos >= ref_length: st.error(f"Position {variant_pos:,} outside sequence range (1-{ref_length:,})."); error_occurred = True
            if not error_occurred:
                try:
                    original_base = ref_cleaned[zero_based_pos]
                    if original_base == variant_base: st.info(f"No change: Base at position **{variant_pos:,}** ('**{seq_id}**') is already '**{original_base}**'."); st.write("**Reference:**"); st.code(ref_cleaned, language='text')
                    else:
                        alt_list = list(ref_cleaned); alt_list[zero_based_pos] = variant_base; alt_sequence = "".join(alt_list)
                        st.subheader("Variant Applied")
                        st.write(f"**Sequence:** '**{seq_id}**' | **Variant:** Pos **{variant_pos:,}**")
                        st.markdown(f"**Change:** Ref (`{original_base}`) → Alt (`{variant_base}`)", unsafe_allow_html=True)
                        s_col1, s_col2 = st.columns(2)
                        with s_col1: st.write("**Reference:**"); st.code(ref_cleaned, language='text')
                        with s_col2: st.write("**Altered:**"); st.code(alt_sequence, language='text')

                        st.divider(); st.write("#### Codon Context & Effect:")
                        codon_start_index = (zero_based_pos // 3) * 3; position_in_codon = (zero_based_pos % 3) + 1
                        if codon_start_index + 3 <= ref_length:
                            ref_codon_str = ref_cleaned[codon_start_index : codon_start_index + 3]
                            alt_codon_str = alt_sequence[codon_start_index : codon_start_index + 3]
                            st.write(f"Variant at pos **{variant_pos:,}** is position **{position_in_codon}** in codon starting at base **{codon_start_index + 1:,}**.")
                            try:
                                ref_codon_obj, alt_codon_obj = Seq(ref_codon_str), Seq(alt_codon_str)
                                ref_aa, alt_aa = str(ref_codon_obj.translate(table=1, cds=False)), str(alt_codon_obj.translate(table=1, cds=False))
                                c_col1, c_col2 = st.columns(2)
                                with c_col1: st.write("**Ref Codon:**"); st.code(f"Codon: {ref_codon_str}\nAA:    {ref_aa}", language='text')
                                with c_col2: st.write("**Alt Codon:**"); st.code(f"Codon: {alt_codon_str}\nAA:    {alt_aa}", language='text')
                                effect = "Unknown"
                                if ref_aa == alt_aa: effect = "✅ **Silent (Synonymous)**"
                                elif alt_aa == '*': effect = "🛑 **Nonsense (Stop Gained)**" if ref_aa != '*' else "❓ **Stop Retained**"
                                elif ref_aa == '*': effect = "➡️ **Stop-Lost**"
                                else: effect = "🔄 **Missense (Non-synonymous)**"
                                st.markdown(f"**Predicted Effect:** {effect}")
                            except Exception as codon_e: st.warning(f"Could not translate codons: {codon_e}")
                        else: st.info("Variant near sequence end; full codon context unavailable.")
                except Exception as e: st.error(f"Error applying variant: {e}")


elif menu == "Codon Usage Analyzer":
    st.header("Codon Usage Analyzer")
    st.markdown("Analyze codon frequency in a Coding DNA Sequence (CDS). Requires A/T/C/G only, length multiple of 3.")

    input_method = st.radio("Input CDS:",("Paste","FASTA"), key="cod_in", horizontal=True)
    cds_sequence, seq_id="","CDS Sequence"
    if input_method=="Paste": cds_sequence=st.text_area("CDS:", height=120, key="cod_seq", placeholder="Paste CDS (multiple of 3, ATCG only)...")
    else: uploaded_file=st.file_uploader("FASTA:", type=['fasta','fa','fna'], key="cod_f"); cds_sequence, seq_id = parse_fasta(uploaded_file) if uploaded_file else (None, None)
    if cds_sequence: st.text_area("Preview:", value=f">**{seq_id}**\n{cds_sequence[:80]}{'...' if len(cds_sequence)>80 else ''}", height=75, key="cod_d", disabled=True)

    st.divider()
    if st.button("Analyze Codon Usage", key="cod_btn"):
        if not cds_sequence: st.warning("Please provide CDS.")
        else:
            cds_cleaned = "".join(re.findall(r'[ATCG]', cds_sequence.upper())); error_occurred = False
            if not cds_cleaned: st.error("Sequence has no valid DNA bases."); error_occurred = True
            elif len(cds_cleaned) == 0: st.error("Sequence empty."); error_occurred = True
            elif len(cds_cleaned) % 3 != 0: st.error(f"Length (**{len(cds_cleaned):,}**) not multiple of 3."); error_occurred = True
            if not error_occurred:
                try:
                    st.subheader(f"Codon Usage for: '**{seq_id}**'")
                    codons_list = [cds_cleaned[i:i+3] for i in range(0, len(cds_cleaned), 3)]
                    total_codons = len(codons_list); codon_counts = Counter(codons_list)
                    st.metric("Total Codons", f"{total_codons:,}")
                    try: standard_table = CodonTable.unambiguous_dna_by_id[1]; all_possible = list(standard_table.forward_table.keys()) + standard_table.stop_codons; aa_map = standard_table.forward_table; stop_codons = set(standard_table.stop_codons)
                    except Exception as table_e: st.error(f"Failed to load codon table: {table_e}"); st.stop()
                    usage_data = []
                    for codon in sorted(all_possible):
                        count = codon_counts.get(codon, 0); frequency = (count / total_codons * 100) if total_codons > 0 else 0
                        amino_acid = aa_map.get(codon, 'Stop' if codon in stop_codons else '?')
                        usage_data.append({"Codon": codon, "Amino Acid (AA)": amino_acid, "Count": count, "Frequency (%)": frequency})
                    usage_df = pd.DataFrame(usage_data)
                    st.write("#### Codon Usage Table:"); st.dataframe(usage_df.style.format({'Count':'{:,}', 'Frequency (%)':'{:.1f}%'}), hide_index=True, use_container_width=True)

                    st.divider(); st.subheader("Visualization")
                    plot_choice = st.selectbox("Plot Type:", ["Frequency of All Codons", "Relative Usage for Specific Amino Acid"], key="cod_plot")
                    present_codons_df = usage_df[usage_df['Count'] > 0].copy()
                    if present_codons_df.empty: st.info("No codons found to visualize.")
                    else:
                        if plot_choice == "Frequency of All Codons":
                            plot_df_all = present_codons_df.sort_values("Codon"); fig, ax = plt.subplots(figsize=(16, 6))
                            unique_aas = sorted(list(plot_df_all['Amino Acid (AA)'].unique())); cmap = plt.get_cmap('tab20', len(unique_aas))
                            aa_colors = {aa: cmap(i) for i, aa in enumerate(unique_aas)}; bar_colors = [aa_colors.get(aa, '#888888') for aa in plot_df_all['Amino Acid (AA)']]
                            ax.bar(plot_df_all['Codon'], plot_df_all['Frequency (%)'], color=bar_colors)
                            style_plot(fig, ax, title=f"Codon Frequency for '{seq_id}'", xlabel="Codon", ylabel="Frequency (%)") # Apply light theme style
                            ax.tick_params(axis='x', rotation=90, labelsize=9)
                            for label in ax.get_xticklabels(): label.set_fontfamily('Times New Roman'); label.set_fontweight('bold')
                            st.pyplot(fig)
                        elif plot_choice == "Relative Usage for Specific Amino Acid":
                             available_aas = sorted([aa for aa in present_codons_df['Amino Acid (AA)'].unique() if aa != '?'])
                             if not available_aas: st.info("No standard AAs/Stop encoded.")
                             else:
                                 selected_aa = st.selectbox("Select AA or 'Stop':", available_aas, key="cod_aa")
                                 aa_specific_df = present_codons_df[present_codons_df['Amino Acid (AA)'] == selected_aa].copy()
                                 total_count_for_aa = aa_specific_df['Count'].sum()
                                 if total_count_for_aa > 0 and not aa_specific_df.empty:
                                     aa_specific_df['Relative Freq (%)'] = (aa_specific_df['Count'] / total_count_for_aa * 100)
                                     aa_specific_df = aa_specific_df.sort_values("Codon")
                                     fig, ax = plt.subplots(figsize=(max(6, len(aa_specific_df)*1.5), 5))
                                     ax.bar(aa_specific_df['Codon'], aa_specific_df['Relative Freq (%)'], color='#ff7f0e') # Orange bars
                                     style_plot(fig, ax, title=f"Relative Usage for '{selected_aa}' in '{seq_id}' ({total_count_for_aa:,} total)", xlabel="Codon", ylabel="Relative Freq (%)") # Apply light theme style
                                     ax.set_ylim(0, 105); st.pyplot(fig)
                                     st.write(f"**Details for '{selected_aa}':**")
                                     st.dataframe(aa_specific_df[['Codon', 'Count', 'Relative Freq (%)']].style.format({'Count':'{:,}', 'Relative Freq (%)':'{:.1f}%'}), hide_index=True, use_container_width=True)
                                 else: st.info(f"No codons for '{selected_aa}' found.")
                except Exception as e: st.error(f"Codon usage analysis error: {e}")


# --- Footer ---
st.markdown("---")
st.markdown('<div class="footer">Copyright © 2024 CCDB Tools | For Educational and Informational Purposes Only</div>', unsafe_allow_html=True)
