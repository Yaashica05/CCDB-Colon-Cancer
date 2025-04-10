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
            uploaded_file.seek(0); content_bytes = uploaded_file.getvalue()
            try: content_str = content_bytes.decode("utf-8")
            except UnicodeDecodeError: content_str = content_bytes.decode("latin-1")
            stringio = io.StringIO(content_str); records = list(SeqIO.parse(stringio, "fasta"))
            if records: first_record = records[0]; sequence = str(first_record.seq).upper(); seq_id = first_record.id; st.success(f"Read FASTA: **{seq_id}** ({len(sequence):,} bases/AAs).")
            else: st.error("FASTA empty/invalid.")
        except Exception as e: st.error(f"FASTA Read Error: {e}")
    return sequence, seq_id

# --- Helper function to create links for the Database Search ---
def format_link(value, column_name):
    """Creates an HTML link based on column name and value."""
    if pd.isna(value): return ""; s_value = str(value).strip();
    if not s_value: return ""
    col_name_std = column_name.strip().title(); link_style = 'style="color:#0056b3 !important;"' # Light theme link color
    try:
        if s_value.startswith("http"): display_url = s_value if len(s_value) < 50 else s_value[:47] + "..."; return f'<a href="{s_value}" target="_blank" {link_style} title="{s_value}">{display_url}</a>'
        elif col_name_std == 'Pubmed Id':
            if s_value.replace('.', '', 1).isdigit(): url = f"https://pubmed.ncbi.nlm.nih.gov/{s_value}"; return f'<a href="{url}" target="_blank" {link_style}>{s_value}</a>'
        elif col_name_std == 'Doi Id':
            if '/' in s_value: encoded_doi = urllib.parse.quote(s_value, safe='/:'); url = f"https://doi.org/{encoded_doi}"; return f'<a href="{url}" target="_blank" {link_style}>{s_value}</a>'
        elif col_name_std == 'Uniprot':
             if re.match(r'^[A-Z0-9_.\-]+$', s_value, re.IGNORECASE): url = f"https://www.uniprot.org/uniprotkb/{s_value}/entry"; return f'<a href="{url}" target="_blank" {link_style}>{s_value}</a>'
        else: return s_value
    except Exception: return s_value

# --- Helper function for Matplotlib Plot Styling ---
def style_plot(fig, ax, title="", xlabel="", ylabel=""):
    """Applies consistent LIGHT theme styling (Times New Roman Bold) to a Matplotlib plot."""
    fig.patch.set_facecolor('#f8f9fa'); fig.patch.set_alpha(1.0); ax.set_facecolor('#ffffff')
    font_props = {'family': 'Times New Roman', 'weight': 'bold'}; text_color = '#1E2A3A'; accent_color = '#004080'
    ax.set_title(title, color=accent_color, fontsize=14, **font_props)
    ax.set_xlabel(xlabel, color=text_color, fontsize=12, **font_props); ax.set_ylabel(ylabel, color=text_color, fontsize=12, **font_props)
    tick_color = '#333333'; ax.tick_params(axis='x', colors=tick_color, labelsize=10); ax.tick_params(axis='y', colors=tick_color, labelsize=10)
    for label in ax.get_xticklabels() + ax.get_yticklabels(): label.set_fontfamily('Times New Roman'); label.set_fontweight('bold'); label.set_color(tick_color)
    ax.grid(True, linestyle=':', alpha=0.6, color='#cccccc'); [spine.set_edgecolor('#dee2e6') for spine in ax.spines.values()]; plt.tight_layout()

# --- Main App Title ---
st.title("🧬 CCDB: Colon Cancer Database and Bioinformatics Tools 🔬")

# --- Sidebar Navigation ---
with st.sidebar:
    st.header("🛠️ Toolbox Menu")
    menu = st.radio("Select Feature:", ["Home", "Colon Cancer Database Search", "DNA Sequence Analysis", "Protein Sequence Analysis", "Primer Design", "Restriction Enzyme Analysis", "Pairwise Sequence Alignment", "Motif Finder Tool", "Bioinformatics Tool (Transcription/Translation)", "Genome Coverage Plotter", "Variant Annotation Tool", "Codon Usage Analyzer"], key="main_menu", help="Choose a tool.")
    st.divider(); st.info("Ensure input files are available.")


# --- Page Content based on Menu Selection ---

if menu == "Home": # Content unchanged
    st.header("Welcome"); st.markdown(
        """
        Access curated colon cancer gene data & bioinformatics tools. User-friendly platform for research, education, & clinical interest.
        ---
        ### **Understanding Colon Cancer (CRC)**
        CRC develops in the large intestine (colon/rectum). Common worldwide. Starts as polyps; some become cancerous. Adenomas are common pre-cancerous polyps. Can spread (metastasize) if not detected early.
        ### **Key Info:**
        *   **Risk Factors:** Age (>50), Personal/Family History (CRC, polyps, IBD), Inherited Syndromes (Lynch, FAP), Lifestyle (diet, obesity, inactivity, smoking, alcohol), Type 2 Diabetes.
        *   **Symptoms (May be absent early):** Bowel changes, rectal bleeding/blood in stool, abdominal discomfort, incomplete emptying feeling, unexplained weight loss/fatigue.
        *   **Screening:** Crucial. Colonoscopy (gold standard), Stool tests (FIT, FOBT, DNA), Sigmoidoscopy, CT Colonography. *Consult provider.*
        *   **Genetics:** *APC, KRAS, BRAF, TP53*, Mismatch Repair (*MLH1, MSH2*, etc.) often involved.
        *   **Treatment:** Surgery, Chemo, Radiation, Targeted Therapy (EGFR, VEGF), Immunotherapy (MSI-High). Depends on stage/factors.
        ---
        ### **Available Tools:**
        Use sidebar: Database Search, DNA Analysis, Protein Analysis, Primer Design, Restriction Analysis, Alignment, Motif Finder, Transcription/Translation, Coverage Plotter, Variant Tool, Codon Usage.
        ---
        *Disclaimer: Educational/informational only. Not medical advice. Consult a healthcare provider.*
        """, unsafe_allow_html=True)
    st.image("https://img.freepik.com/free-vector/colorectal-cancer-crc-infographic-education_1308-47971.jpg", caption="CRC Infographic (Source: Freepik)", use_container_width=True)
    st.markdown("""<a href="https://www.cancer.org/cancer/types/colon-rectal-cancer.html" target="_blank" style="color:#0056b3;">Learn more (American Cancer Society)</a>""", unsafe_allow_html=True)
    st.divider(); st.subheader("Feedback"); feedback=st.text_area("Share feedback:",key="fb_home",placeholder="Thoughts, suggestions...");
    if st.button("Submit",key="sub_fb_home"): st.success("Thank you!") if feedback else st.warning("Feedback empty.")


elif menu == "Colon Cancer Database Search": # Link list updated
    st.header("Colon Cancer Gene DB Search"); st.markdown(f"Search local file (`{DATA_PATH}`). Enter gene symbol.")
    if not os.path.exists(DATA_PATH): st.error(f"`{DATA_PATH}` missing."); st.stop()
    query=st.text_input("Gene Symbol:",key="gene_search",placeholder="e.g., APC, KRAS")
    if st.button("Search",key="gene_search_btn"):
        if not query: st.warning("Enter gene symbol.")
        else:
            try:
                data=pd.read_excel(DATA_PATH); data.columns=data.columns.str.strip().str.title(); gene_col=None; poss_cols=['Gene Symbol','Gene Name','Symbol','Gene']
                for col in poss_cols: if col in data.columns: gene_col=col; break
                if not gene_col: st.error(f"Gene col not found. Found:{list(data.columns)}"); st.stop()
                data[gene_col]=data[gene_col].astype(str).str.strip(); q=query.strip(); results=data[data[gene_col].str.fullmatch(q,case=False,na=False)].copy()
                if results.empty: st.info(f"No exact match for '{q}'. Trying contains..."); results=data[data[gene_col].str.contains(q,case=False,na=False,regex=False)].copy()
                if not results.empty:
                    st.success(f"Found {len(results)} result(s) for '{query}'."); st.write("### Results:")
                    fmt_res=results.copy(); link_cols=['Pubmed Id','Doi Id','Uniprot','Blast','Conserved Domain','Link','Url','Reference','Entry','History','Variant Viewer','Feature Viewer','Genomic Coordinates','Publications','External Links','Sequence']
                    link_cols.extend([c for c in fmt_res.columns if ('link' in c.lower() or 'url' in c.lower()) and c not in link_cols]); link_cols=list(set(link_cols))
                    for col in fmt_res.columns:
                        if col in link_cols:
                            try: fmt_res[col]=fmt_res[col].apply(lambda x: format_link(x, col))
                            except Exception as e: st.warning(f"Link format fail '{col}': {e}.")
                    html=fmt_res.to_html(escape=False,index=False,na_rep='-',justify='left',classes=['st-table']); st.write(html,unsafe_allow_html=True)
                else: st.warning(f"No results found for '{query}' in '{gene_col}'.")
            except ImportError: st.error("`openpyxl` needed. `pip install openpyxl`.")
            except FileNotFoundError: st.error(f"`{DATA_PATH}` not found.")
            except Exception as e: st.error(f"DB search error: {e}"); st.exception(e)

elif menu == "DNA Sequence Analysis": # Logic unchanged, styling updates via CSS/style_plot
    st.header("DNA Analysis"); st.markdown("Basic DNA properties (A/T/C/G/N).")
    in_m=st.radio("Input:",("Paste","FASTA"),key="dna_in",horizontal=True); seq,seq_id="","Seq"
    if in_m=="Paste": seq=st.text_area("DNA:",h=150,key="dna_seq_in",placeholder="Paste DNA...")
    else: up_f=st.file_uploader("FASTA:",type=['fasta','fa','fna'],key="dna_up"); seq,seq_id=parse_fasta(up_f) if up_f else (None,None)
    if seq: st.text_area("Preview:",value=f">**{seq_id}**\n{seq[:100]}{'...' if len(seq)>100 else ''}",h=100,key="dna_disp",disabled=True)
    st.divider()
    if st.button("Analyze",key="dna_analyze_btn"):
        if not seq: st.warning("Provide DNA.")
        else:
            seq_cl="".join(re.findall(r'[ATCGN]',seq.upper()))
            if not seq_cl: st.error("No valid DNA bases.")
            else:
                try:
                    st.subheader(f"Results: '**{seq_id}**'")
                    seq_o=Seq(seq_cl); l=len(seq_o); st.metric("Length",f"{l:,} bp")
                    seq_nn=seq_cl.replace('N',''); l_nn=len(seq_nn)
                    if l_nn>0: gc=gc_fraction(seq_nn)*100; at=100.0-gc; st.metric("GC% (no N)",f"{gc:.1f}%"); st.metric("AT% (no N)",f"{at:.1f}%")
                    else: st.info("GC/AT% N/A.")
                    st.divider(); st.write("#### Composition:"); comp=Counter(seq_cl); data=[{"Base":b,"Count":c,"%":(c/l*100) if l>0 else 0} for b,c in sorted(comp.items())]
                    df=pd.DataFrame(data); st.dataframe(df.style.format({"Count":"{:,}","%":"{:.1f}%"}),hide_index=True,use_container_width=True)
                    p_comp={k:v for k,v in comp.items() if k in STANDARD_DNA_BASES}
                    if p_comp:
                        fig,ax=plt.subplots(figsize=(6,4)); bases=sorted(p_comp.keys()); counts=[p_comp[b] for b in bases]
                        cols={'A':'#007bff','T':'#ff7f0e','C':'#2ca02c','G':'#d62728'}; b_cols=[cols.get(b,'#888') for b in bases]; ax.bar(bases,counts,color=b_cols)
                        style_plot(fig,ax,title="Counts (no N)",xlabel="Base",ylabel="Count"); ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x,p:format(int(x),','))); st.pyplot(fig)
                    elif 'N' in comp and len(comp)==1: st.info("Plot omitted (N only).")
                    elif not p_comp: st.info("Plot omitted (no std bases).")
                    st.divider(); st.write("#### ORF Finder"); min_aa=st.number_input("Min ORF (AA):",10,500,30,5,key="orf_len"); orfs=[]
                    st.info("Searching ORFs (M to *) fwd...");
                    with st.spinner("Scanning..."):
                         for frame in range(3):
                            try: t=str(seq_o[frame:].translate(1,to_stop=True,cds=False))
                            except: t=""
                            if not t: continue; pos=0
                            while 'M' in t[pos:]:
                                start_i=t.find('M',pos); if start_i==-1: break
                                p_orf=t[start_i:]
                                if len(p_orf)>=min_aa: s_dna=frame+(start_i*3)+1; e_dna=s_dna+(len(p_orf)*3)-1; if e_dna<=l and not any(o['Start (DNA)']==s_dna and o['Frame']==frame+1 for o in orfs): orfs.append({"Frame":frame+1,"Start (DNA)":s_dna,"End (DNA)":e_dna,"Len (AA)":len(p_orf),"Prot (Prev)":p_orf[:40]+("..." if len(p_orf)>40 else "")})
                                pos=start_i+1
                    if orfs: st.success(f"{len(orfs)} ORF(s) ≥ {min_aa} AA."); df=pd.DataFrame(orfs).sort_values(["Frame","Start (DNA)"]); st.dataframe(df.style.format({"Start (DNA)":"{:,}","End (DNA)":"{:,}","Len (AA)":"{:,}"}),hide_index=True,use_container_width=True)
                    else: st.info(f"No ORFs ≥ {min_aa} AA.")
                except Exception as e: st.error(f"DNA analysis error: {e}"); st.exception(e)

elif menu == "Protein Sequence Analysis": # Logic unchanged, styling updates via CSS/style_plot
    st.header("Protein Analysis"); st.markdown("Analyze properties.")
    in_m=st.radio("Input:",("Paste","FASTA"),key="prot_in",horizontal=True); seq,seq_id="","Seq"
    if in_m=="Paste": seq=st.text_area("Protein:",h=150,key="prot_seq",placeholder="Paste protein...")
    else: up_f=st.file_uploader("FASTA:",type=['fasta','fa','faa'],key="prot_fasta"); seq,seq_id=parse_fasta(up_f) if up_f else (None,None)
    if seq: st.text_area("Preview:",value=f">**{seq_id}**\n{seq[:100]}{'...' if len(seq)>100 else ''}",h=100,key="prot_disp",disabled=True)
    st.divider()
    if st.button("Analyze",key="prot_analyze_btn"):
        if not seq: st.warning("Provide protein.")
        else:
            seq_cl="".join(seq.split()).upper(); chars=set(seq_cl); std=set(STANDARD_AA); amb=set(AMBIGUOUS_AA)
            n_std=chars-std-amb; amb_p=chars.intersection(amb); seq_c="".join(c for c in seq_cl if c in std)
            if n_std: st.warning(f"Ignoring unknown: `{'`, `'.join(sorted(list(n_std)))}`")
            if amb_p: st.warning(f"Ignoring ambiguous: `{'`, `'.join(sorted(list(amb_p)))}`")
            if not seq_c: st.error("No standard AAs.")
            else:
                try:
                    st.subheader(f"Results: '**{seq_id}**'")
                    pa=ProteinAnalysis(seq_c); c1,c2=st.columns(2)
                    with c1: st.metric("In Len",f"{len(seq_cl):,} aa"); st.metric("Calc Len",f"{len(seq_c):,} aa"); st.metric("Mol Wt",f"{pa.molecular_weight():,.1f} Da"); st.metric("pI",f"{pa.isoelectric_point():.2f}")
                    with c2: gr=pa.gravy(); hy="Hydrophobic" if gr>0 else ("Hydrophilic" if gr<0 else "Neutral"); st.metric("GRAVY",f"{gr:.3f}",delta=hy,delta_color="off"); inst=pa.instability_index(); stab="Stable" if inst<40 else "Unstable"; st.metric("Instability",f"{inst:.2f}",delta=stab,delta_color="normal" if stab=="Stable" else "inverse")
                    st.divider(); st.write("#### AA Comp:"); comp=Counter(seq_cl); data=[{"AA":a,"Count":c,"%":(c/len(seq_cl)*100) if seq_cl else 0} for a,c in sorted(comp.items())]
                    df=pd.DataFrame(data); st.dataframe(df.style.format({'Count':'{:,}','%':'{:.1f}%'}),hide_index=True,use_container_width=True)
                    st.divider(); st.write("#### Sec Structure:")
                    try:
                        h,t,s=pa.secondary_structure_fraction(); labs='H','T','S'; sizes=[h*100,t*100,s*100]
                        if sum(sizes)>0.1: cols=['#007bff','#ff7f0e','#2ca02c']; fig,ax=plt.subplots(figsize=(6,4)); w,tx,au=ax.pie(sizes,labels=labs,startangle=90,colors=cols,autopct='%1.1f%%',wedgeprops={'edgecolor':'#fff','lw':1.5},textprops={'color':'#333'}); [a.set_color('#fff')or a.set_fontfamily('TN R')or a.set_fontweight('bold')or a.set_fontsize(10) for a in au]; [t.set_fontfamily('TN R')or t.set_fontweight('bold') for t in tx]; ax.axis('equal'); style_plot(fig,ax,title="Sec Structure"); ax.set_xlabel(""); ax.set_ylabel(""); st.pyplot(fig)
                        else: st.info("Structure fractions near zero.")
                    except Exception as ss_e: st.warning(f"Structure fail: {ss_e}")
                except Exception as e: st.error(f"Protein analysis error: {e}"); st.exception(e)

elif menu == "Primer Design": # Indentation fix applied
    st.header("Basic Primer Design"); st.markdown("Simple Fwd/Rev primers from ends (A/T/C/G) & Tm estimate."); st.warning("⚠️ **Edu Tool:** Basic, no validation.", icon="🔬")
    in_m=st.radio("Input:",("Paste","FASTA"),key="p_in",horizontal=True); seq,seq_id="","Seq"
    if in_m=="Paste": seq=st.text_area("DNA Template:",h=150,key="p_seq",placeholder="Paste DNA...")
    else: up_f=st.file_uploader("FASTA:",type=['fasta','fa','fna'],key="p_fasta"); seq,seq_id=parse_fasta(up_f) if up_f else (None,None)
    if seq: st.text_area("Preview:",value=f">**{seq_id}**\n{seq[:100]}{'...' if len(seq)>100 else ''}",h=100,key="p_disp",disabled=True)
    st.divider(); st.subheader("Params"); c1,c2=st.columns(2)
    with c1: p_len=st.slider("Length (bp):",15,35,20,1,key="p_len")
    with c2: tm_m=st.selectbox("Tm Method:",["Tm_NN","Tm_GC","Tm_Wallace"],0,key="p_tm_m"); dna,na,mg,dntp=50.0,50.0,0.0,0.0 # Defaults
        # *** CORRECTED INDENTATION ***
        if tm_m=="Tm_NN":
             with st.expander("Tm_NN Params"): dna=st.number_input("Primer(nM):",1.,1000.,50.,1.,key="p_dna",fmt="%.1f"); na=st.number_input("Na+(mM):",1.,200.,50.,1.,key="p_na",fmt="%.1f"); mg=st.number_input("Mg++(mM):",0.,50.,0.,0.1,key="p_mg",fmt="%.1f"); dntp=st.number_input("dNTPs(mM):",0.,10.,0.,0.1,key="p_dntp",fmt="%.1f")
    st.divider()
    if st.button("Design",key="p_design_btn"):
        if not seq: st.warning("Provide DNA.")
        else:
            seq_cl="".join(re.findall(r'[ATCG]',seq.upper()))
            if not seq_cl: st.error("No valid DNA bases.")
            elif len(seq_cl)<p_len*2: st.error(f"Seq too short ({len(seq_cl):,} bp) for {p_len} bp primers.")
            else:
                try:
                    fw=Seq(seq_cl[:p_len]); rv=Seq(seq_cl[-p_len:]).reverse_complement(); fw_gc=gc_fraction(fw)*100; rv_gc=gc_fraction(rv)*100; fw_tm,rv_tm="N/A","N/A"; tm_disp=f"Method: {tm_m}"
                    try: args={'strict':False};
                        if tm_m=="Tm_NN": args.update({'Na':na,'Mg':mg,'dNTPs':dntp,'dnac1':dna,'nn_table':MeltingTemp.DNA_NN4}); tm_disp+=f", P={dna:.1f}nM, Na={na:.1f}mM, Mg={mg:.1f}mM, dNTP={dntp:.1f}mM"
                        if hasattr(MeltingTemp,tm_m): func=getattr(MeltingTemp,tm_m); fw_tm=func(fw,**args); rv_tm=func(rv,**args)
                        else: st.error(f"Tm method '{tm_m}' missing.")
                    except Exception as tm_e: st.warning(f"Tm calc fail: {tm_e}")
                    st.subheader(f"Primers ('**{seq_id}**')"); p1,p2=st.columns(2)
                    with p1: st.markdown("#### Fwd (5'→3')"); st.code(str(fw)); st.metric("Len",f"{len(fw)} bp"); st.metric("GC%",f"{fw_gc:.1f}%"); st.metric("Tm",f"{fw_tm:.1f}°C" if isinstance(fw_tm,(float,int)) else "N/A")
                    with p2: st.markdown("#### Rev (5'→3')"); st.code(str(rv)); st.metric("Len",f"{len(rv)} bp"); st.metric("GC%",f"{rv_gc:.1f}%"); st.metric("Tm",f"{rv_tm:.1f}°C" if isinstance(rv_tm,(float,int)) else "N/A")
                    st.divider(); st.metric(f"Amplicon ('**{seq_id}**')",f"{len(seq_cl):,} bp")
                    if isinstance(fw_tm,(float,int)): st.caption(f"Tm calc: {tm_disp}")
                except Exception as e: st.error(f"Design error: {e}")

elif menu == "Restriction Enzyme Analysis": # Syntax fix applied
    st.header("Restriction Enzyme Analysis"); st.markdown("Identify cut sites (A/T/C/G only)."); st.info("Uses Biopython REBASE.")
    in_m=st.radio("Input:",("Paste","FASTA"),key="re_in",horizontal=True); seq,seq_id="","Seq"
    if in_m=="Paste": seq=st.text_area("DNA:",h=150,key="re_seq",placeholder="Paste DNA...")
    else: up_f=st.file_uploader("FASTA:",type=['fasta','fa','fna'],key="re_fasta"); seq,seq_id=parse_fasta(up_f) if up_f else (None,None)
    if seq: st.text_area("Preview:",value=f">**{seq_id}**\n{seq[:100]}{'...' if len(seq)>100 else ''}",h=100,key="re_disp",disabled=True)
    st.divider(); st.subheader("Enzymes")
    common=sorted(['EcoRI','BamHI','HindIII','NotI','XhoI','PstI','SacI','SalI','SpeI','KpnI','SmaI','EcoRV','HaeIII','AluI','BglII','NcoI','NdeI','XbaI'])
    sel_com=st.multiselect("Common:",common,default=['EcoRI','BamHI','HindIII'],key="re_sel")
    cust_str=st.text_input("Custom (comma-sep):",key="re_cust",placeholder="e.g., BsaI")
    # *** SYNTAX FIX APPLIED HERE ***
    all_sel=set(sel_com)
    if cust_str: all_sel.update({e.strip() for e in cust_str.split(',') if e.strip()})
    # *** END FIX ***
    final=sorted(list(all_sel)); is_lin=st.checkbox("Linear",value=True,key="re_lin",help="Uncheck for circular.")
    st.divider()
    if st.button("Analyze",key="re_analyze_btn"):
        if not seq: st.warning("Provide DNA."); elif not final: st.warning("Select enzyme(s).")
        else:
            seq_cl="".join(re.findall(r'[ATCG]',seq.upper()))
            if not seq_cl: st.error("No valid DNA bases.")
            else:
                try:
                    seq_o=Seq(seq_cl); st.info(f"Analyzing '**{seq_id}**' ({len(seq_o):,} bp, {'Lin' if is_lin else 'Circ'}).")
                    valid,invalid=[],[]
                    try: from Bio.Restriction import AllEnzymes
                    except ImportError: st.error("Cannot import RE data."); st.stop()
                    with st.spinner("Validating..."):
                        for n in final:
                            try: enz=AllEnzymes.get(n); valid.append(enz) if enz else invalid.append(n)
                            except ValueError: invalid.append(n)
                    if invalid: st.warning(f"Ignored unknown: `{'`, `'.join(invalid)}`")
                    if not valid: st.error("No valid enzymes."); st.stop()
                    st.write(f"**Using:** {', '.join(map(str,valid))}")
                    with st.spinner("Searching..."): rb=RestrictionBatch(valid); analysis=Analysis(rb,seq_o,linear=is_lin)
                    st.subheader("Results"); st.write(f"**Seq:** '**{seq_id}**'|**Len:** {len(seq_o):,}|**Type:** {'Lin' if is_lin else 'Circ'}"); st.divider()
                    res=analysis.full(); summary,cuts,found=[],0,False
                    for enz,sites in res.items(): c=len(sites); cuts+=c; found=found or (c>0); summary.append({"Enzyme":str(enz),"Site":str(enz.site),"Cuts":c,"Pos(1-based)":", ".join(f"{s:,}" for s in sites) if sites else "None"})
                    if not found: st.success("✅ No cut sites.")
                    else:
                        st.write("#### Cuts:"); df=pd.DataFrame(summary).sort_values(["Cuts","Enzyme"],ascending=[False,True]); st.dataframe(df.style.format({"Cuts":"{:,}"}),hide_index=True,use_container_width=True)
                        st.divider(); st.write("#### Fragments:")
                        try: f_lens=sorted(analysis.fragments(),reverse=True); fc1,fc2=st.columns(2); fc1.metric("Cuts",f"{cuts:,}"); fc2.metric("Frags",f"{len(f_lens):,}"); st.write("**Lengths(bp):**");
                            if f_lens: st.text(", ".join(f"{l:,}" for l in f_lens[:15])+ (f", ...({len(f_lens)-15} more)" if len(f_lens)>30 else ""))
                            else: st.text("N/A.")
                        except Exception as fe: st.warning(f"Frag fail: {fe}")
                        st.divider();
                        if st.checkbox("Show Map",False,key="re_map"):
                            st.write("#### Map:")
                            try: mo=io.StringIO(); mw=min(100,max(60,len(seq_o)//10)); analysis.print_that(out=mo,title=f"Map:{seq_id}",top=True,rc=True,nc=mw); st.code(mo.getvalue()); mo.close()
                            except Exception as me: st.error(f"Map fail: {me}")
                except Exception as e: st.error(f"Analysis error: {e}")

elif menu == "Pairwise Sequence Alignment": # Syntax fix applied
    st.header("Pairwise Alignment"); st.markdown("Align DNA/protein (global/local)."); st.info("Scores affect results.")
    c1,c2=st.columns(2)
    with c1: st.markdown("#### Seq 1"); i_m1=st.radio("In1:",("Paste","FASTA"),key="al_in1",horizontal=True,label_visibility="collapsed"); s1,id1="","S1";
        if i_m1=="Paste": s1=st.text_area("S1:",h=120,key="al_s1",placeholder="Paste...")
        else: up1=st.file_uploader("F1:",type=['fasta','fa','fna','faa'],key="al_f1"); s1,id1=parse_fasta(up1) if up1 else (None,None)
        if s1: st.text_area("P1:",value=f">**{id1}**\n{s1[:80]}{'...' if len(s1)>80 else ''}",h=75,key="al_d1",disabled=True)
    with c2: st.markdown("#### Seq 2"); i_m2=st.radio("In2:",("Paste","FASTA"),key="al_in2",horizontal=True,label_visibility="collapsed"); s2,id2="","S2";
        if i_m2=="Paste": s2=st.text_area("S2:",h=120,key="al_s2",placeholder="Paste...")
        else: up2=st.file_uploader("F2:",type=['fasta','fa','fna','faa'],key="al_f2"); s2,id2=parse_fasta(up2) if up2 else (None,None)
        if s2: st.text_area("P2:",value=f">**{id2}**\n{s2[:80]}{'...' if len(s2)>80 else ''}",h=75,key="al_d2",disabled=True)
    st.divider(); st.subheader("Params"); pc1,pc2,pc3=st.columns(3)
    with pc1: mode=st.selectbox("Mode:",["Global","Local"],key="al_m"); seq_t=st.radio("Type:",("DNA","Protein"),key="al_t",horizontal=True)
    with pc2: st.write("**Score:**"); mat,m_s,mm_p=None,None,None;
        if seq_t=="DNA": m_s=st.number_input("Match:",value=2.0,step=0.5,key="al_mt"); mm_p=st.number_input("Mismatch:",value=-1.0,step=-0.5,key="al_ms"); mm_p=min(0.0,mm_p); st.caption(f"M={m_s}, MM={mm_p}")
        else: mats=sorted(substitution_matrices.list_matrices()); def_idx=mats.index('BLOSUM62') if 'BLOSUM62' in mats else 0; sel_mat=st.selectbox("Matrix:",mats,index=def_idx,key="al_mat");
            try: mat=substitution_matrices.load(sel_mat); st.caption(f"**{sel_mat}**")
            except Exception as e: st.error(f"Matrix load error: {e}")
    with pc3: st.write("**Gaps:**"); g_o=st.number_input("Open:",value=-10.0,step=-0.5,key="al_go"); g_e=st.number_input("Extend:",value=-0.5,step=-0.1,key="al_ge"); g_o=min(0.0,g_o); g_e=min(0.0,g_e); st.caption(f"O={g_o}, E={g_e}")
    st.divider()
    if st.button("Align",key="al_btn"):
        valid=True; if not s1: st.warning("S1 missing.");valid=False; if not s2: st.warning("S2 missing.");valid=False; if seq_t=="Protein" and not mat: st.error("Matrix needed.");valid=False; if not valid: st.stop()
        s1_c,s2_c="".join(s1.split()).upper(),"".join(s2.split()).upper(); if not s1_c: st.error("S1 empty.");valid=False; if not s2_c: st.error("S2 empty.");valid=False; if not valid: st.stop()
        try:
            st.subheader("Results"); prefix=mode.lower(); gaps={'open':g_o,'extend':g_e}; func=None; params={};
            if mat: fname=f"{prefix}dx"; params={'matrix':mat,**gaps}
            else: fname=f"{prefix}mx"; params={'match':m_s,'mismatch':mm_p,**gaps}
            if hasattr(pairwise2.align,fname): func=getattr(pairwise2.align,fname)
            else: st.error(f"Func '{fname}' not found."); st.stop()
            st.info(f"Aligning '**{id1}**' vs '**{id2}**' ({mode} {seq_t})...");
            with st.spinner("Aligning..."): aligns=func(s1_c,s2_c,**params,one_alignment_only=True)
            if not aligns: st.warning("No alignment.")
            else:
                a1,a2,scr,beg,end=aligns[0]; st.metric("Score",f"{scr:.2f}")
                if prefix=="local": st.write(f"**Region (0-based):** Start={beg}, End={end}")
                st.divider(); st.write("#### Alignment:")
                fmt_str=pairwise2.format_alignment(a1,a2,scr,beg,end,full_sequences=(prefix=='global'))
                disp_txt=f"S1: {id1}\nS2: {id2}\n\n{fmt_str}"; st.code(disp_txt)
                st.caption("Key: | Match, . Mismatch(+),   Mismatch(-), - Gap.")
                ident,g1,g2,pairs=0,0,0,0; align_len=len(a1)
                for i in range(align_len):
                    r1,r2=a1[i],a2[i]; gap1,gap2=(r1=='-'),(r2=='-')
                    # *** SYNTAX FIX APPLIED HERE ***
                    if gap1: g1+=1
                    if gap2: g2+=1
                    # *** END FIX ***
                    if not gap1 and not gap2: pairs+=1; if r1==r2: ident+=1
                tot_g=g1+g2; id_pct=(ident/pairs*100) if pairs>0 else 0; gap_pct=(tot_g/(align_len*2)*100) if align_len>0 else 0
                st.divider(); st.write("#### Stats:")
                sc1,sc2,sc3=st.columns(3); sc1.metric("Len",f"{align_len:,}"); sc2.metric("Ident",f"{id_pct:.1f}%",f"{ident:,}/{pairs:,}"); sc3.metric("Gaps",f"{gap_pct:.1f}%",f"{tot_g:,}/{align_len*2:,}")
                st.caption("Ident=Matches/Aligned Cols. Gaps=Total Gap Chars/(2*Len).")
        except Exception as e: st.error(f"Align error: {e}"); st.exception(e)

elif menu == "Motif Finder Tool": # Logic unchanged, styling updates via CSS
    st.header("Motif Finder"); st.markdown("Search motifs (exact/regex)."); st.info("Regex: `GAATTC`, `G[AT]TC`, `^ATG`, `TAG$`. Escape with `\\`.")
    in_m=st.radio("Input:",("Paste","FASTA"),key="mo_in",horizontal=True); seq,seq_id="","Seq"
    if in_m=="Paste": seq=st.text_area("Seq:",h=150,key="mo_seq",placeholder="Paste DNA/protein...")
    else: up_f=st.file_uploader("FASTA:",type=['fasta','fa','fna','faa'],key="mo_f"); seq,seq_id=parse_fasta(up_f) if up_f else (None,None)
    if seq: st.text_area("Preview:",value=f">**{seq_id}**\n{seq[:80]}{'...' if len(seq)>80 else ''}",h=75,key="mo_d",disabled=True)
    st.divider(); motif=st.text_input("Motif/Regex:",key="mo_pat",placeholder="e.g., GAATTC"); m1,m2=st.columns(2); with m1: case=st.checkbox("Case Sens",False,key="mo_case"); with m2: overlap=st.checkbox("Overlaps",True,key="mo_ov")
    st.divider()
    if st.button("Find",key="mo_find"):
        if not seq: st.warning("Provide seq."); elif not motif: st.warning("Enter motif.")
        else:
            seq_cl="".join(seq.split());
            if not seq_cl: st.warning("Seq empty.")
            else:
                try:
                    flags=0 if case else re.IGNORECASE; matches=[]; st.info(f"Searching '**{motif}**' in '**{seq_id}**'...")
                    with st.spinner("Searching..."): pat=re.compile(motif,flags);
                        if overlap: matches=list(pat.finditer(seq_cl))
                        else: pos=0; while pos<len(seq_cl): m=pat.search(seq_cl,pos); if m: matches.append(m); pos=m.end(); else: break
                    st.subheader(f"Results for '**{motif}**' in '**{seq_id}**'")
                    if matches:
                        st.success(f"Found **{len(matches):,}** match(es)."); m_data=[{"#":i+1,"Start":m.start()+1,"End":m.end(),"Len":m.end()-m.start(),"Match":m.group()} for i,m in enumerate(matches)]
                        m_df=pd.DataFrame(m_data); st.dataframe(m_df.style.format({"Start":"{:,}","End":"{:,}","Len":"{:,}"}),hide_index=True,use_container_width=True)
                        st.divider()
                        if st.checkbox("Highlight (2k chars)",False,key="mo_hi"):
                             lim=2000; s_hl=seq_cl[:lim]; html=""; last=0; style="background-color:#fff3cd;color:#856404;font-weight:bold;"; op=f"<mark style='{style}'>"; cl="</mark>"; sort_m=sorted(matches,key=lambda m:m.start())
                             for m in sort_m:
                                 if m.start()>=lim: break; s,e=m.start(),min(m.end(),lim)
                                 if s>=last: txt=s_hl[last:s].replace("&","&").replace("<","<").replace(">",">"); html+=txt; m_txt=s_hl[s:e].replace("&","&").replace("<","<").replace(">",">"); html+=op+m_txt+cl; last=e
                                 elif e>last: ov=s_hl[last:e].replace("&","&").replace("<","<").replace(">",">"); html+=op+ov+cl; last=e
                             rem=s_hl[last:].replace("&","&").replace("<","<").replace(">",">"); html+=rem
                             st.markdown(f"**Highlighted ({lim:,} chars):**"); st.markdown(f"<div style='font-family:monospace;word-wrap:break-word;border:1px solid #dee2e6;padding:10px;border-radius:5px;background-color:#f8f9fa;color:#333;'>{html}{'...' if len(seq_cl)>lim else ''}</div>",unsafe_allow_html=True)
                    else: st.info("Motif not found.")
                except re.error as r_e: st.error(f"Invalid Regex: {r_e}")
                except Exception as e: st.error(f"Motif search error: {e}"); st.exception(e)

elif menu == "Bioinformatics Tool (Transcription/Translation)": # Logic unchanged, styling updates via CSS
    st.header("DNA Tx & Tl"); st.markdown("Transcribe DNA->RNA, translate RNA->Protein.")
    in_m=st.radio("Input:",("Paste","FASTA"),key="tr_in",horizontal=True); dna_seq,seq_id="","Seq"
    if in_m=="Paste": dna_seq=st.text_area("DNA (Coding):",h=120,key="tr_seq",placeholder="Paste DNA...")
    else: up_f=st.file_uploader("FASTA:",type=['fasta','fa','fna'],key="tr_f"); dna_seq,seq_id=parse_fasta(up_f) if up_f else (None,None)
    if dna_seq: st.text_area("Preview:",value=f">**{seq_id}**\n{dna_seq[:80]}{'...' if len(dna_seq)>80 else ''}",h=75,key="tr_d",disabled=True)
    st.divider()
    if st.button("Run",key="tr_btn"):
        if not dna_seq: st.warning("Provide DNA.")
        else:
            dna_cl="".join(re.findall(r'[ATCG]',dna_seq.upper()));
            if not dna_cl: st.error("No valid DNA.")
            else:
                try:
                    st.subheader(f"Results: '**{seq_id}**'")
                    dna_o=Seq(dna_cl); dna_l=len(dna_o); c1,c2=st.columns(2)
                    with c1: st.metric("DNA Len",f"{dna_l:,} bp"); st.write("**In DNA:**"); st.code(str(dna_o))
                    with c2: gc=gc_fraction(dna_o)*100 if dna_l>0 else 0; st.metric("GC%",f"{gc:.1f}%"); st.write("**Rev Comp:**"); st.code(str(dna_o.reverse_complement()))
                    st.divider(); st.write("#### Tx"); rna_o=dna_o.transcribe(); st.write("**RNA:**"); st.code(str(rna_o)); st.metric("RNA Len",f"{len(rna_o):,}")
                    st.divider(); st.write("#### Tl"); st.caption("Std code. '*' STOP.")
                    rem=len(rna_o)%3; rna_t=rna_o; warn=False
                    if len(rna_o)<3: st.error("RNA < 3."); rna_t=None
                    elif rem!=0: warn=True; rna_t=rna_o[:-rem]; if not rna_t: st.error("Eff RNA len 0."); rna_t=None
                    if rna_t:
                        p_stop=rna_t.translate(1,to_stop=True,cds=False); st.write("**Prot (to STOP):**"); st.code(str(p_stop)); st.metric("Len",f"{len(p_stop):,} aa")
                        if warn: st.warning(f"RNA len not /3; last {rem} ignored.")
                        st.divider();
                        if st.checkbox("Full Tl",False,key="tr_full"): p_full=rna_t.translate(1,to_stop=False,cds=False); st.write("**Full Prot:**"); st.code(str(p_full)); st.metric("Full Len",f"{len(p_full):,} aa")
                except Exception as e: st.error(f"Error: {e}")

elif menu == "Genome Coverage Plotter": # Logic unchanged, styling updates via CSS/style_plot
    st.header("Coverage Plotter"); st.markdown("Visualize coverage depth."); st.info("Needs pos/cov cols.")
    up_f=st.file_uploader("Coverage Data:",type=['csv','tsv','txt','bed','bedgraph','depth'],key="cov_up")
    if up_f:
        st.divider(); st.subheader("Parsing")
        ext=os.path.splitext(up_f.name)[-1].lower(); def_s='\t' if ext!='.csv' else ','
        c1,c2=st.columns(2)
        with c1: s_opts={',':'Comma', '\t':'Tab', ' ':'Space'}; s_keys=list(s_opts.keys()); def_i=s_keys.index(def_s) if def_s in s_keys else 1; sel_k=st.selectbox("Sep:",s_keys,index=def_i,format_func=lambda x:s_opts[x],key="cov_sep"); s_re=r'\s+' if sel_k==' ' else sel_k
        with c2: def_h=ext in ['.csv','.tsv','.txt']; head=st.checkbox("Header",value=def_h,key="cov_head"); h_arg=0 if head else None; comm=st.text_input("Comment:",'#',max_chars=1,key="cov_comm")
        st.divider(); st.subheader("Loading & Cols")
        try:
            df_c=pd.read_csv(io.BytesIO(up_f.getvalue()),sep=s_re,header=h_arg,engine='python',comment=comm if comm else None,low_memory=False)
            if df_c.empty: st.error("File empty."); st.stop()
            if h_arg is None and all(isinstance(c,int) for c in df_c.columns): n_c=df_c.shape[1]; df_c.columns=[f'C{i+1}' for i in range(n_c)] if n_c>=2 else df_c.columns; st.info(f"Generic names: {list(df_c.columns)}")
            if df_c.shape[1]<2: st.error("Need >= 2 cols."); st.stop()
            st.write("Preview:"); st.dataframe(df_c.head(),height=200,use_container_width=True)
            cols=list(df_c.columns); norm={str(c).strip().lower():c for c in cols}; p_g,c_g,ch_g=None,None,None; ch_k=['chr','chrom','#chr']; p_k=['pos','position','start']; c_k=['depth','cov','coverage','score']
            for k in ch_k:   if k in norm: ch_g=norm[k]; break
            for k in p_k:   if k in norm: p_g=norm[k]; break
            for k in c_k:   if k in norm: c_g=norm[k]; break
            ch_i=cols.index(ch_g) if ch_g else 0; p_i=cols.index(p_g) if p_g else 1; c_i=cols.index(c_g) if c_g else (3 if len(cols)>3 and ext in ['.bed','.bedgraph'] else 2); ch_i,p_i,c_i=min(ch_i,len(cols)-1),min(p_i,len(cols)-1),min(c_i,len(cols)-1)
            cs1,cs2,cs3=st.columns(3);
            with cs1: chr_c=st.selectbox("Chrom (Opt):",[None]+cols,index=(cols.index(ch_g)+1 if ch_g else 0),key="cov_chr")
            with cs2: pos_c=st.selectbox("Pos:",cols,index=p_i,key="cov_pos")
            with cs3: cov_c=st.selectbox("Cov:",cols,index=c_i,key="cov_cov")
            if pos_c==cov_c or (chr_c and (chr_c==pos_c or chr_c==cov_c)): st.error("Cols must differ.")
            else:
                st.success(f"Using '{pos_c}'(Pos), '{cov_c}'(Cov)" + (f", '{chr_c}'(Chr)" if chr_c else ".")); st.divider(); st.subheader("Filter & Plot")
                try:
                    use=[pos_c,cov_c]; if chr_c: use.insert(0,chr_c); df_full=df_c[use].copy(); new=['Pos','Cov']; if chr_c: new.insert(0,'Chr'); df_full.columns=new; sel_chr=None
                    if chr_c: u_chrs=sorted(df_full['Chr'].astype(str).unique()); opts=['All']+u_chrs; sel_chr=st.selectbox("Filter Chr:",opts,0,key="cov_filter_chr"); df_plot=df_full[df_full['Chr'].astype(str)==str(sel_chr)].copy() if sel_chr!='All' else df_full.copy(); st.info(f"Filter: **{sel_chr}**.") if sel_chr!='All' else None
                    else: df_plot=df_full.copy()
                    df_plot['Pos']=pd.to_numeric(df_plot['Position'],errors='coerce'); df_plot['Cov']=pd.to_numeric(df_plot['Coverage'],errors='coerce'); r_b=len(df_plot); df_plot.dropna(subset=['Pos','Cov'],inplace=True); r_a=len(df_plot);
                    if r_a<r_b: st.warning(f"Removed **{r_b-r_a:,}** invalid rows.")
                    if df_plot.empty: st.error("No data left."); st.stop()
                    df_plot=df_plot.sort_values('Pos').reset_index(drop=True); mn,mx=df_plot['Pos'].min(),df_plot['Pos'].max(); sfx=f" ({sel_chr})" if sel_chr and sel_chr!='All' else f" (All)" if chr_c else ""; st.write(f"Plotting **{mn:,.0f}**-**{mx:,.0f}**{sfx}.")
                    cp1,cp2=st.columns(2); with cp1: fill=st.checkbox("Fill",True,key="cov_fill"); log_y=st.checkbox("Log Y",False,key="cov_log")
                    with cp2: mx_sm=min(201,max(1,len(df_plot)//2)); sm_w=st.slider("Smooth:",0,mx_sm,0,2,key="cov_smooth"); if sm_w>0 and sm_w%2==0: sm_w+=1; st.caption(f"Smooth:{sm_w}")
                    fig,ax=plt.subplots(figsize=(14,5)); ax.plot(df_plot['Pos'],df_plot['Cov'],label='Raw',lw=0.7,color='#007bff',alpha=0.8); if fill: ax.fill_between(df_plot['Pos'],df_plot['Cov'],alpha=0.2,color='#80ccff')
                    if sm_w>1: df_plot['Smooth']=df_plot['Cov'].rolling(sm_w,center=True,min_periods=1).mean(); ax.plot(df_plot['Pos'],df_plot['Smooth'],label=f'Smooth ({sm_w}bp)',lw=1.8,color='#ff7f0e',ls='-')
                    p_t=f"Coverage: {up_f.name}{sfx}"; p_x=f"Pos ({pos_c})"; p_y=f"Depth ({cov_c}){ ' (Log)' if log_y else ''}"; style_plot(fig,ax,title=p_t,xlabel=p_x,ylabel=p_y)
                    if log_y: ax.set_yscale('log'); mn_pc=df_plot[df_plot['Cov']>0]['Cov'].min(); ax.set_ylim(bottom=max(0.1,mn_pc*0.5) if pd.notna(mn_pc) and mn_pc>0 else 0.1)
                    else: ax.set_ylim(bottom=0)
                    ax.ticklabel_format(style='plain',axis='x',useOffset=False); ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x,p:format(int(x),','))); leg=ax.legend(facecolor='#fff',labelcolor='#333'); [t.set_fontfamily('Times New Roman')or t.set_fontweight('bold') for t in leg.get_texts()]; st.pyplot(fig)
                    st.divider(); st.subheader(f"Stats (Raw{sfx})"); stats=df_plot['Cov'].describe(); sc1,sc2,sc3,sc4=st.columns(4); sc1.metric("Mean",f"{stats.get('mean',0):,.1f}x"); sc2.metric("Median",f"{stats.get('50%',0):,.1f}x"); sc3.metric("Min",f"{stats.get('min',0):,.0f}x"); sc4.metric("Max",f"{stats.get('max',0):,.0f}x")
                    mx_c=int(stats.get('max',1000)); df_th=min(10,max(1,int(stats.get('mean',10)))); mx_th=max(1,mx_c); cov_th=st.number_input("Breadth ≥:",0,mx_th,df_th,1,key="cov_b_thresh")
                    if cov_th>=0: above=(df_plot['Cov']>=cov_th).sum(); total=len(df_plot); if total>0: brd=(above/total*100); st.metric(f"Breadth (≥{cov_th}x)",f"{brd:.1f}%",f"{above:,}/{total:,} pos"); else: st.metric(f"Breadth (≥{cov_th}x)","N/A")
                except Exception as plot_e: st.error(f"Plot error: {plot_e}")
        except Exception as read_e: st.error(f"Read error: {read_e}"); st.exception(read_e)

elif menu == "Variant Annotation Tool": # Logic unchanged, styling updates via CSS
    st.header("Simple Variant Substitution"); st.markdown("Introduce SNP into DNA & see codon effect.")
    in_m=st.radio("Input Ref:",("Paste","FASTA"),key="var_in",horizontal=True); ref_seq,seq_id="","Ref"
    if in_m=="Paste": ref_seq=st.text_area("Ref DNA:",h=100,key="var_ref",placeholder="Paste ref...")
    else: up_f=st.file_uploader("Ref FASTA:",type=['fasta','fa','fna'],key="var_f"); ref_seq,seq_id=parse_fasta(up_f) if up_f else (None,None)
    if ref_seq: st.text_area("Ref Preview:",value=f">**{seq_id}**\n{ref_seq[:80]}{'...' if len(ref_seq)>80 else ''}",h=75,key="var_d",disabled=True)
    st.divider(); st.subheader("Variant"); v1,v2=st.columns(2); with v1: var_pos=st.number_input("Pos (1-based):",1,step=1,value=1,key="var_p"); with v2: var_base=st.selectbox("New Base:",["A","T","C","G"],key="var_b")
    st.divider()
    if st.button("Apply & Annotate",key="var_btn"):
        if not ref_seq: st.warning("Provide ref.")
        else:
            ref_cl="".join(re.findall(r'[ATCG]',ref_seq.upper())); err=False
            if not ref_cl: st.error("Ref invalid."); err=True
            else: ref_len=len(ref_cl); z_pos=var_pos-1; if z_pos<0 or z_pos>=ref_len: st.error(f"Pos {var_pos:,} out of range (1-{ref_len:,})."); err=True
            if not err:
                try:
                    orig=ref_cl[z_pos]
                    if orig==var_base: st.info(f"No change: Base at **{var_pos:,}** ('**{seq_id}**') already '{orig}'."); st.write("**Ref:**"); st.code(ref_cl)
                    else:
                        alt_l=list(ref_cl); alt_l[z_pos]=var_base; alt_seq="".join(alt_l)
                        st.subheader("Variant Applied"); st.write(f"**Seq:** '**{seq_id}**'|**Var:** Pos **{var_pos:,}**"); st.markdown(f"**Change:** Ref (`{orig}`)→Alt (`{var_base}`)")
                        s1,s2=st.columns(2); with s1: st.write("**Ref:**"); st.code(ref_cl); with s2: st.write("**Alt:**"); st.code(alt_seq)
                        st.divider(); st.write("#### Codon & Effect:")
                        c_s=(z_pos//3)*3; p_in_c=(z_pos%3)+1
                        if c_s+3<=ref_len:
                            ref_c,alt_c=ref_cl[c_s:c_s+3],alt_seq[c_s:c_s+3]; st.write(f"Var at **{var_pos:,}** is pos **{p_in_c}** in codon at **{c_s+1:,}**.")
                            try:
                                ref_a,alt_a=str(Seq(ref_c).translate(1,cds=False)),str(Seq(alt_c).translate(1,cds=False));
                                cc1,cc2=st.columns(2); with cc1: st.write("**Ref Codon:**"); st.code(f"Codon:{ref_c}\nAA:   {ref_a}"); with cc2: st.write("**Alt Codon:**"); st.code(f"Codon:{alt_c}\nAA:   {alt_a}")
                                eff="?"; if ref_a==alt_a: eff="✅ Silent"; elif alt_a=='*': eff="🛑 Nonsense" if ref_a!='*' else "❓ Stop Retained"; elif ref_a=='*': eff="➡️ Stop-Lost"; else: eff="🔄 Missense"; st.markdown(f"**Effect:** {eff}")
                            except Exception as ce: st.warning(f"Translate fail: {ce}")
                        else: st.info("Variant near end.")
                except Exception as e: st.error(f"Error: {e}")

elif menu == "Codon Usage Analyzer": # Logic unchanged, styling updates via CSS/style_plot
    st.header("Codon Usage"); st.markdown("Analyze codon freq in CDS (A/T/C/G, len / 3).")
    in_m=st.radio("Input CDS:",("Paste","FASTA"),key="cod_in",horizontal=True); cds_seq,seq_id="","CDS"
    if in_m=="Paste": cds_seq=st.text_area("CDS:",h=120,key="cod_seq",placeholder="Paste CDS...")
    else: up_f=st.file_uploader("FASTA:",type=['fasta','fa','fna'],key="cod_f"); cds_seq,seq_id=parse_fasta(up_f) if up_f else (None,None)
    if cds_seq: st.text_area("Preview:",value=f">**{seq_id}**\n{cds_seq[:80]}{'...' if len(cds_seq)>80 else ''}",h=75,key="cod_d",disabled=True)
    st.divider()
    if st.button("Analyze",key="cod_btn"):
        if not cds_seq: st.warning("Provide CDS.")
        else:
            cds_cl="".join(re.findall(r'[ATCG]',cds_seq.upper())); err=False
            if not cds_cl: st.error("No valid DNA."); err=True; elif len(cds_cl)==0: st.error("Empty."); err=True; elif len(cds_cl)%3!=0: st.error(f"Length ({len(cds_cl):,}) not /3."); err=True
            if not err:
                try:
                    st.subheader(f"Usage: '**{seq_id}**'")
                    codons=[cds_cl[i:i+3] for i in range(0,len(cds_cl),3)]; tot_c=len(codons); counts=Counter(codons); st.metric("Codons",f"{tot_c:,}")
                    try: std_t=CodonTable.unambiguous_dna_by_id[1]; all_p=list(std_t.forward_table.keys())+std_t.stop_codons; aa_map=std_t.forward_table; stops=set(std_t.stop_codons)
                    except Exception as te: st.error(f"Table fail: {te}"); st.stop()
                    usage=[];
                    for cdn in sorted(all_p): cnt=counts.get(cdn,0); freq=(cnt/tot_c*100) if tot_c>0 else 0; aa=aa_map.get(cdn,'Stop' if cdn in stops else '?'); usage.append({"Codon":cdn,"AA":aa,"Count":cnt,"Freq(%)":freq})
                    df=pd.DataFrame(usage); st.write("#### Table:"); st.dataframe(df.style.format({'Count':'{:,}','Freq(%)':'{:.1f}%'}),hide_index=True,use_container_width=True)
                    st.divider(); st.subheader("Plot"); plot_c=st.selectbox("Plot:",["Freq (All)","Rel Usage (AA)"],key="cod_plot")
                    pres=df[df['Count']>0].copy()
                    if pres.empty: st.info("No codons.")
                    else:
                        if plot_c=="Freq (All)":
                            df_a=pres.sort_values("Codon"); fig,ax=plt.subplots(figsize=(16,6)); u_aas=sorted(list(df_a['AA'].unique())); cmap=plt.get_cmap('tab20',len(u_aas)); aa_c={aa:cmap(i) for i,aa in enumerate(u_aas)}; b_c=[aa_c.get(aa,'#888') for aa in df_a['AA']]; ax.bar(df_a['Codon'],df_a['Freq(%)'],color=b_c)
                            style_plot(fig,ax,title=f"Freq '{seq_id}'",xlabel="Codon",ylabel="Freq (%)"); ax.tick_params(axis='x',rotation=90,labelsize=9); [l.set_fontfamily('Times New Roman')or l.set_fontweight('bold') for l in ax.get_xticklabels()]; st.pyplot(fig)
                        elif plot_c=="Rel Usage (AA)":
                            avail=[aa for aa in pres['AA'].unique() if aa!='?'];
                            if not avail: st.info("No std AAs/Stop.")
                            else:
                                sel=st.selectbox("Select:",sorted(avail),key="cod_aa"); aa_df=pres[pres['AA']==sel].copy(); tot_aa=aa_df['Count'].sum()
                                if tot_aa>0 and not aa_df.empty:
                                    aa_df['Rel(%)']=(aa_df['Count']/tot_aa*100); aa_df=aa_df.sort_values("Codon"); fig,ax=plt.subplots(figsize=(max(6,len(aa_df)*1.5),5)); ax.bar(aa_df['Codon'],aa_df['Rel(%)'],color='#ff7f0e')
                                    style_plot(fig,ax,title=f"Rel Usage '{sel}' ({tot_aa:,} total)",xlabel="Codon",ylabel="Rel (%)"); ax.set_ylim(0,105); st.pyplot(fig)
                                    st.write(f"**Details '{sel}':**"); st.dataframe(aa_df[['Codon','Count','Rel(%)']].style.format({'Count':'{:,}','Rel(%)':'{:.1f}%'}),hide_index=True,use_container_width=True)
                                else: st.info(f"No codons for '{sel}'.")
                except Exception as e: st.error(f"Usage error: {e}"); st.exception(e)


# --- Footer ---
st.markdown("---")
st.markdown('<div class="footer">Copyright © 2024 CCDB Tools | For Educational and Informational Purposes Only</div>', unsafe_allow_html=True)
