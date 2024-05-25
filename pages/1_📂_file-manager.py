#!/usr/bin/env python3
__description__ =\
"""
Purpose: Streamlit wrapper for sanger-sequence-trim.
"""
__author__ = "Erick Samera"
__version__ = "1.5.0"
__comments__ = "stable enough"
# --------------------------------------------------
import streamlit as st
import streamlit_ext as ste
# --------------------------------------------------
import plotly.graph_objects as go
from Bio import SeqIO
import zipfile
import io
import copy
import time
# --------------------------------------------------
class App:
    def __init__(self):
        """
        """
    def _init_page(self) -> None:
        """
        Function instantiates the main page.
        """
        self.title = "sanger-sequence-trim"
        self.emoji = ':scissors:'
        st.set_page_config(
            page_title=f"abi-sauce | {self.title}",
            page_icon=self.emoji,
            layout='wide',
            initial_sidebar_state='expanded')
        st.markdown(
            """
            <style>
            [data-testid="stSidebar"][aria-expanded="true"]{
                min-width: 450px;
                max-width: 450px;
            }""",
            unsafe_allow_html=True,
            )   
        self._init_sidebar()
        self._init_file_uploader()

    def _init_sidebar(self) -> None:
        """
        Instantiate the sidebar.
        """
    
        with st.sidebar:
            st.title(f"{self.emoji} abi-sauce | {self.title}")
            st.markdown('This script is intended for processing `.ab1` files into Mott algorithm-trimmed FASTAs.')
            st.markdown('Check out the better-maintained command-line interface on [GitHub](https://github.com/KPU-AGC/general-resources/blob/main/sanger-processing/sanger-sequence-trim.py)!')

            with st.expander('MORE INFO'):
                st.markdown(
                    'In Sanger sequencing, the beginning and end of the ' 
                    'electropherogram generally end up a litle messy due to the '
                    'inherent randomness of the chain-termination method.')
            
                st.markdown(
                    'The low-quality basecalls at the beginning and end of the '
                    'electropherogram are likely not real -- and therefore not '
                    'useful for BLAST or variant identification.')

                st.markdown(
                    'These are usually trimmed off manually, but you can use this '
                    'tool to do it too!')
                st.caption('')
                st.caption(f'[@{__author__}](https://github.com/ericksamera)\t|\tv{__version__}\t|\t{__comments__}')
    def _init_file_uploader(self) -> None:
        """
        """
        st.header('Select trace files (`.ab1`)')
        uploaded_files: list = st.file_uploader(
                        'Upload `.ab1` files which come right off the SeqStudio.',
                        type=['ab1'],
                        accept_multiple_files=True)
    
        
        with st.form("sanger upload form", clear_on_submit=False):
            submit_button = st.form_submit_button("Process files", on_click=self._upload_files, args=(uploaded_files,))
            if submit_button and not uploaded_files:
                st.error('Select some files!')
    def _upload_files(self, _st_uploaded_files: list) -> None:
        if not _st_uploaded_files: return None
        st.session_state.UPLOADED_FILES = _st_uploaded_files
        st.session_state.PROCESSED_FILES = self._process_files(_st_uploaded_files)
    def _process_files(self, _st_uploaded_files) -> None:
        """
        Function creates a new instance of the copied file and processes it.
        """
        processed_files = {}

        trace_scoring = {
            0: "🟩",
            1: "🟨",
            2: "🟥",
        }

        for file in sorted(_st_uploaded_files, key=lambda x: x.name):
            _trimmed_file_instance = copy.deepcopy(file)
            seq_object_raw = SeqIO.read(file, 'abi')
            seq_object_trimmed = SeqIO.read(_trimmed_file_instance, 'abi-trim')

            seq_object_dict = {
                'name': seq_object_raw.name,
                'well': seq_object_raw.annotations['abif_raw']['TUBE1'].decode(),
                '_raw': seq_object_raw,
                '_trimmed': seq_object_trimmed,
                'pup_score': seq_object_raw.annotations['abif_raw']['PuSc1'] if 'PuSc1' in seq_object_raw.annotations['abif_raw'] else -1,
                'trace_score': seq_object_raw.annotations['abif_raw']['TrSc1'] if 'TrSc1' in seq_object_raw.annotations['abif_raw'] else -1,
                'crl_score': seq_object_raw.annotations['abif_raw']['CRLn1'] if 'CRLn1' in seq_object_raw.annotations['abif_raw'] else -1,
                'left_trim': seq_object_raw.seq.find(seq_object_trimmed.seq[0:5]),
                'right_trim': len(seq_object_raw.seq) - len(seq_object_trimmed) - seq_object_raw.seq.find(seq_object_trimmed.seq[0:5])-1,
                'ambig_count': sum([seq_object_trimmed.seq.count(nuc) for nuc in 'NRYKMSWBDHV']) if st.session_state.TRIM else sum([seq_object_raw.seq.count(nuc) for nuc in 'NRYKMSWBDHV'])
                }
            seq_object_dict['scoring_str'] = f" ({seq_object_dict['pup_score']}/{seq_object_dict['trace_score']}/{seq_object_dict['crl_score']}) "
            if all([
                seq_object_dict['trace_score'] > 25,
                seq_object_dict['pup_score'] > 20,
                seq_object_dict['crl_score'] > 100,
                ]):
                seq_object_dict['color_code'] = trace_scoring[0]
            elif all([
                seq_object_dict['trace_score'] > 0,
                seq_object_dict['pup_score'] > 0,
                seq_object_dict['crl_score'] > 0,
                ]):
                seq_object_dict['color_code'] = trace_scoring[1]
            else:
                seq_object_dict['color_code'] = trace_scoring[2]

            processed_files[file.name] = seq_object_dict
        return processed_files
# --------------------------------------------------
if __name__=="__main__":
    streamlit_app = App()
    streamlit_app._init_page()
