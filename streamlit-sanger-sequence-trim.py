#!/usr/bin/env python3
__description__ =\
"""
Purpose: Streamlit wrapper for sanger-sequence-trim.
"""
__author__ = "Erick Samera"
__version__ = "1.0.0"
__comments__ = "stable enough"
# --------------------------------------------------
import streamlit as st
import pandas as pd
# --------------------------------------------------
from Bio import SeqIO
from datetime import datetime
import zipfile
import io
import copy
import time
# --------------------------------------------------
class App:
    def __init__(self):
        """
        """
        self._loaded_abi_files = {}
    def _init_page(self) -> None:
        """
        Function instantiates the main page.
        """
        st.set_page_config(page_title='sanger-sequence-trim', initial_sidebar_state='collapsed')
        st.title('sanger-sequence-trim')
        st.markdown('This script is intended for processing a `.ab1` files into Mott algorithm-trimmed FASTAs.')
        st.markdown('Check out the better-maintained command-line interface on [GitHub](https://github.com/KPU-AGC/general-resources/blob/main/sanger-processing/sanger-sequence-trim.py)!')
        self._init_file_uploader()
        self._init_sidebar()
    def _init_file_uploader(self) -> None:
        """
        Function instantiates file uploading and processing for the files.
        """

        form_instance = st.empty()

        # lots of context-managing -- disgusting.
        with io.BytesIO() as buffer:
            with zipfile.ZipFile(buffer, "w") as zip:

                # Electropherogram upload form
                with form_instance.form("sanger upload form", clear_on_submit=True):

                    ## ========================================================
                    _uploaded_files: list = st.file_uploader(
                        'Upload `.ab1` files which come right off the SeqStudio.',
                        type=['ab1'],
                        accept_multiple_files=True)
                    concatenate_checkbox: bool = st.checkbox("Concatenate entries into single fasta.")
                    with st.expander('Advanced options'):

                        self.trim_flag: bool = st.checkbox("Trim sequences using Mott's trimming algorithm.", value=True)
                        self.trim_option = '_trimmed' if self.trim_flag else '_raw'

                    _submitted = st.form_submit_button("Process files")
                    ## ========================================================

                    if _submitted:
                        filename = datetime.now().strftime("%m%d%Y_%H%M%S")
                        if not concatenate_checkbox:
                            for file in _uploaded_files:
                                processed_file = self._process_seq_object(file)
                                zip.writestr(f"{file.name.split('.')[0]}.fasta", processed_file[self.trim_option].format("fasta"))
                        else:
                            fasta_string: str = ''
                            for file in _uploaded_files:
                                processed_file = self._process_seq_object(file)
                                fasta_string += processed_file[self.trim_option].format("fasta")
                            zip.writestr(f"sanger-sequence-trim_{filename}.fasta", fasta_string)
            buffer.seek(0)

            if _submitted and _uploaded_files:
                with st.spinner('Trimming...'):
                    time.sleep(0.5)
                form_instance.empty()
                self._init_dataframe(self._loaded_abi_files)
                button = st.download_button(
                    label="Download ZIP",
                    data=buffer,
                    file_name=f"sanger-sequence-trim_{filename}.zip" 
                )
            if _submitted and not _uploaded_files:
                st.error('Attach some files first!')
    def _process_seq_object(self, _file) -> None:
        """
        Function creates a new instance of the copied file and processes it.
        """
        _trimmed_file_instance = copy.deepcopy(_file)
        seq_object_raw = SeqIO.read(_file, 'abi')
        seq_object_trimmed = SeqIO.read(_trimmed_file_instance, 'abi-trim')

        seq_object_dict = {
            'name': seq_object_raw.name,
            '_raw': seq_object_raw,
            '_trimmed': seq_object_trimmed,
            'trace_score': seq_object_raw.annotations['abif_raw']['TrSc1'] if 'TrSc1' in seq_object_raw.annotations['abif_raw'] else -1,
            'pup_score': seq_object_raw.annotations['abif_raw']['PuSc1'] if 'PuSc1' in seq_object_raw.annotations['abif_raw'] else -1,
            'left_trim': seq_object_raw.seq.find(seq_object_trimmed.seq[0:5])-1,
            'right_trim': len(seq_object_raw.seq) - len(seq_object_trimmed) - seq_object_raw.seq.find(seq_object_trimmed.seq[0:5])-1
            }
        self._loaded_abi_files[_file.name] = seq_object_dict

        return seq_object_dict
    def _init_dataframe(self, _abi_data) -> None:
        """
        Function instantiates the DataFrame for a given set of input sequences.
        """

        dataframe_header = [key for key in _abi_data[list(_abi_data.keys())[0]].keys() if not key.startswith('_')]
        if not self.trim_flag: dataframe_header = [key for key in dataframe_header if 'trim' not in key]

        dataframe_values = [item for item in _abi_data.values()]
        dataframe = pd.DataFrame(dataframe_values, columns=dataframe_header)
        st.dataframe(dataframe)
    def _init_sidebar(self) -> None:
        """
        Instantiate the sidebar.
        """
        
        with st.sidebar:
            st.header('About')
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
# --------------------------------------------------
if __name__=="__main__":
    ct = App()
    ct._init_page()
