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
        if 'UPLOADED_FILES' not in st.session_state: self._init_file_uploader()
        else: self._update_plot_window()

    def _init_sidebar(self) -> None:
        """
        Instantiate the sidebar.
        """
    
        with st.sidebar:
            st.title(f"{self.emoji} abi-sauce | {self.title}")
            st.markdown('This script is intended for processing `.ab1` files into Mott algorithm-trimmed FASTAs.')
            st.markdown('Check out the better-maintained command-line interface on [GitHub](https://github.com/KPU-AGC/general-resources/blob/main/sanger-processing/sanger-sequence-trim.py)!')

            if 'UPLOADED_FILES' in st.session_state:
                st.divider()
                with st.expander('**TRACE FILES:**', expanded=True):

                    st.session_state.FLAGS = st.multiselect(
                        'Select property to display with name.',
                        [
                            'Qualitative Score',
                            'Trace Score',
                            'PUP Score',
                            'Contiguous Read Length (CRL)',
                            'Ambiguous Nucleotide Count'
                            ],
                        ['Qualitative Score'],
                    )
                    _translate_labels = {
                        'Sequence ID': 'name',
                        'Well Position': 'well',
                        'Trace Score': 'trace_score',
                        'PUP Score': 'pup_score',
                        'Contiguous Read Length (CRL)': 'crl_score',
                        'Ambiguous Nucleotide Count' : 'ambig_count',
                    }

                    st.session_state.REVERSE_SORT = st.checkbox(
                        'Reverse sort'
                    )

                    st.session_state.TRACE_LIST_SORT = st.selectbox(
                        'Sort by:',
                        list(_translate_labels.keys()))

                    def show_flags_with_name(_seq_name, _flags_list) -> str:

                        _flags_dict = {
                            'Qualitative Score': st.session_state.PROCESSED_FILES[f"{_seq_name}.ab1"]['color_code'],
                            'Trace Score': str(st.session_state.PROCESSED_FILES[f"{_seq_name}.ab1"]['trace_score']).zfill(2),
                            'PUP Score': str(st.session_state.PROCESSED_FILES[f"{_seq_name}.ab1"]['pup_score']).zfill(2),
                            'Contiguous Read Length (CRL)': str(st.session_state.PROCESSED_FILES[f"{_seq_name}.ab1"]['crl_score']).zfill(2),
                            'Ambiguous Nucleotide Count': str(st.session_state.PROCESSED_FILES[f"{_seq_name}.ab1"]['ambig_count']).zfill(2),
                        }

                        flags_and_name = '|'.join([str(_flags_dict[i]) for i in st.session_state.FLAGS] + ['']) + st.session_state.PROCESSED_FILES[f'{_seq_name}.ab1']['name']

                        return flags_and_name                
                    
                    st.session_state.SORTED_LIST = sorted(st.session_state.PROCESSED_FILES.values(), reverse=st.session_state.REVERSE_SORT, key=lambda x: x[_translate_labels[st.session_state.TRACE_LIST_SORT]])
                    st.session_state.SELECTED_TRACE = st.radio(
                        'Select trace file to view:', 
                        options=[seq_object['name'] for seq_object in st.session_state.SORTED_LIST],
                        format_func=lambda x: show_flags_with_name(x, st.session_state.FLAGS))
                with st.expander('**DOWNLOAD OPTIONS:**'):
                    st.session_state.DEFAULT_FILENAME = 'abi-sauce-trim'
                    st.session_state.USER_FILENAME = st.text_input('File name', placeholder='abi-sauce-trim')
                    st.session_state.FILENAME =  st.session_state.USER_FILENAME if st.session_state.USER_FILENAME else st.session_state.DEFAULT_FILENAME
                    st.session_state.CONCATENATE: bool = st.checkbox("Concatenate entries into single fasta.", value=True)
                    if not st.session_state.CONCATENATE: st.caption('Individual FASTAs will be compiled into a single ZIP file.')
                    file_type, file_buffer = self._prepare_download()
                    download_button = ste.download_button(
                        label="Download Sequences",
                        data=file_buffer,
                        file_name=f"{st.session_state.FILENAME}.{file_type}",
                    )
            st.divider()

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

            with st.expander('Advanced options'):
                st.session_state.TRIM: bool = st.checkbox("Trim sequences using Mott's trimming algorithm.", value=True)
            submit_button = st.form_submit_button("Process files", on_click=self._upload_files, args=(uploaded_files,))
            if submit_button and not uploaded_files:
                st.error('Select some files!')
    def _upload_files(self, _st_uploaded_files: list) -> None:
        if not _st_uploaded_files: return None
        st.session_state.UPLOADED_FILES = _st_uploaded_files
        st.session_state.TRIM_STR = '_trimmed' if st.session_state.TRIM else '_raw'
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
    def _prepare_download(self):
        """
        """
        if not st.session_state.CONCATENATE:
            with io.BytesIO() as file_buffer:
                with zipfile.ZipFile(file_buffer, "w") as zip:
                    for seq_object in st.session_state.PROCESSED_FILES.values():
                        fasta_string = seq_object[st.session_state.TRIM_STR].format("fasta")
                        zip.writestr(f"{seq_object['name']}.fasta", fasta_string)
                file_buffer.seek(0)
                return ('zip', file_buffer.getbuffer().tobytes())
        else:
            fasta_string: str = ''
            for seq_object in st.session_state.PROCESSED_FILES.values():
                fasta_string += seq_object[st.session_state.TRIM_STR].format("fasta")
            return ('txt', str.encode(fasta_string))
    def _plot_electropherogram(self, _seq_object_dict):
        """
        """
        raw_annotations = _seq_object_dict['_raw'].annotations['abif_raw']
        #trim_annotations = seq_object_dict['_trimmed'].annotations['abif_raw']
        left_trim = raw_annotations['PLOC2'][_seq_object_dict['left_trim']]
        try: right_trim = raw_annotations['PLOC2'][-_seq_object_dict['right_trim']-2] + raw_annotations['SPAC3']/2
        except IndexError: right_trim = raw_annotations['PLOC2'][-_seq_object_dict['right_trim']] + raw_annotations['SPAC3']/2
        
        phred_scores = _seq_object_dict['_raw'].letter_annotations['phred_quality']
        nucleotide_plots = {
            'A': {
                'complement': 'T',
                'peaks': raw_annotations['DATA10'],
                'color': 'green'},
            'T': {
                'complement': 'A',
                'peaks': raw_annotations['DATA11'],
                'color': 'red'},
            'C': {
                'complement': 'G',
                'peaks': raw_annotations['DATA12'],
                'color': 'blue'},
            'G': {
                'complement': 'C',
                'peaks': raw_annotations['DATA9'],
                'color': 'black'}
            }

        if st.session_state.TRIM_STR == '_trimmed':
            max_peak_height = max([max(nucleotide_plots[nucleotide]['peaks'][int(left_trim):int(right_trim)]) for nucleotide in nucleotide_plots])
        else:
            max_peak_height = max([max(nucleotide_plots[nucleotide]['peaks']) for nucleotide in nucleotide_plots])
        relative_heights = {
            'screen_height': max_peak_height+200,
            'basepos_height': max_peak_height+130,
            'basecall_height': max_peak_height+75,
            'highlight_height': max_peak_height+50
        }

        fig = go.Figure()
        fig.add_trace(
            go.Bar(
            #x=[math.ceil(positions-(widths/2)) for widths, positions in zip(adjusted_dists,raw_annotations['PLOC2'])],
            x=raw_annotations['PLOC2'],
            y=[(5/6)*(relative_heights['highlight_height'])*(score/60) for score in phred_scores],
            hoverinfo='skip',
            name='Phred scores',
            width=abs(raw_annotations['SPAC3']),
            marker=dict(
                color='#DFF0FA',
                opacity=1,
                ),
            ))
        fig.add_trace(
            go.Scattergl(
                x=raw_annotations['PLOC2'],
                y=[relative_heights['basecall_height'] for i in raw_annotations['PLOC2']],
                hoverinfo='skip',
                name='Nucleotides',
                text=list(char for char in raw_annotations['PBAS2'].decode()),
                mode="text",
                textfont={
                    'color': [nucleotide_plots[char]['color'] if char in nucleotide_plots else '#ff3aff' for char in raw_annotations['PBAS2'].decode()]
                },
            ))
    
        if str(_seq_object_dict['_raw'].seq) != 'NNNNN':
            fig.add_vline(
                x=left_trim - raw_annotations['SPAC3']/2,
                line_width=2,
                fillcolor='#88ccee',
                opacity=0.5)
            fig.add_vline(
                x=right_trim + raw_annotations['SPAC3']/2,
                line_width=2,
                fillcolor='#88ccee',
                opacity=0.5)

        for nuc, values in nucleotide_plots.items():
            fig.add_trace(
                go.Scattergl(
                    y=values['peaks'],
                    hoverinfo='skip',
                    line=dict(width=1),
                    name=nuc,
                    marker=dict(
                        size=20,
                        color=values['color'])))

        fig.update_layout(
            margin=dict(l=0, r=0, t=50, b=0),
            dragmode='pan',
            xaxis=dict(rangeslider=dict(visible=True, thickness=0.1), tickvals=[None], range=[left_trim-(10*raw_annotations['SPAC3']), right_trim+(10*raw_annotations['SPAC3'])], constrain='domain'),
            yaxis=dict(fixedrange=True, tickvals=[None], range=[0, relative_heights['screen_height']]))

        return fig
    def _update_plot_window(self):
        """
        """
        st.header(f"{st.session_state.SELECTED_TRACE}")
        st.text(f"{st.session_state.PROCESSED_FILES[f'{st.session_state.SELECTED_TRACE}.ab1']['_raw'].id}")
        st.plotly_chart(
            self._plot_electropherogram(st.session_state.PROCESSED_FILES[f"{st.session_state.SELECTED_TRACE}.ab1"]),
            use_container_width=True)
        help_str = "Trimmed with Mott's trimming algorithm!" if st.session_state.TRIM_STR == '_trimmed' else "Untrimmed!"
        st.caption(f'FASTA sequence: ({len(st.session_state.PROCESSED_FILES[f"{st.session_state.SELECTED_TRACE}.ab1"][st.session_state.TRIM_STR])} bp)', help=help_str)
        st.code(st.session_state.PROCESSED_FILES[f"{st.session_state.SELECTED_TRACE}.ab1"][st.session_state.TRIM_STR].format("fasta"))
        if st.session_state.TRIM_STR == '_trimmed':
            with st.expander(f'Untrimmed FASTA ({len(st.session_state.PROCESSED_FILES[f"{st.session_state.SELECTED_TRACE}.ab1"]["_raw"])} bp):'):
                st.code(st.session_state.PROCESSED_FILES[f"{st.session_state.SELECTED_TRACE}.ab1"]["_raw"].format("fasta"))
            st.subheader('Trimming statistics')
            left_column, right_column = st.columns(2)
            with left_column:
                l_trim_count: int = st.session_state.PROCESSED_FILES[f"{st.session_state.SELECTED_TRACE}.ab1"]["left_trim"]
                r_trim_count: int = st.session_state.PROCESSED_FILES[f"{st.session_state.SELECTED_TRACE}.ab1"]["right_trim"]+1
                total_trim_count = l_trim_count+r_trim_count

                phred_scores_raw = st.session_state.PROCESSED_FILES[f"{st.session_state.SELECTED_TRACE}.ab1"]['_raw'].letter_annotations['phred_quality']
                phred_scores_trimmed = st.session_state.PROCESSED_FILES[f"{st.session_state.SELECTED_TRACE}.ab1"]['_trimmed'].letter_annotations['phred_quality']

                average_phred_raw = sum(phred_scores_raw)/len(phred_scores_raw)
                average_phred_trimmed = sum(phred_scores_trimmed)/len(phred_scores_trimmed)

                st.markdown(f'Left trim: `{l_trim_count}` | Right trim:`{r_trim_count}` (Total: `{total_trim_count}`)', help='The right side is usually worse.')
                st.markdown(f'Average PHRED score (before): `{average_phred_raw:.1f}` | (after) `{average_phred_trimmed:.1f}`')

            with right_column:
                """"""
                #st.markdown(f'Right trim,: {st.session_state.PROCESSED_FILES[f"{st.session_state.SELECTED_TRACE}.ab1"]["right_trim"]+1} nucleotides')
# --------------------------------------------------
if __name__=="__main__":
    streamlit_app = App()
    streamlit_app._init_page()
