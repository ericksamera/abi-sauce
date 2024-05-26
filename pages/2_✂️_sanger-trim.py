import streamlit as st
import _file_manager
import _trace_visualizer
import io
import zipfile

class App:
    def __init__(self) -> None:
        """
        """
    def _init_page(self) -> None:
        st.set_page_config(layout='wide')
        self._init_sidebar()
        if 'UPLOADED_FILES' not in st.session_state: self._init_file_uploader()
        else: self._plot_electropherogram()
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
            if submit_button and not uploaded_files: st.error('Select some files!')
    def _init_sidebar(self) -> None:
        """
        Instantiate the sidebar.
        """
    
        with st.sidebar:
            #st.title(f"{self.emoji} abi-sauce | {self.title}")
            st.markdown('This script is intended for processing `.ab1` files into Mott algorithm-trimmed FASTAs.')
            st.markdown('Check out the better-maintained command-line interface on [GitHub](https://github.com/KPU-AGC/general-resources/blob/main/sanger-processing/sanger-sequence-trim.py)!')

            if 'UPLOADED_FILES' in st.session_state:
                st.divider()
                with st.expander('**TRACE FILES:**', expanded=True):
                    st.session_state.SELECTED_TRACE = st.selectbox(
                        'Select trace file to view:',
                        options=[seq_object_name for seq_object_name, seq_object in st.session_state.PROCESSED_FILES.items()],
                        format_func=lambda x:'* ' + x if st.session_state.PROCESSED_FILES[x].modified is True else x)
                    st.session_state.SELECTED_TRACE = st.session_state.PROCESSED_FILES[st.session_state.SELECTED_TRACE]
    
                with st.expander('**DOWNLOAD OPTIONS:**'):
                    st.session_state.DEFAULT_FILENAME = 'abi-sauce-trim'
                    st.session_state.USER_FILENAME = st.text_input('File name', placeholder='abi-sauce-trim')
                    st.session_state.FILENAME =  st.session_state.USER_FILENAME if st.session_state.USER_FILENAME else st.session_state.DEFAULT_FILENAME
                    st.session_state.CONCATENATE: bool = st.checkbox("Concatenate entries into single fasta.", value=True)
                    if not st.session_state.CONCATENATE: st.caption('Individual FASTAs will be compiled into a single ZIP file.')
                    file_type, file_buffer = self._prepare_download()
                    download_button = st.download_button(
                        label="Download Sequences",
                        data=file_buffer,
                        file_name=f"{st.session_state.FILENAME}.{file_type}")
                    st.divider()
    def _upload_files(self, _st_uploaded_files: list) -> None:
        if not _st_uploaded_files: return None
        st.session_state.UPLOADED_FILES = _st_uploaded_files
        st.session_state.PROCESSED_FILES = self._process_files(_st_uploaded_files)
        self._set_defaults()
    def _process_files(self, _st_uploaded_files: list) -> None:
        """
        """
        processed_files = {}
        for file in sorted(_st_uploaded_files, key=lambda x: x.name):
            processed_files[file.name] = _file_manager.ABIFile(file)
        return processed_files
    def _set_defaults(self) -> None:
        """
        """
        if 'SELECTED_TRACE' not in st.session_state:
            st.session_state.SELECTED_TRACE = st.session_state.PROCESSED_FILES[list(st.session_state.PROCESSED_FILES.keys())[0]]
            st.session_state.TRIMMING_ACTIVE = False
            st.session_state.TRIMMING_WHICH = None
    def _update_global(self, key: str, action: str, value) -> None:
        """
        """
        if action == 'toggle': setattr(st.session_state, key, not getattr(st.session_state, key))
        elif action == 'set': setattr(st.session_state.PROCESSED_FILES[_selected_trace], key, value)
    def _update_selected_property(self, _selected_trace: str, key: str, action: str, value) -> None:
        """
        """
        if action == 'toggle': setattr(st.session_state.PROCESSED_FILES[_selected_trace], key, not getattr(st.session_state.PROCESSED_FILES[_selected_trace], key))
        elif action == 'set': setattr(st.session_state.PROCESSED_FILES[_selected_trace], key, value)
    def _prepare_download(self):
        """
        """
        if not st.session_state.CONCATENATE:
            with io.BytesIO() as file_buffer:
                with zipfile.ZipFile(file_buffer, "w") as zip:
                    for seq_object in st.session_state.PROCESSED_FILES.values():
                        fasta_str = seq_object.fasta_str
                        zip.writestr(f"{seq_object.stem}.fasta", fasta_str)
                file_buffer.seek(0)
                return ('zip', file_buffer.getbuffer().tobytes())
        else:
            fasta_str: str = ''
            for seq_object in st.session_state.PROCESSED_FILES.values():
                fasta_str += seq_object.fasta_str
            return ('txt', str.encode(fasta_str))
    def _plot_electropherogram(self) -> None:
        """
        """
        
        def _activate_trim(which: str) -> None:
            """
            """

            if st.session_state.TRIMMING_ACTIVE:
                if st.session_state.TRIMMING_WHICH == which:
                    setattr(st.session_state, 'TRIMMING_ACTIVE', False)
                    setattr(st.session_state, 'TRIMMING_WHICH', None)
                else:
                    setattr(st.session_state, 'TRIMMING_WHICH', which)
            else:
                setattr(st.session_state, 'TRIMMING_ACTIVE', True)
                setattr(st.session_state, 'TRIMMING_WHICH', which)

        if 'SELECTED_TRACE' in st.session_state:

            self.plotter = _trace_visualizer.TraceVisualizer(st.session_state.SELECTED_TRACE, st.session_state.TRIMMING_ACTIVE)
            events = self.plotter.plot_electropherogram()
            st.toggle(
                "Reverse Complement (EXPERIMENTAL)",
                value=st.session_state.PROCESSED_FILES[st.session_state.SELECTED_TRACE.filename].reverse_complement,
                on_change=self._update_selected_property,
                args=(st.session_state.SELECTED_TRACE.filename, "reverse_complement", 'toggle', None))
            
            trim_l_col, trim_r_col = st.columns(2)
            with trim_l_col: st.button('⬅️ Edit left trim', type="primary" if st.session_state.TRIMMING_WHICH=='left' else "secondary", on_click=_activate_trim, args=("left", ), use_container_width=True)
            with trim_r_col: st.button('Edit right trim ➡️', type="primary" if st.session_state.TRIMMING_WHICH=='right' else "secondary", on_click=_activate_trim, args=("right", ), use_container_width=True)

            if events and st.session_state.TRIMMING_ACTIVE:
                trim_index = events[0]['pointNumber'] if not st.session_state.TRIMMING_WHICH == 'right' else len(st.session_state.SELECTED_TRACE.seq) - events[0]['pointNumber']
                self._update_selected_property(st.session_state.SELECTED_TRACE.filename, f'trim_{st.session_state.TRIMMING_WHICH}', 'set', trim_index)
                st.session_state.PROCESSED_FILES[st.session_state.SELECTED_TRACE.filename].modified = True

            #st.caption(f'FASTA sequence: ({len(st.session_state.PROCESSED_FILES[f"{st.session_state.SELECTED_TRACE}.ab1"][st.session_state.TRIM_STR])} bp)', help=help_str)
            window_width: int = 60
            fasta_seq: str = st.session_state.PROCESSED_FILES[st.session_state.SELECTED_TRACE.filename].seq[st.session_state.PROCESSED_FILES[st.session_state.SELECTED_TRACE.filename].trim_left:-st.session_state.PROCESSED_FILES[st.session_state.SELECTED_TRACE.filename].trim_right]
            fasta_header: str = f"{st.session_state.SELECTED_TRACE.name} trimmed: {st.session_state.SELECTED_TRACE.trim_left}L/{st.session_state.SELECTED_TRACE.trim_right}R"
            fasta_output_str: str = f'>{fasta_header}\n' + '\n'.join([fasta_seq[i:i + window_width] for i in range(0, len(fasta_seq), 60)])
            fasta_output_twoline: str = f'>{fasta_header}\n' + fasta_seq
            self._update_selected_property(st.session_state.SELECTED_TRACE.filename, 'fasta_header', 'set', fasta_header)

            
            print_two_line = st.toggle("Two-line FASTA")
            st.code(fasta_output_str if not print_two_line else fasta_output_twoline)
# --------------------------------------------------
if __name__=="__main__":
    streamlit_app = App()
    streamlit_app._init_page()
