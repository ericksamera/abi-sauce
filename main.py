import streamlit as st
import _file_manager
import _trace_visualizer

class App:
    def __init__(self) -> None:
        """
        """
    def _init_page(self) -> None:
        st.set_page_config(layout='wide')
        self._init_sidebar()
        self._init_file_uploader()
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
                        options=[seq_object_name for seq_object_name, seq_object in st.session_state.PROCESSED_FILES.items()])
                    st.session_state.SELECTED_TRACE = st.session_state.PROCESSED_FILES[st.session_state.SELECTED_TRACE]
    def _upload_files(self, _st_uploaded_files: list) -> None:
        if not _st_uploaded_files: return None
        st.session_state.UPLOADED_FILES = _st_uploaded_files
        st.session_state.setdefault('PROCESSED_FILES', {})
        st.session_state.PROCESSED_FILES.update(self._process_files(_st_uploaded_files))
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
# --------------------------------------------------
if __name__=="__main__":
    streamlit_app = App()
    streamlit_app._init_page()
