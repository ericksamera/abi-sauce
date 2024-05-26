import streamlit as st
import _file_manager
import io
import zipfile

class App:
    def __init__(self) -> None:
        """
        """
    def _init_page(self) -> None:
        st.set_page_config(layout='wide')
        self._init_sidebar()
        self._init_main()
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
    def _init_main(self) -> None:
        """
        """
        processed_files_df = []
        for seq_object_name, seq_object in st.session_state.PROCESSED_FILES.items():
            dict_to_add: dict ={
                'Quality': seq_object.quality_color,
                'Filename': seq_object_name,
                'Well': seq_object.well,
                'PUP score': seq_object.quality_metrics['pup_score'],
                'Trace score': seq_object.quality_metrics['trace_score'],
                'CRL score': seq_object.quality_metrics['crl_score'],
            }
            processed_files_df.append(dict_to_add)
        st.dataframe(processed_files_df, hide_index=False, use_container_width=True, selection_mode="single-row")

        self._init_file_uploader()
# --------------------------------------------------
if __name__=="__main__":
    streamlit_app = App()
    streamlit_app._init_page()
