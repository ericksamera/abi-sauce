#!/usr/bin/env python3
__description__ =\
"""
Purpose: Streamlit wrapper for sanger-sequence-trim.
"""
__author__ = "Erick Samera"
__version__ = "1.1.0"
__comments__ = "stable enough"
# --------------------------------------------------
import streamlit as st
import streamlit_ext as ste
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
        st.set_page_config(
            page_title=f"abi-sauce | {self.title}",
            page_icon=':apple:',
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
    def _init_sidebar(self) -> None:
        """
        Instantiate the sidebar.
        """
    
        with st.sidebar:
            st.title(f":apple: abi-sauce | {self.title}")
            st.markdown('This script is intended for processing a `.ab1` files into Mott algorithm-trimmed FASTAs.')
            st.markdown('Check out the better-maintained command-line interface on [GitHub](https://github.com/KPU-AGC/general-resources/blob/main/sanger-processing/sanger-sequence-trim.py)!')
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
# --------------------------------------------------
if __name__=="__main__":
    streamlit_app = App()
    streamlit_app._init_page()