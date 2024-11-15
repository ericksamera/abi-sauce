#!/usr/bin/env python3
__description__ =\
"""
Purpose: Streamlit wrapper for sanger-sequence-trim.
"""
__author__ = "Erick Samera"
__version__ = "1.1.1"
__comments__ = "stable enough"
# --------------------------------------------------
import streamlit as st
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
        self._init_main_window()
    def _init_sidebar(self) -> None:
        """
        Instantiate the sidebar.
        """
    
        with st.sidebar:
            st.title(f":apple: abi-sauce")
            st.success('Select one of the apps above to get started!')
    def _init_main_window(self) -> None:
        """
        """
        st.title('Welcome to :apple: `abi-sauce`!')
        st.markdown('These are a set of open-source tools to make handling electropherogram data easier.')
        st.markdown('Check out my GitHub with the link below for some of my other projects.')
        st.caption(f'[@{__author__}](https://github.com/ericksamera)\t|\tv{__version__}\t|\t{__comments__}')
        return None
# --------------------------------------------------
if __name__=="__main__":
    streamlit_app = App()
    streamlit_app._init_page()