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
import pandas as pd
# --------------------------------------------------
class App:
    def __init__(self):
        """
        """
    def _init_page(self) -> None:
        """
        Function instantiates the main page.
        """
        self.title = "dilution-calculator"
        self.page_icon = "💧"
        st.set_page_config(
            page_title=f"abi-sauce | {self.title}",
            page_icon=f'{self.page_icon}',
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
        self._init_dilution_calc()
    def _init_sidebar(self) -> None:
        """
        Instantiate the sidebar.
        """
    
        with st.sidebar:
            st.title(f"{self.page_icon} abi-sauce | {self.title}")
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
    
    def _init_dilution_calc(self) -> None:
        """
        """
        col1, col2 = st.columns(2)

        data = [
            {"Sample concentration (ng/uL)": 0,"sample volume (mL)": 2,"Target concentration (ng/uL)": 12,},
            {"Sample concentration (ng/uL)": 0,"sample volume (mL)": 2,"Target concentration (ng/uL)": 12,},    
        ]
        df = pd.DataFrame(data).reindex(index=range(1, len(data)+1))
        with col1: x = st.experimental_data_editor(df, num_rows='dynamic', use_container_width=True)
        with col2: st.dataframe(x, use_container_width=True)

# --------------------------------------------------
if __name__=="__main__":
    streamlit_app = App()
    streamlit_app._init_page()