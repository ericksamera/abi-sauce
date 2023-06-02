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
import _calc_ta
# --------------------------------------------------
class App:
    def __init__(self):
        """
        """
    def _init_page(self) -> None:
        """
        Function instantiates the main page.
        """
        self.title = "tm-calc"
        self.page_icon = "🔥"
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
            st.markdown('This script is intended for calculating melting temperatures and annealing temperatures for primers.')
            st.markdown('Check out the better-maintained command-line interface on [GitHub](https://github.com/KPU-AGC/general-resources/blob/main/miscellaneous/calc-ta.py)!')
            st.divider()

            with st.expander('MORE INFO'):
                # st.markdown(
                #     'In Sanger sequencing, the beginning and end of the ' 
                #     'electropherogram generally end up a litle messy due to the '
                #     'inherent randomness of the chain-termination method.')
            
                # st.markdown(
                #     'The low-quality basecalls at the beginning and end of the '
                #     'electropherogram are likely not real -- and therefore not '
                #     'useful for BLAST or variant identification.')

                # st.markdown(
                #     'These are usually trimmed off manually, but you can use this '
                #     'tool to do it too!')
                st.caption('')
                st.caption(f'[@{__author__}](https://github.com/ericksamera)\t|\tv{__version__}\t|\t{__comments__}')
    
    def _init_dilution_calc(self) -> None:
        """
        """
    
        st.warning('This page is experimental and may crash!')

        pol_arg = st.selectbox(
            'Choose a polymerase:',
            ('SuperFi', 'Phusion', 'DreamTaq')
            )
        
        col1, col2 = st.columns(2)

        data = [
            {"primer 1 (5'-3')": 'GATCGTACAGTCGATCG', "primer 2 (5'-3')": 'CTAGCATGCTAGCATG'},
        ]

        
        df = pd.DataFrame(data)
        with col1: primers_table: pd.DataFrame = st.experimental_data_editor(df, num_rows='dynamic', use_container_width=True)
        with col2:
            results = []
            for row_i, value in primers_table.iterrows():
                primer_1, primer_2 = list(value.to_numpy())
                primer_1_tm = None
                primer_2_tm = None

                result = {
                    'primer 1 tm (°C)': '',
                    'primer 2 tm (°C)': '',
                    'predicted Ta (°C)': 'Enter a valid sequence!',
                    'notes': ''
                    }

                if primer_1: 
                    if not [char for char in primer_1 if char.upper() not in 'ACGTUWSMKRYBDHVN']:
                        primer_1_obj  = _calc_ta.Primer('', primer_1.upper(), 0.5)
                        primer_1_tm = primer_1_obj.return_Tm(pol_arg)
                        result['primer 1 tm (°C)'] = primer_1_tm
                if primer_2: 
                    if not [char for char in primer_2 if char.upper() not in 'ACGTUWSMKRYBDHVN']:
                        primer_2_obj  = _calc_ta.Primer('', primer_2.upper(), 0.5)
                        primer_2_tm = primer_2_obj.return_Tm(pol_arg)
                        result['primer 2 tm (°C)'] = primer_2_tm

                if primer_1_tm and primer_2_tm: 
                    ta_result = _calc_ta.calculate_Ta(primer_1_obj, primer_2_obj, pol_arg)
                    result['predicted Ta (°C)'], result['notes'] = ta_result['Ta'], ta_result['note']
                results.append(result)
            results_df = pd.DataFrame(results)
            st.table(results)
        warnings = ''.join(list(results_df['notes'].to_numpy()))
        if not warnings: st.success('Primers look good to me!')
        else: st.warning('Warning! One of your primer sets has an issue. Check the notes.')

        st.text('Todo: will add option to download text file containing batch of primer info')            

# --------------------------------------------------
if __name__=="__main__":
    streamlit_app = App()
    streamlit_app._init_page()