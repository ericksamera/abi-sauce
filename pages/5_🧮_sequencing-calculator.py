#!/usr/bin/env python3
__description__ =\
"""
Purpose: Streamlit calculator for illumina coverage.
"""
__author__ = "Erick Samera"
__version__ = "0.1.0"
__date__ = "20240507"
__comments__ = "minimally functional"
# --------------------------------------------------
import streamlit as st
import streamlit_ext as ste
# --------------------------------------------------
class App:
    def __init__(self):
        """
        """
        self.title = "sequencing-calculator"
        self.emoji = '🧮'
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

    def _parse_region_size(self, input_str: str) -> int:
        """
        """
        def _remove_numbers(input_str: str) -> int:
            input_str = '.'.join(input_str.split('.')[:2])
            return float(''.join([i for i in input_str if i in '1234567890.']))

        if input_str.lower().endswith('bp'): input_str = input_str.lower().split('bp')[0]
        elif input_str.lower().endswith('b'): input_str = input_str.lower().split('b')[0]
        else: input_str = input_str

        if input_str.upper()[-1] in 'MBG':
            conversion_factor = {
                'G': 1_000_000_000,
                'M': 1_000_000,
                'K': 1_000
            }
            return int(_remove_numbers(input_str[:-1]) * conversion_factor[input_str.upper()[-1]])
        else: return int(_remove_numbers(input_str))

    def _init_page(self) -> None:
        """
        Function instantiates the main page.
        """
        self._init_sidebar()

        #x_samples_per, x_depth = st.tabs(["Samples per flow cell", "Dog",])
        output_per_unit_dict = {
            "MiSeq v3 (2x300): 15000000000": 15000000000,
            "MiSeq v2 (2x250): 15000000000": 15000000000,
            "MiSeq v2 Nano (2x250): 600000000": 600000000,
            "Miseq v2 Micro (2x150): 2400000000": 2400000000,
            "MinION: 15000000000": 15000000000
        }

        variable_of_interest_list = [
            "samples per flow cell", 
            "depth", 
            "genome size"
        ]
        variable_of_interest = st.radio("Variable of interest:", variable_of_interest_list)

        region_size_col1, region_size_col2 = st.columns(2)
        with region_size_col1:
            region_size = st.text_input("Region/Genome size (bp)", value="3.3 GB", key="region_size_samples", disabled=True if variable_of_interest=="genome size" else False)
            if not region_size: region_size = "1 bp"
            region_size_int = self._parse_region_size(region_size)
        with region_size_col2:
            st.caption(f"Interpreted as:")
            st.caption(f"{region_size_int} bp / {region_size_int  / 1_000} Kbp / {region_size_int  / 1_000_000} Mbp / {region_size_int / 1_000_000_000} Gbp")

        samples_per_unit = st.number_input("Samples per flow cell", value=1, disabled=True if variable_of_interest=="samples per flow cell" else False)
        depth = st.number_input("Sequencing Depth (X)", value=20, disabled=True if variable_of_interest=="depth" else False)
        dup_col1, dup_col2 = st.columns(2)
        with dup_col1:
            duplication = st.number_input("Percent Read Duplication (%)", value=3, min_value=0, max_value=100)
        with dup_col2:
            st.caption("Duplications occur when multiple copies of the same original molecule are sequenced. The frequency of duplication is protocol dependent.")
            st.caption("The default (2%) is reasonable, according to experimental data.")
        on_target = st.number_input("Percent Reads On-target (%)", value=85, min_value=0, max_value=100)

        output_per_unit = st.selectbox(
            "Sequencing Output (bp)",
            options=(
            "MiSeq v3 (2x300): 15000000000",
            "MiSeq v2 (2x250): 15000000000",
            "MiSeq v2 Nano (2x250): 600000000",
            "Miseq v2 Micro (2x150): 2400000000",
            "MinION: 15000000000"),
            index=0,
            placeholder="Select sequencing chemistry",)
        
        if variable_of_interest != "depth":
            total_output_required = region_size_int * depth / ((1-(duplication/100)) * on_target/100)

        if variable_of_interest == "samples per flow cell":
                output = output_per_unit_dict[output_per_unit] / total_output_required
                output_unit = "samples per flow cell."
        elif variable_of_interest == "samples per flow cell":
                output = ((output_per_unit_dict[output_per_unit]/samples_per_unit) * ((1-(duplication/100)) * on_target/100) / region_size_int)
                output_unit = f"X sequencing depth across {samples_per_unit} samples."
        elif variable_of_interest == "samples per flow cell":
                reverse_conversion = st.radio("Unit", ["bp", "Kbp", "Mbp", "Gbp"])
                reverse_conversion_dict = {
                    'bp': 1,
                    'Kbp': 1_000,
                    'Mbp': 1_000_000,
                    'Gbp': 1_000_000_000
                }
                output = ((output_per_unit_dict[output_per_unit]/samples_per_unit) * ((1-(duplication/100)) * (on_target/100) / depth)) / reverse_conversion_dict[reverse_conversion]
                output_unit = f"{reverse_conversion} region/genome assuming {depth} depth across {samples_per_unit} samples."
        else:
                output = 0
                output_unit = ""

        st.header(f"{output:.1f} {output_unit}")

    def _init_sidebar(self) -> None:
        """
        Instantiate the sidebar.
        """
    
        with st.sidebar:
            st.title(f"{self.emoji} abi-sauce | {self.title}")
            st.markdown('For calculating the coverage required.')
            st.caption(f"@{__author__} | {__version__} ({__date__}): {__comments__}")

    def _reset_state(self) -> None:
        """
        """
        for key in st.session_state.keys():
            del st.session_state[key]
# --------------------------------------------------
if __name__=="__main__":
    streamlit_app = App()
    streamlit_app._init_page()