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

        if input_str.upper()[-1] in 'GMK':
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
            "MiSeq v3 (2x300)": 15000000000,
            "MiSeq v2 (2x250)": 15000000000,
            "MiSeq v2 Nano (2x250)": 600000000,
            "Miseq v2 Micro (2x150)": 2400000000,
            "MinION, FLO-MIN114": 15000000000
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
            duplication = st.number_input("Percent Read Duplication (%)", value=float(2.5), min_value=float(0), max_value=float(99.9), step=float(1.0), format="%f")
        with dup_col2:
            st.caption("Duplications occur when multiple copies of the same original molecule are sequenced. The frequency of duplication is protocol dependent.")
            st.caption("The default (2.5%) is reasonable for Illumina sequencing, according to experimental data.")
        on_target = st.number_input("Percent Reads On-target (%)", value=85, min_value=0, max_value=100)

        seq_out_col1, seq_out_col2 = st.columns(2)
        with seq_out_col1:
            output_per_unit = st.selectbox(
                "Sequencing Output (bp)",
                options=(
                    "MiSeq v3 (2x300)",
                    "MiSeq v2 (2x250)",
                    "MiSeq v2 Nano (2x250)",
                    "Miseq v2 Micro (2x150)",
                    "MinION, FLO-MIN114"),
                index=0,
                placeholder="Select sequencing chemistry")
        
        unit_output = output_per_unit_dict[output_per_unit]
        with seq_out_col2:
            if output_per_unit in ("MinION, FLO-MIN114",):
                
                exp_runtime = st.slider(
                    "Runtime (hrs)",
                    value=72,
                    min_value=0,
                    max_value=72,
                    step=1)
                
                model_option = st.selectbox(
                    "Model",
                    ("Worst so far", "Probable", "Optimistic") if exp_runtime != 72 else ("Worst so far", "Worst WGS", "Probable", "Optimistic", "Theoretical"),
                    help="Base output is calculated by a third-order polynomial model from experimental data.",)

                model_option_to_eq = {
                    'Optimistic': lambda x: (-114458843.02911) + (7555412.38270*x) + (-279.48619*(x**2)) + (-0.09880*(x**3)),
                    'Probable': lambda x: 0.8 * ((-114458843.02911) + (7555412.38270*x) + (-279.48619*(x**2)) + (-0.09880*(x**3))),
                    'Worst WGS': lambda x: 0.85 * ((-198480166.79748) + (6999111.75559*x) + (-1270.41713*(x**2)) + (0.08046*(x**3))),
                    'Worst so far': lambda x: (-240161649.78617) + (6337746.57404*x) + (-2252.15270*(x**2)) + (0.26408*(x**3)),
                    'Theoretical': lambda x: 3472222.22222 * x,
                }

                unit_output = model_option_to_eq[model_option](exp_runtime * 60)
                st.caption(f"Sequencing output: {unit_output / 1_000_000_000:.1f} Gbp with a runtime of {exp_runtime} hrs.")
            elif output_per_unit in ("MiSeq v3 (2x300)", "MiSeq v2 (2x250)"):
                model_option = st.selectbox(
                    "Model",
                    ("Average", "Best so far", "Theoretical"),
                    help="Base output is calculated by a third-order polynomial model from experimental data.",)
                model_option_to_eq = {
                    'Average': 12.14e9,
                    'Best so far': 22.1e9,
                    'Theoretical': 15000000000,
                }
                unit_output = model_option_to_eq[model_option]

        if variable_of_interest != "depth":
            total_output_required = region_size_int * depth / ((1-(duplication/100)) * on_target/100)

        if variable_of_interest == variable_of_interest_list[0]: 
                output = unit_output / total_output_required
                output_unit = "samples per flow cell."
        elif variable_of_interest == variable_of_interest_list[1]:
                output = ((unit_output/samples_per_unit) * ((1-(duplication/100)) * on_target/100) / region_size_int)
                output_unit = f"X sequencing depth across {samples_per_unit} samples."
        elif variable_of_interest == variable_of_interest_list[2]:
                region_size_col1, region_size_col2 = st.columns(2)
                with region_size_col2:
                    reverse_conversion = st.radio("Unit", ["bp", "Kbp", "Mbp", "Gbp"])
                    reverse_conversion_dict = {
                        'bp': 1,
                        'Kbp': 1_000,
                        'Mbp': 1_000_000,
                        'Gbp': 1_000_000_000
                    }
                with region_size_col1:
                    output = ((unit_output/samples_per_unit) * ((1-(duplication/100)) * (on_target/100) / depth)) / reverse_conversion_dict[reverse_conversion]
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