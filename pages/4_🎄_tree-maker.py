#!/usr/bin/env python3
__description__ =\
"""
Purpose: Streamlit wrapper for fla-viewer-advanced.
"""
__author__ = "Erick Samera"
__version__ = "1.1.0"
__comments__ = "stable enough"
# --------------------------------------------------
import streamlit as st
import streamlit_ext as ste
# --------------------------------------------------
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo
import io
import re
from itertools import combinations_with_replacement
# --------------------------------------------------
class App:
    def __init__(self):
        """
        """
        self.title = "tree-maker"
        self.emoji = '🎄'
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
    def _init_page(self) -> None:
        """
        Function instantiates the main page.
        """
        self._init_sidebar()
        self._init_file()
        if 'TREE' in st.session_state:
            st.caption("Copy this into the 'Tree text' box in  <a href=https://itol.embl.de/upload.cgi>Interactive Tree of Life (ITOL)</a>", unsafe_allow_html=True)
            st.code(st.session_state.TREE)
            st.caption("This is for Erick to debug.")
            st.code(st.session_state.PHYLIP)
        return None
    def _add_tree_text(self, newick_tree_text: str) -> None: st.session_state.TREE = newick_tree_text.replace('-', '')
    def _add_phylip_out(self, out):
        parsed_out = '\n'.join(out.split('\n')[:-1])
        st.session_state.PHYLIP = re.sub(' {1,}', '\t', parsed_out)
    def _parse_input_table(self, input_data) -> None:
        keys = ['Name', '339-1', '339-2', '339-3', '339-4', '339-5', '407-1', '407-2', '407-3', '407-4', '407-5', '121-1', '121-2', '121-3', '121-4', '121-5','73-1', '73-2', '73-3', '73-4', '73-5','285-1', '285-2', '285-3', '285-4', '285-5', '223-1', '223-2', '223-3', '223-4', '223-5', '441-1', '441-2', '441-3', '441-4', '441-5', '157-1', '157-2', '157-3', '157-4', '157-5', '95-1', '95-2', '95-3', '95-4', '95-5', '311-1', '311-2', '311-3', '311-4', '311-5']
        if not input_data: return None
        input_list = []
        samples_list = []
        for sample_line in input_data.split('\n')[1:]:
            dict_to_add = {key: value for key, value in zip(keys, sample_line.split('\t'))}
            dict_to_add_processed = {key: value if key == "Name" else float(value) for key, value in dict_to_add.items()}
            input_list.append(dict_to_add_processed)
            if dict_to_add['Name'] != '': samples_list.append(dict_to_add['Name'])
        
        if samples_list:
            self._update_genetic_distance(input_list, samples_list)


    def _update_genetic_distance(self, input_list, samples_list) -> None:
        """
        """
        
        def find_common_difference(numbers):
            if len(numbers) < 2:
                return None

            # Sort the list
            sorted_numbers = sorted(numbers)

            # Calculate the initial difference
            common_diff = sorted_numbers[1] - sorted_numbers[0]

            # Verify that all differences are the same
            for i in range(2, len(sorted_numbers)):
                if sorted_numbers[i] - sorted_numbers[i - 1] != common_diff:
                    return "No common difference"

            return common_diff
        def myround(x, base=5):return base * round(x/base)

        def _get_bruvo_distance(input_dict, marker):
            
            def _generate_matrix_per_allele(alleles_1: list, alleles_2: list):
                full_allele_matrix = np.zeros(array_size, np.float64)
                for allele in range(len(alleles_1)):
                    allele_matrix = np.zeros(array_size, np.float64)
                    sample_1_allele = alleles_1[allele]
                    sample_2_allele = alleles_2[allele]
                    repeat_units = int(myround(abs((sample_1_allele - sample_2_allele)), marker_repeat_adj[marker]) / marker_repeat_adj[marker])
                    distance = float(1 - (2 ** (-repeat_units)))
                    allele_matrix[i][ii] = distance
                    full_allele_matrix += allele_matrix
                return full_allele_matrix
            
            genome_addition = True
            genome_loss = True
            
            marker_matrix = np.zeros(array_size, np.float64)
            for i, sample_1 in enumerate(samples_list):
                for ii, sample_2 in enumerate(samples_list):
                    sample_1_alleles = [allele for allele in input_dict[marker][sample_1] if allele > 1]
                    sample_2_alleles = [allele for allele in input_dict[marker][sample_2] if allele > 1]
                    
                    if len(sample_1_alleles) == len(sample_2_alleles): 
                        smaller_genotype = sample_1_alleles
                        larger_genotype = sample_2_alleles
                        marker_matrix += _generate_matrix_per_allele(smaller_genotype, larger_genotype)
                    else:
                        smaller_genotype = sample_1_alleles if len(sample_1_alleles) < len(sample_2_alleles) else sample_2_alleles
                        larger_genotype = sample_2_alleles if len(sample_1_alleles) < len(sample_2_alleles) else sample_1_alleles

                        genome_adjust_matrices = []
                        if genome_addition:
                            genome_addition_combinations = [smaller_genotype + list(i) for i in set(sorted(combinations_with_replacement(smaller_genotype, len(larger_genotype) - len(smaller_genotype))))]
                            for smaller_genotype_combination in genome_addition_combinations:
                                genome_adjust_matrices.append(_generate_matrix_per_allele(smaller_genotype_combination, larger_genotype))
                        if genome_loss:
                            genome_loss_combinations = [smaller_genotype + list(i) for i in set(sorted(combinations_with_replacement(larger_genotype, len(larger_genotype) - len(smaller_genotype))))]
                            for smaller_genotype_combination in genome_loss_combinations:
                                genome_adjust_matrices.append(_generate_matrix_per_allele(smaller_genotype_combination, larger_genotype))
                        if (not genome_addition) and (not genome_loss):
                            smaller_genotype = smaller_genotype + [1 * (len(larger_genotype) - len(smaller_genotype))]
                            genome_adjust_matrices.append(_generate_matrix_per_allele(smaller_genotype, larger_genotype))
                    
                        genome_adjust_matrices_stack = np.stack(genome_adjust_matrices)
                        marker_matrix += np.mean(genome_adjust_matrices_stack, axis=0)

            return marker_matrix

        parsed_data: dict = {}
        marker_allele_calc = {}
        for sample_dict in input_list:
            marker_list = sorted(set([key.split('-')[0] for key in sample_dict.keys() if key != "Name"]))
            for marker in marker_list:
                if marker not in parsed_data: 
                    parsed_data[marker] = {}
                    marker_allele_calc[marker] = []

                sample_name = sample_dict['Name']
                marker_keys = [key for key in sample_dict.keys() if key.startswith(marker)]
                alleles = sorted([float(sample_dict[key]) if float(sample_dict.get(key, 0)) > 0 else 1.0 for key in marker_keys], reverse=True)
                marker_allele_calc[marker].append(find_common_difference([allele for allele in alleles if allele > 1.0]))
                parsed_data[marker].update({sample_name: alleles})
        
        # marker repeat calculator
        marker_repeat_adj = {}
        for marker_a in marker_allele_calc: marker_repeat_adj[marker_a] = sorted([m for m in marker_allele_calc[marker_a] if isinstance(m, float)])[0]

        samples_list = sorted(samples_list)

        if samples_list:
            array_size = (len(samples_list), len(samples_list))
            starting_matrix = np.zeros(array_size, np.float64)
            for marker in sorted(parsed_data.keys()):
                starting_matrix += _get_bruvo_distance(parsed_data, marker)

            distance_matrix = [row[:i+1] for i, row in enumerate(np.tril(starting_matrix).tolist())]

            x = DistanceMatrix([sample.replace(' ', '_') for sample in samples_list], distance_matrix)
            output = io.StringIO("")
            if x:
                Phylo.write(DistanceTreeConstructor().nj(x), output, 'newick')
                newick_tree = output.getvalue()
                self._add_phylip_out(str(len(samples_list)) +"\n" + str(x))
                self._add_tree_text(newick_tree)

    def _init_file(self) -> None:
        """
        """
    
        with st.form("upload form", clear_on_submit=True):
            raw_input_text = st.text_area('Raw Input Table')
            submit_button = st.form_submit_button("Upload and process", on_click=self._parse_input_table, args=(raw_input_text,), help="If you entered your data and clicked submit but it didn't output anything, click it again.")

    def _init_sidebar(self) -> None:
        """
        Instantiate the sidebar.
        """
    
        with st.sidebar:
            st.title(f"{self.emoji} abi-sauce | {self.title}")
            st.markdown('This script should be used to process the table of Hydrangea alleles for the lab.')

    def _reset_state(self) -> None:
        """
        """
        for key in st.session_state.keys():
            del st.session_state[key]
# --------------------------------------------------
if __name__=="__main__":
    streamlit_app = App()
    streamlit_app._init_page()