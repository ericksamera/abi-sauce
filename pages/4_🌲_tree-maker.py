#!/usr/bin/env python3
__description__ =\
"""
Purpose: Streamlit wrapper for fla-viewer-advanced.
"""
__author__ = "Erick Samera"
__version__ = "1.3.0"
__comments__ = "stable enough"
# --------------------------------------------------
import streamlit as st
# --------------------------------------------------
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo
import io
import re
from itertools import combinations_with_replacement
from collections import Counter
# --------------------------------------------------
class App:
    def __init__(self):
        """
        """
        self.title = "TREE-maker"
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
            st.link_button("Open tree in ITOL", f"https://itol.embl.de/upload.cgi?ttext={st.session_state.TREE}", help=None, type="secondary", disabled=False, use_container_width=False)
            with st.expander("PHYLIP distance matrix output"):
                st.caption("This is for Erick to debug.")
                st.code(st.session_state.PHYLIP)
        return None
    def _add_tree_text(self, newick_tree_text: str) -> None: st.session_state.TREE = newick_tree_text.replace('-', '')
    def _add_phylip_out(self, out) -> None: st.session_state.PHYLIP = re.sub(' {1,}', '\t', '\n'.join(out.split('\n')[:-1]))
    def _parse_input_table(self, input_data: str) -> None:
        """
        Function parses the input table pasted into the text box.

        Parameters:
            input_data (str): input data in tabular format.
        
        Returns:
            (None)
        """
        hardcoded_alleles: list = ['Name', '339-1', '339-2', '339-3', '339-4', '339-5', '407-1', '407-2', '407-3', '407-4', '407-5', '121-1', '121-2', '121-3', '121-4', '121-5','73-1', '73-2', '73-3', '73-4', '73-5','285-1', '285-2', '285-3', '285-4', '285-5', '223-1', '223-2', '223-3', '223-4', '223-5', '441-1', '441-2', '441-3', '441-4', '441-5', '157-1', '157-2', '157-3', '157-4', '157-5', '95-1', '95-2', '95-3', '95-4', '95-5', '311-1', '311-2', '311-3', '311-4', '311-5']
        if not input_data: return None
        input_list = []
        samples_list = []
        for sample_line in input_data.split('\n')[1:]:
            dict_to_add = {key: value for key, value in zip(hardcoded_alleles, sample_line.split('\t'))}
            dict_to_add_processed = {key: value if key == "Name" else float(value) for key, value in dict_to_add.items()}
            input_list.append(dict_to_add_processed)
            if dict_to_add['Name'] != '': samples_list.append(dict_to_add['Name'])
        
        if samples_list: self._update_genetic_distance(input_list, samples_list)

    def _update_genetic_distance(self, input_list: list[dict], samples_list: list[str]) -> None:
        """
        Function generates a genetic distance matrix with a given input list of sample dictionaries and sample list.

        Parameters:
            input_list (list[dict]):
                list of sample dictionaries containing alleles per marker
            samples_list (list[str]):
                list of sample names, could probably be parsed from the input_list tbh
        
        Returns:
            (None)
        """
        
        def find_common_difference(numbers: list[float]):
            """
            Function finds the common difference in a list of numbers, i.e., [136, 139, 142] = 3.

            Parameters:
                numbers (list[float]):
                    list of input numbers
            
            Returns:
                (None) if invalid input else (Int) common difference
            """
            if len(numbers) < 2: return None
            sorted_numbers = sorted(numbers)
            common_diff = sorted_numbers[1] - sorted_numbers[0]

            for i in range(2, len(sorted_numbers)):
                if sorted_numbers[i] - sorted_numbers[i - 1] != common_diff: return None
            return common_diff
        
        def myround(x, base=5):return base * round(x/base)

        def _infer_marker_nuc_repeat_num(input_alleles_per_marker: dict) -> dict:
            """
            Function takes an input dictionary of alleles present in a given marker and infers the number of nucleotides in the repeated motif.

            Parameters:
                input_alleles_per_marker (dict[str]):
                    dictionary containing a list of alleles per marker.
            
            Returns:
                dict[str] = int
                    dictionary containing inferred number of nucleotides in the repeated motif for each given marker
            """
            inferred_marker_repeat_dict: dict = {}
            for marker in input_alleles_per_marker:
                marker_alleles = [m for m in input_alleles_per_marker[marker] if isinstance(m, float)]
                inferred_marker_repeat_dict[marker] = sorted(Counter(marker_alleles).items(), key=lambda x: x[1], reverse=True)[0][0]
            return inferred_marker_repeat_dict

        def _parse_data() -> tuple[dict, dict]:
            """
            Function parses the input
            """
            parsed_data: dict = {}
            nuc_repeat_n_raw_per_marker_dict: dict = {}
            for sample_dict in input_list:
                marker_list: list = sorted(set([key.split('-')[0] for key in sample_dict.keys() if key != "Name"]))
                for marker in marker_list:
                    if marker not in parsed_data: 
                        parsed_data[marker]: dict = {}
                        nuc_repeat_n_raw_per_marker_dict[marker]: list = []

                    sample_name = sample_dict['Name']
                    marker_keys = [key for key in sample_dict.keys() if key.startswith(marker)]
                    alleles = sorted([float(sample_dict[key]) if float(sample_dict.get(key, 0)) > 0 else 1.0 for key in marker_keys], reverse=True)
                    nuc_repeat_n_raw_per_marker_dict[marker].append(find_common_difference([allele for allele in alleles if allele > 1.0]))
                    parsed_data[marker].update({sample_name: alleles})
            return parsed_data, nuc_repeat_n_raw_per_marker_dict

        def _get_marker_bruvo_dist(input_dict: dict[str], marker: str) -> np.array:
            """
            Function generates bruvo distance for a given marker.

            Parameters:
                input_dict (dict):
                    dictionary of markers, samples per marker, alleles per sample
                input_marker (str):
                    marker
            
            Returns:
                (np.array) matrix distance matrix.
            """

            def __generate_matrix_per_allele(alleles_1: list, alleles_2: list) -> np.array:
                """
                Function generates a distance matrix per allele between 2 individuals (for a given marker).

                Parameters:
                    alleles_1 (list):
                        list of alleles for individual 1
                    alleles_2 (list):
                        list of alleles for individual 2
                
                Returns:
                    (np.array) containing distance matrix of all alleles stacked.
                """
                full_allele_matrix = np.zeros(array_size, np.float64)
                for allele in range(len(smaller_genotype)):
                    allele_matrix = np.zeros(array_size, np.float64)
                    sample_1_allele = alleles_1[allele]
                    sample_2_allele = alleles_2[allele]
                    repeat_units = int(myround(abs((sample_1_allele - sample_2_allele)), marker_nuc_repeat_dict[marker]) / marker_nuc_repeat_dict[marker])
                    distance = float(1 - (2 ** (-repeat_units)))
                    allele_matrix[i][ii] = distance
                    full_allele_matrix += allele_matrix
                return full_allele_matrix
            
            genome_addition: bool = True
            genome_loss: bool = True
            
            marker_matrix = np.zeros(array_size, np.float64)
            for i, sample_1 in enumerate(samples_list):
                for ii, sample_2 in enumerate(samples_list):
                    sample_1_alleles: list = [allele for allele in input_dict[marker][sample_1] if allele > 1]
                    sample_2_alleles: list = [allele for allele in input_dict[marker][sample_2] if allele > 1]

                    if len(sample_1_alleles) == len(sample_2_alleles): 
                        smaller_genotype = sample_1_alleles
                        larger_genotype = sample_2_alleles
                        marker_matrix += __generate_matrix_per_allele(smaller_genotype, larger_genotype)
                    else:
                        smaller_genotype = sample_1_alleles if len(sample_1_alleles) < len(sample_2_alleles) else sample_2_alleles
                        larger_genotype = sample_2_alleles if len(sample_1_alleles) < len(sample_2_alleles) else sample_1_alleles

                        genome_adjust_matrices: list = []
                        if genome_addition:
                            genome_addition_combinations = [smaller_genotype + list(i) for i in set(sorted(combinations_with_replacement(smaller_genotype, len(larger_genotype) - len(smaller_genotype))))]
                            for smaller_genotype_combination in genome_addition_combinations:
                                genome_adjust_matrices.append(__generate_matrix_per_allele(smaller_genotype_combination, larger_genotype))
                        if genome_loss:
                            genome_loss_combinations = [smaller_genotype + list(i) for i in set(sorted(combinations_with_replacement(larger_genotype, len(larger_genotype) - len(smaller_genotype))))]
                            for smaller_genotype_combination in genome_loss_combinations:
                                genome_adjust_matrices.append(__generate_matrix_per_allele(smaller_genotype_combination, larger_genotype))
                        if (not genome_addition) and (not genome_loss):
                            smaller_genotype = smaller_genotype + [1 * (len(larger_genotype) - len(smaller_genotype))]
                            genome_adjust_matrices.append(__generate_matrix_per_allele(smaller_genotype, larger_genotype))
                    
                        genome_adjust_matrices_stack = np.stack(genome_adjust_matrices)
                        marker_matrix += np.mean(genome_adjust_matrices_stack, axis = 0)

            return marker_matrix

        parsed_data, nuc_repeat_n_raw_per_marker_dict = _parse_data()

        # marker repeat calculator
        marker_nuc_repeat_dict: dict = _infer_marker_nuc_repeat_num(nuc_repeat_n_raw_per_marker_dict)
        samples_list: list = sorted(samples_list)
        array_size: tuple = (len(samples_list), len(samples_list))

        if samples_list:            
            marker_matrices = []
            for marker in sorted(parsed_data.keys()): marker_matrices.append(_get_marker_bruvo_dist(parsed_data, marker))

            # generate lower triangle distance matrix
            distance_matrix = [row[:i+1] for i, row in enumerate(np.tril(np.sum(marker_matrices, axis=0)).tolist())]

            distance_matrix_output = DistanceMatrix([sample.replace(' ', '_') for sample in samples_list], distance_matrix)
            string_io_output = io.StringIO("")
            if distance_matrix_output:
                Phylo.write(DistanceTreeConstructor().nj(distance_matrix_output), string_io_output, 'newick')
                newick_tree = string_io_output.getvalue()
                self._add_phylip_out(str(len(samples_list)) +"\n" + str(distance_matrix_output))
                self._add_tree_text(newick_tree)

    def _init_file(self) -> None:
        """
        Instantiate the text upload area.
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
            st.markdown('This tool should be used to process the table of Hydrangea alleles for the lab.')
            st.divider()
            st.markdown('Genetic distance is calculated using Bruvo\'s genetic distance ([Bruvo et al., 2004](https://pubmed.ncbi.nlm.nih.gov/15189230/)) '
                        'and uses repeat number differences to calculate genetic relatedness.')
            st.markdown("Annoying to implement, but special considerations are made for polyploid individuals of different allele numbers.")

    def _reset_state(self) -> None:
        """
        """
        for key in st.session_state.keys():
            del st.session_state[key]
# --------------------------------------------------
if __name__=="__main__":
    streamlit_app = App()
    streamlit_app._init_page()