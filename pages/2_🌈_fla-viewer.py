#!/usr/bin/env python3
__description__ =\
"""
Purpose: Streamlit wrapper for sanger-sequence-trim.
"""
__author__ = "Erick Samera"
__version__ = "1.5.0"
__comments__ = "stable enough"
# --------------------------------------------------
import streamlit as st
import streamlit_ext as ste
# --------------------------------------------------
import plotly.graph_objects as go
from Bio import SeqIO
import zipfile
import io
import csv
# --------------------------------------------------
class App:
    def __init__(self):
        """
        """
    def _init_page(self) -> None:
        """
        Function instantiates the main page.
        """
        self.title = "FLA-viewer"
        self.emoji = ':rainbow:'
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
        self._init_sidebar()
        if 'UPLOADED_FLA' not in st.session_state: self._init_file_uploader()
        else: self._update_plot_window()

    def _init_sidebar(self) -> None:
        """
        Instantiate the sidebar.
        """
    
        with st.sidebar:
            st.title(f"{self.emoji} abi-sauce | {self.title}")
            st.markdown('This script is intended for processing `.fsa.csv` files for fragment length analysis (FLA).')
            #st.markdown('Check out the better-maintained command-line interface on [GitHub](https://github.com/KPU-AGC/general-resources/blob/main/sanger-processing/sanger-sequence-trim.py)!')

            if 'UPLOADED_FLA' in st.session_state:
                st.divider()
                with st.expander('**TRACE FILES:**', expanded=True):

                    st.session_state.SORTED_LIST = st.session_state.PROCESSED_FLA.values()
                    st.session_state.SELECTED_TRACE = st.radio(
                        'Select trace file to view:', 
                        options=[trace_object['name'] for trace_object in st.session_state.SORTED_LIST])
                # with st.expander('**DOWNLOAD OPTIONS:**'):
                #     st.session_state.DEFAULT_FILENAME = 'abi-sauce-trim'
                #     st.session_state.USER_FILENAME = st.text_input('File name', placeholder='abi-sauce-trim')
                #     st.session_state.FILENAME =  st.session_state.USER_FILENAME if st.session_state.USER_FILENAME else st.session_state.DEFAULT_FILENAME
                #     st.session_state.CONCATENATE: bool = st.checkbox("Concatenate entries into single fasta.", value=True)
                #     if not st.session_state.CONCATENATE: st.caption('Individual FASTAs will be compiled into a single ZIP file.')
            st.divider()

            with st.expander('MORE INFO'):
                st.markdown(
                    'Fragment length analysis (FLA) is a technique for determining '
                    'the size of DNA fragments which is useful for genotyping ' 
                    'microsatellites.')
            
                st.markdown(
                    'In "the old days", agarose gel electrophoresis would be used '
                    'to separate bands out for genotyping. With our fancy SeqStudio, '
                    'capillary electrophoresis with fluorescently tagged fragments '
                    'gives us much better resolution--even down to the nucleotide.')

                st.markdown(
                    'Transducing fluorescence data to digital data allows us to '
                    'process data with specialized software to facilitate genotyping.')
                st.caption('')
                st.caption(f'[@{__author__}](https://github.com/ericksamera)\t|\tv{__version__}\t|\t{__comments__}')
    def _init_file_uploader(self) -> None:
        """
        """
        st.header('Select fragment length analysis (FLA) files (`.fsa.csv`)')
        uploaded_files: list = st.file_uploader(
                        'Upload `.fsa.csv` files which come right off the SeqStudio.',
                        type=['csv'],
                        accept_multiple_files=True)
    
        
        with st.form("upload form", clear_on_submit=False):
            submit_button = st.form_submit_button("Process files", on_click=self._upload_files, args=(uploaded_files,))
            if submit_button and not uploaded_files:
                st.error('Select some files!')
    def _upload_files(self, _st_uploaded_files: list) -> None:
        if not _st_uploaded_files: return None
        st.session_state.UPLOADED_FLA = [file for file in _st_uploaded_files if 'QCReport' not in file.name]
        st.session_state.PROCESSED_FLA = self._process_files(st.session_state.UPLOADED_FLA)
    def _process_files(self, _st_uploaded_files) -> None:
        """
        Function creates a new instance of the copied file and processes it.
        """
        processed_files = {}

        for file in sorted(_st_uploaded_files, key=lambda x: x.name):
            processed_name = str(file.name).strip('.fla.csv')
            processed_files[processed_name] = {
                'name': processed_name.strip('.fla.csv'),
                'channels': {},
                'peaks': {}
            }
            for line_dict in csv.DictReader(io.StringIO(file.getvalue().decode("utf-8"))):
                try: channel_color = line_dict['Dye'].lower()
                except KeyError:
                    try: channel_color = line_dict['Dye Color'].lower()
                    except KeyError: continue
                if channel_color not in processed_files[processed_name]['channels']: 
                    processed_files[processed_name]['channels'][channel_color] = {
                        'name': channel_color,
                        'color_name': channel_color,
                        'x': [],
                        'y': []
                    }
                    processed_files[processed_name]['peaks'][channel_color] = {
                        'x': [],
                        'y': []
                    }
                try:
                    start = float(line_dict['Begin BP'])
                    end = float(line_dict['End BP'])
                    middle = float(line_dict['Size'])

                    height = float(line_dict['Height'])
                except KeyError:
                    try:
                        start = float(line_dict['Begin Point (Base Pairs)'])
                        end = float(line_dict['End Point (Base Pairs)'])
                        middle = float(line_dict['Size'])

                        height = float(line_dict['Height'])
                    except KeyError: continue
                except ValueError: continue
                processed_files[processed_name]['channels'][channel_color]['x'] += [None,start,middle,end,None]
                processed_files[processed_name]['channels'][channel_color]['y'] += [None,0,height,0,None]

                processed_files[processed_name]['peaks'][channel_color]['x'] += [middle]
                processed_files[processed_name]['peaks'][channel_color]['y'] += [height]
            processed_files[processed_name]['predicted_genotypes'] = self._predict_genotypes(processed_files[processed_name])
        return processed_files
    def _predict_genotypes(self, _processed_file_dict: dict) -> None:
        """
        """
        predicted_genotypes_dict = {}
        for color in _processed_file_dict['peaks']:
            if color in ('orange'): continue
            if color not in predicted_genotypes_dict: predicted_genotypes_dict[color] = []
            try: max_height = max([y_value for x_value, y_value in zip(_processed_file_dict['peaks'][color]['x'], _processed_file_dict['peaks'][color]['y']) if x_value>50])
            except ValueError: continue
            for x_value, y_value in zip(_processed_file_dict['peaks'][color]['x'], _processed_file_dict['peaks'][color]['y']):
                if x_value<50: continue
                if y_value > max_height/5:
                    predicted_genotypes_dict[color].append(
                        {'channel': color.upper(), 'position (bp)': f"{x_value:.1f}", 'height': f"{y_value:.0f}"})
        return predicted_genotypes_dict
    def _plot_total_trace(self, _trace_dict):
        """
        """

        fig = go.Figure()

        all_x_vals = []
        for x_vals in _trace_dict['channels'].values():
            all_x_vals += [x for x in x_vals['x'] if x]
        max_x_val = max(all_x_vals)

        maxes_no_ladder = []
        for color in [key for key in _trace_dict['channels'].keys() if key not in ['orange']]:
            try:
                list_to_add = max([i for i in [y_vals if x_vals>50 else None for x_vals, y_vals in zip(_trace_dict['channels'][color]['x'], _trace_dict['channels'][color]['y']) if x_vals] if i])
                maxes_no_ladder.append(list_to_add)
            except:continue
        max_ladder = max([i for i in _trace_dict['channels']['orange']['y'] if i])
        average_height = max(maxes_no_ladder) if max(maxes_no_ladder)>max_ladder else max_ladder

        for color, values in _trace_dict['channels'].items():
            fig.add_trace(
                go.Scatter(
                    mode='lines',
                    x=values['x'],
                    y=values['y'],
                    hoverinfo='skip',
                    marker=dict(
                        size=0.8,
                        color=color.lower() if color.lower() not in ['yellow'] else 'black'),
                    name=color.upper(),
                    fill='tozeroy',
                    legendgroup=color.lower()
                ))
            fig.add_trace(
                go.Scatter(
                    mode='markers',
                    x=_trace_dict['peaks'][color]['x'],
                    y=_trace_dict['peaks'][color]['y'],
                    #hoverinfo='skip',
                    marker=dict(
                        size=0.8,
                        color=color.lower() if color.lower() not in ['yellow'] else 'black'),
                    name=color.upper(),
                    showlegend=False,
                    legendgroup=color.lower()
                ))

        fig.update_layout(
            legend_title_text='Channels',
            xaxis=dict(
                range=(0, max_x_val+20),
                rangeslider=dict(visible=True, thickness=0.1, range=(0, max_x_val+20)), constrain='domain'),
            yaxis=dict(range=(0, average_height))
            )
        return fig
    def _per_channel_bar_plot(self, _per_channel_trace_dict, _per_channel_peaks_dict):
        """
        """

        averages_past_50 = []
        for x_val, y_val in zip(_per_channel_peaks_dict['x'], _per_channel_peaks_dict['y']):
            if x_val<50:continue
            averages_past_50.append(y_val)
        average_height = max(averages_past_50)

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            mode='lines',
            x=_per_channel_trace_dict['x'],
            y=_per_channel_trace_dict['y'],
            name=_per_channel_trace_dict['color_name'].upper(),
            hoverinfo='skip',
            fill='tozeroy',
            marker=dict(
                color=_per_channel_trace_dict['color_name'].lower() if _per_channel_trace_dict['color_name'].lower() not in ['yellow'] else 'black'
                ),
            )
        )
        fig.add_trace(go.Scatter(
            mode='markers',
            x=_per_channel_peaks_dict['x'],
            y=_per_channel_peaks_dict['y'],
            #hoverinfo='skip',
            marker=dict(
                size=0.5,
                color=_per_channel_trace_dict['color_name'].lower() if _per_channel_trace_dict['color_name'].lower() not in ['yellow'] else 'black'),
            name=_per_channel_trace_dict['color_name'].upper(),
            showlegend=False,
            ))

        fig.update_layout(
            margin=dict(l=50, r=150, t=0, b=0),
            height=100,
            legend_title_text='Channels',
            #xaxis=dict(rangeslider=dict(visible=True, thickness=0.1, range=(0,max([val for val in _trace_dict['peaks'][color]['x'] if val]))), constrain='domain'),
            xaxis=dict(range=(0, max([val for val in _per_channel_trace_dict['x'] if val])+50)),
            yaxis=dict(range=(0, average_height)),
            modebar=dict(remove=[ "autoScale2d", "autoscale", "editInChartStudio", "editinchartstudio", "hoverCompareCartesian", "hovercompare", "lasso", "lasso2d", "orbitRotation", "orbitrotation", "pan", "pan2d", "pan3d", "reset", "resetCameraDefault3d", "resetCameraLastSave3d", "resetGeo", "resetSankeyGroup", "resetScale2d", "resetViewMapbox", "resetViews", "resetcameradefault", "resetcameralastsave", "resetsankeygroup", "resetscale", "resetview", "resetviews", "select", "select2d", "sendDataToCloud", "senddatatocloud", "tableRotation", "tablerotation", "toImage", "toggleHover", "toggleSpikelines", "togglehover", "togglespikelines", "toimage", "zoom", "zoom2d", "zoom3d", "zoomIn2d", "zoomInGeo", "zoomInMapbox", "zoomOut2d", "zoomOutGeo", "zoomOutMapbox", "zoomin", "zoomout"]),
            )
        return fig
    def _plot_per_channel_traces(self, _trace_dict):
        """
        """
        with st.expander('Individual channels'):
            for color, values in [(color, value) for color, value in _trace_dict['channels'].items() if color not in ['orange']] + [('orange', _trace_dict['channels']['orange'])]:
                if not any([y_val if x_val>50 else None for x_val, y_val in zip(values['x'], values['y']) if x_val]): continue

                st.plotly_chart(
                    self._per_channel_bar_plot(_trace_dict['channels'][color], _trace_dict['peaks'][color]),
                    use_container_width=True)
                if color not in ('orange'): st.dataframe(_trace_dict['predicted_genotypes'][color], use_container_width=True)
        
        with st.expander('Full predicted genotypes'):
            full_genotype_table = []
            for color in _trace_dict['predicted_genotypes']:
                full_genotype_table += _trace_dict['predicted_genotypes'][color]
            st.dataframe(sorted(full_genotype_table, key=lambda x: float(x['position (bp)']), reverse=False), use_container_width=True)
    def _update_plot_window(self):
        """
        """
        st.header(f"{st.session_state.SELECTED_TRACE}")
        #st.text(f"{st.session_state.PROCESSED_FLA[f'{st.session_state.SELECTED_TRACE}']}")
        st.plotly_chart(
            self._plot_total_trace(st.session_state.PROCESSED_FLA[f"{st.session_state.SELECTED_TRACE}"]),
            use_container_width=True)
        self._plot_per_channel_traces(st.session_state.PROCESSED_FLA[f"{st.session_state.SELECTED_TRACE}"])
    #     help_str = "Trimmed with Mott's trimming algorithm!" if st.session_state.TRIM_STR == '_trimmed' else "Untrimmed!"
    #     st.caption(f'FASTA sequence: ({len(st.session_state.PROCESSED_FLA[f"{st.session_state.SELECTED_TRACE}.ab1"][st.session_state.TRIM_STR])} bp)', help=help_str)
    #     st.code(st.session_state.PROCESSED_FLA[f"{st.session_state.SELECTED_TRACE}.ab1"][st.session_state.TRIM_STR].format("fasta"))
    #     if st.session_state.TRIM_STR == '_trimmed':
    #         with st.expander(f'Untrimmed FASTA ({len(st.session_state.PROCESSED_FLA[f"{st.session_state.SELECTED_TRACE}.ab1"]["_raw"])} bp):'):
    #             st.code(st.session_state.PROCESSED_FLA[f"{st.session_state.SELECTED_TRACE}.ab1"]["_raw"].format("fasta"))
    #         st.subheader('Trimming statistics')
    #         left_column, right_column = st.columns(2)
    #         with left_column:
    #             l_trim_count: int = st.session_state.PROCESSED_FLA[f"{st.session_state.SELECTED_TRACE}.ab1"]["left_trim"]
    #             r_trim_count: int = st.session_state.PROCESSED_FLA[f"{st.session_state.SELECTED_TRACE}.ab1"]["right_trim"]+1
    #             total_trim_count = l_trim_count+r_trim_count

    #             phred_scores_raw = st.session_state.PROCESSED_FLA[f"{st.session_state.SELECTED_TRACE}.ab1"]['_raw'].letter_annotations['phred_quality']
    #             phred_scores_trimmed = st.session_state.PROCESSED_FLA[f"{st.session_state.SELECTED_TRACE}.ab1"]['_trimmed'].letter_annotations['phred_quality']

    #             average_phred_raw = sum(phred_scores_raw)/len(phred_scores_raw)
    #             average_phred_trimmed = sum(phred_scores_trimmed)/len(phred_scores_trimmed)

    #             st.markdown(f'Left trim: `{l_trim_count}` | Right trim:`{r_trim_count}` (Total: `{total_trim_count}`)', help='The right side is usually worse.')
    #             st.markdown(f'Average PHRED score (before): `{average_phred_raw:.1f}` | (after) `{average_phred_trimmed:.1f}`')

    #         with right_column:
    #             """"""
# --------------------------------------------------
if __name__=="__main__":
    streamlit_app = App()
    streamlit_app._init_page()
