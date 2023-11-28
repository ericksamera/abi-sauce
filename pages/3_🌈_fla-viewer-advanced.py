#!/usr/bin/env python3
__description__ =\
"""
Purpose: Streamlit wrapper for fla-viewer-advanced.
"""
__author__ = "Erick Samera"
__version__ = "1.2.0"
__comments__ = "stable enough"
# --------------------------------------------------
import streamlit as st
import streamlit_ext as ste
# --------------------------------------------------
from scipy.signal import find_peaks
import plotly.graph_objects as go
from Bio import SeqIO
import csv
# --------------------------------------------------
class App:
    def __init__(self):
        """
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
    def _init_page(self) -> None:
        """
        Function instantiates the main page.
        """
        self._init_sidebar()
        if 'UPLOADED_FSA' not in st.session_state: self._init_file_uploader()
        else: self._update_plot_window()
        return None
    def _init_sidebar(self) -> None:
        """
        Instantiate the sidebar.
        """
    
        with st.sidebar:
            st.title(f"{self.emoji} abi-sauce | {self.title}")
            st.markdown('This script is intended for processing `.fsa` files for fragment length analysis (FLA).')
            #st.markdown('Check out the better-maintained command-line interface on [GitHub](https://github.com/KPU-AGC/general-resources/blob/main/sanger-processing/sanger-sequence-trim.py)!')

            if 'UPLOADED_FSA' in st.session_state:
                st.divider()

                if 'FAILED_FILES' in st.session_state:
                        with st.expander('**⚠ WARNINGS:**', expanded=True):
                            for failed_file_name in st.session_state.FAILED_FILES:
                                st.error(f"{failed_file_name} failed to process!", icon="⚠")

                with st.expander('**TRACE FILES:**', expanded=True):
                    st.session_state.SORTED_LIST = st.session_state.PROCESSED_FLA.values()
                    st.session_state.SELECTED_TRACE = st.selectbox(
                        'Select trace file to view:', 
                        options=[trace_object['name'] for trace_object in st.session_state.SORTED_LIST])
                st.button('Reset & Upload New', type='primary', on_click=self._reset_state, use_container_width=True)

            st.divider()

            with st.expander('MORE INFO'):
                st.markdown(
                    'Fragment length analysis (FLA) is a technique for determining '
                    'the size of DNA fragments which is useful for genotyping ' 
                    'microsatellites.')
            
                st.markdown(
                    'In "the good old days", gel electrophoresis would be used to separate'
                    ' bands out for genotyping. With our fancy SeqStudio, capillary '
                    'electrophoresis with fluorescently tagged fragments '
                    'gives us much better resolution--even down to the nucleotide.')

                st.markdown(
                    'Transducing fluorescence data to digital data allows us to '
                    'process data with specialized software to facilitate genotyping.')
                st.caption('')
                st.caption(f'[@{__author__}](https://github.com/ericksamera)\t|\tv{__version__}\t|\t{__comments__}')
    def _init_file_uploader(self) -> None:
        """
        Function initializes file upload handling.
        """
        st.header('Select fragment length analysis (FLA) files (`.fsa`)')
        uploaded_files: list = st.file_uploader(
                        'Upload `.fsa` files which come right off the SeqStudio.',
                        type=['fsa'],
                        accept_multiple_files=True)
    
        
        with st.form("upload form", clear_on_submit=True):
            submit_button = st.form_submit_button("Process files", on_click=self._upload_files, args=(uploaded_files,))
            if submit_button and not uploaded_files:
                st.error('Select some files!')
    def _upload_files(self, _st_uploaded_files: list) -> None:
        """
        Function handles file upload and adds files to cache.
        """
        if not _st_uploaded_files: return None
        st.session_state.UPLOADED_FSA = [file for file in _st_uploaded_files if 'QCReport' not in file.name]
        st.session_state.PROCESSED_FLA = self._process_files(st.session_state.UPLOADED_FSA)
    def _process_files(self, _st_uploaded_files) -> None:
        """
        Function creates a new instance of the copied file and processes it.
        """
        processed_files = {}

        for file in sorted(_st_uploaded_files, key=lambda x: x.name):
            fsa_object = SeqIO.read(file, 'abi')
            try: smap2_list = fsa_object.annotations['abif_raw']['SMap2']
            except KeyError: 
                if 'FAILED_FILES' not in st.session_state: st.session_state.FAILED_FILES = []
                st.session_state.FAILED_FILES.append(file.name)

            channels_dict = {
                '6-FAM': fsa_object.annotations['abif_raw']['DATA9'],
                'VIC': fsa_object.annotations['abif_raw']['DATA10'],
                'NED': fsa_object.annotations['abif_raw']['DATA11'],
                'PET': fsa_object.annotations['abif_raw']['DATA12'],
                'LIZ': fsa_object.annotations['abif_raw']['DATA205'],
                }

            seq_object_dict = {
                'name': fsa_object.name,
                'channels_peaks': channels_dict,
                'colors': {
                    '6-FAM': "blue",
                    'VIC': "green",
                    'NED': "black",
                    'PET': "red",
                    'LIZ': "orange"},
                'smap': smap2_list,
                'detected_peaks_x': {dye_name: self._detect_channels(channels, smap2_list)[0] for dye_name, channels in channels_dict.items()},
                'detected_peaks_y': {dye_name: self._detect_channels(channels, smap2_list)[1] for dye_name, channels in channels_dict.items()}
                }
            processed_files[fsa_object.name] = seq_object_dict
        return processed_files
    def _detect_channels(self, channel_list: list, smap2_list: list, threshold_peak_height: int = 75) -> tuple:
        """
        """
        indices = find_peaks(channel_list)[0]

        filtered_y = [channel_list[i] for i in indices if channel_list[i] > threshold_peak_height]
        filtered_x = [smap2_list[i] for i in indices if channel_list[i] > threshold_peak_height]
        return filtered_x, filtered_y

    def _plot_total_trace(self, _trace_dict):
        """
        """

        all_x_vals = []
        for dye_color in _trace_dict['channels_peaks'].keys():
            all_x_vals += [x for x in _trace_dict['detected_peaks_x'][dye_color]]
        max_x_val = max(all_x_vals)

        all_y_vals = []
        for dye_color in _trace_dict['channels_peaks'].keys():
            all_y_vals += [y for x, y in zip(_trace_dict['detected_peaks_x'][dye_color], _trace_dict['detected_peaks_y'][dye_color]) if x > 25]
        max_y_val = max(all_y_vals) * 1.05

        max_heights_with_pad = {key: max([y for x, y in zip(_trace_dict['detected_peaks_x'][key], _trace_dict['detected_peaks_y'][key]) if x > 25]) * 1.05 for key in _trace_dict['channels_peaks'].keys() if [y for x, y in zip(_trace_dict['detected_peaks_x'][key], _trace_dict['detected_peaks_y'][key]) if x > 25]}

        fig = go.Figure()

        for dye_color, values in _trace_dict['channels_peaks'].items():
            # plot the peak data, x-axis is nucleotides
            fig.add_trace(
                go.Scatter(
                    mode='lines',
                    x=_trace_dict['smap'],
                    y=values,
                    hoverinfo='skip',
                    marker=dict(
                        size=0.8,
                        color=_trace_dict['colors'][dye_color]),
                    name=dye_color,
                    fill='tozeroy',
                    legendgroup=dye_color
                )
            )
            
            fig.add_trace(
                go.Scatter(
                    mode='markers',
                    x=_trace_dict['detected_peaks_x'][dye_color],
                    y=_trace_dict['detected_peaks_y'][dye_color],
                    hoverinfo='x+y+text+name' if dye_color not in ('LIZ') else 'skip',
                    hovertemplate = "size (bp): %{x}<br>" + "height: %{y}<br>" if dye_color not in ('LIZ') else "",
                    marker=dict(
                        size=0.8,
                        color=_trace_dict['colors'][dye_color]),
                    name=dye_color,
                    #fill='tozeroy',
                    showlegend=False,
                    legendgroup=dye_color
                )
            )

        fig.update_layout(
            legend_title_text='Channels',
            dragmode="zoom",
            xaxis=dict(
                range=(0, max_x_val+20),
                #rangeslider=dict(visible=True, thickness=0.1, range=(0, max_x_val+20))
                ),
            yaxis=dict(range=(0, max_y_val))
            )

        buttons_list = [dict(args=[{"yaxis": dict(autorange=True)}],label="Scale to tallest",method="relayout")]

        for dye_color in _trace_dict['channels_peaks'].keys():
            if not dye_color in max_heights_with_pad: continue
            buttons_list.append(dict(args=[{"yaxis": dict(autorange=False, range=(0, max_heights_with_pad[dye_color]))}],label=f"Scale to {dye_color}",method="relayout"))

        fig.update_layout(
            updatemenus=[
                dict(
                    type = "buttons",
                    #direction = "left",
                    buttons=buttons_list,
                    active=0,
                    showactive=True,
                ),
                ]
            )

        return fig
    def _per_channel_bar_plot(self, _trace_dict: str, _dye_name: str, _average_height: float):
        """
        """

        fig = go.Figure()
        fig.add_trace(go.Scattergl(
            mode='lines',
            x=_trace_dict['smap'],
            y=_trace_dict['channels_peaks'][_dye_name],
            name=_dye_name,
            hoverinfo='skip',
            fill='tozeroy',
            marker=dict(
                color=_trace_dict['colors'][_dye_name]
                ),
            )
        )
        fig.add_trace(go.Scattergl(
            mode='markers',
            x=_trace_dict['detected_peaks_x'][_dye_name],
            y=_trace_dict['detected_peaks_y'][_dye_name],
            #hoverinfo='skip',
            hovertemplate = "size (bp): %{x}<br>" + "height: %{y}<br>",
            marker=dict(
                size=0.5,
                color=_trace_dict['colors'][_dye_name]
                ),
            name=_dye_name,
            showlegend=False,
            ))

        fig.update_layout(
            margin=dict(l=50, r=150, t=0, b=0),
            height=100,
            legend_title_text='Channels',
            #xaxis=dict(rangeslider=dict(visible=True, thickness=0.1, range=(0,max([val for val in _trace_dict['peaks'][color]['x'] if val]))), constrain='domain'),
            xaxis=dict(range=(0, max([val for val in _trace_dict['detected_peaks_x'][_dye_name] if val])*1.05)),
            yaxis=dict(range=(0, int(_average_height*1.25))),
            modebar=dict(remove=[ "autoScale2d", "autoscale", "editInChartStudio", "editinchartstudio", "hoverCompareCartesian", "hovercompare", "lasso", "lasso2d", "orbitRotation", "orbitrotation", "pan", "pan2d", "pan3d", "reset", "resetCameraDefault3d", "resetCameraLastSave3d", "resetGeo", "resetSankeyGroup", "resetScale2d", "resetViewMapbox", "resetViews", "resetcameradefault", "resetcameralastsave", "resetsankeygroup", "resetscale", "resetview", "resetviews", "select", "select2d", "sendDataToCloud", "senddatatocloud", "tableRotation", "tablerotation", "toImage", "toggleHover", "toggleSpikelines", "togglehover", "togglespikelines", "toimage", "zoom", "zoom2d", "zoom3d", "zoomIn2d", "zoomInGeo", "zoomInMapbox", "zoomOut2d", "zoomOutGeo", "zoomOutMapbox", "zoomin", "zoomout"]),
            )
        
        return fig
    def _plot_per_channel_traces(self, _trace_dict) -> None:
        """
        """
        with st.expander('Individual channels'):
            st.session_state.FULL_GENOTYPES = st.checkbox("Show all predicted peaks", value=False)
            st.session_state.SHOW_INDIV_TABLES = st.checkbox("Show table with individual channels", value=False)

            for dye_color in _trace_dict['colors']:
                averages_past_threshold = []
                full_genotype = []
                for x_val, y_val in zip(_trace_dict['detected_peaks_x'][dye_color], _trace_dict['detected_peaks_y'][dye_color]):
                    if x_val < 25: continue
                    averages_past_threshold.append(y_val)
                    full_genotype.append({'size (bp)': x_val, 'height': y_val})
                
                if averages_past_threshold: 
                    average_height = max(averages_past_threshold)
                    st.plotly_chart(
                        self._per_channel_bar_plot(_trace_dict, dye_color, average_height),
                        use_container_width=True)
            
                if full_genotype:
                    average_height = sum([genotype['height'] for genotype in full_genotype])/len(full_genotype)
                    filtered_genotype = [genotype for genotype in full_genotype if genotype['height']>average_height]
                    if st.session_state.SHOW_INDIV_TABLES:
                        if dye_color not in ('LIZ'): st.dataframe(full_genotype if st.session_state.FULL_GENOTYPES else filtered_genotype, use_container_width=True)

        # with st.expander('Predicted genotype'):
        #     full_genotype_table = []
        #     for color in _trace_dict['predicted_genotypes']:
        #         full_genotype_table += _trace_dict['predicted_genotypes'][color]
        #     st.dataframe(sorted(full_genotype_table, key=lambda x: float(x['position (bp)']), reverse=False), use_container_width=True)
        return None
    def _update_plot_window(self) -> None:
        """
        """
        st.header(f"{st.session_state.SELECTED_TRACE}")
        #st.text(f"{st.session_state.PROCESSED_FLA[f'{st.session_state.SELECTED_TRACE}']}")
        st.plotly_chart(
            self._plot_total_trace(st.session_state.PROCESSED_FLA[f"{st.session_state.SELECTED_TRACE}"]),
           use_container_width=True)
        self._plot_per_channel_traces(st.session_state.PROCESSED_FLA[f"{st.session_state.SELECTED_TRACE}"])
        return None
    def _reset_state(self) -> None:
        """
        """
        for key in st.session_state.keys():
            del st.session_state[key]
# --------------------------------------------------
if __name__=="__main__":
    streamlit_app = App()
    streamlit_app._init_page()