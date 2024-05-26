import _file_manager
import plotly.graph_objects as go
from streamlit_plotly_events import plotly_events

class TraceVisualizer:
    """
    """
    def __init__(self, input_trace_data: _file_manager.ABIFile, trimming_active) -> None:
        self.trace_data: _file_manager.ABIFile = input_trace_data
        self.trimming_active: bool = trimming_active
        self.base_bounds: list = self._calclate_base_bounds()

    def _calclate_base_bounds(self):
        positions = self.trace_data.base_positions
        halfway_points = [0] + [(positions[i] + positions[i - 1]) / 2 for i in range(1, len(positions))] + [(positions[-1] + (positions[-1] - positions[-2]) / 2)]
        
        # Pair the halfway points
        pairs = [(halfway_points[i], halfway_points[i + 1]) for i in range(len(halfway_points) - 1)]
        
        return pairs

    def _build_phred_rectangle(self, shape_tuple: tuple, height: int) -> tuple:
        """
        """
        x = [shape_tuple[0], shape_tuple[0], shape_tuple[1], shape_tuple[1], shape_tuple[0], None]
        y = [0, height, height, 0, 0, None]
        return [int(round(x)) if x is not None else None for x in x], [int(round(y)) if y is not None else None for y in y]

    def plot_electropherogram(self):
        self.rev_comp_adjust = 1 if not self.trace_data.reverse_complement else -1
        fig = go.Figure()
        max_peak_height = max([max(self.trace_data.channels[nucleotide]['peaks']) for nucleotide in self.trace_data.channels])

        plot_trim_left = self.trace_data.seq_object_raw.annotations['abif_raw']['PLOC2'][self.trace_data.trim_left] if not self.trace_data.reverse_complement else self.trace_data.seq_object_raw.annotations['abif_raw']['PLOC2'][-self.trace_data.trim_right]
        plot_trim_right = self.trace_data.seq_object_raw.annotations['abif_raw']['PLOC2'][-self.trace_data.trim_right] if not self.trace_data.reverse_complement else self.trace_data.seq_object_raw.annotations['abif_raw']['PLOC2'][self.trace_data.trim_left]    

        relative_heights = {
            'screen_height': max_peak_height + 200,
            'basepos_height': max_peak_height + 130,
            'basecall_height': max_peak_height + 75,
            'highlight_height': max_peak_height + 50
        }

        NUCLEOTIDES_KEY: dict = {
            'A': {'complement': 'T', 'color': 'green'},
            'T': {'complement': 'A', 'color': 'red'},
            'C': {'complement': 'G', 'color': 'blue'},
            'G': {'complement': 'C', 'color': 'black'},
            'Y': {'complement': 'R', 'color': '#ff3aff'},
            'R': {'complement': 'Y', 'color': '#ff3aff'},
            'W': {'complement': 'W', 'color': '#ff3aff'},
            'S': {'complement': 'S', 'color': '#ff3aff'},
            'K': {'complement': 'M', 'color': '#ff3aff'},
            'M': {'complement': 'K', 'color': '#ff3aff'},
            'D': {'complement': 'H', 'color': '#ff3aff'},
            'V': {'complement': 'B', 'color': '#ff3aff'},
            'H': {'complement': 'D', 'color': '#ff3aff'},
            'B': {'complement': 'V', 'color': '#ff3aff'},
            'N': {'complement': 'N', 'color': '#ff3aff'},
        }

        self.phred_x, self.phred_y = [], []
        for i, pair in enumerate(self.base_bounds):
            height = int(self.trace_data.phred_scores[i] / 60 * (5 / 6) * relative_heights['highlight_height'])
            x_add, y_add = self._build_phred_rectangle(pair, height)
            self.phred_x.extend(x_add)
            self.phred_y.extend(y_add)

        fig.add_trace(
            go.Scatter(
            x=self.phred_x[self.trace_data.trim_left * 6:-(self.trace_data.trim_right-1) * 6:self.rev_comp_adjust],
            y=self.phred_y[self.trace_data.trim_left * 6:-(self.trace_data.trim_right-1) * 6:self.rev_comp_adjust],
            hoverinfo='skip',
            name='Phred scores',
            fill="toself",
            marker=dict(
                color='#DFF0FA',
                opacity=0.5,
                ),
            ))
        # not active
        fig.add_trace(
            go.Scatter(
            x=self.phred_x[:self.trace_data.trim_left * 6:self.rev_comp_adjust] + self.phred_x[-(self.trace_data.trim_right-1) * 6::self.rev_comp_adjust],
            y=self.phred_y[:self.trace_data.trim_left * 6:self.rev_comp_adjust] + self.phred_y[-(self.trace_data.trim_right-1) * 6::self.rev_comp_adjust],
            hoverinfo='skip',
            name='Phred scores',
            fill="toself",
            marker=dict(
                color='lightgrey',
                opacity=0.25,
                ),
            ))


        channels_with_peaks = [(nuc, values) for nuc, values in self.trace_data.channels.items() if nuc in 'ATCG']
        for nuc, values in channels_with_peaks:
            nuc = nuc if not self.trace_data.reverse_complement else NUCLEOTIDES_KEY[nuc]['complement']
            fig.add_trace(
                go.Scattergl(
                    y=values['peaks'],
                    line=dict(width=1),
                    hoverinfo='skip',
                    name=nuc,
                    marker=dict(size=20, color=NUCLEOTIDES_KEY[nuc]['color'])))

        fig.add_vline(x=self.base_bounds[self.trace_data.trim_left][0], line_width=2, line_color='#88ccee', opacity=0.5)
        fig.add_vline(x=self.base_bounds[-self.trace_data.trim_right][1], line_width=2, line_color='#88ccee', opacity=0.5)

        fig.add_vline(x=self.base_bounds[self.trace_data.mott_trim_left][0], line_width=2, line_color='purple', opacity=0.5)
        fig.add_vline(x=self.base_bounds[-self.trace_data.mott_trim_right][1], line_width=2, line_color='purple', opacity=0.5)

        if self.trace_data.reverse_complement: basecalls = list(NUCLEOTIDES_KEY[char]['complement'] for char in self.trace_data.seq[::-1])
        else: basecalls = list(char for char in self.trace_data.seq)
            
        fig.add_trace(
            go.Scattergl(
                x=self.trace_data.base_positions[::self.rev_comp_adjust],
                y=[relative_heights['basecall_height'] for i in self.trace_data.seq][::self.rev_comp_adjust],
                name='Nucleotides',
                text=basecalls,
                mode="text",
                textfont={'color': list(NUCLEOTIDES_KEY[char]['color'] for char in basecalls)},
            ))

        fig.update_layout(
            margin=dict(l=0, r=0, t=50, b=0),
            dragmode=False,
            xaxis=dict(rangeslider=dict(visible=True, thickness=0.1), tickvals=[None], 
                       range=(plot_trim_left, plot_trim_right),
                       constrain='domain'),
            yaxis=dict(fixedrange=True, tickvals=[None], range=[0, relative_heights['screen_height']]),
            showlegend=False
        )
        fig.show(config={'displayModeBar': False})
        return plotly_events(fig, click_event=self.trimming_active)