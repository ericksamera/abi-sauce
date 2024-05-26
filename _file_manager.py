# file_manager.py
import streamlit as st
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import copy

class ABIFile:
    def __init__(self, streamlit_file: st.runtime.uploaded_file_manager.UploadedFile):
        self.filename: str = streamlit_file.name
        self.stem: str = '.'.join(streamlit_file.name.split('.')[:-1])
        self.seq_object_trimmed: SeqRecord = SeqIO.read(copy.deepcopy(streamlit_file), 'abi-trim')
        self.seq_object_raw: SeqRecord = SeqIO.read(streamlit_file, 'abi')
        self.modified = False
        self._parse_channels()
        self._parse_trimming()
        self._parse_properties()
        self._assign_color_code()
    def _parse_channels(self) -> None:
        """
        """
        self.reverse_complement: bool = False
        self.channels: dict = {
            'A': {
                'complement': 'T',
                'peaks': self.seq_object_raw.annotations['abif_raw']['DATA10']},
            'T': {
                'complement': 'A',
                'peaks': self.seq_object_raw.annotations['abif_raw']['DATA11']},
            'C': {
                'complement': 'G',
                'peaks': self.seq_object_raw.annotations['abif_raw']['DATA12']},
            'G': {
                'complement': 'C',
                'peaks': self.seq_object_raw.annotations['abif_raw']['DATA9']},
            }
        self.base_positions: list = self.seq_object_raw.annotations['abif_raw']['PLOC1']
    def _parse_trimming(self) -> None:
        """
        """

        self.mott_trim_left: int = self.seq_object_raw.seq.find(self.seq_object_trimmed.seq[0:5])
        self.mott_trim_right: int = len(self.seq_object_raw.seq) - len(self.seq_object_trimmed) - self.seq_object_raw.seq.find(self.seq_object_trimmed.seq[0:5])-1
        self.trim_left: int = self.mott_trim_left
        self.trim_right: int = self.mott_trim_right
    def _parse_properties(self) -> None:
        """
        """
        self.name: str = self.seq_object_raw.name
        self.well: str = self.seq_object_raw.annotations['abif_raw']['TUBE1'].decode(),
        self.quality_metrics: dict = {
            'pup_score': self.seq_object_raw.annotations['abif_raw']['PuSc1'] if 'PuSc1' in self.seq_object_raw.annotations['abif_raw'] else -1,
            'trace_score': self.seq_object_raw.annotations['abif_raw']['TrSc1'] if 'TrSc1' in self.seq_object_raw.annotations['abif_raw'] else -1,
            'crl_score': self.seq_object_raw.annotations['abif_raw']['CRLn1'] if 'CRLn1' in self.seq_object_raw.annotations['abif_raw'] else -1,
        }
        self.quality_metrics_str: str = f" ({self.quality_metrics['pup_score']}/{self.quality_metrics['trace_score']}/{self.quality_metrics['crl_score']}) "
        self.seq = str(self.seq_object_raw.seq)
        window_width = 60
        self.fasta_header: str = '.'.join(self.filename.split('.')[:-1])
        self.fasta_twoline: str = f'>{self.fasta_header}\n{self.seq}\n'
        self.fasta_str: str = f'>{self.fasta_header}\n' + '\n'.join([self.seq[i:i + window_width] for i in range(0, len(self.seq), 60)]) + '\n'
        self.phred_scores: list = self.seq_object_raw.letter_annotations['phred_quality']
    def _assign_color_code(self) -> None:
        """
        """

        trace_scoring: dict = {0: "🟩",1: "🟨",2: "🟥"}

        if all([
            self.quality_metrics['trace_score'] > 25,
            self.quality_metrics['pup_score'] > 20,
            self.quality_metrics['crl_score'] > 100,
            ]): self.quality_color: str = trace_scoring[0]
        elif all([
            self.quality_metrics['trace_score'] > 0,
            self.quality_metrics['pup_score'] > 0,
            self.quality_metrics['crl_score'] > 0,
            ]): self.quality_color: str = trace_scoring[1]
        else: self.quality_color: str = trace_scoring[2] 