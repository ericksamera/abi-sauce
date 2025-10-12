from __future__ import annotations
import streamlit as st

from abi_sauce.services.file_manager import FileManager
from abi_sauce.services.sample_manager import SampleManager
from abi_sauce.services.alignment import pairwise_align

def align_page():
    st.title("🧬 Pairwise alignment")
    st.caption("Pick any two Samples. Scores are adjustable. This lays groundwork for trace alignment & consensus assembly later.")

    fm: FileManager = st.session_state._manager
    sm: SampleManager = st.session_state._samples

    samples = sm.list()
    if len(samples) < 2:
        st.info("Upload at least two samples.")
        return

    left, right = st.columns(2)
    with left:
        a_id = st.selectbox("Sequence A", [s.id for s in samples], format_func=lambda sid: sm.get(sid).name)
    with right:
        b_id = st.selectbox("Sequence B", [s.id for s in samples if s.id != a_id], format_func=lambda sid: sm.get(sid).name)

    seqA = (sm.get(a_id).effective_sequence(fm._assets) or "").upper()
    seqB = (sm.get(b_id).effective_sequence(fm._assets) or "").upper()

    st.caption(f"{sm.get(a_id).name}: {len(seqA)} bp  •  {sm.get(b_id).name}: {len(seqB)} bp")
    if not seqA or not seqB:
        st.warning("One or both samples lack a sequence.")
        return

    mode = st.selectbox("Mode", ["globalms", "localms", "globalxx", "localxx"], index=0)
    c1, c2, c3, c4 = st.columns(4)
    with c1: match = st.number_input("match", value=2.0)
    with c2: mismatch = st.number_input("mismatch", value=-1.0)
    with c3: gap_open = st.number_input("gap_open", value=-10.0)
    with c4: gap_extend = st.number_input("gap_extend", value=-0.5)

    if st.button("Align"):
        results = pairwise_align(seqA, seqB, mode=mode, match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend, top_n=1)
        if not results:
            st.error("No alignment produced.")
            return
        r = results[0]
        st.markdown(f"**Score:** {r.score:.2f} • **Identity:** {100*r.identity:.2f}%")
        mid = "".join("|" if a == b and a != "-" else " " for a, b in zip(r.a_aln, r.b_aln))
        st.code(f"{r.a_aln}\n{mid}\n{r.b_aln}")
