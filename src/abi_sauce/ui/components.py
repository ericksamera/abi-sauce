from __future__ import annotations
from typing import List, Optional
import streamlit as st
import plotly.graph_objects as go

from abi_sauce.models import AssetBase, SequenceAsset, TraceAsset, AssetKind


def asset_table(assets: List[AssetBase]) -> Optional[str]:
    if not assets:
        st.info("No files uploaded yet.")
        return None
    # radio selection
    labels = [f"[{a.kind}] {a.name} • {getattr(a, 'length', '')}" for a in assets]
    ids = [a.id for a in assets]
    selected = st.radio("Assets", ids, format_func=lambda x: labels[ids.index(x)])
    return selected


def asset_detail(asset: AssetBase) -> None:
    if isinstance(asset, SequenceAsset):
        st.subheader("Sequence")
        st.code(asset.sequence[:5000] + ("…" if asset.length > 5000 else ""))
        with st.expander("Features"):
            if asset.features:
                st.json(asset.features)
            else:
                st.caption("(No features)")
        fa = asset.to_fasta()
        st.download_button("Download FASTA", data=fa, file_name=f"{asset.name}.fasta")

    elif isinstance(asset, TraceAsset):
        st.subheader("Trace (AB1)")
        if asset.sequence:
            st.text(f"Sequence length: {len(asset.sequence)} bases")
            fa = asset.to_fasta()
            if fa:
                st.download_button("Download called sequence (FASTA)", data=fa, file_name=f"{asset.name}.fasta")
        if any(asset.channels.values()):
            st.caption("Chromatogram (raw channel intensities)")
            fig = go.Figure()
            for base in "ACGT":
                y = asset.channels.get(base)
                if y:
                    fig.add_trace(go.Scatter(y=y, mode="lines", name=base))
            fig.update_layout(height=250, margin=dict(l=10, r=10, t=10, b=10))
            st.plotly_chart(fig, use_container_width=True)
        with st.expander("Metadata"):
            st.json(asset.meta)

    else:
        st.write(asset)