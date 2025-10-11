#!/usr/bin/env python3
import streamlit as st
from abi_sauce.services.file_manager import FileManager
from abi_sauce.models import AssetKind
from abi_sauce.ui.components import asset_table, asset_detail

# ---- Session-scoped manager ----
if "_manager" not in st.session_state:
    st.session_state._manager = FileManager()
manager: FileManager = st.session_state._manager

st.set_page_config(page_title="abi-sauce", layout="wide")
st.title("🍏 abi-sauce — Uploads manager")
st.caption("Upload FASTA / GenBank / ApE / AB1 files. View sequences, features, and traces.")

# ---- Upload UI ----
uploaded = st.file_uploader(
    "Drop files here",
    type=["fa", "fasta", "fna", "gb", "gbk", "gbff", "ape", "ab1", "abi"],
    accept_multiple_files=True,
)

if uploaded:
    for f in uploaded:
        raw = f.getvalue()
        try:
            manager.add_bytes(name=f.name, raw=raw)
            st.toast(f"Imported {f.name}")
        except Exception as e:
            st.warning(f"{f.name}: {e}")

# ---- Asset list + detail ----
assets = manager.list()
selected_id = asset_table(assets)
if selected_id:
    asset = manager.get(selected_id)
    asset_detail(asset)

# ---- Utilities ----
with st.sidebar:
    st.header("Actions")
    if st.button("Clear all", use_container_width=True):
        manager.clear()
        st.rerun()
    st.markdown(":small: Built with BioPython + Streamlit + Plotly.")