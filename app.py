#!/usr/bin/env python3
import streamlit as st

from abi_sauce.services.file_manager import FileManager
from abi_sauce.services.sample_manager import SampleManager

from abi_sauce.app_pages.uploads import uploads_page
from abi_sauce.app_pages.samples import samples_page
from abi_sauce.app_pages.align import align_page
from abi_sauce.app_pages.viewer import viewer_page
from abi_sauce.app_pages.projects import projects_page

st.set_page_config(page_title="abi-sauce", layout="wide")

if "_manager" not in st.session_state:
    st.session_state._manager = FileManager()
if "_samples" not in st.session_state:
    st.session_state._samples = SampleManager(st.session_state._manager)

pages = [
    st.Page(uploads_page, title="Uploads", icon=":material/upload:", default=True),
    st.Page(viewer_page, title="Viewer", icon=":material/visibility:"),
    st.Page(samples_page, title="Samples", icon=":material/inventory_2:"),
    st.Page(align_page, title="Align", icon=":material/compare_arrows:"),
    st.Page(projects_page, title="Projects", icon=":material/folder_open:"),
]

nav = st.navigation(pages, position="sidebar", expanded=True)
nav.run()
