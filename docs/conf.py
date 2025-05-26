import os

project = "matRad"
copyright = "2025, e0404"
author = "e0404"

version = "3.1.0"
release = "3.1.0"

html_theme = "sphinx_rtd_theme" #pip install sphinx-rtd-theme

extensions = [
    'sphinx.ext.autodoc',
    'sphinxcontrib.matlab', 
    'sphinxcontrib.youtube',
    'sphinx_togglebutton',
    ]
primary_domain = "mat"

matlab_src_dir = os.path.join(os.path.dirname(__file__), '../')
matlab_auto_link = "basic"

matlab_show_property_default_value = True
matlab_show_property_specs = True
matlab_class_signature = True
matlab_short_links = True

autoclass_content = "both"
autodoc_member_order = "bysource"