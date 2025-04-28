import os

project = "matRad"
copyright = "2025, e0404"
author = "e0404"

version = "3.1.0"
release = "3.1.0"

extensions = [
    'sphinx.ext.autodoc',
    'sphinxcontrib.matlab', 
    'sphinxcontrib.youtube',
    ]
primary_domain = "mat"

matlab_src_dir = os.path.join(os.path.dirname(__file__), '../')
matlab_auto_link = "all"

matlab_show_property_default_value = True
matlab_class_signature = True

autoclass_content = "both"
autodoc_member_order = "bysource"