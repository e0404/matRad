import os

# Compatibility patch: sphinxcontrib-matlabdomain 0.22.1 was written for
# Sphinx <8.x. Sphinx 8.x's autodoc base class now emits ':no-index:'
# (hyphenated) in generated RST, but MatObject/MatModule option_spec only
# register the old 'noindex' name, causing "unknown option" errors.
# Remove this block once sphinxcontrib-matlabdomain is updated upstream.
def _patch_mat_modindex() -> None:
    """
    Register 'mat-modindex' as a virtual docname in Sphinx's StandardDomain so
    it can be used directly in toctree directives (like genindex/modindex/search).

    Background: sphinxcontrib-matlabdomain generates mat-modindex.html via the
    domain index machinery, but the StandardDomain._virtual_doc_names dict only
    knows about 'genindex', 'modindex' (→ py-modindex), and 'search'.  Any name
    in that dict bypasses the normal "does this RST file exist?" check in both
    the toctree directive (sphinx/directives/other.py) and the toctree adapter
    (sphinx/environment/adapters/toctree.py) and is resolved directly to the
    given output filename.  Adding 'mat-modindex' here lets you write::

        .. toctree::
           mat-modindex

    and have it link to the generated mat-modindex.html page with the correct
    relative URL regardless of nesting depth.
    """
    from sphinx.domains.std import StandardDomain
    StandardDomain._virtual_doc_names['mat-modindex'] = (
        'mat-modindex', 'MATLAB Module Index'
    )

_patch_mat_modindex()


def _patch_matlabdomain_noindex() -> None:
    from docutils.parsers.rst import directives as _directives
    from sphinxcontrib.matlab import MatObject, MatModule
    for cls in (MatObject, MatModule):
        if "no-index" not in cls.option_spec:
            cls.option_spec["no-index"] = _directives.flag

_patch_matlabdomain_noindex()

project = "matRad"
copyright = "2025, e0404"
author = "e0404"

version = "3.2.2"
release = "3.2.2"

html_theme = "sphinx_rtd_theme" #pip install sphinx-rtd-theme

html_static_path = ["_static"]
html_css_files = ["style.css"]
html_logo = "../matRad/gfx/matrad_logo.png"
html_theme_options = {
    'logo_only': True,
    'navigation_depth': 5,
    'collapse_navigation': True,
}

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx_toolbox.collapse',
    'sphinxcontrib.matlab',
    'sphinxcontrib.youtube',
    'sphinx_togglebutton',
    ]

napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_custom_sections = [
    ('input', 'params_style'),
    ('output', 'params_style'),
    'call'
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
