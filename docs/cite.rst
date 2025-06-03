.. _cite:

===================
How to cite matRad?
===================

Whenever you use matRad, we kindly ask yuo to cite our publications and/or the software (with appropriate version) itself.
This helps us to keep the project alive and to continue the development of matRad.

.. _citepubs:
Cite our publications
---------------------

.. image:: https://img.shields.io/endpoint?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.1002%2Fmp.12251&style=flat&color=blue
   :alt: Citation Badge


The main citation for matRad is `Development of the open-source dose calculation and optimization toolkit matRad <https://www.ncbi.nlm.nih.gov/pubmed/28370020>`_, which details the development process of matRad and provides a detailed evaluation of the dose calculation and optimization engine.

.. collapse:: matRad main publication - BibTeX
    
    .. code-block:: bibtex

        @article{matRad2017,
            title = {Development of the open-source dose calculation and optimization toolkit {{matRad}}},
            author = {Wieser, Hans-Peter and Cisternas, Eduardo and Wahl, Niklas and Ulrich, Silke and Stadler, Alexander and Mescher, Henning and M{\"u}ller, Lucas-Raphael and Klinge, Thomas and Gabrys, Hubert and Burigo, Lucas and Mairani, Andrea and Ecker, Swantje and Ackermann, Benjamin and Ellerbrock, Malte and Parodi, Katia and J{\"a}kel, Oliver and Bangert, Mark},
            year = {2017},
            journal = {Medical Physics},
            volume = {44},
            number = {6},
            pages = {2556--2568},
            issn = {2473-4209},
            doi = {10.1002/mp.12251},
        }

----

The first paper about matRad was submitted to the `2015 World Congress on Medical Physics & Biomedical Engineering <http://wc2015.org/>`_ by Eduardo Cisternas Jimenéz:
`matRad – a multi-modality open source 3D treatment planning toolkit <https://github.com/e0404/matRad/wiki/documents/matRad_Full_Paper_Eduardo.pdf>`_.

.. collapse:: matRad conference paper - BibTeX
    
    .. code-block:: bibtex

        @inproceedings{Cisternas2015a,
            title = {{{matRad}} - a multi-modality open source {{3D}} treatment planning toolkit},
            booktitle = {World {{Congress}} on {{Medical Physics}} and {{Biomedical Engineering}}, {{June}} 7-12, 2015, {{Toronto}}, {{Canada}}},
            author = {Cisternas, Eduardo and Mairani, Andrea and Ziegenhein, Peter and J{\"a}kel, Oliver and Bangert, Mark},
            editor = {Jaffray, David A.},
            year = {2015},
            series = {{{IFMBE Proceedings}}},
            pages = {1608--1611},
            publisher = {Springer International Publishing},
            abstract = {We present matRad, an open source software for three-dimensional radiation treatment planning of intensitymodulated photon, proton, and carbon ion therapy. matRad is developed for educational and research purposes; it is entirely written in MATLAB. A first beta release is available for download. The toolkit features a highly modular design with a set of individual functions modeling the entire treatment planning workflow based on a segmented patient CT. All algorithms, e.g. for ray tracing, photon/proton/carbon dose calculation, fluence optimization, and multileaf collimator sequencing, follow well-established approaches and operate on clinically adequate voxel and bixel resolution. Patient data as well as base data for all required computations is included in matRad. We achieve computation times of 60-100s (60-400s) for realistic patient cases including photon (particle) dose calculation and fluence optimization. Memory consumption ranges between 0.2GB and 2.2GB. Dose distributions of a treatment planning study for a phantom and prostate patient case considering multiple radiation modalities are shown. Both the computational and dosimetric results encourage a future use of matRad in an educational and scientific setting.},
            isbn = {978-3-319-19387-8},
            langid = {english},
            keywords = {matRadGrantOther,open source software,particle therapy,Radiation therapy,treatment planning}
        }


----

In addition, a joint open-source paper about matRad focusing on educational purposes can be found at the following link:
`matRad - An open-source treatment planning toolkit for educational purposes <http://mpijournal.org/pdf/2018-01/MPI-2018-01-p119.pdf>`_.

.. collapse:: matRad educational publication - BibTeX

    .. code-block:: bibtex

        @article{Wieser2018b,
            title = {{{matRad}} - an open-source treatment planning toolkit for educational purposes},
            author = {Wieser, Hans-Peter and Wahl, Niklas and Gabry{\'s}, Hubert S. and M{\"u}ller, Lucas-Raphael and Pezzano, Giuseppe and Winter, Johanna and Ulrich, Silke and Burigo, Lucas Noberto and J{\"a}kel, Oliver and Bangert, Mark},
            year = {2018},
            journal = {Medical Physics International Journal},
            volume = {6},
            number = {1},
            pages = {119--127},
            issn = {2306-4609},
            langid = {english}
        }

----

.. _citesoftware:
Direct Software Citation (zenodo)
---------------------------------

You can also directly cite the software through its Publication on `Zenodo with a citable DOI <https://doi.org/10.5281/zenodo.3879615>`_.

You can either cite the generic DOI which will always resolve to the latest release of matRad: `10.5281/zenodo.3879615 <https://doi.org/10.5281/zenodo.3879615>`_:

.. image:: https://zenodo.org/badge/doi/10.5281/zenodo.3879615.svg
   :target: https://doi.org/10.5281/zenodo.3879615
   :alt: DOI

Alternatively, you can cite the specific version you are using, e.g. 3.1.0 with DOI: `10.5281/zenodo.14181851 <https://doi.org/10.5281/zenodo.14181851>`_:

.. image:: https://zenodo.org/badge/29671667.svg
    :target: https://zenodo.org/badge/latestdoi/29671667
    :alt: DOI

The respective Zenodo entry gives you an overview over all available versions to cite. Below you find example BibTeX and BibLaTex (with a few more configuration options than vanilla BibTeX) entries for the latest release 3.1.0.

.. collapse:: matRad software publication (Release 3.1.0) - BibTeX

    .. code-block:: bibtex

        @misc{matRad310,
            title = {{matRad (v3.1.0)}},
            author = {Abbani, Nelly and {Al-Hasnawi}, Nabe and Ackermann, Benjamin and Bangert, Mark and Becher, Tobias and Bennan, Amit Ben Antony and Burigo, Lucas and Cabal, Gonzalo and Cisternas, Eduardo and Charton, Louis and Christiansen, Eric and Cristoforetti, Remo and Dallas, Marios and Doerner, Edgardo and Ecker, Swantje and Ellerbrock, Malte and Facchiano, Simona and Gabry{\'s}, Hubert and Handrack, Josefine and Hardt, Jennifer and Heath, Emily and Hermann, Cindy and Homolka, Noa and Ibragim, Raed and J{\"a}ger, Fabian and J{\"a}kel, Oliver and {Hueso-Gonz{\'a}lez}, Fernando and Khaledi, Navid and Klinge, Thomas and Kunz, Jeremias and Mairani, Andrea and Meder, Paul Anton and Mescher, Henning and M{\"u}ller, Lucas-Raphael and Neishabouri, Ahmad and Palkowitsch, Martina and Parodi, Katia and Pezzano, Giuseppe and Ramirez, Daniel and Sarnighausen, Claus and Scholz, Carsten and Sevilla, Camilo and Stadler, Alexander and Ulrich, Silke and Titt, Uwe and Wahl, Niklas and Welsch, Jona and Wieser, Hans-Peter and Winter, Johanna and Xu, Tong},
            year = {2024},
            month = nov,
            address = {Heidelberg},
            doi = {10.5281/zenodo.14181851},
            howpublished = {Deutsches Krebsforschungszentrum}
        }

|

.. collapse:: matRad software publication (Release 3.1.0) - BibLaTeX

    .. code-block:: bibtex

        @software{matRad310,
            title = {{matRad}},
            author = {Abbani, Nelly and Al-Hasnawi, Nabe and Ackermann, Benjamin and Bangert, Mark and Becher, Tobias and Bennan, Amit Ben Antony and Burigo, Lucas and Cabal, Gonzalo and Cisternas, Eduardo and Charton, Louis and Christiansen, Eric and Cristoforetti, Remo and Dallas, Marios and Doerner, Edgardo and Ecker, Swantje and Ellerbrock, Malte and Facchiano, Simona and Gabryś, Hubert and Handrack, Josefine and Hardt, Jennifer and Heath, Emily and Hermann, Cindy and Homolka, Noa and Ibragim, Raed and Jäger, Fabian and Jäkel, Oliver and Hueso-González, Fernando and Khaledi, Navid and Klinge, Thomas and Kunz, Jeremias and Mairani, Andrea and Meder, Paul Anton and Mescher, Henning and Müller, Lucas-Raphael and Neishabouri, Ahmad and Palkowitsch, Martina and Parodi, Katia and Pezzano, Giuseppe and Ramirez, Daniel and Sarnighausen, Claus and Scholz, Carsten and Sevilla, Camilo and Stadler, Alexander and Ulrich, Silke and Titt, Uwe and Wahl, Niklas and Welsch, Jona and Wieser, Hans-Peter and Winter, Johanna and Xu, Tong},
            date = {2024-11},
            location = {Heidelberg},
            doi = {10.5281/zenodo.14181851},
            url = {https://doi.org/10.5281/zenodo.14181851},
            organization = {Deutsches Krebsforschungszentrum},
            version = {3.1.0}
        }

|

