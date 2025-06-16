.. include:: ../includes/logo.rst

.. _get_code:
########################################
Get a local copy of |matRad_logo_header|
########################################

To get a local copy of `matRad <https://github.com/e0404/matRad>`_ you have two options:

.. contents::
    :depth: 1
    :local:

Cloning the |matRad_logo_header2| Repository
--------------------------------------------

We encourage everybody working with the code to get familiar with ``git``. Using the git workflow does not allow more sophisticated code development and collaboration with the |matRad_logo| development team, but also facilitates keeping a history of changes and integrating updates regularly. If you have never worked with a git repository before, you might want to have a look at the `github guides <https://guides.github.com/>`_ first and rely on a software with a graphical user interface (such as GitHub Desktop) to manage your local |matRad_logo| copy.

If you follow the instructions on cloning the repository, you will not only get a local copy of the |matRad_logo| source ode, but get a copy (a "clone") of the entire repository, including all its history. This is part of the decentralized philosophy of git, where each copy is its own repository with which you can work independently, and which you can later synchronize with the original repository (or others). The originial remote repository, in this case `<https://github.com/e0404/matRad.git>`_ (HTTPS) or `<git@github.com:e0404/matRad.git>`_ (SSH), is usually refered to as "origin" or "upstream" repository. git provides the necessary syncing capabilities as "pulling" and "pushing" from/to a remote repository. 

If you plan to contribute and/or publish your changes on GitHub, follow below instructions to create a fork. If you only plan to pull some updates from the main repository, you can skip the forking step and directly clone the repository. Nothing will be lost then: If you decide later to contribute, you can always delay creation of a fork to contribute to the main repository.


Forking |matRad_logo_header3|
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A fork serves as a second personal online copy of the repository within your GitHub Account. Forks are mainly used to propose changes to the original repository, or to use it as a starting point for a derived project. Forking a repository allows you to freely experiment with changes without affecting the original project, and it is a common practice in open source software development.

.. note:: 

    GitHub Forks are usually the only way to contribute to a GitHub project within normal git workflow practices, as most upstream repositories (including the |matRad_logo| repository) will not allow you to push changes directly to the main repository, and there are no other mechanisms to propose changes from other remote repositories. 

Forks can be created directly from the repository main page on GitHub. To create a fork of the |matRad_logo| repository, follow these steps: 
Afterwards your fork of |matRad_logo| will be available under ``https://github.com/YOURGITHUBUSERNAME/matRad``.

Cloning |matRad_logo_header3|
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can `clone <https://help.github.com/articles/cloning-a-repository/>`_ either the main repository or your fork of matRad. With the command line, you can clone the repository with::
    
    git clone --recurse-submodules https://github.com/e0404/matRad.git matRad

This will clone the main repository in the folder ``matRad``. If you want to clone your fork, replace ``e0404`` in the URL with your GitHub Username. The ``--recurse-submodules`` option is not necessary, but will also clone other dependencies of |matRad_logo| (like our unit testing framework).

The repository you clone with the command above will be referred to as "origin" by git by default. You can add further remotes (like your fork) with::
    
    git remote add upstream https://github.com/myaccount/matRad.git


Downloading the matRad.zip folder
---------------------------------

To get a copy of |matRad_logo| you can either download the *.zip folder from the `matRad repository on Github <https://github.com/e0404/matRad>`_ or use `this link <https://github.com/e0404/matRad/archive/master.zip>`_.

Once you downloaded the *matRad.zip* folder you can unzip the files to your desired location.

