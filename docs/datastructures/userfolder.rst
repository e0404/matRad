.. _userfolder:

====================
The Userfolder Setup
====================

matRad uses one or more *userfolders* to store user-provided data and generated files outside the version-controlled source tree. The folders are managed centrally by the global configuration object :class:`MatRad_Config` and are available through ``matRad_cfg.userfolders`` if ``matRad_cfg = MatRad_Config.instance();`` or ``matRad_rc`` has been called before.

The userfolder mechanism is intended for resources such as custom machine files, patients, HLUTs, scripts, and temporary working directories created by external dose engines. It complements the repository content by providing a writable location for user-specific data.

The source tree of matRad already contains a default userfolder at ``<matRadRoot>/userdata``. This folder is included in the default search path and can be used immediately for custom data. However, users are free to configure additional folders in their home directory or via environment variable.

Default behavior
----------------

The configured userfolders are collected during startup in :meth:`MatRad_Config.updateUserfolders`.

When matRad is run from source, the default search order is:

1. ``<matRadRoot>/userdata/``
2. ``<home>/matRad`` if this folder already exists
3. Any additional folders listed in the ``MATRAD_USERDATA`` environment variable

When matRad is deployed, the default is:

1. ``<home>/matRad``
2. Any additional folders listed in the ``MATRAD_USERDATA`` environment variable

The ``MATRAD_USERDATA`` variable is interpreted as a semicolon-separated list of existing directories.

.. note::

    The first entry in ``matRad_cfg.userfolders`` is exposed separately as ``matRad_cfg.primaryUserFolder``. This is the default writable location used by several matRad components.

Automatically created subfolders
--------------------------------

If a configured userfolder does not yet exist, matRad creates it on demand. It also tries to create the following standard subfolders:

* ``hluts``
* ``machines``
* ``patients``
* ``scripts``

These subfolders are used throughout matRad for user-supplied resources. For example, machine data are searched in ``machines`` and HLUT files in ``hluts``.

Search and path handling
------------------------

Each configured userfolder is added to the MATLAB path together with its subfolders. This makes helper scripts placed in the userfolder tree directly accessible from MATLAB.

If the list of userfolders is changed at runtime, ``MatRad_Config`` updates the MATLAB path accordingly:

* newly added folders are added to the path
* removed folders are removed from the path

If an empty list is assigned, matRad falls back to the default ``userdata`` folder in the matRad installation.

How matRad uses the userfolders
-------------------------------

The userfolder list is used as a search path for several kinds of data. Typical examples are:

* custom machines in ``machines``
* patient files in ``patients``
* HLUT definitions in ``hluts``
* user scripts in ``scripts``
* engine-specific working directories below ``primaryUserFolder``

Because ``primaryUserFolder`` always points to the first configured entry, changing the order of ``matRad_cfg.userfolders`` also changes where matRad writes default output for features that rely on this property.

MATLAB examples
---------------

Inspect the currently active configuration:

.. code-block:: matlab

    global matRad_cfg
    matRad_cfg = MatRad_Config.instance();

    matRad_cfg.userfolders
    matRad_cfg.primaryUserFolder

Set a custom list of userfolders at runtime:

.. code-block:: matlab

    global matRad_cfg
    matRad_cfg = MatRad_Config.instance();

    matRad_cfg.userfolders = { ...
        fullfile(matRad_cfg.matRadRoot,'userdata'), ...
        fullfile(getenv('USERPROFILE'),'matRad'), ...
        'D:\matRadData\projectA' ...
    };

Add additional userfolders via the environment before starting MATLAB:

.. code-block:: matlab

    setenv('MATRAD_USERDATA', 'D:\matRadData\projectA;D:\matRadData\shared')
    global matRad_cfg
    matRad_cfg = MatRad_Config.instance();

.. tip::

    For local scripts and experimental data, the repository-local ``userdata`` folder is usually the safest starting point because it is separated from the tracked source files.
