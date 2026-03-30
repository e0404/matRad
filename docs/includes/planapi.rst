.. list-table:: Properties and their corresponding API functions
    :header-rows: 1
    :widths: 20 20 20 20 20

    * - Plan property
      - API function
      - Description
      - ID
      - Folder
    * - ``propStf``
      - :func:`matRad_generateStf`
      - Create beam Geometry
      - generator
      - :mod:`matRad.steering`
    * - ``propDoseCalc``
      - :func:`matRad_calcDoseInfluence`
        :func:`matRad_calcDoseForward`
      - Calculate dose matrix / distribution
      - engine
      - :mod:`matRad.doseCalc`
    * - ``propOpt``
      - :func:`matRad_fluenceOptimization`
      - Optimization of beam fluences
      - problem
      - :mod:`matRad.optimization`
    * - ``propSeq``
      - :func:`matRad_sequencing`
      - Sequencing of beams
      - sequencer
      - :mod:`matRad.sequencing`