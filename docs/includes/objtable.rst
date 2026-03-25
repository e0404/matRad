.. list-table::
    :header-rows: 1

    * - Objective
      - Class
      - Effect on the objective function
    * - **square underdosing**
      - :class:`DoseObjectives.matRad_SquaredUnderdosing`
      - Only dose values lower than the threshold dose for this VOI are considered for the objective function. The deviations are squared, multiplied with the penalty factor and added to the objective function value. The penalty is normalized to the number of voxels per VOI.
    * - **square overdosing** 
      - :class:`DoseObjectives.matRad_SquaredOverdosing`
      - Only dose values larger than the threshold dose for this VOI are considered for the objective function. The deviations are squared, multiplied with the penalty factor and added to the objective function value. The penalty is normalized to the number of voxels per VOI.
    * - **square deviation** 
      - :class:`DoseObjectives.matRad_SquaredDeviation`
      - All deviation from a reference dose for this VOI are considered for the objective function. All deviations are squared, multiplied with the penalty factor and added to the objective function value. The penalty is normalized to the number of voxels per VOI.
    * - **mean** 
      - :class:`DoseObjectives.matRad_MeanDose`
      - All dose values inside this VOI are weighted with the specified, and normalized (see above), penalty factor and added to the objective function value.
    * - **EUD** 
      - :class:`DoseObjectives.matRad_EUD`
      - EUD is the abbreviation for equivalent uniform dose. For this method a weighting factor and an exponent *a* have to be defined. For the calculation of the objective function value the dose in each voxel is taken to the power of *a*. Then the sum of all these values is taken (Σ D\ :sub:`i`\ :sup:`a`) and divided by the number of voxels. The *a*-th root of this value is then multiplied with the weighting factor and added to the objective function value.
    * - **Min/max DVH objective**  
      - :class:`DoseObjectives.matRad_MaxDVHObjective`
        :class:`DoseObjectives.matRad_MinDVHObjective`
      - Only deviations from the reference dose over/under a reference volume are considered according to `Wu & Mohan (2000 Medical Physics) <http://scitation.aip.org/content/aapm/journal/medphys/27/4/10.1118/1.598932>`_.