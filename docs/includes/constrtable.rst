.. list-table::
    :header-rows: 1

    * - Constraint
      - Class
      - Description
    * - **Min/Max Dose**
      - :class:`DoseConstraints.matRad_MinMaxDose`
      - Keeps dose above and below the set minimum and maximum dose. Can use a LogSumExp Approximation or a voxel-wise constraint.
    * - **Min/Max Mean Dose**
      - :class:`DoseConstraints.matRad_MinMaxMeanDose`
      - Keeps the mean dose above and below the set minimum and maximum mean dose.
    * - **Min/Max EUD** 
      - :class:`DoseConstraints.matRad_MinMaxEUD`
      - Keeps the EUD above and below the set minimum and maximum EUD. The EUD is calculated using the same method as for the EUD objective.
    * - **min/max DVH** 
      - :class:`DoseConstraints.matRad_MinMaxDVH`
      - Keeps the dose volume histogram above and below the set minimum and maximum *volunme* for a given dose level.