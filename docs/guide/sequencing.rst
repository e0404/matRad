.. _sequencing_dao:

===============================================
Sequencing, Direct Aperture Optimization & VMAT
===============================================

:ref:`Fluence optimization <plan_opt>` yields a bixel weight per beamlet, i.e. an intensity map per beam.
A multileaf collimator (MLC) cannot deliver such a map directly - it can only deliver a sequence of *apertures*, each irradiating the field for a certain number of monitor units.
Turning the optimized fluence into deliverable apertures, and then re-optimizing those apertures, is the subject of this section.

For photons matRad supports three delivery techniques, in increasing order of sophistication:

step-and-shoot IMRT
   The fluence maps are decomposed into apertures by a leaf sequencing algorithm and delivered as they come out. Set ``pln.propSeq.runSequencing = true``.

direct aperture optimization (DAO)
   The sequenced apertures are used as the starting point of a second optimization that treats the leaf positions themselves as variables. Set ``pln.propOpt.runDAO = true`` (this implies sequencing).

volumetric modulated arc therapy (VMAT)
   Direct aperture optimization on a continuously rotating gantry, subject to the machine's mechanical limits. Set ``pln.propOpt.runVMAT = true``.

Particle plans are not collimated; for them "sequencing" means computing the spot delivery order and timing instead, see `Spot sequencing for particles`_.

.. _leaf_sequencing:

Leaf sequencing
---------------

A leaf sequencer stratifies each fluence map into ``pln.propSeq.numLevels`` discrete intensity levels and decomposes the result into a set of deliverable apertures whose weighted sum reproduces the stratified map.
More levels reproduce the optimized fluence more faithfully at the price of more apertures, and therefore more monitor units and a longer delivery.

matRad ships three algorithms, selected with ``pln.propSeq.sequencer``:

``'siochi'`` (default)
   Rod pushing with collision and tongue-and-groove correction after `Siochi (1999) <https://www.ncbi.nlm.nih.gov/pubmed/10078655>`_, implemented in :class:`matRad_PhotonLeafSequencerSiochi`. The only sequencer that supports VMAT.

``'xia'``
   Reducing level or sliding window after `Xia et al. (1998) <http://online.medphys.org/resource/1/mphya6/v25/i8/p1424_s1>`_, implemented in :class:`matRad_PhotonLeafSequencerXia`.

``'engel'``
   Minimizes the total number of monitor units after `Engel (2005) <http://www.sciencedirect.com/science/article/pii/S0166218X05001411>`_, implemented in :class:`matRad_PhotonLeafSequencerEngel`.

Sequencers are classes derived from :class:`matRad_SequencerBase`, mirroring the design of the :ref:`dose engines <dosecalc_guide>`.
They are discovered automatically, so a new algorithm only has to derive from the photon base class :class:`matRad_PhotonLeafSequencerAbstract`, declare a ``shortName`` and implement ``sequence()`` to become selectable through ``pln.propSeq.sequencer`` and to appear in the GUI.
The top level entry point is :func:`matRad_sequencing`, which returns the resulting apertures in ``resultGUI.apertureInfo``.

.. _dao_guide:

Direct aperture optimization
----------------------------

Sequencing reproduces the optimized fluence, but the result is no longer optimal: the stratification into a finite number of levels, and the requirement that every aperture be a valid MLC shape, both cost dose quality.
Direct aperture optimization recovers part of that loss by optimizing the aperture *shapes* and *weights* directly, so that the optimization can only ever produce a plan that is deliverable by construction.

matRad's implementation largely follows `Wild et al. (2015) <https://doi.org/10.1118/1.4914863>`_, which builds on `Bzdusek et al. (2009) <http://www.ncbi.nlm.nih.gov/pubmed/19610322>`_ and, with some modification, `Unkelbach & Cassioli (2012) <http://iopscience.iop.org/article/10.1088/0031-9155/58/2/301/meta>`_.
It reuses the objective and constraint machinery of the fluence optimization: :class:`matRad_OptimizationProblemDAO` derives from :class:`matRad_OptimizationProblem` and adds the mapping from the aperture vector (leaf positions and shape weights) to bixel weights, together with its Jacobian.
Everything documented under :ref:`fluence optimization <plan_opt>` - objectives, constraints, optimizers - therefore applies unchanged.
The entry point is :func:`matRad_directApertureOptimization`.

Since the optimizer works on leaf positions rather than bixel weights, the dose influence matrix has to be recomputed from the apertures in every iteration, which makes DAO considerably more expensive than fluence optimization.
``pln.propSeq.preconditioner = true`` applies Jacobi preconditioning to the aperture weights, which usually improves convergence.

.. _vmat_guide:

VMAT
----

In VMAT the gantry rotates while the beam is on, so the plan is a continuous trajectory rather than a set of static fields.
matRad models this by three nested sets of gantry angles along each arc, generated by the ``PhotonVMAT`` steering generator (:class:`matRad_StfGeneratorPhotonVMAT`):

fluence optimization (FMO) angles
   The coarsest set. A full fluence map is optimized at each of these angles.

DAO control points
   Each FMO angle is subdivided into several control points, each of which carries exactly one aperture and is a variable of the direct aperture optimization.

dose calculation angles
   The finest set, used to sample the dose along the arc. Every DAO control point is also a dose calculation angle.

The three sets are nested, so the spacings must satisfy ``maxGantryAngleSpacing <= maxDAOGantryAngleSpacing <= maxFMOGantryAngleSpacing``.
A VMAT plan is set up like this:

.. code-block:: matlab

   pln.propStf.generator                = 'PhotonVMAT';
   pln.propStf.gantryAngles             = [-180, 180]; % arc anchor points
   pln.propStf.couchAngles              = [0, 0];
   pln.propStf.maxGantryAngleSpacing    = 8;   % [deg] dose calculation
   pln.propStf.maxDAOGantryAngleSpacing = 8;   % [deg] DAO control points
   pln.propStf.maxFMOGantryAngleSpacing = 24;  % [deg] fluence optimization
   pln.propStf.minAperturesPerFMOBeam   = 3;   % apertures sequenced per FMO beam
   pln.propOpt.runVMAT                  = true;

Several arcs can be defined by giving more than two anchor points and assigning them to arcs with ``pln.propStf.arcIndex``.

The workflow is then the usual one.
:func:`matRad_fluenceOptimization` optimizes only the FMO angle subset - which angles those are travels in the dose influence matrix as ``dij.isFMOBeam``, so no steering information is needed at that point.
:func:`matRad_sequencing` decomposes each optimized fluence map and spreads the resulting apertures over that beam's DAO control points, and :func:`matRad_directApertureOptimization` then optimizes aperture weights, leaf positions and the time spent at each control point.

Machine limits
^^^^^^^^^^^^^^

Because gantry rotation, leaf motion and dose rate all have to be physically achievable, VMAT optimization is constrained by the machine's mechanical limits, read from ``machine.constraints``.
If the machine file does not define them, matRad falls back to:

.. list-table::
   :header-rows: 1
   :widths: 30 25 45

   * - Field
     - Default ``[min max]``
     - Meaning
   * - ``gantryRotationSpeed``
     - ``[0 6]`` deg/s
     - how fast the gantry may rotate between control points
   * - ``leafSpeed``
     - ``[0 60]`` mm/s
     - how fast an individual leaf may travel
   * - ``monitorUnitRate``
     - ``[1.25 10]`` MU/s
     - dose rate the linac can hold

Sequencing options for arcs
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Every DAO control point receives exactly one aperture, so the sequencer must produce at least ``minAperturesPerFMOBeam`` apertures from each fluence map.
Two options control how it gets there:

``pln.propSeq.arcFluenceSmoothing`` (``'gaussian'`` by default, or ``'none'``)
   Blurs each fluence map before stratifying it. A fluence that is flat across its whole field decomposes into a single aperture no matter how many levels are used, so without smoothing the sequencer can fail to produce enough apertures. The price is that the sequenced fluence no longer reproduces the optimized one exactly, which is why this is applied for VMAT only.

``pln.propSeq.apertureSelection`` (``'doseAreaProduct'`` by default, or ``'leastSquares'``)
   The stratification generally yields more apertures than the arc can deliver, so the surplus is discarded. ``'doseAreaProduct'`` keeps the apertures with the largest open area and rescales them by a common factor. ``'leastSquares'`` instead keeps the apertures that best reconstruct the optimized fluence and refits their weights to it, preserving the total dose-area product; it generally starts the DAO closer to the fluence optimization result.

``pln.propSeq.continuousAperture`` (``false`` by default)
   Interpolates the leaf positions between DAO control points, modelling a dynamic delivery in which the leaves move while the gantry rotates, rather than a sequence of discrete shapes.

.. note::

   ``maxFMOGantryAngleSpacing`` deserves more attention than it looks.
   Each fluence map is optimized at one gantry angle but delivered as ``minAperturesPerFMOBeam`` separate apertures spread over its whole FMO sector.
   The sum reproduces the optimized fluence, yet every individual control point only receives a fragment of it, so the wider the sector, the more the delivered dose clusters into a few entrance directions - and direct aperture optimization tends to amplify that rather than repair it.
   Measured on ``matRad_example22_photonsVMAT`` as the peak-to-mean of the out-of-target dose binned by gantry angle, sequencing adds about 2 % angular non-uniformity at 24 deg, 30 % at 45 deg and 38 % at 72 deg.
   Target coverage and maximum dose look healthy in all three cases, so this is only visible in the spatial dose distribution.

Checking the result
^^^^^^^^^^^^^^^^^^^

:func:`matRad_calcDeliveryMetrics` collects the total monitor units and delivery time of the plan, plus the MU rate, gantry speed and leaf speeds per control point, in ``resultGUI.deliveryMetrics``.

:func:`matRad_visApertureInfo` visualizes the apertures. Besides the per-beam and grid views it offers a ``'trajectory'`` view of leaf position against gantry angle, a ``'metrics'`` view plotting the delivery quantities against their machine limits, and an ``'animate'`` view that plays the arc back in delivery time.

The plan is optimized at the DAO control points, which are coarser than what is actually delivered.
:func:`matRad_refineApertureArc` resamples the optimized apertures onto a finer gantry angle grid so that the delivered dose can be forward-calculated and compared against the optimized one.

A complete, commented walkthrough is given in ``matRad_example22_photonsVMAT``.

.. _particle_sequencing:

Spot sequencing for particles
-----------------------------

Scanned particle delivery needs no collimator, so nothing has to be decomposed into apertures.
For particles the sequencer instead computes the order in which the spots are delivered and how long each takes, using a synchrotron spill and lateral scanning model, and leaves the fluence untouched.
This is implemented in :class:`matRad_ParticleScanningSequencerSpill` and is what :ref:`4D dose calculation <4d>` uses to map spots onto breathing phases.
