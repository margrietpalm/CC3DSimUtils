=================
AnalysisUtils
=================
Set of function that do a single analysis.

------------------------------
Compactness
------------------------------
.. autofunction:: AnalysisUtils.getCompactness
.. autofunction:: AnalysisUtils.getLCC

------------------------------
Order parameter
------------------------------
The `order parameter <http://en.wikipedia.org/wiki/Liquid_crystal#Order_parameter>`_ describes the orientational order of a liquid crystal : :math:`s = \left \langle \frac{ 3\text{cos}^2(\theta)}{2}\quad\right \rangle`; :math:`s = 0\;` for a random sample and :math:`s = 1\;` for an isotropic sample. :math:`\theta\;` is the angle between the cell direction and the `director <http://en.wikipedia.org/wiki/Liquid_crystal#Director>`_. The director is a dimensionless unit vector. It represents the direction of preferred orientation of cells in the neighborhood of any point. Because there is no physical polarity along the director axis, n and -n are fully equivalent. Here the neighborhood is defined as a circle with center com and radius r. 

.. autofunction:: AnalysisUtils.getDirector
.. autofunction:: AnalysisUtils.getOrderParameterFromGrid
.. autofunction:: AnalysisUtils.getLocalOrderParameterFromGrid
.. autofunction:: AnalysisUtils.getGlobalOrderParameterFromGrid
.. autofunction:: AnalysisUtils.getRelativeDirField

------------------------------
Cell clusters
------------------------------
.. autofunction:: AnalysisUtils.getCellClusters
.. autoclass:: AnalysisUtils.Cluster
	:members:
.. autoclass:: AnalysisUtils.ClusterCellTC
	:members:
	
------------------------------
Cell angles
------------------------------
The angle of a cell is calculated from the inertia tensor of a cell. From the intertia tensor we calucate the eigenvalues and eigenvectors; the eigenvector that corresponds with the largest eigenvalue represents the direction of the long axis of a cell. The angle between the long axis and the x-axis is the cell angle.

.. autofunction:: AnalysisUtils.getCellAngle
.. autofunction:: AnalysisUtils.getAngleField

------------------------------
Mean squared displacement
------------------------------
The `mean squared displacement <http://en.wikipedia.org/wiki/Mean_squared_displacement>`_ describes the displacement of a cell over time with respect to the initial position : :math:`MSD = \left \langle (x(t) - x(0))^2 \right \rangle\;`. In a similar manner the mean squared angular displacement can be calculated : :math:`MSD = \left \langle (\theta(t) - \theta(0))^2 \right \rangle\;`

.. autofunction:: AnalysisUtils.calcMSDTransForCellTC
.. autofunction:: AnalysisUtils.calcMSDRotForCellTC


