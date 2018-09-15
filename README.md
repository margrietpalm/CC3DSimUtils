Python library for setting up, running and analyzing parameter sweeps with cell based models. This package was specifically designed for [CompuCell3D](http://compucell3d.org) but can also be used with other tools, as shown in this [bookchapter](http://dx.doi.org/10.1007/978-1-4939-1164-6_20).


# Installation
Please see the [documentation](#documentation) for more extensive installation instructions.

1. Install python 2.x
2. Install prerequisites
	* [numpy](http://www.numpy.org/)
	* [scipy](https://www.scipy.org/)
	* [python imaging library (PIL)](http://www.pythonware.com/products/pil/)
		* The PIL fork [Pillow](https://pillow.readthedocs.io/en/5.2.x/) may work as well, but this is not tested.
	* [mahotas](https://mahotas.readthedocs.io/en/latest/)
	* [pymorph](https://mahotas.readthedocs.io/en/latest/)
3. Move CC3DSimUtils to your preferred location.
4. Make sure python can find CC3DSimUtils (use a OR b).
	a. Add the path to CC3DSimUtils to the global pythonpath.
	b. Add the path to CC3DSimUtils to the pythonpath in your python code:
		```
		import sys
		sys.path.append('PATHTOCC3DSIMUTILS')
		```
		
# Documentation
* HTML documentation: doc/html/CC3DSimUtils.html
* PDF documentation: [doc/latex/CC3DSimUtils.pdf](https://github.com/margrietpalm/CC3DSimUtils/blob/master/doc/latex/CC3DSimUtils.pdf)
