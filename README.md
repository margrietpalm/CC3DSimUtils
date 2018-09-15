Python library for setting up, running and analyzing parameter sweeps with cell based models. This package was specifically designed for [CompuCell3D]() but can also be used with other tools, as shown in this [bookchapter](http://dx.doi.org/10.1007/978-1-4939-1164-6_20).


# Installation
Please see the [documentation](#Documentation) for more extensive installation instructions.

1. Install python 2.x
2. Install prerequisites
	- numpy
	- scipy
	- python imaging library
	- mahotas
	- pymorph
3. Move CC3DSimUtils to your preferred location.
4. Make sure python can find CC3DSimUtils (use a OR b).
	a. Add the path to CC3DSimUtils to the global pythonpath.
	b. Add the path to CC3DSimUtils to the pythonpath in your python code:
		import sys
		sys.path.append('PATHTOCC3DSIMUTILS')

# Documentation
HTML documentation: doc/html/CC3DSimUtils.html
PDF documentation: [doc/latex/CC3DSimUtils.pdf](https://github.com/margrietpalm/CC3DSimUtils/blob/master/doc/latex/CC3DSimUtils.pdf)