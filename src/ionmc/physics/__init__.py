"""Physics models.

Modules in this package are *shared-source*: they are written against the
math namespace of :mod:`ionmc.backend.mathlib` and can be executed as Warp
functions (CPU and CUDA), as pure Python in float64, or vectorised with numpy.
See ``decisions/0005`` for the execution model.
"""
