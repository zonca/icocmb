icocmb: Project HEALPix maps to a Truncated Icosahedron
=======================================================

Requirements
------------

* `matplotlib` (tested only 1.3.0)
* `healpy`, needs this patch to be applied to `projaxes.py`:

    diff --git a/healpy/projaxes.py b/healpy/projaxes.py
    index 354fd99..9818de9 100644
    --- a/healpy/projaxes.py
    +++ b/healpy/projaxes.py
    @@ -317,6 +317,8 @@ class SphericalProjAxes(axes.Axes):
             vec = R.dir2vec(theta,phi,lonlat=lonlat)
             vec = (R.Rotator(rot=rot,coord=coord,eulertype='Y')).I(vec)
             x,y = self.proj.vec2xy(vec,direct=kwds.pop('direct',False))
    +        if kwds.pop("return_projected_points", False):
    +            return x, y
             s = self.scatter(x, y, *args, **kwds)
             if save_input_data:
                 if not hasattr(self, '_scatter_data'):

or get the `icocmb` branch from the `healpy` repository:

https://github.com/healpy/healpy/tree/icocmb

How to run
----------

It is important always to use the non-interactive `matplotlib` backend,
so better not run inside `ipython --pylab`.

Steps:

* the script first creates a grid of all the hexagons and pentagons
* reads the `HEALPix` Planck map (you can modify to read other maps)
* makes gnomviews of all the pentagons and rescales them to fit in the figure
* makes gnomviews of all the hexagons and rescales them to fit in the figure
* saves as pdf

It is a quite hacky script (sorry), I could not find a way to compute automatically
the scale factor, so I estimated it by eye, it should be good enough even
for printing at high resolution, but it might be different on other
platforms/matplotlib versions.

Output
------

All Planck frequency maps created with this script:

http://figshare.com/articles/Planck_maps_projected_on_a_Truncated_Icosahedron/799823
