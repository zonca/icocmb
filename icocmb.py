import matplotlib
### DO NOT USE AN INTERACTIVE BACKEND, BECAUSE IT RESCALES EVERYTHING
matplotlib.use("agg")
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import sys

import truncated_icosahedron

#################################
## Scaling factors, this might depend on matplotlib/machine configuration
hexa_dimfac = 1.89
penta_dimfac = 2.00 

TI = truncated_icosahedron.TruncatedIcosahedron()

print 50*"*" + "Create the canvas"

from matplotlib.patches import RegularPolygon

polygons = [[]] * len(TI.faces) 
centers  = [[]] * len(TI.faces) 

fig = plt.figure(figsize=(25, 25))

ax = fig.add_subplot(111, zorder=100)
#ax = fig.add_subplot(111)
ax.set_xlim([0, 1])
ax.set_ylim([0, 1])
ax.set_autoscale_on(False)
print ax.axis

edge = .05
penta_radius = edge / 2. / np.sin(np.pi / 5)
penta_apo    = penta_radius * np.cos(np.pi / 5)
hexa_radius = edge / 2. / np.sin(np.pi / 6)
hexa_apo    = hexa_radius * np.cos(np.pi / 6)

margin = .1

def create_pentagon(center, top=True):
    orientation = 0 if top else np.pi
    return RegularPolygon(center, 5, radius = penta_radius, orientation=orientation, fc="none", ec="black")

def create_hexagon(center):
    return RegularPolygon(center, 6, orientation=np.pi/6., radius = hexa_radius, fc="none", ec="black")

# north pole pentagon
centers[2]  = np.array([penta_radius+margin, 1-penta_radius-margin])
polygons[2] = create_pentagon(centers[2]) 

ref_edge = TI.faces[2][:2]
ii = TI.find_adjacent(2, ref_edge) 
centers[ii]  = centers[2] + np.array([0, - penta_apo - hexa_apo])
polygons[ii] = create_hexagon(centers[ii])

i2 = TI.find_face_south(ii) 
centers[i2]  = centers[ii] + np.array([0, -2*hexa_apo])
polygons[i2] = create_hexagon(centers[i2])

i3 = TI.find_face_south(i2) 
centers[i3]  = centers[i2] + np.array([0, -hexa_apo-penta_apo])
polygons[i3] = create_pentagon(centers[i3], False)

prev_face = i2

for cent_col in range(9):
#for cent_col in range(1):
    print "COLUMN", cent_col
    #cent_col = 0
    if cent_col % 2 == 0:
        central_face = TI.find_face_se(prev_face)
        centers[central_face] = centers[prev_face] + np.array([hexa_radius+edge/2., -hexa_apo])
        polygons[central_face] = create_hexagon(centers[central_face])
        bot_hexa = TI.find_face_south(central_face)
        centers[bot_hexa] = centers[central_face]  + np.array([0, -2*hexa_apo])
        polygons[bot_hexa] = create_hexagon(centers[bot_hexa])
        top_penta = TI.find_face_north(central_face)
        centers[top_penta] = centers[central_face]  + np.array([0, hexa_apo+penta_apo])
        polygons[top_penta] = create_pentagon(centers[top_penta])
        if cent_col == 8:
            top_penta = TI.find_face_south(bot_hexa)
            centers[top_penta] = centers[bot_hexa]  - np.array([0, hexa_apo+penta_apo])
            polygons[top_penta] = create_pentagon(centers[top_penta], top=False)
    else:
        central_face = TI.find_face_ne(prev_face)
        centers[central_face] = centers[prev_face] + np.array([hexa_radius+edge/2., +hexa_apo])
        polygons[central_face] = create_hexagon(centers[central_face])
        top_penta = TI.find_face_south(central_face)
        centers[top_penta] = centers[central_face]  - np.array([0, hexa_apo+penta_apo])
        polygons[top_penta] = create_pentagon(centers[top_penta], top=False)
        bot_hexa = TI.find_face_north(central_face)
        centers[bot_hexa] = centers[central_face]  + np.array([0, +2*hexa_apo])
        polygons[bot_hexa] = create_hexagon(centers[bot_hexa])
    prev_face = central_face

for ii in range(len(TI.faces)):
    if polygons[ii]:
        ax.add_patch(polygons[ii])
        #plt.text(centers[ii][0], centers[ii][1], str(ii))

plt.axis("off")
#plt.show()
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from glob import glob

def GetMapName(freq):
    indir = '/data/NEVERCOMMIT/'
    f = np.int(freq)
    filename = glob('/data/NEVERCOMMIT/*_%d_*fits' % f)[0]
    return filename

############### Universal colormap
# setup linear colormap
from matplotlib.colors import ListedColormap
planck_freqmap_cmap = ListedColormap(np.loadtxt("Planck_FreqMap_RGB.txt")/255.)
planck_freqmap_cmap.set_bad("gray") # color of missing pixels
planck_freqmap_cmap.set_under("white") # color of background, necessary if you want to use

# Graticule spacing in degrees
grat_spacing = 10.0
imgsize = 750
imgresolution = 5.

faces_to_do = np.arange(32)
freq = sys.argv[1]
# load a map
mapfname = GetMapName(freq)

I = hp.ma(hp.read_map(mapfname,field=0,verbose=True))
if np.int(freq)<500:
    I = I * 1e6
else:
    I = I * 1e3
D = np.log10(0.5*(I+np.sqrt(I**2+4)))

print 50*"*" + "CMB pentagons" + freq

for f in np.arange(32):
 if len(TI.faces[f]) == 5:
    print f

    c = fig.transFigure.inverted().transform(ax.transData.transform(centers[f]))
    #c = ax.transAxes.inverted().transform(ax.transData.transform(centers[f]))

    center_offset = penta_dimfac*edge/2
    local_ax = plt.axes([c[0]-center_offset, c[1]-center_offset, penta_dimfac*edge, penta_dimfac*edge], frameon=False)

    imglim = np.radians(imgresolution * imgsize /60)/2.
    local_ax.axis("off")
    #local_ax.set_xlim((-imglim, imglim))
    #local_ax.set_ylim((-imglim, imglim))

    theta,phi = hp.vec2ang(TI.centroids[f])
    face_theta = theta[0]
    face_phi = phi[0]
    rotang = [phi[0]*180.0/np.pi,90.0-theta[0]*180.0/np.pi,0.0]
    print "Face", f, "Face centroid: " , np.degrees([face_theta, face_phi])

    planck_freqmap_cmap.set_bad("none")
    gnomD = hp.gnomview(D,cmap=planck_freqmap_cmap,rot=rotang,min=-3,max=7,title='%s GHz face %d'%(freq,f),reso=imgresolution,xsize=imgsize,ysize=imgsize,coord='G',notext=True,cbar=False, return_projected_map=True)
    hp.graticule(dpar=grat_spacing,dmer=grat_spacing,local=False,force=True)

    gnomD.mask = np.zeros(gnomD.shape, dtype=np.bool)
    #gnomD.mask[300:400] = 1
    cut_y = 110
    gnomD.mask[:cut_y, :] = 1
    gnomD.mask[-cut_y:, :] = 1
    cut_x = 67
    gnomD.mask[:, :cut_x] = 1
    gnomD.mask[:, -cut_x:] = 1
    gnom_ax = plt.gca()
    thetas = []
    phis = []

    for v in TI.v[TI.faces[f]]:
        theta, phi = hp.vec2ang(v/np.linalg.norm(v))
        thetas.append(theta[0]); phis.append(phi[0])
        vec = (hp.rotator.Rotator(rot=rotang,eulertype='Y')).I(v/np.linalg.norm(v))

    x, y = gnom_ax.projscatter(thetas,phis,color='k',s=10, return_projected_points=True)
    plt.close(plt.gcf())
    #local_ax.scatter(x, y, color="k", s=2)

    local_v = np.hstack([x[:,None], y[:,None]])

    # cut
    # bottom
    lower_side = np.min(local_v[:,1])
    upper_side = np.max(local_v[:,1])
    #gnomD.mask[:-np.floor((imglim-lower_side)*(imgsize/(2*imglim))), :] = 1
    # bottom right
    right = np.max(local_v[:,0])

    local_ax.imshow(gnomD[:,::-1], extent=[-imglim, imglim, -imglim, imglim],
            cmap=planck_freqmap_cmap, vmin=-3, vmax=7, alpha=1, origin='lower', interpolation='nearest')

print 50*"*" + "CMB hexagons" + freq

f=24
for f in np.arange(32):
 if len(TI.faces[f]) == 6:
    print f

    c = fig.transFigure.inverted().transform(ax.transData.transform(centers[f]))
    #c = ax.transAxes.inverted().transform(ax.transData.transform(centers[f]))
    local_ax = plt.axes([c[0]- hexa_dimfac*edge/2, c[1]-hexa_dimfac*edge/2, hexa_dimfac*edge, hexa_dimfac*edge], frameon=False)

    imglim = np.radians(imgresolution * imgsize /60)/2.
    local_ax.axis("off")
    local_ax.set_xlim((-imglim, imglim))
    local_ax.set_ylim((-imglim, imglim))
    local_ax.set_autoscale_on(False)

    theta,phi = hp.vec2ang(TI.centroids[f])
    face_theta = theta[0]
    face_phi = phi[0]
    rotang = [phi[0]*180.0/np.pi,90.0-theta[0]*180.0/np.pi,0.0]
    print "Face", f, "Face centroid: " , np.degrees([face_theta, face_phi])

    planck_freqmap_cmap.set_bad("none")
    gnomD = hp.gnomview(D,cmap=planck_freqmap_cmap,rot=rotang,min=-3,max=7,title='%s GHz face %d'%(freq,f),reso=imgresolution,xsize=imgsize,ysize=imgsize,coord='G',notext=True,cbar=False, return_projected_map=True)
    hp.graticule(dpar=grat_spacing,dmer=grat_spacing,local=False,force=True)

    gnomD.mask = np.zeros(gnomD.shape, dtype=np.bool)
    cut_y = 110
    gnomD.mask[:cut_y, :] = 1
    gnomD.mask[-cut_y:, :] = 1
    cut_x = 67
    gnomD.mask[:, :cut_x] = 1
    gnomD.mask[:, -cut_x:] = 1
    grid_x = np.arange(len(gnomD))
    grid_y = np.arange(len(gnomD))[:, None]
    line = np.tan(np.radians(60)) * grid_x + len(gnomD)/2 - cut_x * np.tan(np.radians(60))
    matline = line < grid_y
    line = -np.tan(np.radians(60)) * grid_x + len(gnomD)/2 + cut_x * np.tan(np.radians(60))
    matline |= line > grid_y
    line = np.tan(np.radians(60)) * grid_x - (len(gnomD) - cut_x) * np.tan(np.radians(60)) + len(gnomD)/2.
    matline |= line > grid_y
    line = -np.tan(np.radians(60)) * grid_x + (len(gnomD) - cut_x) * np.tan(np.radians(60)) + len(gnomD)/2.
    matline |= line < grid_y
    gnom_ax = plt.gca()
    gnomD.mask |= matline
    thetas = []
    phis = []

    for v in TI.v[TI.faces[f]]:
        theta, phi = hp.vec2ang(v/np.linalg.norm(v))
        thetas.append(theta[0]); phis.append(phi[0])
        vec = (hp.rotator.Rotator(rot=rotang,eulertype='Y')).I(v/np.linalg.norm(v))

    x, y = gnom_ax.projscatter(thetas,phis,color='k',s=10, return_projected_points=True)
    plt.close(plt.gcf())
    #local_ax.scatter(x, y, color="k", s=2)

    local_v = np.hstack([x[:,None], y[:,None]])

    # cut
    # bottom
    lower_side = np.min(local_v[:,1])
    upper_side = np.max(local_v[:,1])
    # bottom right
    right = np.max(local_v[:,0])

    local_ax.imshow(gnomD[:, ::-1], extent=[-imglim, imglim, -imglim, imglim],
            cmap=planck_freqmap_cmap, vmin=-3, vmax=7, alpha=1, origin='lower', interpolation='nearest')

print 50*"*" + "Save" + freq
plt.savefig("icocmb_planck_%s.pdf" % freq, dpi=300)
