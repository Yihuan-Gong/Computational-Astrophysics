import yt
import cv2
import matplotlib.pyplot as plt
import numpy as np
from multiprocessing import Pool

'''
Add necessary quantities
'''
# 
# Freqently used physical quantities
# 

def electron_density(field, data):
    return data[('flash', 'dens')]/(yt.units.atomic_mass_unit_cgs*1.18)
yt.add_field(
    ('gas', 'electron_density'), 
    function=electron_density, 
    sampling_type='cell',
    units='cm**(-3)',
    force_override=True
)

def temp_in_keV(field, data): 
    return yt.units.boltzmann_constant_cgs * data[('flash', 'temp')]
yt.add_field(
    ("gas", "temp_in_keV"), 
    function = temp_in_keV, 
    sampling_type='cell', 
    units='keV',
    force_override=True
)

def entropy(field, data):
    return data[("gas", "temp_in_keV")] * data[('gas', 'electron_density')]**(-2/3)
yt.add_field(
    ("gas", "entropy"), 
    function = entropy, 
    sampling_type='cell', 
    units='keV*cm**2',
    force_override=True
)

def potential(field, data):
    return data[('flash', 'gpot')] * (-1)
yt.add_field(
    ("gas", "potential"), 
    function = potential, 
    sampling_type='cell', 
    units='erg/g',
    force_override=True
)


# 
# Modified "code unit" to real unit
# 

def density(field, data):
    return data[('flash', 'dens')]
yt.add_field(
    ('gas', 'density'), 
    function=density, 
    sampling_type='cell',
    units='g*cm**(-3)',
    force_override=True
)


def pressure(field, data):
    return data[('flash', 'pres')]
yt.add_field(
    ('gas', 'pressure'), 
    function=pressure, 
    sampling_type='cell',
    units='Ba',
    force_override=True
)

'''
Plot temperature, pressure, density and entropy
'''
def ploter(i, L=8):
    ds = yt.load('/data/yhgong/galaxy_cluster_merger/mass_scalar/perseus_merger_hdf5_plt_cnt_%04d'%i)
    
    # # Temperature
    # plot = yt.SlicePlot(ds, "z", ("gas", "temp_in_keV"), width=(L, 'Mpc'))
    # plot.set_zlim(("gas", "temp_in_keV"), 0.8, 15)
    # plot.annotate_text([0.05,0.95],'time = %.2f Gyr' %ds.current_time.in_units('Gyr'),coord_system='axis',text_args={'color':'black'})
    # plot.set_cmap(('gas', 'temp_in_keV'), 'Blue-Red')
    # plot.save('picture/temp/%1dMpc/%04d.png'%(L, i))

    # # Density
    # plot = yt.SlicePlot(ds, 'z', [('gas', 'density')],width = (L, 'Mpc'))
    # plot.set_cmap([('gas', 'density')], 'Blue-Red')
    # plot.set_zlim(('gas', 'density'), 1e-30, 1e-25)
    # plot.annotate_text([0.05,0.95],'time = %.2f Gyr' %ds.current_time.in_units('Gyr'),coord_system='axis',text_args={'color':'black'})
    # plot.save('picture/dens/%1dMpc/%04d.png'%(L, i))

    # # Pressure
    # plot = yt.SlicePlot(ds, 'z', [('gas', 'pressure')],width = (L, 'Mpc'))
    # plot.set_cmap([('gas', 'pressure')], 'Blue-Red')
    # plot.set_zlim(('gas', 'pressure'), 1e-14, 1e-9)
    # plot.annotate_text([0.05,0.95],'time = %.2f Gyr' %ds.current_time.in_units('Gyr'),coord_system='axis',text_args={'color':'black'})
    # plot.save('picture/pres/%1dMpc/%04d.png'%(L, i))

    # # Entropy
    # plot = yt.SlicePlot(ds, 'z', [('gas', 'entropy')],width = (L, 'Mpc'))
    # plot.set_cmap([('gas', 'entropy')], 'Blue-Red')
    # plot.set_zlim(('gas', 'entropy'), 1e1, 1e4)
    # plot.annotate_text([0.05,0.95],'time = %.2f Gyr' %ds.current_time.in_units('Gyr'),coord_system='axis',text_args={'color':'black'})
    # plot.save('picture/enpy/%1dMpc/%04d.png'%(L, i))

    # Mass scalar
    plot = yt.SlicePlot(ds, 'z', [('flash','sub2')],width = (8, 'Mpc'))
    plot.set_cmap([('flash', 'sub2')], 'Blue-Red')
    plot.set_zlim(('flash', 'sub2'), 0, 1)
    plot.set_log(('flash', 'sub2'), False)
    plot.annotate_text([0.05,0.95],'time = %.2f Gyr' %ds.current_time.in_units('Gyr'),coord_system='axis',text_args={'color':'black'})
    plot.annotate_sphere([0,0,0], radius=(300,  "kpc"), circle_args={"color": "black"})
    plot.annotate_sphere([0,0,0], radius=(2500, "kpc"), circle_args={"color": "black"})
    plot.save('picture/sub2/%1dMpc/%04d.png'%(L, i))

# Multiprocessing in 8 cores
pool = Pool(16)
output = pool.map(ploter, range(180,365), chunksize=10)






'''
Make it a video
'''
# img = cv2.imread('picture/temp/0180.png')
# fps = 12
# size = (img.shape[1], img.shape[0])  # (width, height)
# fourcc = cv2.VideoWriter_fourcc(*'mp4v')
# video = cv2.VideoWriter('video/temp_3Mpc.mp4', fourcc, fps, size)

# for j in range(180,365):
#     img = cv2.imread('picture/temp/%04d.png'%j)
#     video.write(img)

# cv2.destroyAllWindows()
# video.release()

# plot.annotate_contour(('gas','potential'), ncont=20, plot_args={"colors": "gray", "linewidths": 0.5})