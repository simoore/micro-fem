import numpy as np
import microfem


# The topology of the cantilever is defined by a binary matrix. `1` indicates
# the element is solid, and `0` indicates the element is void. The cantilever
# is fixed along the left edge of the matrix. The following example has a
# large base of size 50 x 25 elements and and skinny tip of 10 x 25 elements.
topology_tip = np.vstack((np.zeros((20, 25)), 
                          np.ones((10, 25)), 
                          np.zeros((20, 25))))
topology_base = np.ones((50,25))
topology = np.hstack((topology_base, topology_tip))

# Along with the topology, the size of each element needs to be defined. `a`
# is half the element width in the x-direction (rows) and `b` is half the 
# element length in the y-direction (cols). The units are in um.
a = 5
b = 5

# One option for normalising the mode shapes is to normalize with respect to
# the deflection at a point on the cantilever. The coordinates are in um and
# should be on a solid element.  
xtip = 250
ytip = 495

# Combine the above parameters into a Cantilever object in the microfem 
# package. Information about the cantilever can be displayed and the topology
# can be plotted.
cantilever = microfem.Cantilever(topology, a, b, xtip, ytip)
cantilever.to_console()
microfem.plot_topology(cantilever)

material = microfem.SoiMumpsMaterial()
fem = microfem.PlateFEM(material, cantilever)

# With the FEM model buitl, modal analysis can be performed. This returns the
# resonance frequencies and mode shapes of the cantilever. These can then be
# plotted.
vals = fem.modal_analysis(n_modes=3)
#for i, w in enumerate(ws):
#    print('Mode %d resonance frequency is (Hz): %g' % (i, np.sqrt(w)/2/np.pi))
#microfem.plot_mode(fem, vall[:, 0])
#microfem.plot_mode(fem, vall[:, 1])
#microfem.plot_mode(fem, vall[:, 2])