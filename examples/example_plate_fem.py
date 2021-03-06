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

# The material properties for the SOI Mumps manufacturing process has been 
# stored in an object. The plate finite element model requires both the 
# cantilever object and material object to initialize. 
material = microfem.SoiMumpsMaterial()
fem = microfem.PlateFEM(material, cantilever)

# With the FEM model built, modal analysis can be performed. This returns the
# resonance frequencies and mode shapes of the cantilever. These can then be
# plotted.
n_modes = 3
ws, _, vall = fem.modal_analysis(n_modes=n_modes)

# Using the results of the modal analysis, the frequency, stiffness, 
# tip displacement, and mode type (flexural, torsional), can be determined.
coords = (cantilever.xtip, cantilever.ytip)
opr = microfem.PlateDisplacement(fem, coords).get_operator()
mode_ident = microfem.ModeIdentification(fem, cantilever, type_='plate')

freq = np.sqrt(ws) / (2*np.pi)
kuu = fem.get_stiffness_matrix(free=False)
phis = [vall[:, [i]] for i in range(n_modes)]
wtips = [opr @ p for p in phis]
kfunc = lambda p, w: np.asscalar(p.T @ kuu @ p / w ** 2)
ks = [kfunc(p, w) for p, w in zip(phis, wtips)]
types = [mode_ident.is_mode_flexural(p) for p in phis]
            
tup = ('Disp', 'Freq (Hz)', 'Stiffness', 'Flexural')
print('\n    %-15s %-15s %-15s %-10s' % tup)
for i in range(n_modes):
    tup = (i+1, wtips[i], freq[i], ks[i], str(types[i]))
    print('%-2d: %-15g %-15g %-15g %-10s' % tup)
    
microfem.plot_mode(fem, vall[:, 0])
microfem.plot_mode(fem, vall[:, 1])
microfem.plot_mode(fem, vall[:, 2])