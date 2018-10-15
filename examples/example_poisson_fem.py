import numpy as np
import microfem


def define_parameters(domain):

    k = np.zeros(domain.shape)
    q = np.zeros(domain.shape)
    
    for i, j in np.ndindex(domain.shape):
        k[i, j] = 1e-2 if (3 < j < 16) or (j > 15 and i > 5) else 1
        q[i, j] = 1e-3 if k[i, j] == 1 else 0
        
    return k, q
    
    
domain = np.ones((20, 25))
domain[15:20, 20:25] = 0
conductivity, source = define_parameters(domain)
a, b = 5, 5

poisson_domain = microfem.PoissonDomain(domain, conductivity, source, a, b)
fem = microfem.PoissonFEM(poisson_domain)
uall, _ = fem.solve()
microfem.plot_poisson_solution(fem, uall, show='domain')
microfem.plot_poisson_solution(fem, uall, show='conductivity')
microfem.plot_poisson_solution(fem, uall, show='source')

