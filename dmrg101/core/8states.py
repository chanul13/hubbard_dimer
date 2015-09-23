# 
# File: sites.py
# Author: Ivan Gonzalez
#
""" A module for single sites.
"""
import numpy as np 

from dmrg_exceptions import DMRGException
from numpy import sqrt

class Site(object):
    """A general single site
    
    You use this class to create a single site. The site comes empty (i.e.
    with no operators included), but for th identity operator. You should
    add operators you need to make you site up.
    
    Parameters
    ----------
    dim : an int 
	Size of the Hilbert space. The dimension must be at least 1. A site of
        dim = 1  represents the vaccum (or something strange like that, it's
        used for demo purposes mostly.)
    operators : a dictionary of string and numpy array (with ndim = 2).
	Operators for the site.

    Examples
    --------
    >>> from dmrg_solution.core.sites import Site
    >>> brand_new_site = Site(2)
    >>> # the Hilbert space has dimension 2
    >>> print brand_new_site.dim
    2
    >>> # the only operator is the identity
    >>> print brand_new_site.operators
    {'id': array([[ 1.,  0.],
           [ 0.,  1.]])}
    """
    def __init__(self, dim):
    	"""Creates an empty site of dimension dim.
    
	Raises
	------
	DMRGException
	    if `dim` < 1.

	Notes	
	-----
	Postcond : The identity operator (ones in the diagonal, zeros elsewhere)
	is added to the `self.operators` dictionary.
    	"""
    	if dim < 1:
    	    raise DMRGException("Site dim must be at least 1")
    	super(Site, self).__init__()
    	self.dim = dim
	self.operators = { "id" : np.eye(self.dim, self.dim) }
    
    def add_operator(self, operator_name):
    	"""Adds an operator to the site.
    
    	Parameters
	----------
    	operator_name : string
	    The operator name.

	Raises
	------
	DMRGException 
	    if `operator_name` is already in the dict.
	    
	Notes
	-----
	Postcond:

        - `self.operators` has one item more, and
        - the newly created operator is a (`self.dim`, `self.dim`)
          matrix of full of zeros.
	
	Examples
	--------
	>>> from dmrg_solution.core.sites import Site
	>>> new_site = Site(2)
	>>> print new_site.operators.keys()
	['id']
	>>> new_site.add_operator('s_z')
	>>> print new_site.operators.keys()
	['s_z', 'id']
	>>> # note that the newly created op has all zeros
	>>> print new_site.operators['s_z']
	[[ 0.  0.]
	 [ 0.  0.]]
        """
	if str(operator_name) in self.operators.keys():
    	    raise DMRGException("Operator name exists already")
    	else:
    	    self.operators[str(operator_name)] = np.zeros((self.dim, self.dim))

class PauliSite(Site):
    """A site for spin 1/2 models.
    
    You use this site for models where the single sites are spin
    one-half sites. The Hilbert space is ordered such as the first state
    is the spin dn, and the second state is the spin up. Therefore e.g.
    you have the following relation between operator matrix elements:

    .. math::
        \langle \dnarrow \left| A \\right|\uparrow \\rangle = A_{0,1}

    Notes
    -----
    Postcond : The site has already built-in the spin operators for s_z, s_p, s_m.

    Examples
    --------
    >>> from dmrg_solution.core.sites import PauliSite
    >>> pauli_site = PauliSite()
    >>> # check all it's what you expected
    >>> print pauli_site.dim
    2
    >>> print pauli_site.operators.keys()
    ['s_p', 's_z', 's_m', 'id']
    >>> print pauli_site.operators['s_z']
    [[-1.  0.]
     [ 0.  1.]]
    >>> print pauli_site.operators['s_x']
    [[ 0.  1.]
     [ 1.  0.]]
    """
    def __init__(self):
	"""Creates the spin one-half site with Pauli matrices.

	Notes
	-----
	Postcond : the dimension is set to 2, and the Pauli matrices
	are added as operators.
	"""
        super(PauliSite, self).__init__(2)
	# add the operators
        self.add_operator("s_z")
        self.add_operator("s_x")
	# for clarity
        s_z = self.operators["s_z"]
        s_x = self.operators["s_x"]
	# set the matrix elements different from zero to the right values
        s_z[0, 0] = -1.0
        s_z[1, 1] = 1.0
        s_x[0, 1] = 1.0
        s_x[1, 0] = 1.0

class SpinOneHalfSite(Site):
    """A site for spin 1/2 models.
    
    You use this site for models where the single sites are spin
    one-half sites. The Hilbert space is ordered such as the first state
    is the spin dn, and the second state is the spin up. Therefore e.g.
    you have the following relation between operator matrix elements:

    .. math::
        \langle \dnarrow \left| A \\right|\uparrow \\rangle = A_{0,1}

    Notes
    -----
    Postcond : The site has already built-in the spin operators for s_z, s_p, s_m.

    Examples
    --------
    >>> from dmrg_solution.core.sites import SpinOneHalfSite
    >>> spin_one_half_site = SpinOneHalfSite()
    >>> # check all it's what you expected
    >>> print spin_one_half_site.dim
    2
    >>> print spin_one_half_site.operators.keys()
    ['s_p', 's_z', 's_m', 'id']
    >>> print spin_one_half_site.operators['s_z']
    [[-0.5  0. ]
     [ 0.   0.5]]
    >>> print spin_one_half_site.operators['s_p']
    [[ 0.  0.]
     [ 1.  0.]]
    >>> print spin_one_half_site.operators['s_m']
    [[ 0.  1.]
     [ 0.  0.]]
    """
    def __init__(self):
	"""Creates the spin one-half site.

	Notes
	-----
	Postcond : the dimension is set to 2, and the Pauli matrices
	are added as operators.
	"""
        super(SpinOneHalfSite, self).__init__(2)
	# add the operators
        self.add_operator("s_z")
        self.add_operator("s_p")
        self.add_operator("s_m")
        self.add_operator("s_x")
	# for clarity
        s_z = self.operators["s_z"]
        s_p = self.operators["s_p"]
        s_m = self.operators["s_m"]
        s_x = self.operators["s_x"]
	# set the matrix elements different from zero to the right values
        s_z[0, 0] = -0.5
        s_z[1, 1] = 0.5
        s_p[1, 0] = 1.0
        s_m[0, 1] = 1.0
        s_x[0, 1] = 0.5
        s_x[1, 0] = 0.5


class ElectronicSite(Site):
    """A site for electronic models
    
    You use this site for models where the single sites are electron
    sites. The Hilbert space is ordered such as:

    - the first state, labelled 0,  is the empty site,
    - the second, labelled 1, is spin dn, 
    - the third, labelled 2, is spin up, and 
    - the fourth, labelled 3, is double occupancy.
    
    Notes
    -----
    Postcond: The site has already built-in the spin operators for: 

    - c_up : destroys an spin up electron,
    - c_up_dag, creates an spin up electron,
    - c_dn, destroys an spin down electron,
    - c_dn_dag, creates an spin down electron,
    - s_z, component z of spin,
    - s_p, raises the component z of spin,
    - s_m, lowers the component z of spin,
    - n_up, number of electrons with spin up,
    - n_dn, number of electrons with spin down,
    - n, number of electrons, i.e. n_up+n_dn, and
    - u, number of double occupancies, i.e. n_up*n_dn.

    Examples
    --------
    >>> from dmrg_solution.core.sites import ElectronicSite
    >>> hubbard_site = ElectronicSite()
    >>> # check all it's what you expected
    >>> print hubbard_site.dim
    4
    >>> print hubbard_site.operators.keys() # doctest: +ELLIPSIS
    ['s_p', ...]
    >>> print hubbard_site.operators['n_dn']
    [[ 0.  0.  0.  0.]
     [ 0.  1.  0.  0.]
     [ 0.  0.  0.  0.]
     [ 0.  0.  0.  1.]]
    >>> print hubbard_site.operators['n_up']
    [[ 0.  0.  0.  0.]
     [ 0.  0.  0.  0.]
     [ 0.  0.  1.  0.]
     [ 0.  0.  0.  1.]]
    >>> print hubbard_site.operators['u']
    [[ 0.  0.  0.  0.]
     [ 0.  0.  0.  0.]
     [ 0.  0.  0.  0.]
     [ 0.  0.  0.  1.]]
    """
    def __init__(self):
        super(ElectronicSite, self).__init__(8)
	# add the operators
        self.add_operator("rp_up")
        self.add_operator("rp_down")
        self.add_operator("rm_up")
        self.add_operator("rm_down")
        self.add_operator("rprm_up_plus")
        self.add_operator("rprm_down_plus")
        self.add_operator("rprm_up_minus")
        self.add_operator("rprm_down_minus")
        self.add_operator("rprm_up_plus_dag")
        self.add_operator("rprm_down_plus_dag")
        self.add_operator("rprm_up_minus_dag")
        self.add_operator("rprm_down_minus_dag")
        self.add_operator("dimer")
        self.add_operator("n_up")
        self.add_operator("n_down")
        self.add_operator("n")
        self.add_operator("u")
	# for clarity
        rm_up = self.operators["rm_up"]
        rm_down = self.operators["rm_down"]
        rp_up = self.operators["rp_up"]
        rp_down = self.operators["rp_down"]
        rprm_up_plus = self.operators["rprm_up_plus"]
        rprm_down_plus = self.operators["rprm_down_plus"]
        rprm_up_minus = self.operators["rprm_up_minus"]
        rprm_down_minus = self.operators["rprm_down_minus"]
        rprm_up_plus_dag = self.operators["rprm_up_plus_dag"]
        rprm_down_plus_dag = self.operators["rprm_down_plus_dag"]
        rprm_up_minus_dag = self.operators["rprm_up_minus_dag"]
        rprm_down_minus_dag = self.operators["rprm_down_minus_dag"]
        dimer = self.operators["dimer"]
        n_up = self.operators["n_up"]
        n_down = self.operators["n_down"]
        n = self.operators["n"]
        u = self.operators["u"]
	# set the matrix elements different from zero to the right values
	# TODO: missing s_p, s_m
        rp_up[4,0] = 1./sqrt(2.)
        rp_up[5,1] = -1./sqrt(2.)
        rp_up[4,2] = 1./sqrt(2.)
        rp_up[5,3] = -1./sqrt(2.)
        rp_up[0,4] = 1./sqrt(2.)
        rp_up[2,4] = 1./sqrt(2.)
        rp_up[1,5] = 1./sqrt(2.)
        rp_up[2,5] = 1./sqrt(2.)
        rp_up[4,6] = 1.
        rp_up[4,7] = 1.
        
        rm_up[5,0] = -1./sqrt(2.)
        rm_up[4,1] = 1./sqrt(2.)
        rm_up[5,2] = 1./sqrt(2.)
        rm_up[4,3] = -1./sqrt(2.)
        rm_up[1,4] = -1./sqrt(2.)
        rm_up[3,4] = 1./sqrt(2.)
        rm_up[0,5] = -1./sqrt(2.)
        rm_up[3,5] = -1./sqrt(2.)
        rm_up[5,6] = -1.
        rm_up[5,7] = 1.
        
#        rp_down[4,0] = 1./sqrt(2.)
#        rp_down[5,1] = 1./sqrt(2.)
#        rp_down[4,2] = 1./sqrt(2.)
#        rp_down[5,3] = -1./sqrt(2.)
#        rp_down[6,4] = 1.
#        rp_down[4,7] = 1.
#        
#        rm_down[5,0] = -1./sqrt(2.)
#        rm_down[4,1] = -1./sqrt(2.)
#        rm_down[5,2] = 1./sqrt(2.)
#        rm_down[4,3] = -1./sqrt(2.)
#        rm_down[6,5] = -1.
#        rm_down[5,7] = 1.
        

        rprm_up_plus =rp_up + rm_up
        rprm_up_minus = rp_up - rm_up
        rprm_down_plus = rp_down + rm_down
        rprm_down_minus = rp_down - rm_down
        
        rprm_up_plus_dag = rp_up.T + rm_up.T
        rprm_up_minus_dag = rp_up.T - rm_up.T
        rprm_down_plus_dag = rp_down.T + rm_down.T
        rprm_down_minus_dag = rp_down.T - rm_down.T

        n_even_up = 0.5*rprm_up_plus_dag.dot(rprm_up_plus)
        n_even_down = 0.5*rprm_down_plus_dag.dot(rprm_down_plus)
        n_odd_up = 0.5*rprm_up_minus_dag.dot(rprm_up_minus)
        n_odd_down = 0.5*rprm_down_minus_dag.dot(rprm_down_minus)
#        
#        n_up = n_even_up + n_odd_up
#        n_down = n_even_down + n_odd_down
#        n = n_up+n_down
#       
        dimer = 2.*(rp_up.T.dot(rp_up) - rm_up.T.dot(rm_up))
        #dimer = rp_up.T.dot(rp_up) + rp_down.T.dot(rp_down) - rm_up.T.dot(rm_up) - rm_down.T.dot(rm_down) 
#        u= (n_even_up - 0.5 *np.eye(8,8)).dot(n_even_down - 0.5 *np.eye(8,8))
#        u= (n_odd_up - 0.5 *np.eye(8,8)).dot(n_odd_down - 0.5 *np.eye(8,8))
        u = np.add((n_even_up - 0.5 *np.eye(8,8)).dot(n_even_down - 0.5 *np.eye(8,8)), (n_odd_up - 0.5 *np.eye(8,8)).dot(n_odd_down - 0.5 *np.eye(8,8)))
        np.set_printoptions(precision=3, suppress=True)

        #u=n_odd_down +n_even_down
#        u = dimer
        for i in range(len(u)):
            for ii in range(len(u)):
                if u[i,ii] != 0.0:
                    print '\t' + '[%s,%s] = %s' %( i, ii, u[i,ii])
#        for ii in range(len(w)):
#            print w[ii], v[:, ii]
        w, v = np.linalg.eigh(rp_up)
        print w
        print v
        print w[0], v[:,0]
        w, v = np.linalg.eigh(rm_up)
        print w
        print v
        #print w[0], v[:,0]
        w, v = np.linalg.eigh(dimer)
        print w
        print v
        #print w[0], v[:,0]
        
        
        
#print rp_up.T.dot(rp_up)
#
#        names = self.operators.keys()
#        names.remove('id')
#        names.remove('rp_up')
#        names.remove('rp_down')
#        names.remove('rm_up')
#        names.remove('rm_down')
#        for name in names: 
#            for i in range(len(eval(name))):
#                for ii in range(len(eval(name))):
#                    if eval(name)[i,ii] != 0.0:
#                        print '\t' + name + '[%s,%s] = %s' %( i, ii, eval(name)[i,ii])
#            print '\n'


