"""A few convenience functions to setup the Hubbard model.

.. math::
    H=\sum_{i}\vec{S}_{i}\cdot\vec{S}_{i+1}=
    \sum_{i}\left[S^{z}_{i}S^{z}_{i+1}+
    \frac{1}{2}\left(S^{\dagger}_{i}S^{-}_{i+1}+
    S^{-}_{i}S^{\dagger}_{i+1}\right)\right]

"""
class HubbardModel(object):
    """Implements a few convenience functions for Hubbard model.
    
    Does exactly that.
    """
    def __init__(self):
        super(HubbardModel, self).__init__()
	self.U = 0.
		
    def set_hamiltonian(self, system):
        """Sets a system Hamiltonian to the Hubbard Hamiltonian.
    
        Does exactly this. If the system hamiltonian has some other terms on
        it, there are not touched. So be sure to use this function only in
        newly created `System` objects.
    
        Parameters
        ----------
        system : a System.
            The System you want to set the Hamiltonian for.
        """
        system.clear_hamiltonian()
        if 'bh' in system.left_block.operators.keys():
            system.add_to_hamiltonian(left_block_op='bh')
        if 'bh' in system.right_block.operators.keys():
            system.add_to_hamiltonian(right_block_op='bh')
        system.add_to_hamiltonian('dimer', 'id', 'id', 'id', -(1. - self.U))
        system.add_to_hamiltonian('id', 'dimer', 'id', 'id', -(1. - self.U))
        system.add_to_hamiltonian('id', 'id', 'dimer', 'id', -(1. - self.U))
        system.add_to_hamiltonian('id', 'id', 'id', 'dimer', -(1. - self.U))
        
#        system.add_to_hamiltonian('dimer', 'id', 'id', 'id', self.U)
#        system.add_to_hamiltonian('id', 'dimer', 'id', 'id', self.U)
#        system.add_to_hamiltonian('id', 'id', 'dimer', 'id', self.U)
#        system.add_to_hamiltonian('id', 'id', 'id', 'dimer', self.U)

        system.add_to_hamiltonian('rprm_up_minus_dag', 'rprm_up_plus', 'id', 'id', -(1. + self.U)/2.)
        system.add_to_hamiltonian('rprm_down_minus_dag', 'rprm_down_plus', 'id', 'id', -(1. + self.U)/2.)
        system.add_to_hamiltonian('rprm_up_minus', 'rprm_up_plus_dag', 'id',  'id', (1. + self.U)/2.)
        system.add_to_hamiltonian('rprm_down_minus', 'rprm_down_plus_dag', 'id',  'id', (1. + self.U)/2.)
   
        system.add_to_hamiltonian('id', 'rprm_up_minus_dag', 'rprm_up_plus', 'id', -(1.+self.U)/2.)
        system.add_to_hamiltonian('id', 'rprm_down_minus_dag', 'rprm_down_plus', 'id', -(1.+self.U)/2.)
        system.add_to_hamiltonian('id', 'rprm_up_minus', 'rprm_up_plus_dag', 'id', (1.+self.U)/2.)
        system.add_to_hamiltonian('id', 'rprm_down_minus', 'rprm_down_plus_dag', 'id', (1.+self.U)/2.)

        system.add_to_hamiltonian('id','id',  'rprm_up_minus_dag', 'rprm_up_plus', -(1.+self.U)/2.)
        system.add_to_hamiltonian('id','id',  'rprm_down_minus_dag', 'rprm_down_plus', -(1.+self.U)/2.)
        system.add_to_hamiltonian('id','id',  'rprm_up_minus', 'rprm_up_plus_dag', (1.+self.U)/2.)
        system.add_to_hamiltonian('id','id',  'rprm_down_minus', 'rprm_down_plus_dag', (1.+self.U)/2.)

    def set_block_hamiltonian(self, tmp_matrix_for_bh, system):
        """Sets the block Hamiltonian to the Hubbard model block Hamiltonian.
    
        Parameters
        ----------
	tmp_matrix_for_bh : a numpy array of ndim = 2.
	    An auxiliary matrix to keep track of the result.
        system : a System.
            The System you want to set the Hamiltonian for.
        """
        # If you have a block hamiltonian in your block, add it
        if 'bh' in system.growing_block.operators.keys():
            system.add_to_block_hamiltonian(tmp_matrix_for_bh, 'bh', 'id')
        system.add_to_block_hamiltonian(tmp_matrix_for_bh, 'id', 'dimer', -(1. - self.U))
        system.add_to_block_hamiltonian(tmp_matrix_for_bh, 'dimer', 'id', -(1. - self.U))
#        system.add_to_block_hamiltonian(tmp_matrix_for_bh, 'id', 'dimer', self.U)
#        system.add_to_block_hamiltonian(tmp_matrix_for_bh, 'dimer', 'id', self.U)
        system.add_to_block_hamiltonian(tmp_matrix_for_bh, 'rprm_up_minus_dag', 'rprm_up_plus', -(1.+self.U)/2.)
        system.add_to_block_hamiltonian(tmp_matrix_for_bh, 'rprm_down_minus_dag', 'rprm_down_plus', -(1.+self.U)/2.)
        system.add_to_block_hamiltonian(tmp_matrix_for_bh, 'rprm_up_minus', 'rprm_up_plus_dag', (1.+self.U)/2.)
        system.add_to_block_hamiltonian(tmp_matrix_for_bh, 'rprm_down_minus', 'rprm_down_plus_dag', (1.+self.U)/2.)
#        system.add_to_block_hamiltonian(tmp_matrix_for_bh, 'id', 'u', self.U)
#        system.add_to_block_hamiltonian(tmp_matrix_for_bh, 'u', 'id', self.U)
    
    def set_operators_to_update(self, system):
        """Sets the operators to update to the ones for the Hubbard model.
    
        Parameters
        ----------
        system : a System.
            The System you want to set the Hamiltonian for.

	Notes
	-----
	The block Hamiltonian, althought needs to be updated, is treated
	separately by the very functions in the `System` class.
        """
        system.add_to_operators_to_update('rprm_up_plus_dag', site_op='rprm_up_plus_dag')
        system.add_to_operators_to_update('rprm_down_plus_dag', site_op='rprm_down_plus_dag')
        system.add_to_operators_to_update('rprm_up_minus_dag', site_op='rprm_up_minus_dag')
        system.add_to_operators_to_update('rprm_down_minus_dag', site_op='rprm_down_minus_dag')
        system.add_to_operators_to_update('rprm_up_plus', site_op='rprm_up_plus')
        system.add_to_operators_to_update('rprm_down_plus', site_op='rprm_down_plus')
        system.add_to_operators_to_update('rprm_up_minus', site_op='rprm_up_minus')
        system.add_to_operators_to_update('rprm_down_minus', site_op='rprm_down_minus')
        system.add_to_operators_to_update('dimer', site_op='dimer')
        #system.add_to_operators_to_update('u', site_op='u')
