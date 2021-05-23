from o3seespy.base_model import OpenSeesObject


class SystemBase(OpenSeesObject):
    op_base_type = "system"


class BandGen(SystemBase):
    """
    The BandGen System Class
    
    This command is used to construct a BandGeneralSOE linear system of equation object. As the name implies, this class
    is used for matrix systems which have a banded profile. The matrix is stored as shown below in a 1dimensional array of
    size equal to the bandwidth times the number of unknowns. When a solution is required, the Lapack routines DGBSV and
    SGBTRS are used.
    """
    op_type = 'BandGen'

    def __init__(self, osi):
        """
        Initial method for BandGen

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.system.BandGen(osi)
        """
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


class BandSPD(SystemBase):
    """
    The BandSPD System Class
    
    This command is used to construct a BandSPDSOE linear system of equation object. As the name implies, this class is
    used for symmetric positive definite matrix systems which have a banded profile. The matrix is stored as shown below
    in a 1 dimensional array of size equal to the (bandwidth/2) times the number of unknowns. When a solution is
    required, the Lapack routines DPBSV and DPBTRS are used.
    """
    op_type = 'BandSPD'

    def __init__(self, osi):
        """
        Initial method for BandSPD

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.system.BandSPD(osi)
        """
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


class ProfileSPD(SystemBase):
    """
    The ProfileSPD System Class
    
    This command is used to construct a profileSPDSOE linear system of equation object. As the name implies, this class
    is used for symmetric positive definite matrix systems. The matrix is stored as shown below in a 1 dimensional array
    with only those values below the first non-zero row in any column being stored. This is sometimes also referred to
    as a skyline storage scheme.
    """
    op_type = 'ProfileSPD'

    def __init__(self, osi):
        """
        Initial method for ProfileSPD

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.system.ProfileSPD(osi)
        """
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


class ParallelProfileSPD(SystemBase):
    """
    The ParallelProfileSPD System Class

    This command is used to construct a parallel version of the profileSPDSOE linear system of equation object.
    As the name implies, this class
    is used for symmetric positive definite matrix systems. The matrix is stored as shown below in a 1 dimensional array
    with only those values below the first non-zero row in any column being stored. This is sometimes also referred to
    as a skyline storage scheme.
    """
    op_type = 'ParallelProfileSPD'

    def __init__(self, osi):
        """
        Initial method for ParallelProfileSPD

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.system.ParallelProfileSPD(osi)
        """
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


class SuperLU(SystemBase):
    """
    The SuperLU System Class
    
    This command is used to construct a SparseGEN linear system of equation object. As the name implies, this class is
    used for sparse matrix systems. The solution of the sparse matrix is carried out using `SuperLU`_.
    """
    op_type = 'SuperLU'

    def __init__(self, osi):
        """
        Initial method for SuperLU

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.system.SuperLU(osi)
        """
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


class UmfPack(SystemBase):
    """
    The UmfPack System Class
    
    This command is used to construct a sparse system of equations which uses the `UmfPack`_ solver.
    """
    op_type = 'UmfPack'

    def __init__(self, osi):
        """
        Initial method for UmfPack

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.system.UmfPack(osi)
        """
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


class FullGeneral(SystemBase):
    """
    The FullGeneral System Class
    
    This command is used to construct a Full General linear system of equation object. As the name implies, the class
    utilizes NO space saving techniques to cut down on the amount of memory used. If the matrix is of size, nxn, then
    storage for an nxn array is sought from memory when the program runs. When a solution is required, the Lapack
    routines DGESV and DGETRS are used... note::This type of system should almost never be used! This is because
    it requires a lot more memory than every other solver and takes more time in the actal solving operation
    than any other solver. It is required if the user is interested in looking at the global system matrix.
    """
    op_type = 'FullGeneral'

    def __init__(self, osi):
        """
        Initial method for FullGeneral

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.system.FullGeneral(osi)
        """
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


class SparseSYM(SystemBase):
    """
    The SparseSYM System Class
    
    This command is used to construct a sparse symmetric system of equations which uses a row-oriented solution method
    in the solution phase.
    """
    op_type = 'SparseSYM'

    def __init__(self, osi):
        """
        Initial method for SparseSYM

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.system.SparseSYM(osi)
        """
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


class Mumps(SystemBase):
    """
    The Mumps System Class
    
    Create a system of equations using the Mumps solver
    """
    op_type = 'Mumps'

    def __init__(self, osi, icntl14=None, icntl7=None):
        """
        Initial method for Mumps

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        icntl14: None, optional
            
        icntl7: None, optional
            
        """
        self.osi = osi
        self.icntl14 = icntl14
        self.icntl7 = icntl7
        self._parameters = [self.op_type]
        if getattr(self, 'icntl14') is not None:
            self._parameters += ['-ICNTL14', self.icntl14]
        if getattr(self, 'icntl7') is not None:
            self._parameters += ['-ICNTL7', self.icntl7]
        self.to_process(osi)


class BandGeneral(SystemBase):
    op_type = "BandGeneral"

    def __init__(self, osi):
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


class ParallelBandGeneral(SystemBase):
    op_type = "ParallelBandGeneral"

    def __init__(self, osi):
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


class SparseGeneral(SystemBase):
    op_type = "SparseGeneral"

    def __init__(self, osi):
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


def apply_mumps_or_sparse_general(osi, **kwargs):
    if osi.mp:
        Mumps(osi, **kwargs)
    else:
        SparseGeneral(osi)

def apply_mumps_or(osi, alt_system, **kwargs):
    if osi.mp:
        Mumps(osi, **kwargs)
    else:
        alt_system(osi)
