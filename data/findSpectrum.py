import numpy as np
from numpy.linalg import eigh

class Spectrum(object):
    def __init__(self, nLat, lenLat, lamb, frac):
        # Simulation parameters
        self.N = nLat
        self.L = lenLat
        self.lamb = lamb
        self.frac = frac

        # Derivatives
        self.Ninv = self.N**-1.
        self.knyq = self.N//2+1
        self.dx = self.L/self.N
        self.dk = 2.*np.pi/self.L

    def getLattice(self):
        xL = self.xlist()
        kL = self.klist()
        dx = self.dx
        dk = self.dk
        knyq = self.knyq
        invFTcoeffs = self.inv_phases()
        dirFTcoeffs = self.dir_phases()
        return xL, kL, knyq, dx, dk, invFTcoeffs, dirFTcoeffs

    def lattice(self):
        return np.arange(self.N)
    def xlist(self):
        return self.lattice()*self.dx
    def klist(self):
        return np.roll(((self.lattice() - self.N//2)*self.dk), self.N//2)

    # Fourrier coefficients for direct and inverse transformation
    # Format: phases(x, y) i.e. axes are phases[x, k]
    def inv_phases(self):
        return np.exp(1j*np.outer(self.xlist(), self.klist()))
    def dir_phases(self):
        return self.Ninv * np.exp(-1j*np.outer(self.xlist(), self.klist()))

    def confOmega(self, x, lamb):
        if x>0.:
            return (1.+np.exp(2.*lamb*(x - self.N/self.frac//2)))**-1.
        else:
            return (1.+np.exp(-2.*lamb*(x + self.N/self.frac//2)))**-1.

    def Ux(self, lamb):
        return np.asarray([self.confOmega(x/self.dx-self.N//2, lamb) for x in self.xlist()])
    def Uk(self):
        return np.tensordot(self.Ux(self.lamb), self.dir_phases(), axes=(0,0))

    # The Laplacian is diagonal in k-space
    def specLaplacian(self):
        return np.diag(self.Ninv * self.klist()**2.)

    # This is the formal inverse fourier transform
    # This is the most time consuming computation in this notebook
    def kineticTerm(self):
        firstFT = np.tensordot(self.specLaplacian(), self.inv_phases(), axes=(0,1))
        doubleFT = np.tensordot(firstFT, self.inv_phases(), axes=(0,1))
        return doubleFT

    # The potential term is diagonal in configuration space
    def potentialTerm(self):
        return np.diag(self.Ux(self.lamb))

    def Hamiltonian(self):
        return np.real(self.kineticTerm() + self.potentialTerm())

    # Use python linalg package function to decompose matrix assumed Hermitian
    def eigenspectrum(self):
        H = self.Hamiltonian()
        En, phin = eigh(H)
        phin = np.transpose(phin)
        phik = np.tensordot(phin, self.dir_phases(), axes=(1,0))
        return En, phin, phik

    # Select only real eigenvectors and corresponding eigenvalues and x-space eigenvectors
    def getNonDegenerateSpectrum(self):
        En, phin, phik = self.eigenspectrum()
        reParts = np.linalg.norm(np.real(phik), axis=1)
        sensitivity = np.sort(reParts)[self.knyq-2]
        indices = [i for i in range(len(En)) if reParts[i] >= sensitivity]
        return sensitivity, En[indices], phin[indices], phik[indices]

    # Transformation matrix is the matrix of valid eigenvectors, but normalized to 1.
    def transfMatrix(self, phikV):
        tM = np.real(phikV)
        norms = np.linalg.norm(tM, axis=1)
        norms = np.transpose(np.tile(norms, (len(tM[0]), 1)))
        return tM/norms

    def writeFrequencies(self):
        sensitivity, EnV, phinV, phikV = self.getNonDegenerateSpectrum()        
        print('For parameters ', self.N, self.L, self.lamb, self.frac, ', sensitivity is:', sensitivity)
        tM = self.transfMatrix(phikV)

        with open('./frequencies_N'+str(self.N)+'_L'+str(int(self.L))+'_lambda'+str('%.3f'%self.lamb)+'.txt', 'w') as f:
            f.write(str(len(EnV))+'\n')
            for i in EnV:
                f.write(str(i)+'\n')
        with open('./transfMatrix_N'+str(self.N)+'_L'+str(int(self.L))+'_lambda'+str('%.3f'%self.lamb)+'.txt', 'w') as f:
            f.write(str(np.shape(tM)[0])+'\n')
            f.write(str(np.shape(tM)[1])+'\n')
            for vector in tM:
                for i in vector:
                    f.write(str(i)+'\n')
        return 'Done'

