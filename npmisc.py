import numpy as np

# ind : indices
# ws : weights
# bins : bin edge array
# Note re digitize : it sets elements < min(bins) to 0 and elements >= max(bins) to len(kbins)
# Therefore, it generates a maximum of len(kbins)+1 bins
# We therefore trim these in what follows
def mybincount(ind, ws, bins):
    nbins = len(bins) + 1
    count=np.bincount(ind,weights=ws,minlength=nbins)
    return count[1:-1] # strip padding bins



# Generate the k-grid
# Make this corresponding to a full FFT, instead of using a real FFT.
# That way, we don't need to worry about the count case, even though it is a little inefficient
def make_k_grid(N, L):
    fac = 2*np.pi*N/L;
    associated_kx=np.fft.fftfreq(N,d=1)*fac # kx vector. Python gives (0,1,...)/N. k=2*pi*q/L, so multiply it by 2*pi*N/L
    associated_ky=np.fft.fftfreq(N,d=1)*fac
    associated_kz=np.fft.fftfreq(N,d=1)*fac
    ksize = (len(associated_kx), len(associated_ky), len(associated_kz))
    kxyz = np.meshgrid(associated_kx, associated_ky, associated_kz, indexing='ij')
    associated_k=np.sqrt(kxyz[0]**2 + kxyz[1]**2 + kxyz[2]**2)
    associated_k[0,0,0] = 1.0e-10

    return ksize, kxyz, associated_k
