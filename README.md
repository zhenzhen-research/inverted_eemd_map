# inverted_eemd_map

# 1) import this package
import inverted_eemd_map

# 2) prepare data cube
file = './THOR_HI_38.3deg_C-D_ft_VGPS.fits'
hdu = fits.open(file)
data = hdu[0].data

# 3) using EEMD algorithm to decompose data
list_IMFs = inverted_eemd_map.eemd_decomposition(cube=data, num_processor=48)
# 4） reconstruct 3D data cube
recons_cube = inverted_eemd_map.reconstruct_inverted_cube(list_IMFs=list_IMFs, cube_shape=data.shape, imfn=0, n_imf_min=4)

# 5） add header and store results
hdu[0].data = recons_cube
hdu.writeto(file + '_invertedIMF' + '.fits', overwrite=True)
