import numpy as np
from netCDF4 import Dataset
from astropy.io import fits


def cumulative_distribution_era5(site, path, month_list, out_file, npix, pix_size):
    sites_dict = {
                    'tenerife': {
                                'lllong': -17.0,
                                'lllat': -27.5,
                                'urlong': -12.9,
                                'urlat': 29.5
                                },
                    'strip':    {
                                'lllong': -16.5110782 - pix_size * npix,
                                'lllat': 28.3003912 - pix_size * npix,
                                'urlong': -16.5110782 + pix_size * npix,
                                'urlat': 28.3003912 + pix_size * npix},
                 }

    lllong = sites_dict[site]['lllong']
    lllat = sites_dict[site]['lllat']
    urlong = sites_dict[site]['urlong']
    urlat = sites_dict[site]['urlat']

    month_hdul = fits.HDUList()
    mese = 0
    for month in month_list:
        print("Mese:", month)

        data = Dataset(month, mode='r')
        lons = data.variables['longitude'][:]
        lats = data.variables['latitude'][:]

        lon, lat = np.meshgrid(lons, lats)

        idx = np.where((lon[1, :] < urlong) & (lon[1, :] > lllong))
        idy = np.where((lat[:, 1] < urlat) & (lat[:, 1] > lllat))

        try:
            WV_l_tqv = data.variables["tcwv"][:, 0, :, :]
            print(np.shape(WV_l_tqv))
        except ValueError as verr:
            WV_l_tqv = data.variables["tcwv"][:, :, :]
            print(np.shape(WV_l_tqv))

        TQV = WV_l_tqv[:, np.amin(idy):np.amax(idy)+1, np.amin(idx):np.amax(idx)+1]
        print("After BIN: ", np.shape(TQV))
        TQV = np.sum(np.sum(TQV, axis=1), axis=1)/(np.shape(idx)[1]*np.shape(idy)[1])
        print("Before BIN: ", np.shape(TQV))

        TQV = np.reshape(TQV, (24, int(len(TQV)/24)))
        cdf_mat_tqv = np.array([])

        for h in range(0, 24):
            cdf_ext_tqv = np.sort(TQV[h, :])
            idx_tqv = np.round(np.linspace(0, len(cdf_ext_tqv)-1, 101)).astype(int)
            cdf_tqv = cdf_ext_tqv[idx_tqv]
            cdf_mat_tqv = np.append(cdf_mat_tqv, cdf_tqv)

        cdf_mat_tqv = np.reshape(cdf_mat_tqv, (24, 101))

        print("Fine creazione distribuzioni")

        hdu = fits.BinTableHDU.from_columns([
              fits.Column(name='TQV', format='101D', array=cdf_mat_tqv),
              ])

        hdu.header['MONTH'] = mese
        hdu.header['PROBSTRT'] = 0.0
        hdu.header['PROBSTOP'] = 1.0
        hdu.header['PROBSTEP'] = 0.01
        hdu.header['NSTEP'] = 101
        hdu.header['SOURCE'] = "ERA-5 data from 1980 to 2019"

        month_hdul.append(hdu)
        mese = mese + 1

    print("scrivo le distribuzioni su file")

    print(month_hdul[1].header["probstrt"])

    month_hdul.writeto(out_file, overwrite=True)

    return 0
