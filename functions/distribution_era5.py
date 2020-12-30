import os
import numpy as np
import pandas as pd
import matplotlib.cm as cm
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from astropy.io import fits

def cumulative_distribution_era5(site, path, month_list, out_file, npix):
    sites_dict = {'strip': {'lllong':-16.5,
                            'lllat':28.0,
                            'urlong':-16.25,
                            'urlat':28.25}
                  }

    lllong = sites_dict[site]['lllong']
    lllat = sites_dict[site]['lllat']
    urlong = sites_dict[site]['urlong']
    urlat = sites_dict[site]['urlat']

    month_hdul = fits.HDUList()
    mese = 0
    for month in month_list:
        print("Month:", month)
        df_tqi = pd.DataFrame()
        df_tql = pd.DataFrame()
        df_tqv = pd.DataFrame()
        # df_qv10m = pd.DataFrame()
        df_ps = pd.DataFrame()
        df_ts = pd.DataFrame()
        df_t10m = pd.DataFrame()
        df_u10m = pd.DataFrame()
        df_v10m = pd.DataFrame()

        file_list = os.listdir(month)

        for day in file_list:
            data = Dataset(month+"/"+day, mode='r')
            # print(month+"/"+day)
            lons = data.variables['longitude'][:]
            lats = data.variables['latitude'][:]
            
            lon, lat = np.meshgrid(lons, lats)
            
            # Single pixel anaylis for Pico del Teide.
            idx = np.where((lon[1, :] < urlong) & (lon[1, :] >= lllong))
            idy = np.where((lat[:, 1] <= urlat) & (lat[:, 1] > lllat))
            
            idx = idx[0]
            idy = idy[0] 
            try:
                WV_l_tqi = data.variables["tciw"][:, idy, idx]
            except:
                WV_l_tqi = np.zeros((24, 1, 1))
            
            try: 
                WV_l_tql = data.variables["tclw"][:, idy, idx]
            except:
                WV_l_tql = np.zeros((24, 1, 1))
                                    
            WV_l_tqv = data.variables["tcwv"][:, idy, idx]
            # There are not data available for QV10M
            # WV_l_qv10m = data.variables["QV10M"][:, np.amin(idy):np.amax(idy)+1, np.amin(idx):np.amax(idx)+1]
            # ....
            WV_l_ps = data.variables["sp"][:, idy, idx]
            WV_l_ts = data.variables["skt"][:, idy, idx]
            WV_l_t2m = data.variables["t2m"][:, idy,idx]
            WV_l_t10m = ((WV_l_t2m - WV_l_ts) / 2.0)*10.0 + WV_l_ts
            WV_l_u10m = data.variables["u10"][:, idy, idx]
            WV_l_v10m = data.variables["v10"][:, idy, idx]
                       
            tmp_tqi = np.sum(np.sum(WV_l_tqi, axis=1), axis=1)/(np.shape(idx)[0]*np.shape(idy)[0])
            tmp_tql = np.sum(np.sum(WV_l_tql, axis=1), axis=1)/(np.shape(idx)[0]*np.shape(idy)[0])
            tmp_tqv = np.sum(np.sum(WV_l_tqv, axis=1), axis=1)/(np.shape(idx)[0]*np.shape(idy)[0])
            # tmp_qv10m = np.sum(np.sum(WV_l_qv10m, axis=1), axis=1)/(np.shape(idx)[0]*np.shape(idy)[0])
            tmp_ps = np.sum(np.sum(WV_l_ps, axis=1), axis=1)/(np.shape(idx)[0]*np.shape(idy)[0])
            tmp_ts = np.sum(np.sum(WV_l_ts, axis=1), axis=1)/(np.shape(idx)[0]*np.shape(idy)[0])
            tmp_t10m = np.sum(np.sum(WV_l_t10m, axis=1), axis=1)/(np.shape(idx)[0]*np.shape(idy)[0])
            tmp_u10m = np.sum(np.sum(WV_l_u10m, axis=1), axis=1)/(np.shape(idx)[0]*np.shape(idy)[0])
            tmp_v10m = np.sum(np.sum(WV_l_v10m, axis=1), axis=1)/(np.shape(idx)[0]*np.shape(idy)[0])
            
            to_app_tqi = pd.Series(tmp_tqi)
            to_app_tql = pd.Series(tmp_tql)
            to_app_tqv = pd.Series(tmp_tqv)
            # to_app_qv10m = pd.Series(tmp_qv10m)
            to_app_ps = pd.Series(tmp_ps)
            to_app_ts = pd.Series(tmp_ts)
            to_app_t10m = pd.Series(tmp_t10m)
            to_app_u10m = pd.Series(tmp_u10m)
            to_app_v10m = pd.Series(tmp_v10m)
            
            if len(to_app_tqv.values[to_app_tqv.values > 12.0]) != 0:
                print("Day skipped for bad weather, TQV: {}", np.amax(to_app_tqv.values))
            else:
                df_tqi = df_tqi.append(to_app_tqi, ignore_index=True)
                df_tql = df_tql.append(to_app_tql, ignore_index=True)
                df_tqv = df_tqv.append(to_app_tqv, ignore_index=True)
                # df_qv10m = df_qv10m.append(to_app_qv10m, ignore_index=True)
                df_ps = df_ps.append(to_app_ps, ignore_index=True)
                df_ts = df_ts.append(to_app_ts, ignore_index=True)
                df_t10m = df_t10m.append(to_app_t10m, ignore_index=True)
                df_u10m = df_u10m.append(to_app_u10m, ignore_index=True)
                df_v10m = df_v10m.append(to_app_v10m, ignore_index=True)

        cdf_mat_tqi = np.array([])
        cdf_mat_tql = np.array([])
        cdf_mat_tqv = np.array([])
        # cdf_mat_qv10m = np.array([])
        cdf_mat_ps = np.array([])
        cdf_mat_ts = np.array([])
        cdf_mat_t10m = np.array([])
        cdf_mat_u10m = np.array([])
        cdf_mat_v10m = np.array([])

        print("Building thr CDFs...")
        # print(np.shape(df_tqi), np.shape(df_tql), np.shape(df_tqv), np.shape(df_v10m))
        for i in range(0, 24):
            out_tqi = np.sort(df_tqi[i])
            index_tqi = np.round(np.linspace(0, np.shape(df_tqi)[0]-1, 101)).astype(int)
            out_tqi = out_tqi[index_tqi]

            out_tql = np.sort(df_tql[i])
            index_tql = np.round(np.linspace(0, np.shape(df_tql)[0]-1, 101)).astype(int)
            out_tql = out_tql[index_tql]

            out_tqv = np.sort(df_tqv[i])
            index_tqv = np.round(np.linspace(0, np.shape(df_tqv)[0]-1, 101)).astype(int)
            out_tqv = out_tqv[index_tqv]

            # out_qv10m = np.sort(df_qv10m[i])
            # index_qv10m = np.round(np.linspace(0, np.shape(df_qv10m)[0]-1, 101)).astype(int)
            # out_qv10m = out_qv10m[index_qv10m]

            out_ps = np.sort(df_ps[i])
            index_ps = np.round(np.linspace(0, np.shape(df_ps)[0]-1, 101)).astype(int)
            out_ps = out_ps[index_ps]

            out_ts = np.sort(df_ts[i])
            index_ts = np.round(np.linspace(0, np.shape(df_ts)[0]-1, 101)).astype(int)
            out_ts = out_ts[index_ts]

            out_t10m = np.sort(df_t10m[i])
            index_t10m = np.round(np.linspace(0, np.shape(df_t10m)[0]-1, 101)).astype(int)
            out_t10m = out_t10m[index_t10m]

            out_u10m = np.sort(df_u10m[i])
            index_u10m = np.round(np.linspace(0, np.shape(df_u10m)[0]-1, 101)).astype(int)
            out_u10m = out_u10m[index_u10m]

            out_v10m = np.sort(df_v10m[i])
            index_v10m = np.round(np.linspace(0, np.shape(df_v10m)[0]-1, 101)).astype(int)
            out_v10m = out_v10m[index_v10m]

            cdf_mat_tqi = np.append(cdf_mat_tqi, out_tqi)
            cdf_mat_tql = np.append(cdf_mat_tql, out_tql)
            cdf_mat_tqv = np.append(cdf_mat_tqv, out_tqv)
            # cdf_mat_qv10m = np.append(cdf_mat_qv10m, out_qv10m)
            cdf_mat_ps = np.append(cdf_mat_ps, out_ps)
            cdf_mat_ts = np.append(cdf_mat_ts, out_ts)
            cdf_mat_t10m = np.append(cdf_mat_t10m, out_t10m)
            cdf_mat_u10m = np.append(cdf_mat_u10m, out_u10m)
            cdf_mat_v10m = np.append(cdf_mat_v10m, out_v10m)

        
        cdf_mat_tqi = np.reshape(cdf_mat_tqi, (24, 101))
        cdf_mat_tql = np.reshape(cdf_mat_tql, (24, 101))
        cdf_mat_tqv = np.reshape(cdf_mat_tqv, (24, 101))
        # cdf_mat_qv10m = np.reshape(cdf_mat_qv10m, (24, 101))
        cdf_mat_ps = np.reshape(cdf_mat_ps, (24, 101))
        cdf_mat_ts = np.reshape(cdf_mat_ts, (24, 101))
        cdf_mat_t10m = np.reshape(cdf_mat_t10m, (24, 101))
        cdf_mat_u10m = np.reshape(cdf_mat_u10m, (24, 101))
        cdf_mat_v10m = np.reshape(cdf_mat_v10m, (24, 101))

        print("The CDFs was built")

        hdu = fits.BinTableHDU.from_columns([
        fits.Column(name='TQI', format='101D', array=cdf_mat_tqi), fits.Column(name='TQL', format='101D', array=cdf_mat_tql), fits.Column(name='TQV', format='101D', array=cdf_mat_tqv), fits.Column(name='PS', format='101D', array=cdf_mat_ps), fits.Column(name='TS', format='101D', array=cdf_mat_ts), fits.Column(name='T10M', format='101D', array=cdf_mat_t10m), fits.Column(name='U10M', format='101D', array=cdf_mat_u10m), fits.Column(name='V10M', format='101D', array=cdf_mat_v10m)])

        hdu.header['MONTH'] = mese
        hdu.header['PROBSTRT'] = 0.0
        hdu.header['PROBSTOP'] = 1.0
        hdu.header['PROBSTEP'] = 0.01
        hdu.header['NSTEP'] = 101
        hdu.header['SOURCE'] = "ERA-5 data from 1980 to 2020"

        month_hdul.append(hdu)
        mese = mese + 1

    print("Writting the CDF on file...")
    print(month_hdul[1].header["probstrt"])
    month_hdul.writeto(out_file, overwrite=True)

    return df_tqv
