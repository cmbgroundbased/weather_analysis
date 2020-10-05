import os
import numpy as np
import pandas as pd
import matplotlib.cm as cm
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from astropy.io import fits

def cumulative_distribution_merra2(site, path, month_list, out_file, npix):
    sites_dict = {'atacama':{'lllong':-67.9774138-0.35*npix,'lllat':-23.0755446-0.35*npix,'urlong':-67.9774138+0.35*npix,'urlat':-23.0755446+0.35*npix},
                  'tenerife':{'lllong':-17.0,'lllat':-27.5,'urlong':-12.9,'urlat':29.5},
                  'test_atc':{'lllong':-63.9774138-0.35,'lllat':-23.0755446-0.35,'urlong':-63.9774138+0.35,'urlat':-23.0755446+0.35},
                  'qubic': {'lllong':-66.474650-0.35*npix,'lllat':-24.191996 -0.35*npix,'urlong':-66.474650+0.35*npix,'urlat':-24.191996+0.35*npix},
                  'strip': {'lllong':-16.5110782-0.35*npix,'lllat':28.3003912-0.35*npix,'urlong':-16.5110782+0.35*npix,'urlat':28.3003912+0.35*npix}}

    lllong = sites_dict[site]['lllong']
    lllat = sites_dict[site]['lllat']
    urlong = sites_dict[site]['urlong']
    urlat = sites_dict[site]['urlat']

    month_hdul = fits.HDUList()
    mese = 0
    for month in month_list:
        print("Mese:", month)
        df_tqi = pd.DataFrame()
        df_tql = pd.DataFrame()
        df_tqv = pd.DataFrame()
        df_qv10m = pd.DataFrame()
        df_ps = pd.DataFrame()
        df_ts = pd.DataFrame()
        df_t10m = pd.DataFrame()
        df_u10m = pd.DataFrame()
        df_v10m = pd.DataFrame()

        file_list = os.listdir(month)

        for day in file_list:
            # print("GIorno: ", day)
            data = Dataset(month+"/"+day, mode='r')
            lons = data.variables['lon'][:]
            lats = data.variables['lat'][:]
            
            lon, lat = np.meshgrid(lons, lats)
            idx = np.where((lon[1, :] < urlong) & (lon[1, :] > lllong))
            idy = np.where((lat[:, 1] < urlat) & (lat[:, 1] > lllat))
            WV_l_tqi = data.variables["TQI"][:, np.amin(idy):np.amax(idy)+1, np.amin(idx):np.amax(idx)+1]
            WV_l_tql = data.variables["TQL"][:, np.amin(idy):np.amax(idy)+1, np.amin(idx):np.amax(idx)+1]
            WV_l_tqv = data.variables["TQV"][:, np.amin(idy):np.amax(idy)+1, np.amin(idx):np.amax(idx)+1]
            WV_l_qv10m = data.variables["QV10M"][:, np.amin(idy):np.amax(idy)+1, np.amin(idx):np.amax(idx)+1]
            WV_l_ps = data.variables["PS"][:, np.amin(idy):np.amax(idy)+1, np.amin(idx):np.amax(idx)+1]
            WV_l_ts = data.variables["TS"][:, np.amin(idy):np.amax(idy)+1, np.amin(idx):np.amax(idx)+1]
            WV_l_t10m = data.variables["T10M"][:, np.amin(idy):np.amax(idy)+1, np.amin(idx):np.amax(idx)+1]
            WV_l_u10m = data.variables["U10M"][:, np.amin(idy):np.amax(idy)+1, np.amin(idx):np.amax(idx)+1]
            WV_l_v10m = data.variables["V10M"][:, np.amin(idy):np.amax(idy)+1, np.amin(idx):np.amax(idx)+1]

            # spatial average (if nedeed)
            tmp_tqi = np.sum(np.sum(WV_l_tqi, axis=1), axis=1)/(np.shape(idx)[1]*np.shape(idy)[1])
            tmp_tql = np.sum(np.sum(WV_l_tql, axis=1), axis=1)/(np.shape(idx)[1]*np.shape(idy)[1])
            tmp_tqv = np.sum(np.sum(WV_l_tqv, axis=1), axis=1)/(np.shape(idx)[1]*np.shape(idy)[1])
            tmp_qv10m = np.sum(np.sum(WV_l_qv10m, axis=1), axis=1)/(np.shape(idx)[1]*np.shape(idy)[1])
            tmp_ps = np.sum(np.sum(WV_l_ps, axis=1), axis=1)/(np.shape(idx)[1]*np.shape(idy)[1])
            tmp_ts = np.sum(np.sum(WV_l_ts, axis=1), axis=1)/(np.shape(idx)[1]*np.shape(idy)[1])
            tmp_t10m = np.sum(np.sum(WV_l_t10m, axis=1), axis=1)/(np.shape(idx)[1]*np.shape(idy)[1])
            tmp_u10m = np.sum(np.sum(WV_l_u10m, axis=1), axis=1)/(np.shape(idx)[1]*np.shape(idy)[1])
            tmp_v10m = np.sum(np.sum(WV_l_v10m, axis=1), axis=1)/(np.shape(idx)[1]*np.shape(idy)[1])

            to_app_tqi = pd.Series(tmp_tqi)
            to_app_tql = pd.Series(tmp_tql)
            to_app_tqv = pd.Series(tmp_tqv)
            to_app_qv10m = pd.Series(tmp_qv10m)
            to_app_ps = pd.Series(tmp_ps)
            to_app_ts = pd.Series(tmp_ts)
            to_app_t10m = pd.Series(tmp_t10m)
            to_app_u10m = pd.Series(tmp_u10m)
            to_app_v10m = pd.Series(tmp_v10m)

            if len(to_app_tqv.values[to_app_tqv.values > 3.5]) != 0:
                print("Day skipped for bad weather")
            else:
                print("Day ok")
                print(to_app_tqv)
                df_tqi = df_tqi.append(to_app_tqi, ignore_index=True)
                df_tql = df_tql.append(to_app_tql, ignore_index=True)
                df_tqv = df_tqv.append(to_app_tqv, ignore_index=True)
                df_qv10m = df_qv10m.append(to_app_qv10m, ignore_index=True)
                df_ps = df_ps.append(to_app_ps, ignore_index=True)
                df_ts = df_ts.append(to_app_ts, ignore_index=True)
                df_t10m = df_t10m.append(to_app_t10m, ignore_index=True)
                df_u10m = df_u10m.append(to_app_u10m, ignore_index=True)
                df_v10m = df_v10m.append(to_app_v10m, ignore_index=True)
                
            # print("Fine giorno: ", day)

        cdf_mat_tqi = np.array([])
        cdf_mat_tql = np.array([])
        cdf_mat_tqv = np.array([])
        cdf_mat_qv10m = np.array([])
        cdf_mat_ps = np.array([])
        cdf_mat_ts = np.array([])
        cdf_mat_t10m = np.array([])
        cdf_mat_u10m = np.array([])
        cdf_mat_v10m = np.array([])

        print("Inizio a creare le distribuzioni")
        print(np.shape(df_tqi), np.shape(df_tql), np.shape(df_tqv), np.shape(df_v10m))
        for i in range(0, 24):
            # print("Ora: ", i)
            out_tqi = np.sort(df_tqi[i])
            index_tqi = np.round(np.linspace(0, np.shape(df_tqi)[0]-1, 101)).astype(int)
            out_tqi = out_tqi[index_tqi]

            out_tql = np.sort(df_tql[i])
            index_tql = np.round(np.linspace(0, np.shape(df_tql)[0]-1, 101)).astype(int)
            out_tql = out_tql[index_tql]

            out_tqv = np.sort(df_tqv[i])
            index_tqv = np.round(np.linspace(0, np.shape(df_tqv)[0]-1, 101)).astype(int)
            out_tqv = out_tqv[index_tqv]

            out_qv10m = np.sort(df_qv10m[i])
            index_qv10m = np.round(np.linspace(0, np.shape(df_qv10m)[0]-1, 101)).astype(int)
            out_qv10m = out_qv10m[index_qv10m]

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
            cdf_mat_qv10m = np.append(cdf_mat_qv10m, out_qv10m)
            cdf_mat_ps = np.append(cdf_mat_ps, out_ps)
            cdf_mat_ts = np.append(cdf_mat_ts, out_ts)
            cdf_mat_t10m = np.append(cdf_mat_t10m, out_t10m)
            cdf_mat_u10m = np.append(cdf_mat_u10m, out_u10m)
            cdf_mat_v10m = np.append(cdf_mat_v10m, out_v10m)

        cdf_mat_tqi = np.reshape(cdf_mat_tqi, (24, 101))
        cdf_mat_tql = np.reshape(cdf_mat_tql, (24, 101))
        cdf_mat_tqv = np.reshape(cdf_mat_tqv, (24, 101))
        cdf_mat_qv10m = np.reshape(cdf_mat_qv10m, (24, 101))
        cdf_mat_ps = np.reshape(cdf_mat_ps, (24, 101))
        cdf_mat_ts = np.reshape(cdf_mat_ts, (24, 101))
        cdf_mat_t10m = np.reshape(cdf_mat_t10m, (24, 101))
        cdf_mat_u10m = np.reshape(cdf_mat_u10m, (24, 101))
        cdf_mat_v10m = np.reshape(cdf_mat_v10m, (24, 101))

        print("Fine creazione distribuzioni")

        hdu = fits.BinTableHDU.from_columns([
        fits.Column(name='TQI', format='101D', array=cdf_mat_tqi), fits.Column(name='TQL', format='101D', array=cdf_mat_tql), fits.Column(name='TQV', format='101D', array=cdf_mat_tqv), fits.Column(name='QV10M', format='101D', array=cdf_mat_qv10m), fits.Column(name='PS', format='101D', array=cdf_mat_ps), fits.Column(name='TS', format='101D', array=cdf_mat_ts), fits.Column(name='T10M', format='101D', array=cdf_mat_t10m), fits.Column(name='U10M', format='101D', array=cdf_mat_u10m), fits.Column(name='V10M', format='101D', array=cdf_mat_v10m)])

        hdu.header['MONTH'] = mese
        hdu.header['PROBSTRT'] = 0.0
        hdu.header['PROBSTOP'] = 1.0
        hdu.header['PROBSTEP'] = 0.01
        hdu.header['NSTEP'] = 101
        hdu.header['SOURCE'] = "MERRA-2 data from 2000 to 2020"



        month_hdul.append(hdu)
        mese = mese + 1

    print("scrivo le distribuzioni su file")

    print(month_hdul[1].header["probstrt"])

    month_hdul.writeto(out_file, overwrite=True)

    return df_tqv
