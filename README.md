# Weather tools 

The `wehater_tools` is a verry dummy set of `python >=3` scripts coincieved to create a statistical representation of the atmospheric condiction above a given observe site. 

# Usage

Help message of `./analyse_data.py -h`

```
usage: analyse_data.py [-h] [-p PATH] [-o OUT_FILE] [-ps PIX_SIZE] [-npix NUM_PIX]
                       [-m] [-var VARIABLE] [-ra REANALYSIS] [-d]
                       site

positional arguments:
  site                  Site name

optional arguments:
  -h, --help            show this help message and exit
  -p PATH, --path PATH  Data Directory
  -o OUT_FILE, --out_file OUT_FILE
                        Output File
  -ps PIX_SIZE, --pix_size PIX_SIZE
                        Select the pixel size of your reanalysis
  -npix NUM_PIX, --num_pix NUM_PIX
                        Number of pixel
  -m, --seasonal_maps   Return the seasonal maps
  -var VARIABLE, --variable VARIABLE
                        Atmospheric Variable
  -ra REANALYSIS, --reanalysis REANALYSIS
                        Select the reanalysis type of your data
  -d, --debug
```

## Site Name 

Available coordinate site:

```
sites_dict = {
  'atacama':{'lllong':-67.9774138-0.35*npix,'lllat':-23.0755446-0.35*npix,'urlong':-67.9774138+0.35*npix,'urlat':-23.0755446+0.35*npix},
  'tenerife':{'lllong':-17.0,'lllat':-27.5,'urlong':-12.9,'urlat':29.5},
  'test_atc':{'lllong':-63.9774138-0.35,'lllat':-23.0755446-0.35,'urlong':-63.9774138+0.35,'urlat':-23.0755446+0.35},
  'qubic': {'lllong':-66.474650-0.35*npix,'lllat':-24.191996 -0.35*npix,'urlong':-66.474650+0.35*npix,'urlat':-24.191996+0.35*npix},
  'strip': {'lllong':-16.5110782-0.35*npix,'lllat':28.3003912-0.35*npix,'urlong':-16.5110782+0.35*npix,'urlat':28.3003912+0.35*npix}
}
```




## Directory format and data organization - PATH

You have to specify the working PATH with all the data in this way `./analyse_data qubic --path ./data_atacama`. The working path design has to be as follow:

```
.
├── data_atacama
│   ├── m_01
│   │   ├── MERRA2_200.inst1_2d_asm_Nx.20000101.SUB.nc
│   │   ├── MERRA2_200.inst1_2d_asm_Nx.20000102.SUB.nc
│   │   └── ...
│   ├── m_02
│   │   ├── MERRA2_200.inst1_2d_asm_Nx.20000201.SUB.nc
│   │   ├── MERRA2_200.inst1_2d_asm_Nx.20000202.SUB.nc
│   │   └── ...
│   ├── m_03
│   │   ├── MERRA2_200.inst1_2d_asm_Nx.20000301.SUB.nc
│   │   ├── MERRA2_200.inst1_2d_asm_Nx.20000302.SUB.nc
│   │   └── ...
...
```
