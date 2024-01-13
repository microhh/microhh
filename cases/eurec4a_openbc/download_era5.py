from datetime import datetime
import ls2d

settings = {
    'central_lon' : -57.7,
    'central_lat' : 13.3,
    'start_date'  : datetime(year=2020, month=2, day=1, hour=0),
    'end_date'    :datetime(year=2020, month=2, day=2, hour=0),
    'area_size'   : 3,
    'case_name'   : 'eurec4a_openbc',
    'era5_path'   : '/home/scratch1/bart/LS2D/',
    'era5_expver' : 1,
    'write_log'   : False,
    'data_source' : 'CDS',
    'ntasks'      : 1
    }

ls2d.download_era5(settings)
