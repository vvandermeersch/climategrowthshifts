import cdsapi

dataset = "reanalysis-era5-land"
client = cdsapi.Client(timeout=3600,quiet=True,debug=False)


def download_era5land(year):
  
  request = {
    "variable": ["2m_temperature"],
    "year": year,
    "month": [
      "01", "02", "03",
      "04", "05", "06",
      "07", "08", "09",
      "10", "11", "12"
    ],
    "day": [
      "01", "02", "03",
      "04", "05", "06",
      "07", "08", "09",
      "10", "11", "12",
      "13", "14", "15",
      "16", "17", "18",
      "19", "20", "21",
      "22", "23", "24",
      "25", "26", "27",
      "28", "29", "30",
      "31"
    ],
    "time": [
        "00:00", "01:00", "02:00",
        "03:00", "04:00", "05:00",
        "06:00", "07:00", "08:00",
        "09:00", "10:00", "11:00",
        "12:00", "13:00", "14:00",
        "15:00", "16:00", "17:00",
        "18:00", "19:00", "20:00",
        "21:00", "22:00", "23:00"
    ],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": [47, -122, 46, -121.5]
  }
  file_name = "hourly_tmp" + "_" + year + ".nc"
  client.retrieve(dataset, request).download(file_name)

