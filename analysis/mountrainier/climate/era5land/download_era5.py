import cdsapi

dataset = "derived-era5-land-daily-statistics"
client = cdsapi.Client(timeout=3600,quiet=True,debug=False)


def download_era5land(year):
  
  request = {
    "variable": ["snow_depth_water_equivalent"],
    "year": year,
    "month": [
      "01", "02", "03",
      "04", "05", "06",
      "07", "08", "09",
      "10", "11", "12"],
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
    "daily_statistic": "daily_mean",
    "time_zone": "utc+00:00",
    "frequency": "1_hourly",
    "area": [47, -122, 46, -121.5]
  }
  file_name = "swe" + "_" + year + ".nc"
  client.retrieve(dataset, request).download(file_name)
