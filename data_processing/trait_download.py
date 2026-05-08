import earthaccess
from datetime import datetime, timedelta
from pathlib import Path

AVIRIS_DATES = [datetime.strptime(d, "%Y-%m-%d") for d in [
    "2022-02-24","2022-02-28","2022-03-08","2022-03-16","2022-03-22",
    "2022-04-05","2022-04-12","2022-04-20","2022-04-29","2022-05-03",
    "2022-05-11","2022-05-17","2022-05-29",
]]
BBOX = (-120.51, 34.43, -119.53, 34.76)
DL_DIR = Path("/home/renatob/data/FluoData1/aviris_dangermond/review/test_field_data/earthdata_traits")
DL_DIR.mkdir(exist_ok=True)

earthaccess.login(strategy="netrc")

for i, flight_date in enumerate(AVIRIS_DATES):
    d0 = (flight_date - timedelta(days=1)).strftime("%Y-%m-%d")
    d1 = (flight_date + timedelta(days=1)).strftime("%Y-%m-%d")
    gr = earthaccess.search_data(short_name="SHIFT_AVNG_Plant_Trait_Mosaics", version="1", temporal=(d0,d1), bounding_box=BBOX)
    if not gr:
        gr = earthaccess.search_data(doi="10.3334/ORNLDAAC/2453", temporal=(d0,d1), bounding_box=BBOX)
    if gr:
        files = earthaccess.download(gr, str(DL_DIR))
        print(f"time_{i:02d} ({flight_date.date()}): {len(files)} file(s)", flush=True)
    else:
        print(f"time_{i:02d} ({flight_date.date()}): not found", flush=True)

print("Done.")
