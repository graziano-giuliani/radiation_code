#!/usr/bin/env python

import os
import sys
import time
import calendar
import cdsapi
from cdo import Cdo

def ldomon(y,m):
    return calendar.monthrange(y,m)[1]

llist = "/".join([ str(x) for x in range(1,138) ])
gbase = "/home/clima-archive4-b/ggiulian/era5_gaussian_reduced"
gdes = gbase + "/era5_original.grid"

c = cdsapi.Client()
cdo = Cdo( )

varlist = {
            "cc"   : { "code"    : "248",
                       "command" : "", },
            "o3"   : { "code"    : "203",
                       "command" : "", },
            "clwc" : { "code"    : "247",
                       "command" : "", },
            "ciwc" : { "code"    : "246",
                       "command" : "", },
            "q"    : { "code"    : "133",
                       "command" : "", },
            "t"    : { "code"    : "130",
                       "command" : "-remapbil,"+gdes+" -sp2gp ", },
            "lnps" : { "code"    : "152",
                       "command" : "-remapbil,"+gdes+" -sp2gp ", },
          }

year_start = 1950
year_end = 1953

for year in range(year_start,year_end):
    yy = '%04d' % year
    for month in range(1,13):
        mm = '%02d' % month
        lm = '%02d' % ldomon(year,month)
        yymm = yy+"-"+mm+"-"
        for var in varlist.keys( ):
            gname = "_".join([var,yy,mm])+".grib"
            ncname = "_".join([var,yy,mm])+".nc"
            if not os.path.exists(ncname):
                drange = yymm+"01/to/"+yymm+lm
                print('Downloading '+gname)
                try:
                    c.retrieve("reanalysis-era5-complete",
                      {
                        "class": "ea",
                        "date": drange,
                        "expver": "1",
                        "levelist": llist,
                        "levtype": "ml",
                        "param": varlist[var]["code"],
                        "stream": "oper",
                        "time": "00:00:00/06:00:00/12:00:00/18:00:00",
                        "type": "an"
                      }, gname)
                    print(gname+" downloaded.")
                    try:
                        cdo.copy(input = varlist[var]["command"]+" "+gname,
                                 output = ncname,
                                 options="-L -f nc4 -z zip_9 -t ecmwf")
                        os.unlink(gname)
                        print(ncname+" created.")
                    except:
                        print('Error creating '+ncname)
                        time.sleep(2)
                except:
                    print('Error download '+gname)
                    time.sleep(2)
            else:
                print(ncname+" already on disk.")
