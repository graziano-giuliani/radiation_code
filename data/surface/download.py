#!/usr/bin/env python

import os
import sys
import time
import calendar
import cdsapi
from netCDF4 import Dataset
from cdo import Cdo

def ldomon(y,m):
    return calendar.monthrange(y,m)[1]

tlist = "/".join([ "{0:02d}:00:00".format(x) for x in range(0,24) ])
gbase = "../"
gdes = gbase + "/era5_original.grid"

c = cdsapi.Client()
cdo = Cdo( )

varlist = {
            "skt"  : { "code"    : "235.128",
                       "command" : "",
                       "rename"  : "SKT,skt", },
            "aluvp": { "code"    : "15.128",
                       "command" : "",
                       "rename"  : "var15,aluvp", },
            "aluvd": { "code"    : "16.128",
                       "command" : "",
                       "rename"  : "var16,aluvd", },
            "alnip": { "code"    : "17.128",
                       "command" : "",
                       "rename"  : "var17,alnip", },
            "alnid": { "code"    : "18.128",
                       "command" : "",
                       "rename"  : "var18,alnid", },
          }

year_start = 1950
year_end = 1951

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
                        "levtype": "sfc",
                        "param": varlist[var]["code"],
                        "stream": "oper",
                        "time": tlist,
                        "type": "an"
                      }, gname)
                    print(gname+" downloaded.")
                    try:
                        cdo.copy(input = varlist[var]["command"]+" "+gname,
                                 output = ncname,
                                 options="-L -f nc4 -z zip_9 -t ecmwf")
                        os.unlink(gname)
                        try:
                            if varlist[var]["rename"]:
                                vfrom,vto = varlist[var]["rename"].split(',')
                                ds = Dataset(ncname,mode='a')
                                ds.renameVariable(vfrom,vto)
                                ds.close( )
                        except:
                            pass
                        print(ncname+" created.")
                    except:
                        print('Error creating '+ncname)
                        time.sleep(2)
                except:
                    print('Error download '+gname)
                    time.sleep(2)
            else:
                print(ncname+" already on disk.")
