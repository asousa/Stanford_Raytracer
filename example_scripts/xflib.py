import ctypes as ct
import datetime
import numpy as np

class xflib(object):
    ''' A wrapper class for the xform-double coordinate transformation library.
        (The fortran one Forrest used in the raytracer)

        To build this thing: You need to compile libxformd.so for your machine.
            - Move into the xform_double directory
            - 'make shared'
            - Copy the new file (libxformd.so) to wherever you'd like it to live.
    '''


    def __init__(self, lib_path='libxformd.so'):

        self.D2R = 3.141592653589793238462643/180.
        self.R2D = 180./3.141592653589793238462643
        # data types
        self.i2 = ct.c_int*2
        self.d3 = ct.c_double*3
        
        # load shared library
        ct.cdll.LoadLibrary(lib_path)
        self.xf = ct.CDLL(lib_path)

        # methods
        self.geo2sm_l = self.xf.geo_to_sm_d_
        self.sm2geo_l = self.xf.sm_to_geo_d_

        self.geo2mag_l= self.xf.geo_to_mag_d_
        self.mag2geo_l= self.xf.mag_to_geo_d_

        self.s2c_l    = self.xf.pol_to_cart_d_
        self.c2s_l    = self.xf.cart_to_pol_d_
        
        self.gse2sm_l = self.xf.gse_to_sm_d_
        self.sm2gse_l = self.xf.sm_to_gse_d_
    
    def s2c(self, x_in):
        ''' spherical to cartesian (degrees)
            x_in: rad, lat, lon
            x_out: x, y, z
        '''
        
#         print x_in
        
        lat_in = ct.c_double(x_in[1]*self.D2R)
        lon_in = ct.c_double(x_in[2]*self.D2R)
        rad_in = ct.c_double(x_in[0])

        cx_out = self.d3()

        self.s2c_l(ct.byref(lat_in), ct.byref(lon_in), ct.byref(rad_in), cx_out)
        
        return [x for x in cx_out]

        
    def c2s(self, x_in):
        ''' cartesian to spherical (degrees)
            x_in: x, y, z
            x_out: rad, lat, lon
        '''
        
        cx_in = self.d3(*x_in)

        lat = ct.c_double()
        lon = ct.c_double()
        rad = ct.c_double()
        
        self.c2s_l(cx_in, ct.byref(lat), ct.byref(lon), ct.byref(rad))
        
        return [rad.value, lat.value*self.R2D, lon.value*self.R2D]
    
    def geo2sm(self, x_in, time_in):
        ''' Geographic (cartesian) to Solar Magnetic
        '''
        # Construct yearday:
        yearday = int(1000*time_in.year + time_in.timetuple().tm_yday)
        milliseconds_day = int((time_in.second + time_in.minute*60 + time_in.hour*60*60)*1e3 + time_in.microsecond*1e-3)

        ct_in = self.i2()
        
        ct_in[0] = yearday
        ct_in[1] = milliseconds_day
        
        # print yearday
        # print milliseconds_day
        
        cx_in = self.d3(*x_in)
        cx_out = self.d3()
        self.geo2sm_l(ct_in, cx_in, cx_out)
        
        return [x for x in cx_out]
    
    def sm2geo(self, x_in, time_in):
        ''' Solar Magnetic to Geographic (cartesian) '''

        # Construct yearday:
        yearday = int(1000*time_in.year + time_in.timetuple().tm_yday)
        milliseconds_day = int((time_in.second + time_in.minute*60 + time_in.hour*60*60)*1e3 + time_in.microsecond*1e-3)

        ct_in = self.i2()
        
        ct_in[0] = yearday
        ct_in[1] = milliseconds_day
        
#         print yearday
#         print milliseconds_day
        
        cx_in = self.d3(*x_in)
        cx_out = self.d3()
        self.sm2geo_l(ct_in, cx_in, cx_out)

        return [x for x in cx_out]

    def geo2mag(self, x_in, time_in):
        ''' Geographic (cartesian) to magnetic dipole (cartesian) '''
        yearday = int(1000*time_in.year + time_in.timetuple().tm_yday)
        milliseconds_day = int((time_in.second + time_in.minute*60 + time_in.hour*60*60)*1e3 + time_in.microsecond*1e-3)

        ct_in = self.i2()
        
        ct_in[0] = yearday
        ct_in[1] = milliseconds_day
        
        cx_in = self.d3(*x_in)
        cx_out = self.d3()
        self.geo2mag_l(ct_in, cx_in, cx_out)

        return [x for x in cx_out]

    
    def mag2geo(self, x_in, time_in):
        ''' Magnetic dipole (cartesian) to geographic (cartesian) '''
        yearday = int(1000*time_in.year + time_in.timetuple().tm_yday)
        milliseconds_day = int((time_in.second + time_in.minute*60 + time_in.hour*60*60)*1e3 + time_in.microsecond*1e-3)

        ct_in = self.i2()
        
        ct_in[0] = yearday
        ct_in[1] = milliseconds_day
        
        cx_in = self.d3(*x_in)
        cx_out = self.d3()
        self.mag2geo_l(ct_in, cx_in, cx_out)

        return [x for x in cx_out]

    def rllgeo2rllmag(self, x_in, time_in):
        ''' Geographic (r, lat, lon) to Geomagnetic (r, lat, lon) '''
        xtmp = self.s2c(x_in)
        xtmp = self.geo2mag(xtmp, time_in)
        return self.c2s(xtmp)

    def rllgeo2sm(self, x_in, time_in):
        ''' geographic (radius, lat, lon) to Solar Magnetic (cartesian) '''
        xtmp = self.s2c(x_in)
        return self.geo2sm(xtmp, time_in)

    def sm2rllgeo(self, x_in, time_in):
        ''' Solar Magnetic (cartesian) geographic (radius, lat, lon) '''
        xtmp = self.sm2geo(x_in, time_in)
        return self.c2s(xtmp)
    
    def rllmag2sm(self, x_in, time_in):
        ''' magnetic dipole (radius, lat, lon) to Solar Magnetic (cartesian) '''
        xtmp = self.s2c(x_in)
        xtmp = self.mag2geo(xtmp, time_in)
        return self.geo2sm(xtmp, time_in)

    def sm2rllmag(self, x_in, time_in):
        ''' Solar Magnetic (cartesian) to magnetic dipole (radius, lat, lon) '''
        xtmp = self.sm2geo(x_in, time_in)
        xtmp = self.geo2mag(xtmp, time_in)
        return self.c2s(xtmp)
    
    def mag2sm(self, x_in, time_in):
        ''' magnetic dipole (cartesian) to Solar Magnetic (cartesian) '''
        xtmp = self.mag2geo(x_in, time_in)
        return self.geo2sm(xtmp, time_in)

    def sm2mag(self, x_in, time_in):
        ''' Solar Magnetic (cartesian) to magnetic dipole (cartesian) '''
        xtmp = self.sm2geo(x_in, time_in)
        return self.geo2mag(xtmp, time_in)

    def transform_data_sph2car(self, lat, lon, d_in):
        D2R = np.pi/180.

        M = np.zeros([3,3])
        d_out = np.zeros(3)

        theta = D2R*(90. - lat)
        phi   = D2R*lon

        st = np.sin(theta)
        sp = np.sin(phi)
        ct = np.cos(theta)
        cp = np.cos(phi)

        M[0,0] = st*cp;    M[0,1] = ct*cp;   M[0,2] = -sp;
        M[1,0] = st*sp;    M[1,1] = ct*sp;   M[1,2] = cp;
        M[2,0] = ct;       M[2,1] = -st;     M[2,2] = 0;

        d_out = np.dot(M, d_in)

        return d_out


    def gse2sm(self, x_in, time_in):
        yearday = int(1000*time_in.year + time_in.timetuple().tm_yday)
        milliseconds_day = int((time_in.second + time_in.minute*60 + time_in.hour*60*60)*1e3 + time_in.microsecond*1e-3)

        ct_in = self.i2()
        
        ct_in[0] = yearday
        ct_in[1] = milliseconds_day
        
        cx_in = self.d3(*x_in)
        cx_out = self.d3()
        self.gse2sm_l(ct_in, cx_in, cx_out)

        return [x for x in cx_out]

    def sm2gse(self, x_in, time_in):
        yearday = int(1000*time_in.year + time_in.timetuple().tm_yday)
        milliseconds_day = int((time_in.second + time_in.minute*60 + time_in.hour*60*60)*1e3 + time_in.microsecond*1e-3)

        ct_in = self.i2()
        
        ct_in[0] = yearday
        ct_in[1] = milliseconds_day
        
        cx_in = self.d3(*x_in)
        cx_out = self.d3()
        self.sm2gse_l(ct_in, cx_in, cx_out)

        return [x for x in cx_out]

        
    def lon2MLT(self, itime, lon):
        # // Input: itime, lon in geomagnetic dipole coords.
        # // Output: MLT in fractional hours
        # // Ref: "Magnetic Coordinate Systems", Laundal and Richmond
        # // Space Science Review 2016, DOI 10.1007/s11214-016-0275-y

        ut_hr = itime.hour + itime.minute/60 # /1000.0/60.0;  #// Milliseconds to fractional hours (UT)
        A1 = [1, 51.48, 0];         #// Location of Greenwich (for UT reference) 
        B1 = [0, 0, 0]; # B1[3]                         // Location of Greenwich in geomag

        self.s2c(A1);
        self.geo2mag(A1, itime);
        self.c2s(A1);

        return np.mod(ut_hr + (lon - A1[2])/15.0,  24);

    def MLT2lon(self, itime, mlt):
        # // Input: itime, mlt in fractional hours
        # // Output: longitude in geomagnetic coordinates
        # // Ref: "Magnetic Coordinate Systems", Laundal and Richmond
        # // Space Science Review 2016, DOI 10.1007/s11214-016-0275-y

        ut_hr = itime.hour + itime.minute/60 # /1000.0/60.0;  #// Milliseconds to fractional hours (UT)
        A1 = [1, 51.48, 0];         #// Location of Greenwich (for UT reference) 
        B1 = [0, 0, 0]; # B1[3]                         // Location of Greenwich in geomag

        self.s2c(A1);
        self.geo2mag(A1, itime);
        self.c2s(A1);
        
        return 15.*(mlt - ut_hr) + A1[2]

# xf = xflib(lib_path='/shared/users/asousa/WIPP/3dWIPP/python/libxformd.so')

# x_in = [1.,45,17]

# time_in = datetime.datetime(2001, 1, 1, 0, 0, 00);


# print x_in
# x_in = xf.s2c(x_in)
# print x_in

# # x_in = xf.c2s(x_in)
# # print x_in

# x_in = xf.geo2sm(x_in, time_in)
# print x_in
# x_in = xf.sm2geo(x_in, time_in)
# print x_in
# x_in = xf.geo2mag(x_in, time_in)
# print x_in
# x_in = xf.mag2geo(x_in, time_in)
# print x_in


