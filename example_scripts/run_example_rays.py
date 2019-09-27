
import numpy as np
import subprocess
import os
import datetime
import xflib  # Fortran xform-double library (coordinate transforms)
# xf = xflib.xflib(lib_path=os.path.join(os.getcwd(),'libxformd.so'))



# ------------------Constants --------------------------
R_E = 6371

# -------------- Simulation parameters -----------------
t_max = 5     # Maximum duration in seconds

dt0 = 1e-3      # Initial timestep in seconds
dtmax = 0.05    # Maximum allowable timestep in seconds
root = 2        # Which root of the Appleton-Hartree equation
                # (1 = negative, 2 = positive)
                # (2=whistler in magnetosphere)
fixedstep = 0   # Don't use fixed step sizes, that's a bad idea.
maxerr = 1.0e-3 # Error bound for adaptive timestepping
maxsteps = 2e5  # Max number of timesteps (abort if reached)
use_IGRF = 0    # Magnetic field model (1 for IGRF, 0 for dipole)
use_tsyg = 0    # Use the Tsyganenko magnetic field model corrections
minalt   = (R_E + 800)*1e3 # cutoff altitude in meters

mag_dump = False # Use magnetic dipole coordinates for the dump, instead of SM.

# -------------- Starting ray parameters -------------
inp_lat = 45;   # launch latitude (geomagnetic)
inp_lon = 256.6;    # launch longitude (geomagnetic)
freq = 440;     # ray frequency in Hz. How about concert A?
inp_alt = (R_E + 1000)*1e3 # Launch altitude in meters
launch_direction = 'field-aligned'

# -------------- Environmental parameters -------------


yearday = '2010001'		#YYYYDDD
milliseconds_day = 0;	# milliseconds into the day

Kp    = 2 # 0.7
AE    = 1.6
Pdyn  = 4
Dst   = 1.0
ByIMF = 0.0
BzIMF = -5
# Tsyganenko correction parameters
W = [0.132,    0.303,    0.083,    0.070,    0.211,    0.308 ]  # Doesn't matter if we're not using Tsyg

modes_to_do = [1,6]

# -------------- Set up the output directory ----------
project_root = os.getcwd();
ray_out_dir = os.path.join(project_root, "test_outputs");

print("output directory:",ray_out_dir)
if not os.path.exists(ray_out_dir):
	os.mkdir(ray_out_dir)


# GCPM model and damping code needs to be run in the same directory 
# as the binary file (and all the misc data files)
cwd = os.getcwd();
os.chdir('../bin')

for mode in modes_to_do:
    if mode == 1:
        # -------------- Test the Ngo model ------------------

        configfile = os.path.join(project_root,"ngo_kp2.in")
        modelnum = 1
        ray_inpfile  = os.path.join(project_root,"ray_inpfile.txt")
        ray_outfile  = os.path.join(ray_out_dir,'example_ray_ngo.ray')
        damp_outfile = os.path.join(ray_out_dir,'example_ray_ngo.damp')

        # ------ Write the ray input file for the Ngo model ---
        f = open(ray_inpfile,'w')
        # pos_mag = [1, inp_lat, inp_lon]   # alt, lat, lon, in magnetic dipole coords
        # dd= datetime.datetime(2010,6,21,16,0,0)
        # pos_SM = xf.rllmag2sm(pos_mag, dd)

        # print("pos_mag",pos_mag)
        # print("pos_SM",pos_SM)
        # print("and back: ", xf.sm2rllmag(pos_SM, dd))

        # pos_flat = [np.sqrt(2)/2, 0, np.sqrt(2)/2]
        # print("pos_flat",xf.sm2rllmag(pos_flat, dd))
        # pos_SM  = xf.rllmag2sm(pos_mag, datetime.datetime(2010,1,1,0,0,0))
        # print(pos_SM)

        pos0 = np.array([np.sqrt(2)/2, 0, np.sqrt(2)/2])*inp_alt

        if (launch_direction is 'up'):
            dir0 = pos_SM/np.linalg.norm(pos0)    # radial outward
        elif (launch_direction is 'field-aligned'):
            dir0 = np.zeros(3)                # Field aligned (set in raytracer)


        w0   = freq*2.0*np.pi
        f.write('%1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e\n'%(pos0[0], pos0[1], pos0[2], dir0[0], dir0[1], dir0[2], w0))
        f.close()

        ngo_cmd= './raytracer --outputper=%d --dt0=%g --dtmax=%g'%(1, dt0, dtmax) + \
                 ' --tmax=%g --root=%d --fixedstep=%d --maxerr=%g'%(t_max, root, fixedstep, maxerr) + \
                 ' --maxsteps=%d --minalt=%d --inputraysfile="%s" --outputfile="%s"'%( maxsteps, minalt, ray_inpfile, ray_outfile) + \
                 ' --modelnum=%d --yearday=%s --milliseconds_day=%d'%(modelnum, yearday, milliseconds_day) + \
                 ' --use_tsyganenko=%d --use_igrf=%d --tsyganenko_Pdyn=%g'%(use_tsyg, use_IGRF, Pdyn) + \
                 ' --tsyganenko_Dst=%g --tsyganenko_ByIMF=%g --tsyganenko_BzIMF=%g'%( Dst, ByIMF, BzIMF ) + \
                 ' --tsyganenko_W1=%g --tsyganenko_W2=%g --tsyganenko_W3=%g'%(W[0], W[1], W[2]) + \
                 ' --tsyganenko_W4=%g --tsyganenko_W5=%g --tsyganenko_W6=%g'%(W[3], W[4], W[5]) + \
        		 ' --ngo_configfile="%s" --kp=%g'%(configfile, Kp) 

        print("------- Running Ngo model -------")
        print("Command is:")
        print(ngo_cmd)
        print();

        # Run it
        os.system(ngo_cmd)

        print("------- Running Ngo damping -------")
        damp_mode = 0
        damp_cmd =  './damping --inp_file "%s" --out_file "%s" '%(ray_outfile, damp_outfile) + \
                '--Kp %g --AE %g --mode %d'%(Kp, AE, damp_mode) + \
                ' --yearday %s --msec %d'%(yearday, milliseconds_day)

        print(damp_cmd)
        os.system(damp_cmd)

    if mode == 6:
        # -------------- Test the Simplified GCPM model ------------------

        modelnum = 6
        MLT = 0
        fixed_MLT = 0 # Force the raytracer to stay in the meridonal plane?

        ray_outfile  = os.path.join(ray_out_dir,'example_ray_mode6.ray')
        damp_outfile = os.path.join(ray_out_dir,'example_ray_mode6.damp')

        mode6_cmd= './raytracer --outputper=%d --dt0=%g --dtmax=%g'%(1, dt0, dtmax) + \
                 ' --tmax=%g --root=%d --fixedstep=%d --maxerr=%g'%(t_max, root, fixedstep, maxerr) + \
                 ' --maxsteps=%d --minalt=%d --inputraysfile="%s" --outputfile="%s"'%( maxsteps, minalt, ray_inpfile, ray_outfile) + \
                 ' --modelnum=%d --yearday=%s --milliseconds_day=%d'%(modelnum, yearday, milliseconds_day) + \
                 ' --use_tsyganenko=%d --use_igrf=%d --tsyganenko_Pdyn=%g'%(use_tsyg, use_IGRF, Pdyn) + \
                 ' --tsyganenko_Dst=%g --tsyganenko_ByIMF=%g --tsyganenko_BzIMF=%g'%( Dst, ByIMF, BzIMF ) + \
                 ' --tsyganenko_W1=%g --tsyganenko_W2=%g --tsyganenko_W3=%g'%(W[0], W[1], W[2]) + \
                 ' --tsyganenko_W4=%g --tsyganenko_W5=%g --tsyganenko_W6=%g'%(W[3], W[4], W[5]) + \
                 ' --MLT="%g" --fixed_MLT=%g --kp=%g'%(MLT, fixed_MLT, Kp)

        print("------- Running mode 6 model -------")

        print(mode6_cmd)
        os.system(mode6_cmd)

        damp_mode = 1
        damp_cmd =  '../bin/damping --inp_file "%s" --out_file "%s" '%(ray_outfile, damp_outfile) + \
                '--Kp %g --AE %g --mode %d'%(Kp, AE, damp_mode) + \
                ' --yearday %s --msec %d'%(yearday, milliseconds_day)

        print("------- Running mode 6 damping -------")

        print(damp_cmd)
        os.system(damp_cmd)
    if mode == 2:

        # -------------- Test the full GCPM model ------------------

        modelnum = 2



        ray_outfile  = os.path.join(ray_out_dir,'example_ray_gcpm.ray')
        damp_outfile = os.path.join(ray_out_dir,'example_ray_gcpm.damp')

        gcpm_cmd= './raytracer --outputper=%d --dt0=%g --dtmax=%g'%(1, dt0, dtmax) + \
                 ' --tmax=%g --root=%d --fixedstep=%d --maxerr=%g'%(t_max, root, fixedstep, maxerr) + \
                 ' --maxsteps=%d --minalt=%d --inputraysfile="%s" --outputfile="%s"'%( maxsteps, minalt, ray_inpfile, ray_outfile) + \
                 ' --modelnum=%d --yearday=%s --milliseconds_day=%d'%(modelnum, yearday, milliseconds_day) + \
                 ' --use_tsyganenko=%d --use_igrf=%d --tsyganenko_Pdyn=%g'%(use_tsyg, use_IGRF, Pdyn) + \
                 ' --tsyganenko_Dst=%g --tsyganenko_ByIMF=%g --tsyganenko_BzIMF=%g'%( Dst, ByIMF, BzIMF ) + \
                 ' --tsyganenko_W1=%g --tsyganenko_W2=%g --tsyganenko_W3=%g'%(W[0], W[1], W[2]) + \
                 ' --tsyganenko_W4=%g --tsyganenko_W5=%g --tsyganenko_W6=%g'%(W[3], W[4], W[5]) + \
                 ' --kp=%g'%(Kp)


        print("------- Running GCPM model -------")

        print(gcpm_cmd)
        os.system(gcpm_cmd)


        damp_mode = 1
        damp_cmd =  './damping --inp_file "%s" --out_file "%s" '%(ray_outfile, damp_outfile) + \
                '--Kp %g --AE %g --mode %d'%(Kp, AE, damp_mode) + \
                ' --yearday %s --msec %d'%(yearday, milliseconds_day)

        print("------- Running GCPM damping -------")

        print(damp_cmd)
        os.system(damp_cmd)


# Move back to the working directory
os.chdir(cwd)