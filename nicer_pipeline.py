import os
import sys
sys.path.append('/Volumes/Samsung_T5/saturnx')
sys.path.append('/Users/xizg0003/AstroBoy/nicer_backdir')
import logging
import datetime
import pathlib
from astropy.io import fits

import tkinter as tk
from tkinter import filedialog

from nicergof.bkg import bkg_estimator as be

from saturnx.utils.my_functions import *
from saturnx.utils.xray import *
from saturnx.utils.nicer_functions import *
from saturnx.utils.my_logging import *

'''
ARGUMENTS
---------
obs_dir: directory containing all the obs IDs
first_obs_ID: first obs_ID to analyse
last_obs_ID: last obs_ID to analyse
single_obs_ID: obs_ID to analyse
obs_ID_list: obs_IDs separated by coma
drama: flag for stopping the script AS SOON AS 
        something goes wrong (default is False) 

TODO
----
- Add log_in argument, the script will look at the value (full path of
    a log file) and it will automatically read the options from there

HISTORY
-------
2021 10 20, Stefano Rapisarda (Uppsala), creation date
    A rennovated version of the original script writting in Shanghai
    taking into account the more detailed data analysis description
    recently provided by NASA people together with the latest
    version of nicerdas
'''

# Recording date and time
# ---------------------------------------------------------------------
now = datetime.now()
date = ('%d_%d_%d') % (now.day,now.month,now.year)
time = ('%d_%d') % (now.hour,now.minute)
# ---------------------------------------------------------------------


# Reading user arguments
# (arguments should be specified in the form arg=value.
# If no value is specified, the value of that argument will be True.
# There is always an argument called name containing the name of the
# script itseld)
# ---------------------------------------------------------------------
default_args = {'obs_dir':False,
                'first_obs_id':False,
                'last_obs_id':False,
                'single_obs_id':False,
                'obs_id_list':False,
                'bary_corr':True,
                'exc_det':True,
                'det_list':'-11,-14,-20,-22,-34,-60',
                'comp_spec':True,
                'comp_arf_rmf':True,
                'bin_spec':True,
                'nibackgen':True,
                'bkg_est':True,
                'zip_evt':True,
                'drama':False,
                'debug':False,
                'target':False,
                'ra':False,
                'dec':False}
arg_dict = read_args(default_args)
# ---------------------------------------------------------------------


# Defining working directories
# ---------------------------------------------------------------------
home = os.getcwd()
if not default_args['obs_dir']:
    root = tk.Tk()
    root.withdraw()
    obs_dir = filedialog.askdirectory(initialdir=home,
        title='Select the folder containing obs_ID directories')
else:
    obs_dir = default_args['obs_dir']
    if type(obs_dir) == str:
        obs_dir = pathlib.Path(obs_dir)
# ---------------------------------------------------------------------


# Initializing logger
# ---------------------------------------------------------------------
log_name = os.path.basename(__file__).replace('.py','')+\
    '_D{}_T{}'.format(date,time)
make_logger(log_name,outdir=obs_dir)

mylogging = LoggingWrapper() 
# ---------------------------------------------------------------------


# Retrieving and printing some info about the dataset
# ---------------------------------------------------------------------
mylogging.info('*'*72)
mylogging.info('*'*24+' Running nicer_pipeline '+'*'*24)
mylogging.info('*'*72+'\n')

mylogging.info('-'*10+'Running options:'+'-'*10)
for key,item in arg_dict.items():
    mylogging.info(f'{key}: {item}')
mylogging.info('-'*36+'\n')

# Reading all digits folders (these are supposed to be NICER IDs)
obs_list = list_items(obs_dir,digits=True)

# Reading target name and position from unfiltered event file
# *********************************************************************
uf_file = list_items(obs_dir/obs_list[0]/'xti/event_uf',itype='file',
    include_and=['_0mpu0_uf'])
if not uf_file[0].is_file():
    mylogging.error('Could not find a cleaned event file to read the target name and position')
    sys.exit()
with fits.open(uf_file[0]) as hdul:
    arg_dict['target'] = hdul[0].header['OBJECT']
    arg_dict['ra'] = str(hdul[0].header['RA_OBJ'])
    arg_dict['dec'] = str(hdul[0].header['DEC_OBJ'])

# Target name
mylogging.info(f"Target: {arg_dict['target']}")

# Right ascension
div_ar = arg_dict['ra'].split()
if len(div_ar) == 3:
    ra = str( float(div_ar[0]) + float(div_ar[1])/60 + \
        float(div_ar[2])/3600 )
elif len(div_ar) == 1:
    ra = div_ar[0]
else:
    mylogging.error('User specified RA format should be 00 00 00.00')
    sys.exit()
mylogging.info(f"AR: {ra} ({arg_dict['ra']})")

# Declination
div_dec = arg_dict['dec'].split()
if len(div_dec) == 3:
    if div_dec[0][0] in ['-','+']:
        dec = div_dec[0][0] + str( float(div_dec[0][1:]) + \
            float(div_dec[1])/60 + float(div_dec[2])/3600 )
elif len(div_dec) == 1:
    dec = div_dec[0]
else:
    mylogging.error('User specified DEC format should be +00 00 00.00')
    sys.exit()
mylogging.info(f"DEC: {dec} ({arg_dict['dec']})\n")

# Confirming
ans = yesno('Important!!! Is this your target?')
if not ans: sys.exit()
# *********************************************************************

mylogging.info('Working directories\n'+72*'-')
mylogging.info('Script running in: {}'.format(home))
mylogging.info('Data dir: {}\n'.format(obs_dir)+72*'-'+'\n')

# Listing obs ID folders in the data folder
n_obs = len(obs_list)
mylogging.info('Total number of observations: {}\n'.\
    format(n_obs))
# -----------------------------------------------------------------------


# Moving to the data folder 
os.chdir(obs_dir) 


mylogging.info('Start processing observations\n'+'='*72+'\n')
# Loop over observations
# LOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOP - START_OBS
for_all = False
for o in range(n_obs):

    obs_id_dir = obs_list[o]
    obs_id = str(obs_id_dir.name)

    # Recording size at the beginning of the process
    start_size = obs_id_dir.stat().st_size

    # Selecting obs. IDs according to user arguments
    # (This assumes that NICER obsID are purely made of integers)
    # -----------------------------------------------------------------
    if not arg_dict['first_obs_id'] is False:
        if int(obs_id) < int(arg_dict['first_obs_id']): continue

    if not arg_dict['last_obs_id'] is False:
        if int(obs_id) > int(arg_dict['last_obs_id']): continue

    if not arg_dict['single_obs_id'] is False:
        if obs_id != arg_dict['single_obs_id']: continue

    if not arg_dict['obs_id_list'] is False:
        if not obs_id in arg_dict['obs_id_list']: continue
    # -----------------------------------------------------------------


    mylogging.info('Processing obs ID: {}, {}/{}\n'.\
        format(obs_id,o+1,n_obs)+'-'*72)

    event_cl_dir = obs_id_dir / 'xti/event_cl'


    # Moving old files in a directory, deleting them, or moving on 
    # -----------------------------------------------------------------
    # Looking for cleaned event file
    cl_found = False
    cl_file_list = list_items(event_cl_dir,itype='file',digits=False,
        include_and='_cl.evt')
    if len(cl_file_list) != 0: 
        cl_found=True
        file_name = str(cl_file_list[0].name)
        mylogging.info('I found an event file: {}'.format(file_name))
    
    run_nicerl2_flag = True
    if cl_found:
        
        # Asking options
        if not for_all:
            opts = ['Delete the files',
                'Move it to the graveyard',
                'Skip observations with already created cleaned event files',
                'Procede with the task, but do not run nicerl2 again']
            opt = ask_opts(opts)

        if opt == 0:
            mylogging.info('Removing old files')
            os.system('rm -f {}/*'.format(event_cl_dir))
     
        elif opt == 1:
            if not for_all:
                ans = yesno('Do you want to give a specific name to the folder?')
                if ans:
                    folder = input('===> ').strip()
                else:
                    folder = 'D{}_T{}'.format(date,time)

            graveyard = event_cl_dir/'graveyard'
            os.makedirs(graveyard, exist_ok=True)
            new_dir = graveyard/folder
            os.makedirs(new_dir, exist_ok=True)
            mylogging.info('Moving old cl and ufa event files to the graveyard')
            os.system('mv {}/* {}'.format(event_cl_dir,new_dir))

        elif opt ==2:
            mylogging.info('Skipping this observation\n'+'-'*72+'\n')
            continue

        elif opt ==3:
            run_nicerl2_flag = False
            cl_file = cl_file_list[0]
            
        # Keeping the decision for all the obs IDs (or not)
        if not for_all:
            for_all = yesno('Do you want to keep your choice for all the obs_ID?')
    # -----------------------------------------------------------------


    # Checking that nicerl2 has everything it needs
    # -----------------------------------------------------------------
    nicerl2_input_flag = check_nicerl2_io(obs_id_dir)
    if not nicerl2_input_flag:
        mylogging.info('Skipping this observation\n'+'-'*72+'\n')
        continue
    else:
        mylogging.info('nicerl2 has everything it needs')
    # -----------------------------------------------------------------


    # Picking some files that will be useful later on
    # -----------------------------------------------------------------
    # Make filter file (also unzipping if zipped)
    auxil_dir = obs_id_dir/'auxil'
    mkf_files = list_items(auxil_dir,itype='file',include_and=['.mkf'],
        exclude_or=['mkf2','mkf3'])
    if mkf_files[0].suffix == '.gz':
        mylogging.info('Unzipping .mkf file')
        os.system('gunzip {}'.format(mkf_files[0]))
        mkf_file = auxil_dir/mkf_files[0].stem
    else:
        mkf_file = mkf_files[0]

    # Attitude file
    att_file = list_items(auxil_dir,itype='file',include_and=['.att'])[0]
    # -----------------------------------------------------------------


    # Running nicerl2
    # -----------------------------------------------------------------
    if run_nicerl2_flag:
        nicerl2_log_name = 'nicerl2_D{}_T{}.log'.format(date,time)
        cl_file = run_nicerl2(obs_id_dir=obs_id,
            log_file=str(obs_id_dir/'xti/event_cl'/nicerl2_log_name))
        if not cl_file:
            mylogging.info('Skipping this observation\n'+'-'*72+'\n')
            continue 
    # -----------------------------------------------------------------


    # Excluding noisy detectors
    # -----------------------------------------------------------------
    if arg_dict['exc_det']:

        log_file = event_cl_dir/'nifpmsel.log'
        cl_bdc_file = run_nifpmsel(cl_file,arg_dict['det_list'],log_file=log_file)

        # Checking if the file was created
        if not cl_bdc_file and arg_dict['drama']:
            mylogging.info('DRAMA! Exiting system\n'+'-'*72+'\n'+72*'=')
            os.chdir(home)
            sys.exit()
    # -----------------------------------------------------------------


    # Applying barycentric correction
    # -----------------------------------------------------------------
    if arg_dict['bary_corr']:

        # Retrieving orbit file
        orb_file = list_items(auxil_dir,itype='file',include_and=['.orb'])[0]

        log_file1 = event_cl_dir/'barycorr.log'
        bc_cl_file = run_barycorr(cl_file,orb_file,log_file=log_file1)
        if not bc_cl_file and arg_dict['drama']:
            mylogging.info('DRAMA! Exiting system\n'+'-'*72+'\n'+72*'=')
            os.chdir(home)
            sys.exit()

        log_file2 = event_cl_dir/'barycorr_bdc.log'
        bc_cl_file = run_barycorr(cl_bdc_file,orb_file,log_file=log_file2)
        if not bc_cl_file and arg_dict['drama']:
            mylogging.info('DRAMA! Exiting system\n'+'-'*72+'\n'+72*'=')
            os.chdir(home)
            sys.exit()        
    # -----------------------------------------------------------------


    # Computing energy spectrum on bad detector corrected cleaned event
    # file
    # -----------------------------------------------------------------
    spectrum_name = f'{obs_id}_spectrum_bdc.pha'
    spectrum_file = event_cl_dir/spectrum_name
    spectrum_flag = False
    if arg_dict['comp_spec']:

        if cl_bdc_file.exists():
            
            # Cleaning up
            if spectrum_file.exists():
                mylogging.info('A spectrum_bdc.pha file already exists, deleting it')
                os.system('rm -f {}'.format(spectrum_file))
            
            mylogging.info('Running xselect (computing energy spectrum) ...')
            spectrum_flag = run_xselect(cl_bdc_file,outfile=spectrum_name)
            if not spectrum_flag:
                mylogging.error('Something went wrong running xselect, check log file')
                mylogging.info('Skipping this observation\n'+'-'*72+'\n')
                continue 
            else:
                mylogging.info('... and done!\n')

            # Checking output files
            if not spectrum_file.exists():
                mylogging.error('Energy spectrum not created, something went wrong')
                spectrum_flag = False
        else:
            mylogging.info('Cleaned bad-detector-corrected file does not exist\n'+\
            'Cannot compute energy spectrum')
    
    else:
        if spectrum_file.exists():
            mylogging.info('A spectrum_bdc.pha file already exists')
            spectrum_flag = True
    # -----------------------------------------------------------------

    # !!! ARF and RMF do not like obsID name in their name

    # Computing ARF file
    # -----------------------------------------------------------------
    if arg_dict['comp_arf_rmf']:
        arf_flag = True
        if spectrum_flag:

            to_remove = arg_dict['det_list']

            # It seems that nicerarf does not like a name that already exists
            # (no matter the extension)
            arf_file = event_cl_dir/f'arf_bdc.arf'
            lis_file = event_cl_dir/f'arf_bdc.lis'

            # Even if the input cl file is not bdc (bad detector 
            # corrected), I specify the same bad detector list used in
            # the procedure to compute bdc cleaned event files.
            # That is why the name of this arf file has _bdc in it.
            mylogging.info('Computing ARF file ...')
            command = f'nicerarf infile={spectrum_file} ra={ra} dec={dec} \
                attfile={mkf_file} selfile={cl_file} \
                detlist="launch,{to_remove}" \
                outfile={arf_file} outwtfile={lis_file} clobber=yes > nicerarf.log'
            if arg_dict['debug']: print(command)
            os.system (command)
            
            mylogging.info('... and done!\n')

            # Checking if output files exist
            if (not arf_file.exists()) or (not lis_file.exists()):
                mylogging.error('Either the arf or lis files have not been created\n')
                arf_flag = False
        else:
            arf_flag = True
    # -----------------------------------------------------------------


    # Computing RMF file
    # -----------------------------------------------------------------
    if arg_dict['comp_arf_rmf']:
        rmf_flag = True
        if spectrum_flag and arf_flag:

            to_remove = arg_dict['det_list']

            rmf_file = event_cl_dir/f'rmf_bdc.rmf'
            mylogging.info('Computing RMF file ...')
            command = f'nicerrmf infile={spectrum_file} mkfile={mkf_file} \
                outfile={rmf_file} detlist="launch,{to_remove}" clobber=yes > nicerrmf.log'
            os.system(command)
            mylogging.info('... and done!\n')

            # Checking if output files exist
            if not rmf_file.exists():
                mylogging.error('The rmf file has not been created\n')  
                rmf_flag = False   
        else:
            rmf_flag = False        
    # -----------------------------------------------------------------


    # Binning energy spectrum
    # -----------------------------------------------------------------
    if spectrum_flag and rmf_flag and arg_dict['bin_spec']:
        bins = 25
        grp_spectrum_name = '{}_spectrum_bdc_grp{}.pha'.format(obs_id,bins)
        grp_spectrum = event_cl_dir/grp_spectrum_name
        command = f'ftgrouppha {spectrum_file} {grp_spectrum} grouptype=optmin \
            groupscale={bins} respfile={rmf_file} clobber=yes'
        mylogging.info('Binning energy spectrum ({} bins) ...'.format(bins))
        os.system(command)
        mylogging.info('... and done!\n')

        # Checking if file exists
        if not grp_spectrum.exists():
            mylogging.error('I could not group the energy spectrum')
    # -----------------------------------------------------------------


    # Computing background via nibackgen3C50
    # -----------------------------------------------------------------
    flag_nibackgen = True
    if arg_dict['nibackgen']:
        out_file1 = event_cl_dir/f'{obs_id}_spectrum_bdc_3C50_tot'
        out_file2 = event_cl_dir/f'{obs_id}_spectrum_bdc_3C50_bkg'
        bkg_model_dir = '/Users/xizg0003/AstroBoy/nicer_backdir/bg_models_3C50'
        log_name = event_cl_dir/'nibackgen3C50.log'
        command = f'nibackgen3C50 \
            rootdir={obs_dir}/ obsid={obs_id} \
            bkgidxdir={bkg_model_dir} bkglibdir={bkg_model_dir} \
            totspec={out_file1} bkgspec={out_file2} \
            gainepoch=2020 clobber=YES > {log_name}'
        mylogging.info('Running nibackgen3C50 (computing energy spectrum and background) ...')
        os.system(command)
        mylogging.info('... and done!\n')

        # Checking if output was created 
        if not out_file1.with_suffix('.pi').exists():
            mylogging.error('Total spectrum (nibackgen3C50) not created')
            flag_nibackgen = False
        if not out_file2.with_suffix('.pi').exists():
            mylogging.error('Background spectrum (nibackgen3C50) not created')    
            flag_nibackgen = False
    else:
        flag_nibackgen = False
    # -----------------------------------------------------------------


    # Checking that the total spectrum exposure is the same of the sum
    # of the GTIs
    # -----------------------------------------------------------------
    if flag_nibackgen:
        check_file1 = event_cl_dir/'check_nibackgen1.txt'
        cmd1 = f"gtisum {out_file1.with_suffix('.pi')} > {check_file1}"
        os.system(cmd1)
        with open(check_file1,'r') as infile:
            lines = infile.readlines()
        time1 = float(lines[6].strip().split()[6])

        check_file2 = event_cl_dir/'check_nibackgen2.txt'
        cmd2 = f"ftlist {out_file2.with_suffix('.pi')} K | grep EXPOSURE > {check_file2}"
        os.system(cmd2)
        with open(check_file2,'r') as infile:
            lines = infile.readlines()
        time2 = float(lines[0].split()[1])

        mylogging.info(f'Time check1: {time1}, Time check2: {time2}')
        if time1 == time2:
            mylogging.info('nibackgen3C50 output seems ok\n')
        else:
            mylogging.warning('Total spectrum and GTI exposure seem different')
            mylogging.info('Changing exposure\n')
            os.system(f"fthedit {out_file1.with_suffix('.pi')} EXPOSURE a {time1}")
    # -----------------------------------------------------------------

    # Computing background via nicer_bkg_estimator
    # -----------------------------------------------------------------
    if arg_dict['bkg_est']:
        # Running niprefilter2
        mkf2_file = auxil_dir/f'ni{obs_id}.mkf2'
        command = f'niprefilter2 indir=./{obs_id} \
            infile=./{obs_id}/auxil/ni{obs_id}.mkf \
            outfile={mkf2_file} clobber=YES'
        mylogging.info('Running niprefilter2 ...')
        os.system(command)
        mylogging.info('... and done!\n')

        if not mkf2_file.exists():
            mylogging.error('mkf2 file not created (niprefilter2), cannot run nicer_bkg_estimator')
        else:
            mylogging.info('Computing mkf3 file ...')
            mkf3_flag = be.add_kp(str(mkf2_file))
            mylogging.info('... and done!\n')

            mkf3_files = list_items(auxil_dir,itype='file',include_and=['.mkf3'])        
            if len(mkf3_files) != 1:
                mylogging.error('I could not fint an .mkf3 file')
            else:
                bevt = '/Users/xizg0003/AstroBoy/caldb/data/nicer/xti/pcf/30nov18targskc_enhanced.evt'
                mylogging.info('Running mk_bkg_spec_evt ...')
                bkg_chan, bkgspectot, btotexpo = be.mk_bkg_spec_evt(
                    str(event_cl_dir/spectrum_name),
                    mkf3file=str(mkf3_files[0]),bevt=bevt)
                mylogging.info('... and done!\n')

                # Checking if the background was created
                bkg_spectrum_name = spectrum_name.replace('.pha','_bkg.pha')
                bkg_spectrum_file = event_cl_dir/bkg_spectrum_name

                if not bkg_spectrum_file.exists():
                    mylogging.error('Background spectrum (mk_bkg_spec_evt) not created')
    # -----------------------------------------------------------------
    
    # Zipping event files
    # ----------------------------------------------------------------------- 
    # -f means "force", this will override already existing files   
    if arg_dict['zip_evt']:
        mylogging.info('Zipping files ...')
        evt_files = list_items(event_cl_dir,itype='file',ext='evt',
            exclude_and=['gz'])
        mylogging.info(f'I found {len(evt_files)} evt files to zip')
        for evt_file in evt_files: 
            mylogging.info('Zipping {}'.format(evt_file.name))
            os.system('gzip -f {}'.format(evt_file))    
        mylogging.info('... and done!\n')       
    # -----------------------------------------------------------------------

    # Checking size at the end of all the processes and printing stats
    # -----------------------------------------------------------------------
    end_size = obs_id_dir.stat().st_size
    mylogging.info(f'Size at the beginning of operations: {start_size/1000000} Mb')
    mylogging.info(f'Size at the end of operations: {end_size/1000000} Mb')
    mylogging.info(f'"Inflation": {round((end_size-start_size)/start_size * 100)}%\n')
    # -----------------------------------------------------------------------

    # Going home
    os.chdir(obs_dir)

    mylogging.info(f'Everything done for obs. ID {obs_id}!\n'+'-'*72+'\n')
# LOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOP - END_OBS

mylogging.info('All is done! Time to go home\n'+'='*72)
os.chdir(home)