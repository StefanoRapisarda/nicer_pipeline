import os
import sys
sys.path.append('/Volumes/Samsung_T5/kronos')
sys.path.append('/Users/xizg0003/AstroBoy/nicer_backdir')
import logging
import datetime
import pathlib

import tkinter as tk
from tkinter import filedialog

from nicergof.bkg import bkg_estimator as be

from kronos.utils.my_functions import *
from kronos.utils.xray import *
from kronos.utils.nicer_functions import *
from kronos.utils.logging import *

default_args = {'obs_dir':False,
                'first_obs_id':False,
                'last_obs_id':False,
                'single_obs_id':False,
                'obs_id_list':False,
                'drama':False,
                'target':'Cygnus-X1',
                'ar':'19 58 21.6758193269',
                'dec':'+35 12 05.782512305'}
arg_dict = read_args(default_args)

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

obs_list = list_items(obs_dir,digits=True)
n_obs = len(obs_list)

# Moving to the data folder 
os.chdir(obs_dir) 

for o in range(n_obs):

    obs_id_dir = obs_list[o]
    obs_id = str(obs_id_dir.name)

    # Selecting obs IDs according to user arguments
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

    print('Processing obs ID: {}, {}/{}\n'.\
        format(obs_id,o+1,n_obs)+'-'*72)

    event_cl_dir = obs_id_dir / 'xti/event_cl'

    # Checking if energy spectrum and rmp exist
    spectrum_name = 'spectrum_bdc.pha'
    spectrum = event_cl_dir/spectrum_name
    rmf = event_cl_dir/'rmf_bdc.rmf'
    if not spectrum.exists():
        print('Energy spectrum does not exist, moving on')
        continue
    if not rmf.exists():
        print('Rmf file does not exist, moving on')
        continue

    bins = 25
    grp_spectrum_name = 'spectrum_bdc_grp{}.pha'.format(bins)
    grp_spectrum = event_cl_dir/grp_spectrum_name
    command = f'ftgrouppha {spectrum} {grp_spectrum} grouptype=optmin ' + \
        f'groupscale={bins} respfile={rmf}'
    print('Grouping bins in the energy spectrum')
    os.system(command)
    print('Done!\n')

