#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#------------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: July, 1st, 2016
#------------------------------------

import yaml
import os, sys, subprocess

def updateConfig(config):
    loadRef(config)
    config["config_file"] = "config.yaml" # trick to force rules on config change	
    for k in ["RPKM_threshold","min_num_samples_expressing_at_threshold", 
                "numgenes_plots","num_kmeans_clust","filter_mirna"]:
        config[k] = str(config[k])

    config = _addExecPaths(config)
    return config

#------------------------------------------------------------------------------
# CHIPS-like refs file
#-----------------------------------------------------------------------------
def loadRef(config):
    """Adds the static reference paths found in config['ref']
    NOTE: if the elm is already defined, then we DO NOT clobber the value
    """
    f = open(config['ref'])
    ref_info = yaml.safe_load(f)
    f.close()
    #print(ref_info[config['assembly']])
    for (k,v) in ref_info[config['assembly']].items():
        #NO CLOBBERING what is user-defined!
        if k not in config:
            config[k] = v

def _addExecPaths(config):
    conda_root = subprocess.check_output('conda info --root',shell=True).decode('utf-8').strip()
    conda_path = os.path.join(conda_root, 'pkgs')
    #NEED the following when invoking python2 (to set proper PYTHONPATH)
    config["python2_pythonpath"] = os.path.join(conda_root, 'envs', 'viper_py2', 'lib', 'python2.7', 'site-packages')
    
    if not "python2" in config or not config["python2"]:
        config["python2"] = os.path.join(conda_root, 'envs', 'viper_py2', 'bin', 'python2.7')

    if not "rseqc_path" in config or not config["rseqc_path"]:
        config["rseqc_path"] = os.path.join(conda_root, 'envs', 'viper_py2', 'bin')

    if not "picard_path" in config or not config["picard_path"]:
        config["picard_path"] = 'picard'

    if not "varscan_path" in config or not config["varscan_path"]:
        config["varscan_path"] = 'varscan'
   
    if "analysis_token" in config and config["analysis_token"]:
        config["token"] = config["analysis_token"]
    else:
        config["token"] = "summary_reports"

    if not "trust_path" in config or not config["trust_path"]:
        config["trust_path"] = os.path.join(conda_root, 'envs', 'viper_py2', 'bin', 'trust')

    if not "optitype_path" in config or not config["optitype_path"]:
        config["optitype_path"] = os.path.join(conda_root, 'envs', 'viper_py2', 'bin')

    return config

