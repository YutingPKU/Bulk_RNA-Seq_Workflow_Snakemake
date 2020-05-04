#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#-------------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: July, 1st, 2016
#-------------------------------------

import pandas as pd
from collections import defaultdict

def updateMeta(config):
    _sanity_checks(config)
    metadata = pd.read_table(config['metasheet'], index_col=0, sep=',', comment='#')
    config["comparisons"] = [c[5:] for c in metadata.columns if c.startswith("comp_")]
    config["comps"] = _get_comp_info(metadata)
    config["metacols"] = [c for c in metadata.columns if c.lower()[:4] != 'comp']
    config["file_info"] = { sampleName : config["samples"][sampleName] for sampleName in metadata.index }
    config["ordered_sample_list"] = metadata.index
    return config


def _sanity_checks(config):
    #metasheet pre-parser: converts dos2unix, catches invalid chars
    _invalid_map = {'\r':'\n', '(':'.', ')':'.', ' ':'_', '/':'.', '$':''}
    _meta_f = open(config['metasheet'])
    _meta = _meta_f.read()
    _meta_f.close()

    _tmp = _meta.replace('\r\n','\n')
    #check other invalids
    for k in _invalid_map.keys():
        if k in _tmp:
            _tmp = _tmp.replace(k, _invalid_map[k])

    #did the contents change?--rewrite the metafile
    if _meta != _tmp:
        #print('converting')
        _meta_f = open(config['metasheet'], 'w')
        _meta_f.write(_tmp)
        _meta_f.close()


def _get_comp_info(meta_info):
    comps_info = defaultdict(dict)
    for comp in meta_info.columns:
        if comp[:5] == 'comp_':
            comps_info[comp[5:]]['control'] = meta_info[meta_info[comp] == 1].index
            comps_info[comp[5:]]['treat'] = meta_info[meta_info[comp] == 2].index
    return comps_info

