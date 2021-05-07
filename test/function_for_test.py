#!/usr/bin/env python3

"""Provide functions for test."""

import os

def check_exist(t, min_size=0):
    """ Make sure target file (t) exist and not empty """
    assert os.path.exists(t) and os.path.getsize(t) > min_size

def check_exists(targets, min_size=0):
    """ Make sure targets files exist and not empty """
    for t in targets:
        check_exist(t, min_size)
