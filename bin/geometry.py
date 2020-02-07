#!/usr/bin/env python
import math

def d2vv(v1, v2):
    """Squared distance between two vectors"""
    dx = v1[0] - v2[0]
    dy = v1[1] - v2[1]
    dz = v1[2] - v2[2]
    return dx*dx+dy*dy+dz*dz

def dvv(v1, v2):
    """Distance between two vectors"""
    return math.sqrt(d2vv(v1,v2))

def inrad(v1, v2, r):
    """Are two points winthin r?"""
    if (abs(v1[0]-v2[0])) > r:
        return False
    elif (abs(v1[1]-v2[1])) > r:
        return False
    elif (abs(v1[2]-v2[2])) > r:
        return False
    elif dvv(v1,v2) > r:
        return False
    else:
        return True
