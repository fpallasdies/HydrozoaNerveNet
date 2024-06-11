#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 12:25:29 2019

@author: pallasdies
"""


from brian2 import *
import numpy as np
shutdown= 1




height = 20 * um
lengths = 175*um
area = lengths*lengths*2 + lengths*height * 4
crosssection = lengths*height
resistance = 700 * ohm * cm
channelmult = 1.5
radius = np.sqrt(lengths*height/np.pi)

Parameters={
    'C' : 1*uF/cm**2 *area,
    'EK':-90*mV,
    'ECa': 60 *mV,
    'EL' :-83.5*mV,
    'Esyn':60*mV,

    'gL': 1* 0.5 * mS/cm**2 * area *channelmult,
    'gCa': 0.9*10* mS/cm**2 * area * channelmult,
    'gK': 0.7*6 * mS/cm**2 *area * channelmult,
    'ggap':radius/(2 * resistance * lengths**2) * area *0.3  ,
    'ggapparallel': radius/(2 * resistance * lengths**2) * area *0.05 ,

    'caV50': -8 *mV,
    'caSc': 20 * mV,

    'kV50': -15 *mV,
    'kSc': 7 * mV,




    'tauV50': 0 * mV,
    'tauSc':  3* mV,
    'htau':3 *ms,
    'ktau':5 *ms,
    'taunalt': 1*ms,

    'phi': 0.02 * ms**(-1),

    'mBase': 50* ms,
    'mAmp': 92 *ms,
    'msigma': 18 *mV,
    'mMean': 10 *mV


}


eqs = """
    dv/dt =(Igapdown + Igapparallel + Isyn - ILeak - ICa - IK +Igap + stim_condition * stimulus(t))/C  :volt
    dm/dt = (mss- m)/(mBase + mAmp * exp(( -(v - mMean)**2)/(msigma ** 2))) :1
    dh/dt = (synin -h)/htau : siemens
    mss = 1/2 * (1 + tanh((v - kV50)/kSc)) :1
    taun = 1/ (phi * cosh((v - tauV50)/tauSc))  : second
    stim_condition : 1 #which neurons get current injection
    ICa = gCa * 1/2 * (1+tanh((v-caV50)/caSc)) *    (v -ECa) :amp
    IK = gK * m * (v - EK) :amp
    ILeak = gL * (v - EL) :amp
    Isyn = h * (Esyn - v) : amp
    Igap : amp # gap junction current
    Igapparallel : amp
    Igapdown : amp
    synin :siemens
    x : meter
    y : meter
    z : meter
    ind: 1
    is_spiking : boolean
"""

def Equilibrium(n, voltage):
    n.v = voltage
    n.m = 1/2 * (1 + np.tanh((voltage - n.namespace['kV50'])/ n.namespace['kSc']))
    n.h = 0
