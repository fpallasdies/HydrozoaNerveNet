#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 16:30:18 2019

@author: pallasdies
"""

from brian2 import *
import numpy as np
shutdown= 1

urest = -70*mV
uthersh = 20*mV
refract = 5*ms

allmult = 1
Namult = 2#
Kmult =1.5
Leak = 0.05
gapmult = 1
#
radius = 8* um # 8 #22

crosssection = np.pi * radius**2

segment =  175 * um,
Parameters={
    'segment': 175 * um,
    'resistance': 56 * ohm * cm,
    'C' : 1.52 * uF/cm**2,



    'EK':-60.752*mV,
    'ENa' :54.145429285085754*mV,
    'ECa' :90.835*mV, #50

    'EL' :-60*mV,
    'Esyn' : 54*mV,
    #
    'gfast' :4.*1.5690576705781536* mS/cm**2 ,
    'gslow' :4*0.28730864456621796 * mS/cm**2 ,
    'gNa' :  2.4*5.749005743113725* mS/cm**2 ,

    'gL' : 0.038 *mS/cm**2,
    'gNasteady' : 2.4* 1.4584215772255933*mS/cm**2 ,
    'gCa': 3.2*0.02335017281 *mS/cm**2,
    'gCaslow': 3.2*0.01189077076  *mS/cm**2,



    'mV50' : -16.490523976874243*mV,
    'mSc' :4.978256225607337*mV,
    'mBase' : 0.13913354538886055*ms,
    'mAmp' : 0.32986446303792966*ms,
    'mMean' : -33.233693183941185*mV,
    'msigma' : 36.71226668720717*mV ,
    'mpower' : 7.573135290908616,

    'hV50' : -18.579701079553796*mV,
    'hSc' : -9.062744701794449*mV,
    'hBase' : 0.5964265854541019*ms,
    'hAmp' : 1.3501768052366485*ms,
    'hMean' : -12.88067749433803*mV,
    'hsigma' : 55.995330796565945*mV,
    'hpower' : 1.6463313527170331,


    'kV50' : -21.501644031414425*mV,
    'kSc' : 2.762622363669918*mV,
    'kBase' : 0.7170227926056161*ms,
    'kAmp' : 1.0444212864905686*ms ,
    'kMean' : -63.75321265302916*mV,
    'ksigma' : 14.686285264224349*mV,
    'kpower' : 3.9583080072724814,

    'lV50' : -30.65328693472969*mV,
    'lSc' : -6.237090854043124*mV,
    'lBase' : 6.615900358423904*ms,
    'lAmp' : 648.681035686458*ms,
    'lMean' : -31.615159468634353*mV,
    'lsigma' : 5.46427467830178*mV,
    'lpower' : 0.7757827015829779,

    'oV50' : -21.5395*mV,
    'oSc' : 13.22204576221433*mV,
    'oBase' : 0.8470694223564093*ms,
    'oAmp' : 1.93190286076001*ms,
    'oMean' : -17.10007953178716*mV,
    'osigma' : 41.575533022696035*mV,
    'opower' : 2.893706227260207,

    'pV50' :-43.33*mV,
    'pSc' : -9.466167783902282*mV,
    'pBase' : 31.343888669000783*ms,
    'pAmp' : 2102.3654833846886*ms  ,
    'pMean' : 32.56312263816121*mV  ,
    'psigma' : 29.98308676766491*mV,
    'ppower' : 1.9304387744636067,

    'nV50' : 34.543*mV,
    'nSc' :19.277529879847897*mV ,
    'nBase' : 121.72978833862781*ms,
    'nAmp' : 637.1920815835433*ms,
    'nMean' : -29.55819*mV,
    'nsigma' : 19.26587408665036*mV,
    'npower' : 1,


    'qV50' :-4.5*mV,
    'qSc' : 5.49610949666*mV,
    'qBase' : 0.33485032257*ms,
    'qAmp' : 46.8266876906*ms  ,
    'qMean' : 47.9*mV  ,
    'qsigma' : 15.5177035589*mV,
    'qpower' : 1,

    'rV50' : -17.828*mV,
    'rSc' : -1*mV ,
    'rBase' : 24.9446338414*ms,
    'rAmp' : 43.3484673197 *ms,
    'rMean' : 5.675*mV,
    'rsigma' : 5.26443425157 *mV,
    'rpower' : 1.13216278434,

    'sV50' :-20.392*mV,
    'sSc' : 5.14526812194*mV,
    'sBase' : 0.25710490182*ms,
    'sAmp' : 1.50644920548*ms  ,
    'sMean' : -54.505*mV  ,
    'ssigma' : 41.6519820594*mV,
    'spower' : 5.44433215606,


    'syntau' : 25 *ms,#25



}
#    dv/dt =( -INa - IK - ICa - ILeak +Igapparallel + IgapparallelB  +Igap +IgapB  + stim_condition * stimulus(t))/(C*area)  :volt


#(1+trigger)*
eqs = """

    dv/dt =( -INa - IK - ICa - ILeak + Isyn +Igapparallel + IgapparallelB  +Igap +IgapB  + stim_condition * stimulus(t))/(C*area)  :volt
    dm/dt = (1/(1+ exp(clip((mV50 - v)/mSc,-15,15))) - m )    / (mBase + mAmp * exp(( -(v - mMean)**2)/(msigma ** 2))) :1
    dh/dt = (1/(1+ exp(clip((hV50 - v)/hSc,-15,15))) - h )    / (hBase + hAmp * exp((- (v - hMean)**2)/(hsigma ** 2))) :1
    do/dt = (1/(1+ exp(clip((oV50 - v)/oSc,-15,15))) - o )    / (oBase + oAmp * exp((- (v - oMean)**2)/(osigma ** 2))) :1
    dp/dt = (1/(1+ exp(clip((pV50 - v)/pSc,-15,15))) - p )    / (pBase + pAmp * exp((- (v - pMean)**2)/(psigma ** 2))) :1
    dn/dt = (1/(1+ exp(clip((nV50 - v)/nSc,-15,15))) - n )    / (nBase + nAmp * exp((- (v - nMean)**2)/(nsigma ** 2))) :1
    dk/dt = (1/(1+ exp(clip((kV50 - v)/kSc,-15,15))) - k )    / (kBase + kAmp * exp((- (v - kMean)**2)/(ksigma ** 2))) :1
    dl/dt = (1/(1+ exp(clip((lV50 - v)/lSc,-15,15))) - l )    / (lBase + lAmp * exp((- (v - lMean)**2)/(lsigma ** 2))) :1
    dq/dt = (1/(1+ exp(clip((qV50 - v)/qSc,-8,8))) - q )    / (qBase + qAmp * exp((- (v - qMean)**2)/(qsigma ** 2))) :1
    dr/dt = (1/(1+ exp(clip((rV50 - v)/rSc,-8,8))) - r )    / (rBase + rAmp * exp((- (v - rMean)**2)/(rsigma ** 2))) :1
    ds/dt = (1/(1+ exp(clip((sV50 - v)/sSc,-8,8))) - s )    / (sBase + sAmp * exp((- (v - sMean)**2)/(ssigma ** 2))) :1
    dsyn/dt = -syn/syntau :siemens/meter**2
    stim_condition : 1 #which neurons get current injection
    INafast = gNa * (m**mpower) * (h**hpower)* (  v -ENa) *area : amp
    INaslow = gNasteady * (k**kpower) * (l**lpower) * (v- ENa) *area :amp
    NaConductance =area* (gNa * (m**mpower) * (h**hpower) + gNasteady * (k**kpower) * (l**lpower)):siemens
    KConductance =area* (gfast* (o**opower) * (p**ppower) + gslow * (n**npower)) :siemens
    CaConductance = area* ((gCa* (q**qpower) * (r**rpower))+ gCaslow* (s**spower) ) : siemens
    INa = INafast+ INaslow : amp
    IKfast = gfast* (o**opower) * (p**ppower) *(v -EK) *area :amp
    IKslow = gslow * (n**npower)* (v - EK) *area :amp
    IK = IKfast + IKslow :amp
    ICa = (gCa* (q**qpower) * (r**rpower) *(v -ECa)+gCaslow * (s**spower) * (v-ECa) ) *area :amp
    ILeak = gL * (v - EL) *area :amp
    Isyn = syn * clip(Esyn - v,0 *mV, 200 *mV) *area: amp
    trigger : 1
    Igap : amp # gap junction current
    IgapB : amp
    Igapparallel : amp
    IgapparallelB : amp
    is_spiking : boolean
    connectable: boolean
    x : meter
    y : meter
    z : meter
    area = 2 * pi  * 175 * um * radius :meter**2
    crosssection = pi * radius**2 : meter**2
    radius : meter
"""


def Equilibrium(n, voltage):
    n.v = voltage
    n.m = 1/(1+ np.exp((n.namespace['mV50'] - voltage )/n.namespace['mSc']))
    n.h = 1/(1+ np.exp((n.namespace['hV50'] - voltage )/n.namespace['hSc']))
    n.o = 1/(1+ np.exp((n.namespace['oV50'] - voltage )/n.namespace['oSc']))
    n.p = 1/(1+ np.exp((n.namespace['pV50'] - voltage )/n.namespace['pSc']))
    n.n = 1/(1+ np.exp((n.namespace['nV50'] - voltage )/n.namespace['nSc']))
    n.k = 1/(1+ np.exp((n.namespace['kV50'] - voltage )/n.namespace['kSc']))
    n.l = 1/(1+ np.exp((n.namespace['lV50'] - voltage )/n.namespace['lSc']))
    n.q = 1/(1+ np.exp((n.namespace['qV50'] - voltage )/n.namespace['qSc']))
    n.r = 1/(1+ np.exp((n.namespace['rV50'] - voltage )/n.namespace['rSc']))
    n.s = 1/(1+ np.exp((n.namespace['sV50'] - voltage )/n.namespace['sSc']))
