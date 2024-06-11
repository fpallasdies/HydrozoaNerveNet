    #!/usr/bin/env python
# coding: utf-8




from brian2 import *

import SMNNeuronModel as SMN

import MorrisLecarMuscleSquare as Mscle

from matplotlib import animation, rc
from mpl_toolkits.mplot3d import Axes3D

from IPython.display import HTML

belllength= 2.2*cm
mouthradius = 0.25 * cm
winkel = np.pi/2*1.1
shapefactor = 2/2.2
multiplier = 2
numberofrings = 12
connectionpercent = 0.22
PSPstrength = 490
steps = 0
xcords = []
zcords = []
def makesteps():
	global steps, xcords, zcords
	steps = int(np.round(belllength /SMN.Parameters['segment']/4)*4) * numberofrings
	xcords.append(0*cm)
	zcords.append(0*cm)
	for i in range(1,steps):
		newangle =winkel * np.min([np.power((abs(i+0.5)/(steps*shapefactor)),multiplier),1])
		#newangle =winkel * np.min([np.power((abs(i+0.5)/(steps)),1),1])
		xcords.append(xcords[-1]+ np.cos(newangle)* belllength/steps)
		zcords.append(zcords[-1]- np.sin(newangle)*belllength/steps )
	#print(newangle)
	zcords -= zcords[-1]




#True length after running with 2cm can be calculated by 1.5 / (xcords[-1]/cm) *2cm

# In[6]:

# 0.3
def RunSim():
    mouthradius = xcords[int(np.round(steps/4))]
    radialheight = int(np.round(3*steps/4))
    jellyradius = xcords[-1]

    neuronsinring = int( (np.round(jellyradius*2 * np.pi /SMN.Parameters['segment'] /4 )*4)*numberofrings)
    neuronsradial =  radialheight
    neuronsinsmallring = int( (np.round(mouthradius*2 * np.pi /SMN.Parameters['segment'] /4 )*4)*numberofrings)
    totalneurons = int(np.round(neuronsinring+ 4* neuronsradial + neuronsinsmallring))
    maxdist =np.repeat(175, totalneurons)


    sequence = np.array([8 * um])
    areaofstim = len(sequence)
    radii =np.tile(sequence, int(np.ceil(totalneurons/len(sequence))) ) * meter
    areaofstim





    SMN.Parameters['stimulus'] = TimedArray([0*nA, 0*nA],dt=15*ms) #30


    neurons = NeuronGroup(totalneurons, SMN.eqs, namespace = SMN.Parameters, method= 'exponential_euler',
                          threshold='v > 0*mV', refractory = 'v>-20*mV',
                          events={'beginspike': 'v > -20*mV and not is_spiking',
                                  'endspike': 'v < -40*mV and is_spiking'},
                          dt = 0.005 * ms)
    neurons.stim_condition ='int(i < areaofstim) '#or (i < neuronsinring and i > neuronsinring-1*areaofstim))'# or (i < neuronsinring and i >= neuronsinring-10))


    neurons.radius = radii[0: len(neurons)]

    neurons.run_on_event( 'beginspike', 'is_spiking = True')
    neurons.run_on_event( 'endspike', 'is_spiking = False')

    SMN.Equilibrium(neurons, -60 * mV)
    neurons.is_spiking = False


    defaultclock.dt = 0.001 *ms #0.001 *ms
    statemon = StateMonitor(neurons, 'v', dt = 0.1 *ms, record=np.arange(0,totalneurons))

    spikes = SpikeMonitor(neurons)


    start_mon = EventMonitor(neurons, 'beginspike')
    stop_mon = EventMonitor(neurons, 'endspike')


    neurons.z = 0*meter
    for i in range(0, neuronsinring+neuronsinsmallring+neuronsradial*4):
        if i < neuronsinring:
            neurons.y[i] = sin(i/(neuronsinring) * 2 * pi)* jellyradius
            neurons.x[i] = cos(i/(neuronsinring) * 2 * pi) *jellyradius
        elif i>= neuronsinring and i < neuronsinring+neuronsradial:
            neurons.x[i]= xcords[-i-1+neuronsinring]
            neurons.z[i] = zcords[-i-1+neuronsinring]
        elif i>= neuronsinring+neuronsradial and i < neuronsinring+neuronsradial*2:
            neurons.y[i]= xcords[-i-1+neuronsinring+neuronsradial]
            neurons.z[i] = zcords[-i-1+neuronsinring+neuronsradial]
        elif i>= neuronsinring+neuronsradial*2 and i < neuronsinring+neuronsradial*3:
            neurons.x[i]= -xcords[-i-1+neuronsinring+neuronsradial*2]
            neurons.z[i] = zcords[-i-1+neuronsinring+neuronsradial*2]
        elif i>= neuronsinring+neuronsradial*3 and i < neuronsinring+neuronsradial*4:
            neurons.y[i]= - xcords[-i-1+neuronsinring+neuronsradial*3]
            neurons.z[i] = zcords[-i-1+neuronsinring+neuronsradial*3]

        elif i>= neuronsinring+neuronsradial*4 :
            neurons.y[i] = sin((i- neuronsinring- 4* neuronsradial)/(neuronsinsmallring) * 2 * pi)* mouthradius
            neurons.x[i] = cos((i- neuronsinring- 4* neuronsradial)/(neuronsinsmallring) * 2 * pi) *mouthradius
            neurons.z[i] = zcords[int(steps/4)]



    synapsematrix = np.full((len(neurons), len(neurons)), 0.0)
    for i in range(len(neurons)):
        if i%numberofrings < numberofrings -1 and np.random.rand() < connectionpercent * (1**(i%numberofrings)):
            synapsematrix[i][i+1] = 1
            synapsematrix[i+1][i] = 1
        if i < neuronsinring:
            synapsematrix[i][(i+numberofrings)%neuronsinring] = 1
            synapsematrix[i][(i-numberofrings)%neuronsinring] = 1
        elif i >= neuronsinring and i < neuronsinring+ neuronsradial:
            if i+ numberofrings < neuronsinring + neuronsradial:
                synapsematrix[i][i+numberofrings] = 1
            else:
                synapsematrix[i][i%numberofrings+ neuronsinring+ neuronsradial*4]=1

            if i- numberofrings >= neuronsinring:
                synapsematrix[i][i-numberofrings] = 1
            else:
                synapsematrix[i][i%numberofrings] = 1

        elif i >= neuronsinring+ neuronsradial*1 and i < neuronsinring+ neuronsradial*2:
            if i+ numberofrings < neuronsinring + neuronsradial*2:
                synapsematrix[i][i+numberofrings] = 1
            else:
                synapsematrix[i][i%numberofrings+ neuronsinring+ neuronsradial*4 + int(1/4*neuronsinsmallring)] =1
            if i- numberofrings >= neuronsinring + neuronsradial * 1:
                synapsematrix[i][i-numberofrings] = 1
            else:
                synapsematrix[i][i%numberofrings+  int(1/4*neuronsinring)] = 1

        elif i >= neuronsinring+ neuronsradial*2 and i < neuronsinring+ neuronsradial*3:
            if i+ numberofrings < neuronsinring + neuronsradial*3:
                synapsematrix[i][i+numberofrings] = 1
            else:
                synapsematrix[i][i%numberofrings+ neuronsinring+ neuronsradial*4 + int(2/4*neuronsinsmallring)] = 1
            if i- numberofrings >= neuronsinring + neuronsradial * 2:
                synapsematrix[i][i-numberofrings] = 1
            else:
                synapsematrix[i][i%numberofrings+  int(2/4*neuronsinring)] = 1

        elif i >= neuronsinring+ neuronsradial*3 and i < neuronsinring+ neuronsradial*4:
            if i+ numberofrings < neuronsinring + neuronsradial*4:
                synapsematrix[i][i+numberofrings] = 1
            else:
                synapsematrix[i][i%numberofrings+ neuronsinring+ neuronsradial*4 + int(3/4*neuronsinsmallring)] = 1
            if i- numberofrings >= neuronsinring + neuronsradial * 3:
                synapsematrix[i][i-numberofrings] = 1
            else:
                synapsematrix[i][i%numberofrings+  int(3/4*neuronsinring)] = 1


        elif i  >= neuronsinring+ neuronsradial*4:
            synapsematrix[i][(i-(neuronsinring+ neuronsradial*4)+ numberofrings)%neuronsinsmallring + (neuronsinring+ neuronsradial*4)] = 1
            synapsematrix[i][(i- (neuronsinring+ neuronsradial*4)-numberofrings)%neuronsinsmallring + (neuronsinring+ neuronsradial*4)] = 1
        temp = np.maximum(synapsematrix[i], synapsematrix[:][i])
        synapsematrix.transpose()[i] = temp




    #synapsematrix = np.tril(synapsematrix)
    sources, targets = synapsematrix.nonzero()
    S = Synapses(neurons, neurons, '''w : siemens # gap junction conductance
                                      Igap_pre = w*(v_post-v_pre): amp (summed)''')

    S.connect(i = sources, j = targets)
    resistance = SMN.Parameters['resistance']
    neuronlen = SMN.Parameters['segment']

    S.w = ' (0.2+ 0.8* int(i%numberofrings == j%numberofrings) )*1.2*radius_pre*radius_post**2 /(resistance* neuronlen**2  *(radius_pre**2+ radius_post**2)) *area_pre ' #N_outgoing N_incoming

    # # Muscles



    steps

    xmuscle = []
    ymuscle = []
    zmuscle = []
    indmuscle = []
    musclelength = Mscle.lengths
    k = 0
    firstmusclegroup = 0
    tested = False
    alternatesafe = 0
    oldmusclelength = 0
    plottingindex = []

    while k <= int(steps/4*3):
        nom = int(np.round(2 * np.pi * xcords[-1-k]/ musclelength/4)*4)
        for i in range(0,nom):
            xmuscle.append(cos(i/nom *2 * np.pi) * xcords[-1-k])
            ymuscle.append(sin(i/nom *2 * np.pi) * xcords[-1-k])
            zmuscle.append(zcords[-1-k])
            indmuscle.append(len(zcords)-1-k)
        if not tested:
            firstmusclegroup = nom
            tested = True
        k = k + numberofrings
        print(len(xmuscle))

        if alternatesafe == 0:
            plottingindex = np.append(plottingindex, arange(oldmusclelength, len(xmuscle),2))

        alternatesafe = (alternatesafe+1)%3

        #plottingindex = np.append(plottingindex, arange(oldmusclelength, len(xmuscle),1))
        oldmusclelength = len(xmuscle)

    plottingindex=plottingindex.astype(int)
    npxm = np.array(xmuscle)
    npym = np.array(ymuscle)
    npzm = np.array(zmuscle)




    #plottingindex[359]



    Mscle.Parameters['stimulus'] = TimedArray([0*nA, 0*nA],dt=10*ms)
    muscles = NeuronGroup(len(xmuscle), Mscle.eqs, namespace = Mscle.Parameters, method= 'exponential_euler', threshold='v > -20*mV', refractory = 'v>-60*mV', dt= 0.01 *ms)
    muscles.stim_condition = 'int(i==0)'
    musclespikes = SpikeMonitor(muscles)
    muscles.x = xmuscle
    muscles.y = ymuscle
    muscles.z = zmuscle
    muscles.ind = indmuscle
    Mscle.Equilibrium(muscles, -80 * mV)
    musclestate = StateMonitor(muscles, 'v', dt = 0.1 * ms, record=np.arange(0,int(np.round(muscles.N))))





    V = Synapses(muscles, muscles, '''w : siemens # gap junction conductance
                                      Igap_post = w*(v_pre-v_post): amp (summed)''')

    V.connect(condition='sqrt( (x_pre - x_post)**2 + (y_pre - y_post) **2 + (z_pre - z_post) **2) <= musclelength * 1.2') # 1.5
    msclstr =Mscle.Parameters['ggap']
    V.w= '0.5*msclstr- 0.4*msclstr * int( abs(z_pre - z_post)>100 * um )'
    # 0.5


    U = Synapses(neurons, muscles, on_pre = 'h_post += PSPstrength/ numberofrings * nS', delay = 5 * ms) #390

    U.connect(condition = '((sqrt( (x_pre - x_post)**2 + (y_pre - y_post) **2+ (z_pre - z_post) **2) <= musclelength * 1.) or (sqrt( (x_pre - x_post)**2 + (y_pre - y_post) **2+ (z_pre - z_post) **2) <= musclelength * 1.2 and i < neuronsinring+neuronsradial*4 and i > neuronsinring) )and ((rand()< 0.25))') # 12

    # # Run and Plot




    store()




    startingneuron = 0
    def inputsynapses(strength, non = 1):
        non-= 1
        start = startingneuron*numberofrings
        neurons[start].syn += strength* mS/cm**2
        index = 1
        while non > 0 :
            neurons[start+index*numberofrings].syn += strength* mS/cm**2
            non -=1
            if non >0:
                neurons[(start-index*numberofrings)%neuronsinring].syn += strength* mS/cm**2
            index += 1




    activatedneurons= 5
    inputsynapses(0.35,activatedneurons)
    run(60*ms)
    if len(spikes.values("t")[0])  == 0:
        inputsynapses(0.4,activatedneurons)
        run(60*ms)
        if len(spikes.values("t")[0])  == 0:
            inputsynapses(0.45,activatedneurons)
            run(20*ms)
            if len(spikes.values("t")[0])  == 0:
                inputsynapses(0.9,activatedneurons)
                run(15*ms)
                if len(spikes.values("t")[0])  == 0:
                    inputsynapses(1,activatedneurons)
                    run(10*ms)
                    if len(spikes.values("t")[0])  == 0:
                        inputsynapses(2.8,activatedneurons)
                        run(30*ms)
    run(300*ms)
    run(350*ms)




    startingtime = int(np.floor(np.min(musclespikes.t/ms)))
    x =savemuscleactivity(muscles, musclestate, startingtime,startingtime+475,0.1)


# In[ ]:

import csv
nogroups = 64
def savemuscleactivity(muscles, musclemonitor, start, stop, dt, filename = "muscles.csv"):
    ####Grenzen aufgeweicht
    sectionlength= (steps*0.75)*2/16
    muscleactivity = []
    musclegroup = []
    for i in range(0, nogroups):
        musclegroup.append([0])


    for i in range(0, len(muscles)):
        angle = np.arctan2(muscles[i].y/cm, muscles[i].x/cm)[0]+1/8*np.pi
        if angle < 0:
            angle = np.pi* 2 + angle
        distance = muscles[i].ind
        #print(distance)
        index = int(np.floor( 8*angle /(np.pi *2)))
        index += int(8* (8-1-np.min([np.round((distance-steps* 0.25)/(sectionlength)),7])  ))
        #print(index)
        musclegroup[index].append(i)
    for i in range(0, len(musclegroup)):
        musclegroup[i].pop(0)

    print("Phase2")
    for i in range(start, stop):
        print(i)
        timestep = int(i/dt)
        row = []
        currentstate = musclemonitor.v[:,timestep]
        for j in range(0, len(musclegroup)):
            sum = 0
            subset = currentstate[musclegroup[j]]
            #sum = np.sum(subset/mV)
            sum = np.sum(np.maximum(subset/mV+10,0))
            #sum = sum/ len(musclegroup[j])
            #sum = (sum+ 80*len(musclegroup[j]))
            #sum = sum/len(musclegroup[j])
            if sum <= 0:
                row.append(0)
            else:
                row.append(sum)
            #print(row)
        muscleactivity.append(row)

    muscleactivity = muscleactivity/np.max(muscleactivity)
    while (0< 1):
        if (np.max(muscleactivity[0]) > 0):
            break
        else:
            muscleactivity = muscleactivity[1:]
    with open(filename,'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter = ' ', quotechar = '|', quoting = csv.QUOTE_MINIMAL)
        for i in range(0, len(muscleactivity)):
            writer.writerow(np.concatenate((muscleactivity[i][0:64:8],muscleactivity[i][4:64:8]  )))
            #writer.writerow(muscleactivity[i])
    return muscleactivity
