#!/usr/bin/env python
# coding: utf-8

# In[1]:



from brian2 import *
import SMNNeuronModel as SMN
import MorrisLecarMuscleSquare as Mscle
from matplotlib import animation, rc


# In[2]:
#Defining size

jellyradius = 1.5 * cm
mouthradius = 0.35 * cm
crossectionrad = 1.5 * cm


# In[3]:
#Calculating number of neurons

connectionpercent = 0.3


numberofrings = 12
neuronsinring = int( (np.round(jellyradius*2 * np.pi /SMN.Parameters['segment'] /4 )*4)*numberofrings)
neuronsradial =  int(np.round((jellyradius-mouthradius)/SMN.Parameters['segment']) * numberofrings)
neuronsinsmallring = int( (np.round(mouthradius*2 * np.pi /SMN.Parameters['segment']/4)*4 )*numberofrings)
totalneurons = int(np.round(neuronsinring+ 4* neuronsradial + neuronsinsmallring))
maxdist =np.repeat(175, totalneurons)


sequence = np.array([8 * um])
areaofstim = len(sequence)
radii =np.tile(sequence, int(np.ceil(totalneurons/len(sequence))) ) * meter

# In[4]:
#Defining neurons

SMN.Parameters['stimulus'] = TimedArray([0*nA, 0*nA],dt=15*ms) #possibility for direct stim


neurons = NeuronGroup(totalneurons, SMN.eqs, namespace = SMN.Parameters, method= 'exponential_euler',
                      threshold='v > 0*mV', refractory = 'v>-20*mV',
                      events={'beginspike': 'v > -20*mV and not is_spiking',
                              'endspike': 'v < -40*mV and is_spiking'},
                      dt = 0.005 * ms)
neurons.stim_condition ='int(i < areaofstim) '


neurons.radius = radii[0: len(neurons)]

neurons.run_on_event( 'beginspike', 'is_spiking = True')
neurons.run_on_event( 'endspike', 'is_spiking = False')

SMN.Equilibrium(neurons, -60 * mV)
neurons.is_spiking = False


defaultclock.dt = 0.001 *ms
statemon = StateMonitor(neurons, 'v', dt = 0.1 *ms, record=np.arange(0,totalneurons))

spikes = SpikeMonitor(neurons)


start_mon = EventMonitor(neurons, 'beginspike')
stop_mon = EventMonitor(neurons, 'endspike')

neurons.y = 'sin(i/(neuronsinring) * 2 * pi)* jellyradius * int(i < neuronsinring)+        int(i>= neuronsinring+ neuronsradial*3 and i < neuronsinring+ neuronsradial*4) * (- jellyradius + (i - neuronsinring-neuronsradial*3) *  (jellyradius-mouthradius) / neuronsradial) +        int(i>= neuronsinring+neuronsradial*1 and i < neuronsinring+ neuronsradial*2) * ( jellyradius - (i - neuronsinring-neuronsradial*1) * (jellyradius-mouthradius) / neuronsradial)  + sin((i- neuronsinring- 4* neuronsradial)/(neuronsinsmallring) * 2 * pi)* mouthradius * int(i >= neuronsinring+ 4* neuronsradial) '


neurons.x = 'cos(i/(neuronsinring) * 2 * pi) * jellyradius* int(i < neuronsinring)+ int(i>= neuronsinring and i < neuronsinring+ neuronsradial) * (jellyradius - (i - neuronsinring) *  (jellyradius-mouthradius) / neuronsradial)                                           +        int(i>= neuronsinring+neuronsradial*2 and i < neuronsinring+ neuronsradial*3) * (- jellyradius + (i - neuronsinring-neuronsradial*2) * (jellyradius-mouthradius) / neuronsradial)      + cos((i- neuronsinring- 4* neuronsradial)/(neuronsinsmallring) * 2 * pi)* mouthradius * int(i >= neuronsinring+ 4* neuronsradial)'


# In[5]:
# Definition of the connection matrix

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









# In[5]:
#Connecting neurons


sources, targets = synapsematrix.nonzero()
S = Synapses(neurons, neurons, '''w : siemens # gap junction conductance
                                  Igap_pre = w*(v_post-v_pre): amp (summed)''')

S.connect(i = sources, j = targets)
resistance = SMN.Parameters['resistance']
neuronlen = SMN.Parameters['segment']

S.w = ' (0.2+ 0.8* int(i%numberofrings == j%numberofrings) )*1.2*radius_pre*radius_post**2 /(resistance* neuronlen**2  *(radius_pre**2+ radius_post**2)) *area_pre ' #N_outgoing N_incoming

# In[6]:
#placing muscles

xmuscle = []
ymuscle = []
musclelength = Mscle.lengths
k = jellyradius
firstmusclegroup = 0
tested = False
alternatesafe = 0
oldmusclelength = 0
plottingindex = []
while k >= mouthradius- 0.5* musclelength:
    nom = int(np.round(2 * np.pi * k/ musclelength/4)*4)
    for i in range(0,nom):
        xmuscle.append(cos(i/nom *2 * np.pi) * k)
        ymuscle.append(sin(i/nom *2 * np.pi) * k)
    if not tested:
        firstmusclegroup = nom
        tested = True
    k = k - musclelength
    print(len(xmuscle))
    if alternatesafe == 0:
        plottingindex = np.append(plottingindex, arange(oldmusclelength, len(xmuscle),2))

    alternatesafe = (alternatesafe+1)%3
    oldmusclelength = len(xmuscle)


# In[7]:
#Defining muscle model

Mscle.Parameters['stimulus'] = TimedArray([0*nA, 0*nA],dt=10*ms)
muscles = NeuronGroup(len(xmuscle), Mscle.eqs, namespace = Mscle.Parameters, method= 'exponential_euler',
        threshold='v > -20*mV', refractory = 'v>-60*mV', dt= 0.01 *ms)
musclespikes = SpikeMonitor(muscles)
muscles.stim_condition = 'int(i==0)'
muscles.x = xmuscle
muscles.y = ymuscle
Mscle.Equilibrium(muscles, -80 * mV)
musclestate = StateMonitor(muscles, 'v', dt = 0.1 * ms, record=np.arange(0,int(np.round(muscles.N))))



# In[8]:
#Muscle connectivity

V = Synapses(muscles, muscles, '''w : siemens # gap junction conductance
                                  Igap_post = w*(v_pre-v_post): amp (summed)''')

V.connect(condition='sqrt( (x_pre - x_post)**2 + (y_pre - y_post) **2) <= musclelength * 1.2') # 1.5
msclstr =Mscle.Parameters['ggap']
V.w= '0.5*msclstr- 0.4*msclstr * int(( sqrt((sqrt(x_pre**2+y_pre**2)-sqrt(x_post**2+y_post**2))**2) >100 * um ))'


U = Synapses(neurons, muscles, on_pre = 'h_post += 520/ numberofrings * nS', delay = 5 * ms) #500

U.connect(condition = '((sqrt( (x_pre - x_post)**2 + (y_pre - y_post) **2) <= musclelength * 1.) or (sqrt( (x_pre - x_post)**2 + (y_pre - y_post) **2) <= musclelength * 1.2 and i < neuronsinring+neuronsradial*4 and i > neuronsinring) )and ((rand()< 0.25))') # 12



# In[9]:
#Define stimulus

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


# In[10]:
#Run Simulation


activatedneurons= 5
inputsynapses(0.3,activatedneurons)
run(70*ms)
if len(spikes.values("t")[0])  == 0:
    inputsynapses(0.3,activatedneurons)
    run(60*ms)
    if len(spikes.values("t")[0])  == 0:
        inputsynapses(0.4,activatedneurons)
        run(20*ms)
        if len(spikes.values("t")[0])  == 0:
            inputsynapses(0.45,activatedneurons)
            run(15*ms)
            if len(spikes.values("t")[0])  == 0:
                inputsynapses(1.9,activatedneurons)
                run(10*ms)
                if len(spikes.values("t")[0])  == 0:
                    inputsynapses(2.6,activatedneurons)
                    run(30*ms)
run(400*ms)

# In[11]:
#plot results




iterator = 2
for i in range(0,int(neuronsinring/(iterator*2)/numberofrings)):
    plot(statemon.t/ms, statemon.v[i*iterator*numberofrings]/mV,  str('C'+str(i%8)) , linewidth=3.0)

#xlim(200,400)
xlabel('Time (ms)')
ylabel('mV');
plt.show()
# In[12]:
plottingindex=plottingindex.astype(int)
nills = np.full((len(plottingindex)), 0.0)
ones = np.full((len(plottingindex)), 1.0)
neuronwhite = np.repeat(1, len(neurons))
musclewhite = np.repeat(1, len(plottingindex))
sigmoiddiv = 30
msigmoiddiv = 50
npxmuscle = np.array(xmuscle)
npymuscle = np.array(ymuscle)
def PlotFrame(time, plotter = plt):
    plotter.xlim(( -(jellyradius/meter + 0.001), jellyradius/meter + 0.001))
    plotter.ylim((-(jellyradius/meter + 0.001), jellyradius/meter + 0.001))
    #timepoint = time * 10
    timepoint = (time+ int(np.round(spikes.values("t")[0][0]/ms))) * 10 -10
    musclecol = np.clip(np.maximum(np.minimum((-50*mV-musclestate.v[plottingindex,timepoint])/mV/40, ones ), nills),0, 1 )
    #for i in range(len(plottingindex)):
    #    musclecol.append(np.min([-np.min([musclestate.v[plottingindex[i],timepoint]/mV/80,0 ]), 1]))
    musclecol = [musclecol]
    #print("iamhere")
    #musclecol = [1 - 1/(1+ np.exp(-((musclestate.v[:,timepoint]/mV + 30) /msigmoiddiv ) ))]
    col = np.array([ musclecol[0] ,musclecol[0],musclewhite, musclewhite])
    plotter.scatter(npxmuscle[plottingindex], npymuscle[plottingindex], s = 150, color = col.transpose())

    neuroncol = [1 - 1/(1+ np.exp(-(statemon.v[:,timepoint]/mV /sigmoiddiv) ))]
    col = np.array([neuronwhite, neuroncol[0] ,neuroncol[0], neuronwhite])
    plotter.scatter(neurons.x[::numberofrings], neurons.y[::numberofrings], color = col.transpose()[::numberofrings])
    plotter.axis('off')
    plotter.text(-0.015, -0.015, str(time)+ " ms",  fontsize = "x-large")
    if time < 10:
        name = "000"+str(time)
    elif time < 100:
        name = "00"+str(time)
    elif time < 1000:
        name = "0" + str(time)
    else:
        name = str(time)
    plotter.savefig("./Animation/Frame"+ name+".png", bbox_inches= 'tight')
    return plotter


# In[12]:
rcParams['figure.figsize'] = [8, 8]

PlotFrame(20)


# In[13]:


def DoAnimation(Startpoint = 0, Endpoint = 180):
    for i in range(Startpoint, Endpoint):
        PlotFrame(i)
        plt.show()
