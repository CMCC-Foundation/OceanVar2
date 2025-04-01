#!/usr/bin/env python
# coding: utf-8

import itertools
import matplotlib as mpl
import matplotlib.pyplot as plt
#import graphviz
import pandas as pd
import numpy as np
pd.set_option("display.precision", 12)
cmapp = plt.get_cmap('Reds')

exe = ["exe_" + i for i in [@EXE@]]
inp = ["inp_" + i for i in [@INP@]]
nmlid = [x[0] for x in [@NMLID@]]
nmlnm = [x.split(' ') for x in [@NMLNM@]]
nml = [ '_'.join([ nmlid[j] + x[j]  for j in range(len(x)) ]) for x in list(itertools.product(*nmlnm)) ]
proc= [@PROC@]
rerun = @RERUN@

#exe = ['']
#inp = ['']
#nml = ['']
#proc = ['1x1', '2x2']

# Parse the processes, sort them according to full amount of processes
nproc = np.array([ np.prod(np.array(i.split('x'),dtype=np.int32)) for i in proc])
proc =[proc[i] for i in sorted(range(len(nproc)), key=lambda k: nproc[k])]
nproc.sort()
nproc_u, idcs = np.unique(nproc, return_counts=True)

#Make graphviz graph
def makegraph(df):
    #Set lists according to maximum level of nested timer
    maxlvl  = df["TIMER_LVL"].max()
    maxtime = df["MEASURED_TIME"].max()
    prvlvl  = ['']*(maxlvl+1)
    prvtim  = [0]*(maxlvl+1)

    dot = graphviz.Digraph()
    nds=[]
    for _, row in df.iterrows():
        lvl         = row['TIMER_LVL']
        time        = row['MEASURED_TIME']
        timeperloop = row['MEASURED_TIMEperLOOP']
        cname       = row['TIMER_NAME (looped)'].replace('\t','').replace(' ','')+f"\n" + "tot " +str(time) +f"\n" + "scaled " + str(timeperloop)

        prvtim[lvl] = time
        prvlvl[lvl] = cname

        if(lvl > 0):
            reltim = time/prvtim[lvl-1]
        else:
            reltim = time/maxtime

        if(reltim> float(.2)): #Only put the most heavy timers in the

            col = mpl.colors.rgb2hex(cmapp(reltim)) #set color according to relative time
        else:
            col = "green"

        #Add node to graph and edge to list of edges
        dot.node(cname, style = 'filled', fillcolor = col, fontcolor = 'yellow') #color = "#%2"+ str(r) + "%2" + str(g) + "%2" + str(b) + "%2" + str(alpha))
        if(lvl>0):
            nds.append([prvlvl[lvl-1],cname])

    dot.edges(nds)
    return dot



#Function to get standard deviation of union of two sets, given std deviation and average of both sets
def combine_std(avg, std, newavg, newstd, jj):
    if jj == 0:
        return newstd
    n = jj*rerun
    m = rerun
    return (( (n-1)*np.square(std) + (m-1)*np.square(newstd)  + (n*m*np.square(newavg - avg))/(n+m))/(n+m-1))




plotvals = np.zeros((len(idcs), 2))
for cpt in itertools.product(exe,inp,nml):
    #Set path
    cpath= '/'.join(cpt)

    #Set average and standard deviation arrays for profiling
    df = pd.read_csv(cpath+f"/profileavg{proc[0]}.txt",sep=':')
    vals_avg = np.zeros((df.shape[0],len(idcs), 2))
    vals_std = np.zeros((df.shape[0],len(idcs), 2))

    idx = 0
    for i in range(len(idcs)):
        for j in range(idcs[i]):
            cpathh = f"{cpath}/profileavg{proc[idx]}.txt"
            df= pd.read_csv(cpathh,sep=':')
            #df= pd.read_csv(f"profileavg{proc[idx]}.txt",sep=':')

            #Compute new standard deviation
            vals_std[:,i,0] = combine_std(vals_avg[:,i,0],vals_std[:,i,0],df['MEASURED_TIME'].to_numpy(),df['MEASURED_TIMEstd'].to_numpy(),j)
            vals_std[:,i,1] = combine_std(vals_avg[:,i,1],vals_std[:,i,1],df['MEASURED_TIMEperLOOP'].to_numpy(),df['MEASURED_TIMEperLOOPstd'].to_numpy(),j)
            #Add average
            vals_avg[:,i,0] += df['MEASURED_TIME']
            vals_avg[:,i,1] += df['MEASURED_TIMEperLOOP']

            #Make graphviz graph
            #curgraph = makegraph(df)
            #curgraph.render(directory=cpath)

            idx = idx + 1
        vals_avg[:,i,0] = vals_avg[:,i,0]/idcs[i]
        vals_avg[:,i,1] = vals_avg[:,i,1]/idcs[i]
        plotvals[i,:] = [vals_avg[0,i,0], vals_std[0,i,0]]
    plt.errorbar(nproc_u, plotvals[:,0], yerr=plotvals[:,1])
    plt.ylabel("time (s)")
    plt.xlabel("nr processors")
    plt.title(df['TIMER_NAME (looped)'][0].replace('\t','').replace(' ','') + " time " + '_'.join(cpt))
    plt.xticks(nproc_u)
    filename = '_'.join(cpt)
    plt.savefig("globaltime/all.pdf", bbox_inches='tight',format='pdf')
    #Plot all



