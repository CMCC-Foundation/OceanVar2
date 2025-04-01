#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 14:54:31 2024

@author: francescocarere
email : francesco.carere@cmcc.it

This file is part of the OceanVar testing suite.
It plots the output of the tests

"""

print("---------------------------")
print("Running Python script to make plots")
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import math
import itertools
import sys

assert sys.version_info[0] >= 3 and sys.version_info[1] >= 11,"Please use Python 3.11.5 or higher"

#To avoid range(len), shorten code
rl =lambda x:range(len(x))

idsys = 1
nrargs_ncsry=2
if len(sys.argv) != nrargs_ncsry+1:
    print (str(nrargs_ncsry), "arguments should be passed, you passed", str(len(sys.argv)))
    print ("----------------")
    print ("Arguments passed:")
    for i in rl(sys.argv) :
        print("arg ", i, ": ", sys.argv[i])
    print ("----------------")
    print ("Please pass the following arguments")
    print (str(idsys) + ") output path for diff plots (i.e. with \"eta\", \"sal\",...) on the axis");               idsys+= 1
    print (str(idsys) + ") output path for variable context plots  on the  x-axis");               idsys+= 1
    sys.exit(0)




path_plot_1  = sys.argv[idsys];     idsys+= 1
path_plot_2   = sys.argv[idsys];    idsys+= 1



files=[@file@]
err_type=@ERRTYPE@
#files=['eta','sal','tem','uvl','vvl']


inp_names=[@IDS@]
inp_types=[@TYPES@]
#inp_names=[["m1","m2"], ["n1","n2","n3"], ["R1","R2","R3","R4"], ["B1","B2","B3","B4","B5"], ["x1","x2"]]
#inp_types=["fixed", "fixed", "var", "var", "fixed"]

assert any((i == 'var' for i in inp_types)), "Error: at least one variable should be set to 'var'"


assert len(inp_names)==len(inp_types), "AssertionError: len(inp_names)==len(inp_types) Please specify if variables are fixed or variable"

inp_var=[inp_names[i] for i in rl(inp_names) if inp_types[i]=='var' ]
inp_fxd=[inp_names[i] for i in rl(inp_names) if inp_types[i]!='var' ]

n_fxd=math.prod(map(len,inp_fxd))
n_var=math.prod(map(len,inp_var))


var_list= np.array([
    @LIST_HERE@
    , dtype=np.float64)
#var_list=np.random.random(math.prod(map(len,inp_fxd)) * len(files) * math.prod(map(len,inp_var))**2)

#reshape into right size
reshape_var_lst=[]
for j in rl(inp_var):
    reshape_var_lst.append(inp_var[j])
    reshape_var_lst.append(inp_var[j])

print("--------------")
print("Reshaping array")
#var_arr2=np.reshape(var_list,list(map(len,inp_fxd)) + [np.prod(2*list(map(len,inp_var)))*len(files)])
var_arr2=np.reshape(var_list,list(map(len,inp_fxd+2*inp_var+[files])))
transp_vec= [ i for i in rl(inp_fxd)] + list(itertools.chain.from_iterable( (  (len(inp_fxd)+i, len(inp_fxd)+i+len(inp_var)) for i in rl(inp_var) ) )) + [len(inp_fxd)+2*len(inp_var)]
var_arr=np.transpose(var_arr2,transp_vec)
#var_arr=np.reshape(var_arr2,list(map(len,reshape_var_lst+[files])))
print("Finished reshaping array")


print("--------------")
print("Defining functions")



#%% Function definitions for making plots
def Make_plot_testing(plot_loop, plot_type,plot_path, plot_idx):
    #This routine makes the plots in a functional style.
    #Loops are over the plots that need to be made ("plot_loop")
    #plot_loop contains sets of indices which indicate which plot need to be made
    #In particular, plot_loop[0] are the indices for the first plot and:
    #   plot_loop[0][0] gives the indices for the titles of the plot
    #   plot_loop[0][1] gives the indices for the labels of the plot (the lines in the plot)
    #   plot_loop[0][2] gives the indices for the x_ticks of the plot (the x-dimension)
    #Together with plot_type, the indices can be mapped to indices of the input list "var_arr"

    #First some internal functions for helping with getting the indices and the names


    ################################################################
    #### START INTERNAL FUNCTION DEFINITIONS "Make_plot_testing" ###
    ################################################################

    def flatten(container):
        for i in container:
            if isinstance(i, (list,tuple)):
                for j in flatten(i):
                    yield j
            else:
                yield i

    def get_indx(pt_arr, pt_type,i_list):
        #Input is:
        #   pt_arr , set of "array of indices"
        #   pt_type, type of elements
        #   i_list , indexes for "array of indices"
        #Output is:
        #   reordered index suitable for "var_arr"

        var_order=["fxd","var","file"]

        idx=[]
        for i in rl(var_order):
            j=pt_type.index(var_order[i])
            idx.append(pt_arr[j][i_list[j]])

        idx[1] = [i[j] for i in idx[1] for j in rl(i)]
        return [*idx[0], *idx[1], idx[2]]

    def get_name(i_arr, i_type,sepr=" AND "):
        #Produce a name the name of the index i_arr with i_type

        #######################################################
        #### START INTERNAL FUNCTION DEFINITIONS "get_name" ###
        #######################################################
        def get_name_file(nm_arr,i_arri):
                return nm_arr[i_arri]

        def get_name_var(nm_arr, i_arri):
            var_name=""
            for j in rl(i_arri):
                if i_arri[j][0] != i_arri[j][1]:
                    var_name+=nm_arr[j][i_arri[j][0]]+"-"+nm_arr[j][i_arri[j][1]]+"_"
                else:
                    var_name+=nm_arr[j][i_arri[j][0]]+"_"
            if var_name=="":
                sys.exit("Error in getting names: the tuple" + i_arri + "should not be a diagonal index")
            return str(var_name[:-1])

        def get_name_fxd(nm_arr, i_arri):
            fxd_name=""
            for j in rl(i_arri):
                fxd_name+="_"+inp_fxd[j][i_arri[j]]
            return fxd_name[1:]

        #######################################################
        ##### END INTERNAL FUNCTION DEFINITIONS "get_name" ####
        #######################################################

        if   (i_type == 'file') :
            name_arr=files
            cfunc = lambda x,y: get_name_file(x,y)
        elif (i_type == "var"):
            name_arr=inp_var
            cfunc = lambda x,y: get_name_var(x,y)
        elif (i_type == "fxd"):
            name_arr=inp_fxd
            cfunc = lambda x,y: get_name_fxd(x,y)
        else:
            sys.exit("ERROR: please provide a correct type.")

        name=""
        try:
            len(i_arr)
        except:
            i_arr=[i_arr]

        namelst=[""]*len(i_arr)
        for i in rl(i_arr):
            new_name=cfunc(name_arr, i_arr[i])
            name+=sepr+new_name
            namelst[i]=new_name


        return name[len(sepr):], namelst




    #############################################################
    ### END INTERNAL FUNCTION DEFINITIONS "Make_plot_testing" ###
    #############################################################


    ###################
    ### START PLOTS ###
    ###################

    for plot_arr in plot_loop:

        ncols=int((len(plot_arr[1])-1)/14)+1
        scale_factor=len(plot_arr[2])*ncols*len(plot_arr[0])*2/21
        plt.figure(figsize=(max(10,scale_factor*10), 8))
        plt.ylabel(err_type + " L_1 differences")
        plt.yscale("log")

        lines_lst = ["-","--","-.",":"]
        curline_idx=-1

        title_name, pstfx_names = get_name(plot_arr[0],plot_type[0])
        temp, label_names = get_name(plot_arr[1], plot_type[1])
        temp, xtick_names = get_name(plot_arr[2], plot_type[2])

        if( len(plot_arr[0]) > 1 ):
            labels=[]
            handles=[]

        #Cycle over colors in the plots
        #color_iter = plt.cm.rainbow(np.linspace(0, 1, len(plot_arr[1])))
        color_iter = plt.cm.nipy_spectral(np.linspace(0,1,len(plot_arr[1])))

        for i_title in rl(plot_arr[0]):
            full_title=plot_type[0] + ": " + title_name
            plt.title(plot_type[0] + ": " + title_name, fontsize = 550/len(full_title))
            plotmarkers = list(matplotlib.markers.MarkerStyle.markers.keys())

            for i_label in rl(plot_arr[1]):

                plot_vls=[]
                for i_xticks in rl(plot_arr[2]):
                    #cur_idx= get_indx(plot_arr,plot_type,[i_title,i_label,i_xticks])
                    #plot_vls.append(var_arr[*cur_idx])
                    plot_vls.append(var_arr[*get_indx(plot_arr,plot_type,[i_title,i_label,i_xticks])])




                #label_pstfx=""
                #if len(plot_arr[0]) > 1:
                #    label_pstfx="__"
                #    label_pstfx+=pstfx_names[i_title]


                #plt.plot(list(rl(plot_vls)),plot_vls,label=label_names[i_label]+label_pstfx, marker = 'x', linestyle='--')
                if(i_label % len(plot_arr[1])  == 0 ):
                        curline_idx= (curline_idx+1)%len(lines_lst)
                plt.plot(list(rl(plot_vls)),plot_vls,label=label_names[i_label], linestyle=lines_lst[curline_idx], marker =
                        plotmarkers[(i_title + 2) % len(plotmarkers) ], c=color_iter[i_label],
                        alpha=1-(i_title+1)*i_label/(len(plot_arr[1])*len(plot_arr[1]))*.6, 
                        markersize = 20 -i_label/len(plot_arr[1])*8)
                #plt.plot(list(rl(plot_vls)),plot_vls,label=label_names[i_label], marker = cur_plotmarker, linestyle='--')

            #print("len plot_arr[0]] = ", len(plot_arr[0]))
            if len(plot_arr[0]) > 1:
                ax = plt.gca()
                h, l = ax.get_legend_handles_labels()
                ph = [plt.plot([],marker="", ls="")[0]]
                idx0=(i_title  )*len(plot_arr[1])
                idx1=(i_title+1)*len(plot_arr[1])
                addlabel=["{0}: {1}".format(plot_type[0],pstfx_names[i_title])] + (ncols-1)*[""]
                labelslist=[]
                handleslist=[]
                for i in range(ncols):
                    idxl=idx0 + i    *int((idx1 - idx0)/ncols)
                    idxr=idx0 + (i+1)*int((idx1 - idx0)/ncols)
                    labelslist = labelslist  + [addlabel[i]] + l[idxl:idxr]
                    handleslist= handleslist + [ph[0]]     + h[idxl:idxr]
                labels.append(labelslist)
                handles.append(handleslist)



        if(plot_type[2] == 'file'):
            rota = 0
        else :
            rota =-60
        plt.xticks(range(0,len(plot_vls)), xtick_names, rotation=-rota)

        #Set legend right with tiles
        if(len(plot_arr[0]) > 1):
            #This code is copied from a stackexchange post about adding titles to the plotlegends
            leg = plt.legend(handles[0], labels, ncol=ncols*len(plot_arr[0]))
            for vpack in leg._legend_handle_box.get_children():
                for hpack in vpack.get_children()[:1]:
                    hpack.get_children()[0].set_width(0)
            leg=plt.legend(flatten(handles),flatten(labels), loc="center left",bbox_to_anchor=(1.03, .5),
                           prop={'size': 225/max(len(plot_arr[1]),5)},ncol=ncols*len(plot_arr[0]),
                           markerfirst=False)
        else :
            plt.legend(loc="center left",bbox_to_anchor=(1.03, .5),
                       prop={'size': 225/max(len(plot_arr[1]),5)},ncol=ncols*len(plot_arr[0]))
        #filename=plot_path + "/" + "TestPlot" + str(i_title) + ".pdf"
        title_name, temp1 = get_name(plot_arr[0],plot_type[0],sepr="")
        label_name, temp2 = get_name(plot_arr[1], plot_type[1],sepr="")
        xtick_name, temp3 = get_name(plot_arr[2], plot_type[2],sepr="")
        #filename=plot_path + "/" + title_name + "-" + label_name[:10] + "-" + xtick_name + ".pdf"
        filename=plot_path + "/plot" + str(plot_idx) + title_name + "-" + label_name[:10] + "-" + xtick_name + ".pdf"
        plt.savefig(filename, bbox_inches='tight',format='pdf')
        plt.show()
        plt.close()

#%% Functions to get all indices to loop over



def get_iterators() :
    ############################################
    #### START INTERNAL FUNCTION DEFINITIONS ###
    ############################################

    #Internal functions are used to get index sets for the inp_var list


    def get_all_idx(mlist):
        if(type(mlist[0]) == list): #if list of lists
            return list( itertools.product(*[ list(itertools.product(i,i)) for i in map(rl,mlist) ]))
        else:
            return list(itertools.product(rl(mlist),rl(mlist)))

    def get_diag_idx(mlist):
        if(type(mlist[0]) == list): #if list of lists
            return list(itertools.product(*[ list(zip(*itertools.repeat(i,2))) for i in map(rl,mlist) ]))
        else: #if list
            return list(zip(*itertools.repeat(rl(mlist),2)))

    def get_lower_idx(mlist):
        if(type(mlist[0]) == list): #if list of lists
            return list(itertools.product(*[ list(itertools.combinations(i,2)) for i in map(rl,mlist) ]))
        else:
            return list(itertools.combinations(rl(mlist),2))

    def get_lower_idx_diag(mlist):
        if(type(mlist[0]) == list): #if list of lists
            return list(itertools.product(*[ list(itertools.combinations_with_replacement(i,2)) for i in map(rl,mlist) ]))
        else:
            return list(itertools.combinations_with_replacement(rl(mlist),2))
    ############################################
    ##### END INTERNAL FUNCTION DEFINITIONS ####
    ############################################


    # Getting iterators iter_{var,fxd} over the variables array (var_arr) for plotting
    symm=[ (i,) for i in get_diag_idx(inp_var[-1]) ]
    iter_var=[ (i,) for i  in get_lower_idx(inp_var[-1])]
    for i in range(len( inp_var)-1):
        iter_var=[ (i,*j) for i in get_lower_idx( inp_var[-i-2]) for j in symm   ] + [ (i,*j) for i in get_all_idx( inp_var[-i-2]) for j in iter_var ]
        symm=get_diag_idx( inp_var[-i-2:])

    iter_fxd =list(itertools.product(*map(rl,inp_fxd)))

    return iter_var,iter_fxd


#%% Functions to get all plot indices, given a maximum number of fixed and var

def get_plot_args(fxd_pp, var_pp):
    ############################################
    #### START INTERNAL FUNCTION DEFINITIONS ###
    ############################################
    from functools import reduce

    #Get list of factors of a number
    def factors(n):
        return set(reduce(list.__add__,
                    ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))

    #binary search
    def bin_search(arr, nr):
        i=0; j=len(arr)-1;
        while(np.abs(i-j)>1):
            new_ij=int((i+j)/2)
            if(arr[new_ij]<=nr):
                i=new_ij
            else:
                j=new_ij
        return arr[i]


    ############################################
    ##### END INTERNAL FUNCTION DEFINITIONS ####
    ############################################

    #Search for the best way to distribute the variables over plots. Not too many variables per plot
    if(len(iter_var) > var_pp):
        max_fct_var = bin_search(sorted(list(factors(len(iter_var)))),max_var_per_plot)
        iter_var_a = list(itertools.batched( iter_var, n=max_fct_var) )
    else:
        iter_var_a= [iter_var]

    if(len(iter_fxd) > fxd_pp):
        max_fct_fxd = bin_search(sorted(list(factors(len(iter_fxd)))),max_fxd_per_plot)
        iter_fxd_a = list(itertools.batched( iter_fxd, n=max_fct_fxd ) )
    else:
        iter_fxd_a = [iter_fxd]

    return iter_var_a, iter_fxd_a


#%% Defining loops to iterate over for the plot titles, plot labels and plot x_ticks

iter_var, iter_fxd = get_iterators()
print("Finished defining functions")
print("--------------")
print("Making plots")
#%% MAKE PLOTS WITH:
#  TITLE : Files (corr_\(...\))
#  LABELS: fixed context
#  XTICKS: variable context
max_fxd_per_plot=28
max_var_per_plot=10


plot_type_1=["file", "fxd", "var"]
iter_var_arg, iter_fxd_arg = get_plot_args(max_fxd_per_plot, max_var_per_plot)


print("Saving plots to " + path_plot_1)
plot_idx=0
for j in rl(iter_fxd_arg):
    for k in rl(iter_var_arg):
        if(len(iter_var_arg[k]) > 1):
            plot_loop_1=[ [[i], [*iter_fxd_arg[j]], [*iter_var_arg[k]]] for i in rl(files)]
            Make_plot_testing(plot_loop_1, plot_type_1, path_plot_1, plot_idx)
            plot_idx+=1



#%% MAKE PLOTS WITH:
#  TITLE : variable context
#  LABELS: fixed context
#  XTICKS: corr_{eta,sal,...}

#To set:
max_fxd_per_plot=28
max_var_per_plot=6


plot_type_2=["var", "fxd", "file"]
iter_var_arg, iter_fxd_arg = get_plot_args(max_fxd_per_plot, max_var_per_plot)

print("Saving plots to " + path_plot_2)
plot_idx=0
for j in rl(iter_fxd_arg):
    if(len(iter_fxd_arg[j]) > 1):
        for k in rl(iter_var_arg):
            plot_loop_2=[ [[*iter_var_arg[k]], [*iter_fxd_arg[j]], rl(files) ]]
            Make_plot_testing(plot_loop_2, plot_type_2,path_plot_2, plot_idx)
            plot_idx+=1
