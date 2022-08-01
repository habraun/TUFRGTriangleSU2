# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import h5py
from matplotlib import cm
import matplotlib


dirtype="local"
subdir="/Users/hbraun/Nextcloud/TUFRGTriangleSU2-main"

N=4
U=4.0
V1=0.0
V2=0.0
V3=0.0
J=0.0
t1=1.0
t2=0.0
t3=0.0
m=2.0
f1=1
f2=1
shell=2
GammaRes="false"
MRes="false"
KRes="false"





##########DISPERISON fcts#############
def boxcheck(kx,ky):

    b=0
    mtest=-(1.0/np.sqrt(3))
    xlim=2*np.pi/3
    ylim=4*np.pi/(np.sqrt(3)*3)

    if (np.abs(kx)<=xlim)&(np.abs(ky)<=ylim):
        mi= ( ( (2*np.pi/(3*np.sqrt(3))) - np.abs(ky)  )/( (2*np.pi/3) - np.abs(kx) + 1e-18))
        if mi>mtest:
            b=b+1
    return b

def S1(kx,ky,phaseinput):
    d1=np.array([1.5,np.sqrt(3)/2])
    d2=np.array([1.5,-np.sqrt(3)/2])
    d3=d1-d2
    k=np.array([kx,ky])
    phi=phaseinput*3.141/6
    #return np.exp(-1j*np.dot(k,d1)-1j*0.5235)+np.exp(-1j*np.dot(k,d2)-1j*0.5235)+np.exp(-1j*np.dot(k,d3)-1j*0.5235) + np.exp(1j*np.dot(k,d1)+1j*0.5235)+np.exp(1j*np.dot(k,d2)+1j*0.5235)+np.exp(1j*np.dot(k,d3)+1j*0.5235)
    return 2*( np.cos(np.dot(k,-d1) + phi) +np.cos(np.dot(k,d2) + phi) +np.cos(np.dot(k,d3) + phi))

def S2(kx,ky):
    d1=np.array([1.5,np.sqrt(3)/2])
    d2=np.array([1.5,-np.sqrt(3)/2])
    d3=d1-d2

    s1=d1+d2
    s2=d3+d1
    s3=d3-d2
    k=np.array([kx,ky])

    return np.exp(-1j*np.dot(k,s1))+np.exp(-1j*np.dot(k,s2))+np.exp(-1j*np.dot(k,s3)) + np.exp(1j*np.dot(k,s1))+np.exp(1j*np.dot(k,s2))+np.exp(1j*np.dot(k,s3))

def S3(kx,ky):
    d1=np.array([1.5,np.sqrt(3)/2])
    d2=np.array([1.5,-np.sqrt(3)/2])
    d3=d1-d2

    s1=d1*2
    s2=d2*2
    s3=d3*2

    k=np.array([kx,ky])

    return np.exp(-1j*np.dot(k,s1))+np.exp(-1j*np.dot(k,s2))+np.exp(-1j*np.dot(k,s3)) + np.exp(1j*np.dot(k,s1))+np.exp(1j*np.dot(k,s2))+np.exp(1j*np.dot(k,s3))



def dispersion_standard1(kx,ky,t1,t2,t3,mu,phaseinput):
    #return (2*t1*(np.cos(kx)+np.cos(ky))+4*t2*np.cos(kx)*np.cos(ky)+2*t3*(np.cos(2*kx)+np.cos(2*ky))) + mu
    return np.real(-t1*S1(kx,ky,phaseinput)-t2*S2(kx,ky)-t3*S3(kx,ky)-mu)





#%%
#########LOAD DATA FUNCTIONS
def idx2coords(coords,f1,f2,vertex):
    N=int(np.shape(coords)[1])
    #print(N)
    x=np.zeros(N).reshape(N)
    y=np.zeros(N).reshape(N)
    z=np.zeros(N).reshape(N)

    for k in range(0,N):
        x[k]=coords[0,k]
        y[k]=coords[1,k]
        z[k]=vertex[k,f2-1,f1-1]

    return x,y,z



def get_vertex_data(dirtype,subdir,f1,f2,N,U,V1,V2,V3,J,t1,t2,t3,m,GammaRes,MRes,KRes,shell):

    try:
        if dirtype=="cluster":
            filename = subdir+"/"+"triangle"+str(N)+str(U)+str(V1)+str(V2)+str(V3)+str(J)+str(t1)+str(t2)+str(t3)+str(m)+str(GammaRes)+str(MRes)+str(KRes)+str(shell)+"/"+"triangle"+str(N)+str(U)+str(V1)+str(V2)+str(V3)+str(J)+str(t1)+str(t2)+str(t3)+str(m)+str(GammaRes)+str(MRes)+str(KRes)+str(shell)+".h5"
        else:
            filename = subdir+"/"+"triangle"+str(N)+str(U)+str(V1)+str(V2)+str(V3)+str(J)+str(t1)+str(t2)+str(t3)+str(m)+str(GammaRes)+str(MRes)+str(KRes)+str(shell)+".h5"



        with h5py.File(filename, "r") as f:
            # List all groups
            #print("Keys: %s" % f.keys())
            #a_group_key = list(f.keys())[1]

            # Get the data
            idxlist = np.transpose(np.array((f["idxlist"])))
            coords = np.array((f["coords"]))

            pvertexv=  np.array(f["pv"])
            cvertexv=  np.array(f["cv"])
            dvertexv=  np.array(f["dv"])

            pmaxv=  np.array((f["pmaxv"]))
            cmaxv=  np.array((f["cmaxv"]))
            dmaxv=  np.array((f["dmaxv"]))

            pvertexw=  np.array(f["pw"])
            cvertexw=  np.array(f["cw"])
            dvertexw=  np.array(f["dw"])

            pmaxw=  np.array((f["pmaxw"]))
            cmaxw=  np.array((f["cmaxw"]))
            dmaxw=  np.array((f["dmaxw"]))

            Lambda= np.array((f["Lambda"]))
    except OSError:
        idxlist=0
        coords=0
        pvertexv=np.array([0])
        cvertexv=np.array([0])
        dvertexv=np.array([0])
        pmaxv=np.array([0])
        cmaxv=np.array([0])
        dmaxv=np.array([0])
        pvertexw=np.array([0])
        cvertexw=np.array([0])
        dvertexw=np.array([0])
        pmaxw=np.array([0])
        cmaxw=np.array([0])
        dmaxw=np.array([0])
        Lambda=np.array([2.5*10**(-6)])


    return idxlist,coords,pvertexv,cvertexv,dvertexv,pmaxv,cmaxv,dmaxv,pvertexw,cvertexw,dvertexw,pmaxw,cmaxw,dmaxw,Lambda



def get_gap_data(dirtype,subdir,N,U,V1,V2,V3,J,t1,t2,t3,m,GammaRes,MRes,KRes,shell):
    #filename = subdir+"/"+"square"+str(N)+str(U)+str(t1)+str(t2)+str(alpha)+str(m)+str(GammaRes)+str(MRes)+str(KRes)+str(shell)+"/"+"squareSC"+str(N)+str(U)+str(t1)+str(t2)+str(alpha)+str(m)+str(GammaRes)+str(MRes)+str(KRes)+str(shell)+".h5"
    #filename = subdir+"/"+"honeycombSC"+str(N)+str(U)+str(V1)+str(V2)+str(t1)+str(t2)+str(t3)+str(m)+str(GammaRes)+str(MRes)+str(KRes)+str(shell)+".h5"
    if dirtype=="cluster":
        filename = subdir+"/"+"triangle"+str(N)+str(U)+str(V1)+str(V2)+str(V3)+str(J)+str(t1)+str(t2)+str(t3)+str(m)+str(GammaRes)+str(MRes)+str(KRes)+str(shell)+"/"+"triangleSC"+str(N)+str(U)+str(V1)+str(V2)+str(V3)+str(J)+str(t1)+str(t2)+str(t3)+str(m)+str(GammaRes)+str(MRes)+str(KRes)+str(shell)+".h5"
    else:
        filename = subdir+"/"+"triangleSC"+str(N)+str(U)+str(V1)+str(V2)+str(V3)+str(J)+str(t1)+str(t2)+str(t3)+str(m)+str(GammaRes)+str(MRes)+str(KRes)+str(shell)+".h5"


    with h5py.File(filename, "r") as f:
        # List all groups
        #print("Keys: %s" % f.keys())
        #a_group_key = list(f.keys())[1]


        gapv =  np.array((f["leadingvecs1111v"]))
        valsv= np.array((f["leadingvals1111v"]))
        gapw =  np.array((f["leadingvecs1111w"]))
        valsw= np.array((f["leadingvals1111w"]))


        coords = np.array((f["coords"]))

        return gapv,valsv,gapw,valsw,coords

def get_bubble_data(dirtype,subdir,N,U,V1,V2,V3,J,t1,t2,t3,m,GammaRes,MRes,KRes,shell):

    if dirtype=="cluster":
        filename = subdir+"/"+"triangle"+str(N)+str(U)+str(V1)+str(V2)+str(V3)+str(J)+str(t1)+str(t2)+str(t3)+str(m)+str(GammaRes)+str(MRes)+str(KRes)+str(shell)+"/"+"triangle"+str(N)+str(U)+str(V1)+str(V2)+str(V3)+str(J)+str(t1)+str(t2)+str(t3)+str(m)+str(GammaRes)+str(MRes)+str(KRes)+str(shell)+".h5"
    else:
        filename = subdir+"/"+"triangle"+str(N)+str(U)+str(V1)+str(V2)+str(V3)+str(J)+str(t1)+str(t2)+str(t3)+str(m)+str(GammaRes)+str(MRes)+str(KRes)+str(shell)+".h5"


    with h5py.File(filename, "r") as f:
        # List all groups
        #print("Keys: %s" % f.keys())
        #a_group_key = list(f.keys())[1]

        # Get the data

        BubblesGamma=  np.array(f["BubblesGamma"])
        BubblesM=  np.array(f["BubblesM"])

        Lambda= np.array((f["Lambda"]))


    return BubblesGamma,BubblesM,Lambda
def get_dap_data(dirtype,subdir,N,U,V1,V2,V3,J,t1,t2,t3,m,phaseinput,GammaRes,MRes,KRes,shell):
    #filename = subdir+"/"+"square"+str(N)+str(U)+str(t1)+str(t2)+str(alpha)+str(m)+str(GammaRes)+str(MRes)+str(KRes)+str(shell)+"/"+"squareSC"+str(N)+str(U)+str(t1)+str(t2)+str(alpha)+str(m)+str(GammaRes)+str(MRes)+str(KRes)+str(shell)+".h5"
    #filename = subdir+"/"+"honeycombSC"+str(N)+str(U)+str(V1)+str(V2)+str(t1)+str(t2)+str(t3)+str(m)+str(GammaRes)+str(MRes)+str(KRes)+str(shell)+".h5"
    if dirtype=="cluster":
        filename = subdir+"/"+"triangle"+str(N)+str(U)+str(V1)+str(V2)+str(V3)+str(J)+str(t1)+str(t2)+str(t3)+str(m)+str(phaseinput)+str(GammaRes)+str(MRes)+str(KRes)+str(shell)+"/"+"triangleSC"+str(N)+str(U)+str(V1)+str(V2)+str(V3)+str(J)+str(t1)+str(t2)+str(t3)+str(m)+str(phaseinput)+str(GammaRes)+str(MRes)+str(KRes)+str(shell)+".h5"
    else:
        filename = subdir+"/"+"triangleD"+str(N)+str(U)+str(V1)+str(V2)+str(V3)+str(J)+str(t1)+str(t2)+str(t3)+str(m)+str(phaseinput)+str(GammaRes)+str(MRes)+str(KRes)+str(shell)+".h5"


    with h5py.File(filename, "r") as f:
        # List all groups
        #print("Keys: %s" % f.keys())
        #a_group_key = list(f.keys())[1]


        gapv =  np.array((f["leadingvecs1111v"]))
        valsw= np.array((f["leadingvals1111v"]))
        gapv =  np.array((f["leadingvecs1111w"]))
        valsw= np.array((f["leadingvals1111w"]))


        coords = np.array((f["coords"]))

        return gapv,valsv,gapw,valsw,coords
#%%
# set inset options
#mpl_axgrid = PyNULL()
#copy!(mpl_axgrid, pyimport_conda("mpl_toolkits.axes_grid1", "mpl_toolkits"))
import mpl_toolkits
import mpl_toolkits.axes_grid1.inset_locator

mpl_axgrid=mpl_toolkits.axes_grid1


zoomed_inset_axes = mpl_toolkits.axes_grid1.inset_locator.zoomed_inset_axes
inset_axes        = mpl_toolkits.axes_grid1.inset_locator.inset_axes
InsetPosition     = mpl_toolkits.axes_grid1.inset_locator.InsetPosition
AutoMinorLocator  = matplotlib.ticker.AutoMinorLocator
LocateAxes        = mpl_axgrid.make_axes_locatable
#%%
# set tex formatting
#rc("text", usetex = true)
matplotlib.rcParams["text.latex.preamble"] = [r"\usepackage{amsmath}", r"\usepackage{amssymb}", r"\usepackage{amsthm}"]

# set fonts
matplotlib.rcParams["font.family"] ="lmr"
matplotlib.rcParams["font.serif"] ="lmodern"
matplotlib.rcParams["font.size"] ="9"

plt.rcParams['text.usetex'] = True


# set axes
matplotlib.rcParams["axes.linewidth"]= 0.2
matplotlib.rcParams["axes.titlesize"]= 9
matplotlib.rcParams["axes.labelsize"]= 9
# set lines

matplotlib.rcParams["lines.linewidth"]  = 1.25

# set ticks

matplotlib.rcParams["xtick.major.width"]=0.2
matplotlib.rcParams["xtick.minor.width"]=0.1
matplotlib.rcParams["ytick.major.width"]=0.2
matplotlib.rcParams["ytick.minor.width"]=0.1

# set legend

matplotlib.rcParams["legend.fontsize"]=9
matplotlib.rcParams["legend.frameon"]=False
matplotlib.rcParams["legend.borderpad"]=0.1
matplotlib.rcParams["legend.borderaxespad"]=0.1
matplotlib.rcParams["legend.labelspacing"]=0.2
matplotlib.rcParams["legend.handletextpad"]=0.75
matplotlib.rcParams["legend.handlelength"]=0.9
matplotlib.rcParams["legend.columnspacing"]=1.0


# function to get proper figsize for paper plot
# columnwidth = 246.0 / 510.0 for single / double column
def get_figsize(columnwidth, wf = 1.0, hf = 0.5 * (np.sqrt(5.0) - 1.0)):

    # get figure width in pt
    fig_width_pt = columnwidth * wf

    # width in inches
    fig_width = fig_width_pt / 72.27

    # height in inches
    fig_height = fig_width * hf

    return (fig_width, fig_height)
#%% Load Data


idxlist,coords,pvertexv,cvertexv,dvertexv,pmaxv,cmaxv,dmaxv,pvertexw,cvertexw,dvertexw,pmaxw,cmaxw,dmaxw,Lambda=get_vertex_data(dirtype,subdir,f1,f2,N,U,V1,V2,V3,J,t1,t2,t3,m,GammaRes,MRes,KRes,shell)
BubblesGamma,BubblesM,Lambda=get_bubble_data(dirtype,subdir,N,U,V1,V2,V3,J,t1,t2,t3,m,GammaRes,MRes,KRes,shell)
gapv,valsv,gapw,valsw,coords             =get_gap_data(dirtype,subdir,N,U,V1,V2,V3,J,t1,t2,t3,m,GammaRes,MRes,KRes,shell)
#%% Plotting competition
fig2        = plt.figure(constrained_layout = True, figsize = get_figsize(246.0, wf = 1.2, hf = 0.45))
gs2         = matplotlib.gridspec.GridSpec(1, 1, figure = fig2)

ax1 = fig2.add_subplot(gs2[0,0])
ax1.set_xlabel(r"$T / t$")
ax1.set_ylabel(r"$\mathrm{max}(|X^{mn}(q)|)$")
ax1.plot(Lambda/abs(t1),pmaxv,label="$X=Pv$",color="deepskyblue",linestyle="solid")
ax1.plot(Lambda/abs(t1),dmaxv,label="$X=Dv$",color="orange",linestyle="solid")
ax1.plot(Lambda/abs(t1),cmaxv,label="$X=Cv$",color="violet",linestyle="solid")
ax1.plot(Lambda/abs(t1),pmaxw,label="$X=Pw$",color="deepskyblue",linestyle="dotted")
ax1.plot(Lambda/abs(t1),dmaxw,label="$X=Dw$",color="orange",linestyle="dotted")
ax1.plot(Lambda/abs(t1),cmaxw,label="$X=Cw$",color="violet",linestyle="dotted")
ax1.set_xscale("log")
#ax1.set_xticks([1e-3, 1e-2, 1e-1, 1e0, 1e1])
ax1.set_yticks([0.0, 9.0, 18.0, 27.0])
ax1.set_yticklabels([r"$0$", r"$1W$", r"$2W$", r"$3W$"])
ax1.legend(ncol = 1, loc = "upper right", handlelength = 1.25)
gs2.tight_layout(fig2, rect = (0.0, 0.0, 1.1, 1.0), pad = 4.0)

### vertices


for vert in range(0,6):
    fig1        = plt.figure(constrained_layout = True, figsize = get_figsize(246.0, wf = 1.2, hf = 0.45))
    gs1         = matplotlib.gridspec.GridSpec(1, 1, figure = fig1)
    channel=["PV","CV","DV","PW","CW","DW"][vert]
    ax2  = fig1.add_subplot(gs1[0,0])

    xi = np.linspace(-2.7,2.7,500)
    yi = np.linspace(-2.7,2.7,500)
    if channel=="CV":
        inputvertex=cvertexv
    elif channel=="DV":
        inputvertex=dvertexv
    elif channel=="PV":
        inputvertex=pvertexv
    elif channel=="CW":
        inputvertex=cvertexw
    elif channel=="DW":
        inputvertex=dvertexw
    elif channel=="PW":
        inputvertex=pvertexw

    x,y,z=idx2coords(coords,f1,f2,inputvertex)
    zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')

    ax2.title.set_text(r"$"+str(channel)+"^{"+str(f1)+","+str(f2)+"}(q)$")

    #ax2.set_xlabel(r"$\phi_1$")
    #ax2.set_ylabel(r"$\phi_2$")

    ax2.set_xticks([-6.28/3, 0.0 , 6.28/3])
    ax2.set_yticks([-6.28/3, 0.0 , 6.28/3])

    ax2.set_xticklabels([r"$\frac{2\pi}{\sqrt{3}}$", r"$0$", r"$\frac{2\pi}{\sqrt{3}}$"])
    ax2.set_yticklabels([r"$\frac{2\pi}{\sqrt{3}}$", r"$0$", r"$\frac{2\pi}{\sqrt{3}}$"])

    ax2.tick_params(
        axis="both",          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False,
        left=False,
        labelleft=False) # labels along the bottom edge are off


    ax2.contourf(xi,yi,zi,50,cmap="spring")
    #ax2.text(0.92 * np.pi, 1.75 * np.pi, r"$\phi_{3} = 55 \Delta \phi$", color = "black", fontsize = 8.0)
    ax2.set_xlabel(r"$k_x$")
    ax2.set_ylabel(r"$k_y$")
    divider = LocateAxes(ax2)
    cax     = divider.append_axes("right", size = "7%", pad = 0.05)
    cbar = plt.colorbar(ax2.contourf(xi,yi,zi,25,cmap="spring"), cax = cax, orientation = "vertical")
    #cbar.set_ticks([0,-5,-10,-15,-20,-25,30])
    #cbar.set_ticklabels([r"$0$", r"$-5$", r"$-10$",r"$-15$",r"$-25$",r"$-30$"])
    ax2.set_aspect(1.0)
    fig1.savefig(subdir+"/vertices"+str(channel)+".pdf",dpi=300,bbox_inches="tight")

gs1.tight_layout(fig1, rect = (0.0, 0.0, 1.1, 1.0), pad = 4.0)
### Bubbles
fig3        = plt.figure(constrained_layout = True, figsize = get_figsize(246.0, wf = 1.2, hf = 0.45))
gs3         = matplotlib.gridspec.GridSpec(1, 1, figure = fig3)
ax3         = fig3.add_subplot(gs3[0, 0])


phm=BubblesM[:,0]
ppm=BubblesM[:,1]
phg=BubblesGamma[:,0]
ppg=BubblesGamma[:,1]

L=(len(Lambda)-1)

ax3.plot(Lambda[:L],np.multiply(Lambda[:L],phm[:L]), label="ph(M)")
ax3.plot(Lambda[:L],np.multiply(Lambda[:L],-ppm[:L]), label="pp(M)")
ax3.plot(Lambda[:L],np.multiply(Lambda[:L],phg[:L]), label="ph(0)")
ax3.plot(Lambda[:L],np.multiply(Lambda[:L],-ppg[:L]), label="pp(0)")

ax3.set_ylabel(r"$T \cdot$ Bubble Value")
ax3.set_xlabel(r"$T / t$")
ax3.set_xscale("log")
ax3.legend(loc=1)

#gaps
fig4        = plt.figure(constrained_layout = True, figsize = get_figsize(246.0, wf = 1.2, hf = 1.2))
gs4         = matplotlib.gridspec.GridSpec(2, 2, figure = fig4)



for k in range(0,4):
    Nres=np.shape(coords)[1]
    ax4         = fig4.add_subplot(gs4[3-k])

    a=gapv[k,:]
    x=np.zeros(Nres).reshape(Nres)
    y=np.zeros(Nres).reshape(Nres)
    z=np.zeros(Nres).reshape(Nres)

    for i in range(0,Nres):
        x[i]=coords[0][i]
        y[i]=coords[1][i]
        z[i]=a[i]

    xi = np.linspace(-1*np.pi,1*np.pi,800)
    yi = np.linspace(-1*np.pi,1*np.pi,800)
        # grid the data.
    zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='linear')
    ax4.set_title(np.round(valsv[k],4))
    #ax4.contour(xi,yi,zi,levels = [-0.000000000000001,0.000000000000001],colors=('k',),linestyles=('-',),linewidths=(2,))
    ax4.contourf(xi,yi,zi,40,cmap="RdBu")
    ax4.set_aspect(1.0)
    ax4.set_xlabel(r"$k_x$")
    ax4.set_ylabel(r"$k_y$")
    ax4.tick_params(
        axis="both",          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False,
        left=False,
        labelleft=False) # labels along the bottom edge are off
    #cb=plt.colorbar(ax4.contourf(xi,yi,zi,20,cmap=plt.cm.jet),ax=ax4)

gs4.tight_layout(fig4, rect = (0.0, 0.0, 1.1, 1.0), pad = 4.0)


fig1.savefig(subdir+"/vertices.pdf",dpi=300,bbox_inches="tight")
fig2.savefig(subdir+"/flow.pdf",dpi=300,bbox_inches="tight")
fig3.savefig(subdir+"/bubbles.pdf",dpi=300,bbox_inches="tight")
fig4.savefig(subdir+"/gaps.pdf",dpi=300,bbox_inches="tight")
