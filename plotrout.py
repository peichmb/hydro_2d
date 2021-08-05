import grid_common as gr
import params_common as par
import matplotlib.pyplot as plt
import numpy as np

first=True
fig=0
ax=0
im1=0
cnt=0
txt=0
point=0
im2=0
txt2=0
fig2=0
ax2=0

def plotrout(vx,vz,rho,pres,time,it):
    global im1,first,fig,ax,cnt,txt,point,im2,txt2,fig2
    if first:
        first=False        
        plt.ion()        
        fig=plt.figure()
        ax=fig.add_subplot(111)
        #im1=ax.imshow(np.transpose(rho),origin='lower')
        ax.set_ylim([-2*par.amp,2*par.amp])
        im1,=ax.plot(gr.xx,(pres[:,0]-par.pres00)/par.pres00)
        txt=plt.text(3.2,0.000225,'t = '+str(round(time,5)))
        point,=ax.plot(1.0,0.00025,marker='x',c='black')
        plt.savefig('img/img'+str(cnt)+'.png')
        #fig2=plt.figure()
        #ax2=fig2.add_subplot(111)
        #ax2.set_ylim([-0.000035,0.000035])
        #ax2.set_xlim([-800,800])
        #fou=np.fft.fft((pres[1:-1,0]-par.pres00)/par.pres00)
        #fou2=fou.copy()
        #fou2[gr.npx/2-1:]=fou[0:gr.npx/2-1]
        #fou2[0:gr.npx/2-2]=fou[gr.npx/2:]
        #fou2=1./(gr.npx-2)*gr.xlen*fou2
        #kk=np.empty(gr.npx-2,dtype=float)    
        #for i in range(gr.npx-2):
        #    kk[i]=2*np.pi/(gr.xf-gr.x0)*(i-(gr.npx-2)/2)
        #im2,=plt.plot(kk,fou2)
        #im3,=plt.plot(kk,fou2,'r--')
        #txt2=plt.text(500,0.000022,'t = '+str(round(time,5)))    

    else:
        import exp_ctrl
        def cs(x):
            return np.sqrt(par.gam*par.p00/(par.rho00*np.exp((x-par.xpack0)/3.)))
        #im1.set_data(np.transpose(rho))
        cnt=cnt+1
        im1.set_ydata((pres[:,0]-par.pres00)/par.pres00)
        txt.set_text('t = '+str(round(time,5)))
        point.set_xdata(exp_ctrl.x0)
        plt.savefig('img/img'+str(cnt)+'.png')
        fig.canvas.draw()
        
        #fou=np.fft.fft((pres[1:-1,0]-par.pres00)/par.pres00)
        #fou2=fou.copy()
        #fou2[gr.npx/2-1:]=fou[0:gr.npx/2-1]
        #fou2[0:gr.npx/2-2]=fou[gr.npx/2:]
        #fou2=1./(gr.npx-2)*gr.xlen*fou2
        #im2.set_ydata(fou2)
        #txt2.set_text('t = '+str(round(time,5)))
        #plt.savefig('img/fou'+str(cnt)+'.png')
        #fig2.canvas.draw()
        

    
    
   
    

