def compute_dists(rjb, mag, hypo_defined, hypo_depth, rake, dip_defined, dip_input, upper_sd, lower_sd, azimuth):
    import numpy as np

    if(hypo_defined==1):
        Z_hyp=hypo_depth
    else:
        if (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
            Z_hyp=5.63+0.68*mag
        else:
            Z_hyp=11.24-0.2*mag

    if(dip_defined==1):
        dip=dip_input
    else:
        if (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
            dip=90
        elif rake > 0:
            dip=40
        else:
            dip=50

    if (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
         # strike slip
        width= 10.0 ** (-0.76 + 0.27 *mag)
    elif rake > 0:
        # thrust/reverse
        width= 10.0 ** (-1.61 + 0.41 *mag)
    else:
        # normal
        width= 10.0 ** (-1.14 + 0.35 *mag)

    source_vertical_width=width*np.sin(np.radians(dip))
    ztor=max(Z_hyp-0.6*source_vertical_width,upper_sd)
    if((ztor+source_vertical_width)>lower_sd):
        source_vertical_width=lower_sd-ztor
        width=source_vertical_width/np.sin(np.radians(dip))
    if(rjb==0):
        rx=0.5*width*np.cos(np.radians(dip))
    else:
        if(dip==90):
            rx=rjb*np.sin(np.radians(azimuth))
        else:
            if (azimuth>=0 and azimuth<90) or (azimuth>90 and azimuth<=180):
                if(rjb*np.abs(np.tan(np.radians(azimuth)))<=width*np.cos(np.radians(dip))):
                    rx=rjb*np.abs(np.tan(np.radians(azimuth)))
                else:
                    rx=rjb*np.tan(np.radians(azimuth))*np.cos(np.radians(azimuth)-np.arcsin(width*np.cos(np.radians(dip))*np.cos(np.radians(azimuth))/rjb))
            elif (azimuth==90): #we assume that Rjb>0 
                rx=rjb+width*np.cos(np.radians(dip))
            else:
                rx=rjb*np.sin(np.radians(azimuth))

    if(azimuth==90 or azimuth==-90):
        ry=0
    elif(azimuth==0 or azimuth==180 or azimuth==-180):
        ry=rjb
    else:
        ry=np.abs(rx*1./np.tan(np.radians(azimuth)))

    if(dip==90):
        rrup=np.sqrt(np.square(rjb)+np.square(ztor))
    else:
        if(rx<ztor*np.tan(np.radians(dip))):
            rrup1=np.sqrt(np.square(rx)+np.square(ztor))
        if(rx>=ztor*np.tan(np.radians(dip)) and rx<=ztor*np.tan(np.radians(dip))+width*1./np.cos(np.radians(dip))):
            rrup1=rx*np.sin(np.radians(dip))+ztor*np.cos(np.radians(dip))
        if(rx>ztor*np.tan(np.radians(dip))+width*1./np.cos(np.radians(dip))):
            rrup1=np.sqrt(np.square(rx-width*np.cos(np.radians(dip)))+np.square(ztor+width*np.sin(np.radians(dip))))
        rrup=np.sqrt(np.square(rrup1)+np.square(ry))
    return rx,rrup,width,ztor,dip,ry,Z_hyp
