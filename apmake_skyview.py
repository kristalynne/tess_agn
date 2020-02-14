#Program outputs a list of pixel coordinates defining the extraction mask, marked as rows and descending columns. It also outputs
#a matrix of size (15,15) by default, but editable in the xsize,ysize parameters, with zeros indicating pixels not to be
#included in the aperture and ones indicating pixels that are to be included. This can be provided to eleanor (and other pipelines) as a custom aperture.

import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from astropy import units as u
import numpy as np
from astroquery.skyview import SkyView
from astroquery.simbad import Simbad
from astroquery.mast import Tesscut
from astropy.coordinates import SkyCoord
from astropy.coordinates import get_icrs_coordinates
from astropy.wcs import WCS

#Parameters for customization:

input_is_coords = 1     #Set equal to 1 the method you wish to use to specify your input method.
input_is_name = 0

dss_pixel_size = 100    #Choose the total number of pixels in the DSS image (it will be square).

xsize = 15              #X-dimension of the output matrix. This should match whatever your intended pipeline assumes.
ysize = 15              #Y-dimension of the output matrix. This should match whatever your intended pipeline assumes.

##########################################################
#Define function to record the positions of clicks in the pixel array image for the extraction mask.
def onclick(event):
    
    global ix,iy
    ix,iy = int(round(event.xdata)),int(round(event.ydata))
    
    global coords
    coords.append((ix,iy))

    plt.plot(event.xdata,event.ydata,marker=u"$\u2713$",color='goldenrod',markersize=10)
    fig.canvas.draw()
    
    print(ix,iy)

#Define function to obtain the tesscuts for all sectors observed, creating a list of the hdu objects.
def get_multi_cutout_lst(coords, sz=(15,15)):

    lst_hdu_elements = Tesscut.get_cutouts(coordinates=coords, size=sz)
    return(lst_hdu_elements)

#Define function to read the flux array image and the header for each hdu object obtained by get_multi_cutout_lst
def get_tpf_single_sector(hdu_element, num=1000):
    
    tpf_cutout = hdu_element[1].data[num][4]
    header_info = hdu_element[2].header
    hdu_element.close()
    return(tpf_cutout, header_info)

#Define function to turn the list of selected aperture pixels into a matrix of zeros and ones, where the ones denote the pixels included in the extraction mask.
def make_custom_aperture(pixel_list,matrix_size = (15,15)):
    '''
    Turn the output pixel lists into matrices with zeros at pixels outside the extraction aperture and ones otherwise.
    '''
    custom_ap = np.zeros(matrix_size)
    for row in pixel_list:
        xval = int(row[0])
        yval = int(row[1])
        custom_ap[yval,xval] = 1
    return(custom_ap)

##########################################################

#Module for manual entry of object name.

if input_is_name == 1:
    source_common_name = input('Common name of source: ')
    source_coordinates = get_icrs_coordinates(source_common_name)

#Module for manual entry of target coordinates.
if input_is_coords ==1:

    source_coords = input('Source coordinates as ra_deg,dec_deg: ')
    source_ra = float(source_coords.split(',')[0])
    source_dec = float(source_coords.split(',')[1])
    source_coordinates = SkyCoord(source_ra, source_dec, unit = u.deg)




#Execute astroquery for the DSS image with the requested number of pixels.

if input_is_coords ==1:
    dss_image = SkyView.get_images(position=source_coords,survey='DSS',pixels=str(dss_pixel_size)) #Use this line for the dss query if the input is manual coor#Use this line for thedinates in decimal degrees.
if input_is_name ==1:
    dss_image = SkyView.get_images(position=source_coordinates,survey='DSS',pixels=str(dss_pixel_size)) #Use this line for the dss query if the input is the source's common name.

dss_image_data = dss_image[0][0].data
dss_pixmin = np.min(dss_image_data)
dss_pixmax = np.max(dss_image_data)
wcs_dss = WCS(dss_image[0][0].header)


#Create the sector masks by allowing user to click the desired pixels, one sector at a time.

cutlist = get_multi_cutout_lst(source_coordinates, sz = (xsize,ysize))

for i in range(0,len(cutlist)):
    
    tesscut_fluxdata,tesscut_header = get_tpf_single_sector(cutlist[i])
    sector_number = tesscut_header['SECTOR']
    
    wcs_tesscut = WCS(tesscut_header)
    image = tesscut_fluxdata
    
    print('Sector '+str(sector_number))
    print('Min pixel value: '+str(np.min(image)))
    print('Max pixel value: '+str(np.max(image)))
    print('Mean pixel value: '+str(np.mean(image)))
    print('Choose object mask pixels.')
    
    pixmin = np.min(image)
    pixmax = np.max(image)
    pixmean = np.mean(image)
    
    #Plot the pixel array and call the function for click-selecting the masks.
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111,projection=wcs_tesscut)
    
    if pixmax > 40000:                                              #Rescales the image contrast for fields that include or don't include bright stars.
        ax.imshow(image,vmin=pixmin,vmax=0.005*pixmax+pixmean)
    elif pixmax > 20000:
        ax.imshow(image,vmin=pixmin,vmax=0.01*pixmax+pixmean)
    else:
        ax.imshow(image,vmin=pixmin,vmax=0.1*pixmax+pixmean)
    
    x1 = -0.5
    x2 = x1 + xsize
    
    y1 = -0.5
    y2 = y1 + ysize
        
    if dss_pixmax < 20000:      #Plots the DSS contours on top of the TESS image. There are two different levels in case the DSS image contains a very bright source, to keep the contour levels reasonable.
        ax.contour(dss_image[0][0].data,transform=ax.get_transform(wcs_dss),levels=[0.3*dss_pixmax,0.5*dss_pixmax,0.75*dss_pixmax],colors='white',alpha=0.9)
        ax.set_xlim(x1,x2)
        ax.set_ylim(y1,y2)
    else:
        ax.contour(dss_image[0][0].data,transform=ax.get_transform(wcs_dss),levels=[0.3*dss_pixmax,0.5*dss_pixmax,0.75*dss_pixmax],colors='white',alpha=0.9)
        ax.set_xlim(x1,x2)
        ax.set_ylim(y1,y2)
        
    ax.scatter(50.5,50.5,transform = ax.get_transform(wcs_dss),marker='x',color='k') #plots a black X on the center.
    
    plt.title('Define extraction pixels:')
    
    coords = []
        
    cid = fig.canvas.mpl_connect('button_press_event',onclick)
    plt.show()
    plt.close()
    
    np.savetxt('sector'+str(sector_number)+'_extraction_aperture_'+source_common_name+'.dat',coords)
    
    matrix_aperture = make_custom_aperture(coords,matrix_size = (xsize,ysize))
    
    np.savetxt('sector'+str(sector_number)+'_matrix_aperture_'+source_common_name+'.dat',matrix_aperture)
   


    



    
    



