#Author: Eric Shore
#Purpose: To calculate the Strehl Ratio of an AO System
import numpy as np
from scipy.special import j1
import matplotlib.pyplot as plt
from matplotlib import rc

#Figure text parameters
rc('text',usetex=True)
rc('font',family='serif')
rc('font',serif='cm')
rc('font',size=28)

#Function for determining the (normalized) maximum intensity for the actual PSF
def actual():
    data = np.zeros([330,320]) #matrix to hold data
    #integrate over 10 files
    for i in range(10):
        temp = plt.imread("strehl_data/no_turb_2017-03-28-170023-000"+str(i)+".tif") #temporary array to hold data
	data=data+temp
    #set up a grid of X and Y positions based on the shape of the data
    h,w = np.shape(data)
    X,Y = np.meshgrid(np.arange(0,w),np.arange(0,h))
    #Initial guess of centroid
    centre1 = 227
    centre2 = 166
    print centre1,centre2
    #Based on Initial Guess, determine centroid
    tol=7     #tolerance on either side of centroid for ignoring noise
    cx = np.sum(X[centre1-tol:centre1+tol+1,centre2-tol:centre2+tol+1]*data[centre1-tol:centre1+tol+1,centre2-tol:centre2+tol+1])/float(np.sum(data[centre1-tol:centre1+tol+1,centre2-tol:centre2+tol+1]))
    cy = np.sum(Y[centre1-tol:centre1+tol+1,centre2-tol:centre2+tol+1]*data[centre1-tol:centre1+tol+1,centre2-tol:centre2+tol+1])/float(np.sum(data[centre1-tol:centre1+tol+1,centre2-tol:centre2+tol+1]))
    #cx = np.sum(X[centre1-tol:centre1+tol,centre2-tol:centre2+tol]*data[centre1-tol:centre1+tol,centre2-tol:centre2+tol])/np.sum(data[centre1-tol:centre1+tol,centre2-tol:centre2+tol])
    #cy = np.sum(Y[centre1-tol:centre1+tol,centre2-tol:centre2+tol]*data[centre1-tol:centre1+tol,centre2-tol:centre2+tol])/np.sum(data[centre1-tol:centre1+tol,centre2-tol:centre2+tol])
    print cx,cy
    #determine normalized central Brightness
    I_real = data[cx-tol:cx+tol+1,cy-tol:cy+tol+1] #array holds entire PSF (above noise threshold)
    cp =  data[cx,cy]/float(np.sum(I_real)) #intensity of centre pixel normalized to above region
    I_real/=float(np.sum(I_real)) #normalized PSF

    #Plot observed PSF
    plt.figure()
    plt.imshow(data,interpolation="nearest",origin="Lower")
    plt.title("Observed")
    plt.colorbar()
    return [cp,I_real]


#Function for determining the expected PSF if the system were diffraction limited
def ideal():
    #subpixel scale
    xscale = 10
    yscale = 10
    #Noise limited threshlod of PSF
    x_room = 7
    y_room = 7
    #Create arrays for X and Y positions of subpixels
    x = np.arange(-(x_room+0.5)*xscale,(x_room+.5)*xscale)*pixel_size/xscale #in m from central position
    y = np.arange(-(y_room+0.5)*yscale,(y_room+.5)*yscale)*pixel_size/yscale # in m from central poisition
    
    #determine expected PSF
    r = np.zeros([len(x),len(y)])#array to hold distances from centre of PSF
    u = np.zeros(np.shape(r))#array to hold Bessel Function variable
    I = np.zeros(np.shape(r))#array to hold intensities
    #loop through every pixel
    for i in range(len(x)):
        for j in range(len(y)):
            #determine distance to centre
            r[i,j] = (x[i]**2+y[j]**2)**0.5
            #determine intensity at every pixel
            if r[i,j] == 0: #if at centre set I to 1
                I[i,j] = 1.
            else: #otherwise calculate I
                u[i,j] = np.pi*d/wavelength * (r[i,j]/(r[i,j]**2+f**2)**0.5)
                I[i,j] = (2*j1(u[i,j])/(u[i,j]))**2
    #normalize I
    I/=np.sum(I)

    #caclulate Intensity at every pixel (from intensity at every subpixel)
    I2 = np.zeros([len(x)/xscale,len(y)/yscale])
    #loop through every pixel
    for i in range(len(x)/xscale):
        for j in range(len(y)/yscale):
            I2[i,j] = np.sum(I[xscale*i:xscale*(i+1),yscale*j:yscale*(j+1)]) #intensity just the sum of all intensities in every subpixel 


    #Plot PSF of subpixels and pixels
    plt.figure()
    plt.imshow(I,interpolation="nearest",origin="Lower")
    plt.colorbar()
    
    plt.figure()
    plt.imshow(I2,interpolation="nearest",origin="Lower")
    plt.colorbar()
    #determine maximum (Normalized) Intensity
    cp = np.amax(I2)
    print "ideal = ", np.amax(I2)
    return[ cp,I2]


#Defrine System Parameters
f = .300 #focal length in m
d = .013 #aperture diameter in m
wavelength = 6.29e-7 #in m
pixel_size = 5.86e-6 # in m

#determine ideal and actual central Intensities
cp_ideal,I_ideal = ideal()f = .300 #focal length in m
d = .013 #aperture diameter in m
wavelength = 6.29e-7 #in m
pixel_size = 5.86e-6 # in m

cp_actual,I_actual = actual()
print "actual =",cp_actual
#Calculate Strehl Ratio
strehl = cp_actual/cp_ideal
print strehl

#Plot Ideal and Actual PSFs
plt.figure()
plt.subplot(121)
plt.imshow(I_ideal,interpolation="nearest",origin="Lower",vmin=0,vmax=0.131)
plt.ylabel("Y pixel")
plt.xlabel("X pixel")
plt.title("Ideal")
plt.colorbar()
plt.subplot(122)
plt.imshow(I_actual,interpolation="nearest",origin="Lower",vmin=0,vmax=0.131)
plt.xlabel("X pixel")
plt.title("Actual")
plt.colorbar()





plt.show()


