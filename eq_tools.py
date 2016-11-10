import numpy as np
import matplotlib.pylab as plt

################################################################################
def get_earth_data(filename='Earth_density_and_earthquake_velocity_data_SMOOTHED.csv'):
    vals = np.loadtxt(filename,skiprows=1,dtype=float,unpack=True,delimiter=',')
    radii = vals[0]
    velocities = vals[3]

    return radii,velocities
################################################################################

################################################################################
def draw_earth(radii,alpha=0.05):
#colors = ['y','b','r','g']
    for i,r in enumerate(radii):
        circle = plt.Circle((0, 0), radius=r,alpha=alpha)#,fc=colors[i])
        plt.gca().add_patch(circle)
    #plt.xlim(-1.2,1.2)
    #plt.ylim(-1.2,1.2)
################################################################################

def normal(r,x,y):
    return r

def snells(theta0,v0,v1):
    # Check for critical angle
    '''
    theta_c = np.arcsin(v0/v1)
    #print "theta_c: ",theta_c
    if theta0>theta_c:
        print "Critical angle!: ",theta0,v1,v0,theta_c
        return None
    '''

    #else:
    sin_theta1 = np.sin(theta0)*(v1/v0)
    theta1 = np.arcsin(sin_theta1)
    if theta1 != theta1:
        # This is a nan from a critical angle issue
        return None
    #print theta0,v1,v0,sin_theta1,theta1
    return theta1

def mag(x,y):
    r = np.sqrt(x*x + y*y)
    return r

def rel_angle(x0,y0, x1,y1):
    mag0 = mag(x0,y0)
    mag1 = mag(x1,y1)
    
    cos_theta = (x0*x1 + y0*y1)/(mag0*mag1)
    
    theta = np.arccos(cos_theta)
    return theta

def linear(m,b,x):
    y = b + m*x
    return y

def radial_dist(x,y):
    r2 = x*x + y*y
    r = np.sqrt(r2)
    return r

def distance_traversed(xvals,yvals):
    dx = np.abs(xvals[0]-xvals[1])
    dy = np.abs(yvals[0]-yvals[1])
    return np.sqrt(dx*dx + dy*dy)

def vector_representation(x0,y0, x1,y1, norm=False):
    #theta = np.arctan2(y1-y0,x1-x0)
    dx = x1-x0
    dy = y1-y0
    vec_mag = 1.0
    if norm==True:
        vec_mag = mag(dx,dy)
    x = dx/vec_mag
    y = dy/vec_mag
    #y = vec_mag*np.sin(theta)
    return x,y

def new_y(x0,y0,x1,theta):
    y = (x1-x0)*np.tan(theta) + y0
    return y

def sgn(x):
    if x<0:
        return -1
    else:
        return +1 

def intersection(x1,y1, x2,y2, r):
    dx = x2-x1
    dy = y2-y1
    dr = mag(dx,dy)
    D = x1*y2 - x2*y1
    
    pts = []
    radical = r**2 * dr**2 - D**2 
    if radical<0:
        return None
    else:
        x = (D*dy + sgn(dy)*dx*np.sqrt(radical))/dr**2
        y = (-D*dx + np.abs(dy)*np.sqrt(radical))/dr**2
        pts.append(np.array([x,y]))
        if radical==0:
            return pts
        else:
            x = (D*dy + sgn(dy)*dx*-1*np.sqrt(radical))/dr**2
            y = (-D*dx + np.abs(dy)*-1*np.sqrt(radical))/dr**2
            pts.append(np.array([x,y]))
            return pts
        
def radial_pts(x,y,radius=1.0):
    theta = rel_angle(1.0, 0.0, x, y)
    rx = radius*np.cos(theta)
    ry = radius*np.sin(theta)
    
    rx *= x/np.abs(x)
    ry *= y/np.abs(y)
    
    return rx,ry

def trace_to_radius(x0,y0,angle,current_radius,new_radius,current_vel,new_vel):
    # x0 and y0 are the points where the ray starts
    # angle is the angle of that ray
    # radius is the radius that we are looking to see if it intercepts

    # return the point at which it intersects the new circle and 
    # the angle at which it enters or exits that radius
    
    #print "input: ",x0,y0,angle,current_radius,new_radius
    # Extend the line
    x = 1.0
    y = new_y(x0,y0,x,angle)
    #print x,y

    # See if it intersects our radius
    pts = intersection(x0,y0,x,y,new_radius)
    #print pts

    closest = None
    if pts is None:
        return None, None, angle
    
    elif pts is not None:
        #print "intersection pts: ",pts
        closest = None
        if len(pts)==1:
            closest = pts[0]
            return closest[0],closest[1],angle
        if len(pts)>1:
            d0 = mag(x0-pts[0][0],y0-pts[0][1])
            d1 = mag(x0-pts[1][0],y0-pts[1][1])
            if d0<d1:
                if new_radius==current_radius:
                    closest = pts[1]
                else:
                    #print "THIS CLOSEST"
                    closest = pts[0]

            else:
                if new_radius==current_radius:
                    closest = pts[0]
                else:
                    closest = pts[1]
                    
    #print "closest: ",closest
    ray = [[x0,closest[0]],[y0,closest[1]]]
    rd0 = radial_pts(closest[0],closest[1])
    
    # Next layer
    #print "TWO"
    vx,vy = vector_representation(ray[0][1],ray[1][1],x0,y0)
    cx = ray[0][1] # Circle x
    cy = ray[1][1] # Circle y
    t0 = rel_angle(cx,cy,vx,vy)

    t1 = snells(t0,current_vel,new_vel)
    if t1 is None:
        # This means the incoming angle was above the critical angle
        # and there is total internal reflection. So return some appropriate flags
        return -999, -999, -999

    norm_angle = rel_angle(cx,cy,1.0, 0.0)
    #print "norm: ",np.rad2deg(norm_angle)
    #print "t0: ",np.rad2deg(t0)
    #print "t1: ",np.rad2deg(t1)
    if new_radius<current_radius:
        #angle -= (t1-t0) # Change in angle is the difference between t0 and t1
        angle = np.pi - (np.pi - norm_angle - t1)
        #print "HERE!"
    else:
        radial_angle = rel_angle(1.0, 0.0, cx,cy)
        if cx>0 and cy<0:
            radial_angle = -np.abs(radial_angle)
        #print "radial angle: ",np.rad2deg(radial_angle)
        angle = radial_angle - np.abs(t1)
        '''
        if cx>0 and cy<0:
            angle = radial_angle - t1
        else:
            angle = radial_angle + t1
        '''
        #print "there"
    #angle = norm_angle - t1
    #print "new angle: ",np.rad2deg(angle)
    
    return closest[0],closest[1],angle



################################################################################
################################################################################
def propagate_rays(radii,velocities,origin,angles):

    # Convert to radians, move this to function???
    angles = np.deg2rad(angles)

    nradii = len(radii)
    nangles = len(angles)

    allrays = []
    tirs = []

    x0=0
    y0=0
    a=0
    for i in range(nangles):

        tir = False # Total Internal Reflection flag
        rays = []

        a = angles[i]
        #print "ANGLE ----- ",a,np.rad2deg(a)

        innermost = False
        for j in range(nradii):
            
            radius = j

            if j==0:
                x0 = origin[0]
                y0 = origin[1]


            # Innermost layer
            if j!=nradii-1:
                #print "A"
                x,y,angle = trace_to_radius(x0,y0, a, radii[j],radii[j+1], velocities[j], velocities[j+1])
            else:
                #print "B"
                x,y,angle = trace_to_radius(x0,y0, a, radii[j],radii[j], velocities[j], velocities[j-1])
                innermost = True

            if x==y and x==-999:
                # There was total internal reflection somewhere
                tir = True
                break

                
            if x is not None:
                ray = [[x0,x],[y0,y]]
                rays.append([ray,radius])
                a = angle
                x0 = x
                y0 = y
                #print "C ",j,a,np.rad2deg(a)

            else:
                # It missed the inner layer
                x,y,angle = trace_to_radius(x0,y0, a, radii[j],radii[j], velocities[j], velocities[j-1])
                ray = [[x0,x],[y0,y]]
                rays.append([ray,radius]) 
                a = angle
                x0 = x
                y0 = y
                #print "D ",j,a,np.rad2deg(a)
                break

            if innermost:
                #print "E"
                radius = j
                break
                
        #############################
        # Come out of the Earth!
        #############################
        if tir==False:
            #print "COMING OUT OF THE EARTH ",radius, a
            for j in range(radius,0,-1):
                
                radius = j

                #print "G ",j
                if j==0:
                    #print "0000000000000000000"
                    x,y,angle = trace_to_radius(x0,y0, a, radii[j],radii[j], velocities[j], velocities[j])
                else:
                    x,y,angle = trace_to_radius(x0,y0, a, radii[j],radii[j-1], velocities[j], velocities[j-1])
                    #print "H",x0,y0,a,x,y,angle
                    
                    
                if x==y and x==-999:
                    # There was total internal reflection somewhere
                    tir = True
                    break
                    
                #print "IIIIIIIIIIIIIIII"
                ray = [[x0,x],[y0,y]]
                rays.append([ray,radius]) 
                a = angle
                x0 = x
                y0 = y

        tirs.append(tir)
        allrays.append(rays)

    return allrays,tirs
        
################################################################################
################################################################################

def plot_velocities(radii,velocities):

    if len(radii) != len(velocities):
        print "Yikes!"
        print "You have %d entries in your list of radii..." % (len(radii))
        print "and have %d entries in your list of velocities." % (len(velocities))
        print 
        print "You should fix this before you go any further."
        print "You need the same number of entries in each."
        print
        print "radii:"
        print radii
        print "\nvelocities:"
        print velocities
        return

    plt.figure(figsize=(7,6))
    tempv = list(velocities)
    tempr = list(radii)
    tempv.append(tempv[-1])
    tempr.append(0.0)
    #plt.plot(tempv,tempr,'o-')
    npts = len(tempv)
    for i in range(0,npts-1):
        v0 =tempv[i]
        r0 =tempr[i]
        r1 =tempr[i+1]
        plt.plot([v0,v0],[r0,r1],'.-',markersize=1,color='k',linewidth=4)
    plt.ylabel("Radius (km)",fontsize=18)
    plt.xlabel("Speed of sound (km/s)",fontsize=18)
    plt.xlim(0,1.2*max(tempv))
    plt.ylim(0,1.2*max(tempr))

    plt.tight_layout()



################################################################################
################################################################################

def make_earthquake(radii,velocities,nrays=10,filename=None,real_earth=False):

    ############################################################################
    # Need to add on a thin shell for now, just for rays that would just
    # never enter a second layer and would just hit the surface.
    # NEED TO FIX THIS!!!!!!!!!!!!!!!!!
    ############################################################################
    if type(radii)==list and type(velocities)==list:
        r = np.zeros(len(radii)+1)
        v = np.zeros(len(radii)+1)

        r[1:] = radii[:]
        v[1:] = velocities[:]

        r[0] = 1.00001*r[1]
        if v[1]-v[2]>0:
            v[0] = 1.00001*v[1]
        else:
            v[0] = 0.99999*v[1]

        radii = r
        velocities = v


	# Earthquake origin (x/y-coordinates)
    origin = [0.0, radii[0]]

    raylo = -90.
    rayhi = -10.
    ray_intervals = np.abs(raylo - rayhi)/nrays
    #print ray_intervals
    angles = np.arange(raylo,rayhi,ray_intervals)

    # Propagate the rays!
    allrays,tirs = propagate_rays(radii,velocities,origin,angles)
    #print angles

    plt.figure(figsize=(12,6))
    plt.subplot(1,2,1)

    draw_earth(radii,alpha=0.05)

    times = []
    ang_distances = []

    for tir,rays in zip(tirs,allrays):
        if tir==False:
        #if len(rays)>1:
            #print "----"
            time = 0
            for rayinfo in rays:
                ray = rayinfo[0]
                radius = rayinfo[1]
                plt.plot(ray[0],ray[1],'k-',linewidth=0.5,alpha=1.0)
                #print radius,ray
                
                vel = velocities[radius]
                d = distance_traversed(ray[0],ray[1])
                time += d/vel
                #print "tot time, time, d, vel: ",time,d/vel,d,vel

            times.append(time)

            last_ray = rays[-1]
            rayinfo = last_ray[0]
            x = rayinfo[0][1]
            y = rayinfo[1][1]
            ang_distance = np.arctan2(y,x)
            ang_distances.append(ang_distance)
                
    plt.xlim(-1.2*radii[0],1.2*radii[0])
    plt.ylim(-1.2*radii[0],1.2*radii[0])
    plt.xlabel("Radius (km)",fontsize=18)
    plt.ylabel("Radius (km)",fontsize=18)


    plt.subplot(1,2,2)
    ang_distances = np.array(ang_distances)
    ang_distances = np.pi/2.0 - ang_distances

    plt.plot(np.rad2deg(ang_distances),times,'bo',label='Your model')

    xe,ye = None,None
    if real_earth:
        xe,ye = np.loadtxt('earth_500pts.csv',delimiter=',',unpack=True,dtype=float)
        #xe = xe.astype(float)
        #ye = ye.astype(float)
        plt.plot(xe,ye,'r+',markersize=10,label='Realistic Earth model')



    plt.xlabel("Angular distance (degrees)",fontsize=18)
    plt.ylabel("Arrival time (seconds)", fontsize=18)
    plt.ylim(0)
    plt.xlim(0)
    plt.legend(loc='upper left',numpoints=1)

    plt.tight_layout()

    if real_earth:
        #plt.plot(np.rad2deg(ang_distances),times,'bo',label='Your model')
        yinterp = np.interp(np.rad2deg(ang_distances),xe[::-1],ye[::-1])
        dy = times - yinterp
        #print dy
        #print times
        #print yinterp
        plt.figure(figsize=(12,6))
        plt.subplot(1,2,2)
        plt.plot([0.0, 180],[0,0],'k--')
        plt.plot(np.rad2deg(ang_distances),dy,'bo',label='Difference between your model and real Earth')
        plt.xlabel("Angular distance (degrees)",fontsize=18)
        plt.ylabel("Your model - real Earth (seconds)", fontsize=18)
        #plt.ylim(0)
        plt.xlim(0)
        #plt.legend(loc='upper left',numpoints=1)
        plt.tight_layout()





    if filename != None:
        a = np.array((np.rad2deg(ang_distances),times))
        np.savetxt(filename, a.transpose(), delimiter=",")
