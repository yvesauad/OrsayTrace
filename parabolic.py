import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider
import time
from tqdm import tqdm

class photon():
    def __init__(self, pos, normal, intensity=1.):
        self.pos = pos
        self.normal = normal / numpy.linalg.norm(normal)
        self.intensity = 1. / (numpy.exp(pos[0]**2) + numpy.exp(pos[1]**2))
        self.n = 1.00
        self.last_surface = [0, 0, 0]
        self.refraction_count = 0
        self.reflection_count = 0

    def set_attr(self, values):
        self.pos, self.normal, self.intensity, self.n, self.last_surface, self.refraction_count, self.reflection_count = values
    
    def get_attr(self):
        return self.pos, self.normal, self.intensity, self.n, self.last_surface, self.refraction_count, self.reflection_count


    def reflection(self):
        inc = self.normal
        inc = inc / numpy.linalg.norm(inc)
        sur_normal = self.last_surface
        sur_normal = sur_normal / numpy.linalg.norm(sur_normal)
        cos_ang1 = -numpy.dot(inc, sur_normal)
        if cos_ang1<0:
            sur_normal = -sur_normal
            cos_ang1 = -numpy.dot(inc, sur_normal)
 
        refl = inc + 2*cos_ang1*sur_normal
        refl = refl / numpy.linalg.norm(refl)
        self.reflection_count+=1

        self.normal = refl
        return True
    
    def refraction(self, n2):
        r = self.n/n2
        inc = self.normal
        inc = inc / numpy.linalg.norm(inc)
        sur_normal = self.last_surface
        sur_normal = sur_normal / numpy.linalg.norm(sur_normal)
        cos_ang1 = -numpy.dot(inc, sur_normal)
        if cos_ang1<0:
            sur_normal = -sur_normal #adjust surface normal
            cos_ang1 = -numpy.dot(inc, sur_normal)

        sin_ang1 = (1 - cos_ang1**2)**0.5
        sin_ang2 = r * sin_ang1

        if sin_ang2>1:
            #refl = inc + 2*cos_ang1*sur_normal
            #refl = refl / numpy.linalg.norm(refl)
            #self.reflection_count+=1
            #self.normal = refl
            print('***WARNING***: TIR.')
            refr = r*inc + (r*cos_ang1 - numpy.sqrt(1-r**2*(1-cos_ang1**2)))*sur_normal
            refr = refr / numpy.linalg.norm(refr)
            self.refraction_count+=1
            self.normal = refr
        else:
            refr = r*inc + (r*cos_ang1 - numpy.sqrt(1-r**2*(1-cos_ang1**2)))*sur_normal
            refr = refr / numpy.linalg.norm(refr)
            self.refraction_count+=1
            self.normal = refr
        
        self.n = n2 #now photon index of refraction
        return True

    def move(self, value):
        self.pos = self.pos + self.normal*value
        return True

    def update(self, value):
        sur_normal, index = value
        if sur_normal != [0, 0, 0]:
            self.last_surface = sur_normal
        
        if index!=self.n and index>0:
            self.refraction(index)
        elif index==-1:
            self.reflection()

        return True

class Simu:
    def __init__(self, x, y, z, res):
        self.size = [x, y, z]
        self.res = res
        self.ss = self.res / 2 #sub sampling for creating elements
        self.grid = (int(x/res), int(y/res), int(z/res))
        self.index = numpy.ones(self.grid)
        self.normal = numpy.asarray([numpy.zeros(self.grid), numpy.zeros(self.grid), numpy.zeros(self.grid)])
        self.photons = list()
        self.good_photons = list()
        self.plans=[list(), list(), list()]
        self.x = numpy.linspace(0, int(x/res))
        self.y = numpy.linspace(0, int(y/res))
        self.x, self.y = numpy.meshgrid(self.x, self.y)
        self.saved_photons=[list(), list(), list()]

    def assign_normal(self, index, normal):
        self.normal[0][index[0], index[1], index[2]] = normal[0]
        self.normal[1][index[0], index[1], index[2]] = normal[1]
        self.normal[2][index[0], index[1], index[2]] = normal[2]

    def assign_n(self, index, value):
        self.index[index[0], index[1], index[2]] = value

    def pos_to_grid_1D(self, value, index):
        return int((value/(self.size[index]/2.)+1)/2*self.grid[index])
        

    def pos_to_grid(self, pos):
        #this returns relative current photon position from [0:self.grid]

        indexes = [int((pos[i]/(self.size[i]/2.)+1)/2*self.grid[i]) for i in range(3)]
        return indexes

    def index_from_pos(self, pos):
        #this returns current index of refraction of a given photon in a given position

        indexes = [int((pos[i]/(self.size[i]/2.)+1)/2*self.grid[i]) for i in range(3)]
        return self.index[indexes[0], indexes[1], indexes[2]]

    def normal_from_pos(self, pos):
        #this returns surface normal (if exists) of a given photon in a given position

        indexes = [int((pos[i]/(self.size[i]/2.)+1)/2*self.grid[i]) for i in range(3)]
        return [self.normal[0][indexes[0], indexes[1], indexes[2]],
            self.normal[1][indexes[0], indexes[1], indexes[2]], 
            self.normal[2][indexes[0], indexes[1], indexes[2]]]
    
    def normal_and_index_from_pos(self, pos):
        #this returns surface normal (if exists) and the refractive index of a given photon in a given position

        indexes = [int((pos[i]/(self.size[i]/2.)+1)/2*self.grid[i]) for i in range(3)]
        return ([self.normal[0][indexes[0], indexes[1], indexes[2]],
            self.normal[1][indexes[0], indexes[1], indexes[2]], 
            self.normal[2][indexes[0], indexes[1], indexes[2]]], self.index[indexes[0], indexes[1], indexes[2]])

    def save_photon(self, sphoton, index, subindex):
        new_photon = photon([0, 0, 0], [0, 0, 1])
        new_photon.set_attr( sphoton.get_attr() )
        self.good_photons.append(new_photon)
        self.saved_photons[index][subindex].append(new_photon)

    def show_photons2D(self, plan='xy'):
        if 'xz' not in plan and 'xz' not in plan and 'xy' not in plan:
            raise Exception('Please pick either xy, xz or yz as plan.')
        
        
        def unpack_2d_photons(photon_list, sindex, label='Source'):
            print(f'Unpacking 2D {len(photon_list)} photons.')
            phx, phy, phz = list(), list(), list()
            nphx, nphy, nphz = list(), list(), list()
            intensity = list()

            for photon in photon_list:
                ipos = (photon.pos)
                nor = photon.normal
                phx.append(ipos[0]); phy.append(ipos[1]); phz.append(ipos[2])
                nphx.append(nor[0]); nphy.append(nor[1]); nphz.append(nor[2])
                intensity.append(photon.intensity)
            
            if plan=='xy': axes[sindex].hist2d(phx, phy, 31, weights = intensity, label='xy')
            if plan=='xz': axes[sindex].hist2d(phx, phz, 31, weights = intensity, label='xz')
            if plan=='yz': axes[sindex].hist2d(phy, phz, 31, weights = intensity, label='yz')
        

        xlen, ylen, zlen = len(self.saved_photons[0]), len(self.saved_photons[1]), len(self.saved_photons[2])
        if xlen*ylen!=0 and xlen*zlen!=0 and ylen*zlen!=0:
            if max(xlen, ylen, zlen)==0:
                raise Exception('Only a single index is supported by now. Please put your plans in same dimension.')

        maxcols = max(xlen, ylen, zlen)

        fig, axes = plt.subplots(nrows=1, ncols=maxcols, sharex=False, sharey=False)

        for index, planes in enumerate(self.saved_photons):
            [unpack_2d_photons(photons, subindex, '') for subindex, photons in enumerate(planes) if planes]
            
        plt.show()
        


    def show_elements(self, good_photons=False, mode='all'):
        
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        xrefl, yrefl, zrefl = numpy.where(self.index==-1.00)
        xrefr, yrefr, zrefr = numpy.where(self.index>1.00)

        def unpack_photons(photon_list, label='Source'):
            print(f'Unpacking {len(photon_list)} photons.')
            phx, phy, phz = list(), list(), list()
            nphx, nphy, nphz = list(), list(), list()

            for photon in photon_list:
                ipos = self.pos_to_grid(photon.pos)
                nor = photon.normal
                phx.append(ipos[0]); phy.append(ipos[1]); phz.append(ipos[2])
                nphx.append(nor[0]); nphy.append(nor[1]); nphz.append(nor[2])
            
            ax.scatter(phx, phz, phy, c='green', label=label)
            ax.quiver(phx, phz, phy, nphx, nphz, nphy)


        if good_photons:
            for index, planes in enumerate(self.saved_photons):
                [unpack_photons(photons, str(subindex)) for subindex, photons in enumerate(planes) if planes]
        else:
            unpack_photons(self.photons)
        

        if mode=='all':
            if xrefl.any(): ax.scatter(xrefl, zrefl, yrefl, c='red', label='Refractive')
            if xrefr.any(): ax.scatter(xrefr, zrefr, yrefr, c='blue', label='Reflective')
            for index, values in enumerate(self.plans):
                for val in values: #index is if x, y, z. val is the value
                    if index==0:
                        ax.plot_surface(0*self.x+self.pos_to_grid_1D(val, index), self.x, self.y, color='yellow', alpha=0.5)
                    elif index==1:
                        ax.plot_surface(self.x, self.y, 0*self.y+self.pos_to_grid_1D(val, index), color='yellow', alpha=0.5)
                    elif index==2:
                        ax.plot_surface(self.x, 0*self.x+self.pos_to_grid_1D(val, index), self.y, color='yellow', alpha=0.5)

        ax.set_xlabel('X')
        ax.set_ylabel('Z')
        ax.set_zlabel('Y')
        plt.legend()
        plt.show()

    def d2_source(self, r, z=0, normal=[0, 0, 1], na = 0.12, angles=1):
        x = numpy.arange(-r, r, self.res)
        if not x.size>0:
            x=[0]
        naper = numpy.linspace(-na, na, angles)
        for xpos in tqdm(x, desc='Source'):
            for ypos in numpy.arange(-numpy.sqrt(r**2-xpos**2), numpy.sqrt(r**2-xpos**2)+self.res, self.res):
                for nax in naper:
                    for nay in naper:
                        normal2 = numpy.add(normal, numpy.multiply([1, 1, 0], [nax, nay, 0]))
                        self.photons.append(photon([xpos, ypos, z], normal2))


    def convergent_lens(self, c=[0, 0, 3.0], e=0.0, radius=0.4, focus=2.0, n=1.43, dir = 0):
        #Put the index of refraction of glass in a given region. radius must be inferior than focus/2.

        R = focus/2
        ss = 2 #sub sampling
        x = numpy.arange(c[0]-radius, c[0]+radius+self.res, self.res/ss)
        y = numpy.arange(c[1]-radius, c[1]+radius+self.res, self.res/ss)
        if e>R or e<0:
            raise Exception('Lens too thick or negative number. Please reduce it')
        if dir:
            print('Plane-Convex Lens')
            z = numpy.arange(c[2]+e, c[2]+R+self.res, self.res/ss)
        else:
            print('Convex-Plane Lens')
            z = numpy.arange(c[2]-R, c[2]-e, self.res/ss)
        
        for xi, xpos in enumerate(x):
            for ypos in y:
                for zpos in z:
                    d = self.distance([xpos, ypos, zpos], c)
                    d1 = self.distance([xpos, ypos], [c[0], c[1]])
                    d2 = self.distance(zpos, c[2]+e) if dir else self.distance(zpos, c[2]-e) 
                    if d1**2<=radius**2:
                        if d**2<=(R+self.res)**2 or d2<=self.res:
                            ind = self.pos_to_grid([xpos, ypos, zpos])
                            self.assign_n(ind, n)
                            if d2<=self.res:
                                self.assign_normal(ind, [0, 0, -1])
                            else:
                                self.assign_normal(ind, [xpos-c[0], ypos-c[1], zpos-c[2]])


    def create_parabolic_section_element(self, c, n, th, wid, pp):
        ss=2
        x0, y0, z0 = c
        #th is thickness and it is related to y. ysym of 0.5 means a symmetric with respect to Z. Thickness is
        #this whole value in z. If ysym is 0, thickness is related to the semi parabola in negative y's. Same
        #applies to x with xsym / width [wid]
        ymax = y0 + 0.5*th
        ymin = y0 - 0.5*th

        xmax = x0 + 0.5*wid
        xmin = x0 - 0.5*wid
        
        zmax = z0
        zmin = z0 - (1/2*pp)*(max(abs(ymax-y0), abs(ymin-y0))**2+max(abs(xmax-x0), abs(xmin-x0))**2)
        x = numpy.arange(xmin, xmax+self.res, self.res/ss)
        y = numpy.arange(ymin, ymax+self.res, self.res/ss)
        z = numpy.arange(zmin, zmax+self.res, self.res/ss)
        for xpos in tqdm(x, desc='Parabolic'):
            for ypos in y:
                for zpos in z:
                    if (zpos-z0)>(-1/(2*pp))*((ypos-y0)**2+(xpos-x0)**2):
                        ind = self.pos_to_grid([xpos, ypos, zpos])
                        self.assign_n(ind, n)
                        self.assign_normal(ind, [(1/pp)*(xpos-x0), (1/pp)*(ypos-y0), 1])
        

    def create_sphere_element(self, c, r, n):
        
        def assign(pos):
            xpos, ypos, zpos = pos
            ind = self.pos_to_grid([xpos, ypos, zpos])
            self.assign_n(ind, n)
            self.assign_normal(ind, [xpos-c[0], ypos-c[1], zpos-c[2]])

        ss = 2 #sub sampling
        x = numpy.arange(c[0], c[0]+r+self.res, self.res/ss)
        y = numpy.arange(c[1], c[1]+r+self.res, self.res/ss)
        z = numpy.arange(c[2], c[2]+r+self.res, self.res/ss)
        for xpos in tqdm(x, desc='Sphere Section'):
            for ypos in y:
                for zpos in z:
                    d = self.distance([xpos, ypos, zpos], c)
                    if d**2<=(r+self.res)**2:
                        assign([xpos, ypos, zpos])
                        assign([xpos, -ypos+2*c[1], zpos])
                        assign([xpos, -ypos+2*c[1], -zpos+2*c[2]])
                        assign([xpos, ypos, -zpos+2*c[2]])
                        
                        assign([-xpos+2*c[0], ypos, zpos])
                        assign([-xpos+2*c[0], -ypos+2*c[1], zpos])
                        assign([-xpos+2*c[0], -ypos+2*c[1], -zpos+2*c[2]])
                        assign([-xpos+2*c[0], ypos, -zpos+2*c[2]])


    def create_rectangle_element(self, val, n, normal):
        xmin, xmax, ymin, ymax, zmin, zmax = val
        ss = 2 #sub sampling
        x = numpy.arange(xmin, xmax+self.res, self.res/ss)
        y = numpy.arange(ymin, ymax+self.res, self.res/ss)
        z = numpy.arange(zmin, zmax+self.res, self.res/ss)
        for xpos in tqdm(x, desc='Rectangle'):
            for ypos in y:
                for zpos in z:
                    ind = self.pos_to_grid([xpos, ypos, zpos])
                    self.assign_n(ind, n)
                    self.assign_normal(ind, normal)


    def distance(self, vec1, vec2):
        return numpy.sum(numpy.power(numpy.subtract(vec1, vec2), 2))**0.5

    def check_position(self, photon):
        pos = photon.pos
        if -self.size[0]/2.<pos[0]<self.size[0]/2. and -self.size[1]/2.<pos[1]<self.size[1]/2. and -self.size[2]/2.<pos[2]<self.size[2]/2.:
           return True
        else:
            return False

    def check_plane(self, photon, plane):
        xp, yp, zp = plane
        pos = photon.pos
        #if yp-self.res/2.<pos[1]<yp+self.res/2.:
        if zp-self.res/2.<pos[2]<zp+self.res/2.:
            return True
        else:
            return False

    def check_analysis_plan(self, photon):
        pos = photon.pos
        for index, values in enumerate(self.plans):
            for subindex, val in enumerate(values):
                if abs(pos[index]-val)<=self.res/2.:
                    self.save_photon(photon, index, subindex)
                    return True

    def create_analysis_plan(self, **kargs):
        if 'z' in kargs:
            self.plans[2].append(kargs['z'])
            self.saved_photons[2].append(list())
        if 'y' in kargs:
            self.plans[1].append(kargs['y'])
            self.saved_photons[1].append(list())
        if 'x' in kargs:
            self.plans[0].append(kargs['x'])
            self.saved_photons[0].append(list())


    def run(self):
        i = 0
        while self.photons:
            i = i + 1
            #print(f'Interaction {i}. Number of photons in the simu is {len(self.photons)}')
            for photon in tqdm(self.photons, desc='Running', ascii=True, ncols=100):
                
                photon.update(self.normal_and_index_from_pos(photon.pos))
                photon.move(self.res)
                
                if not self.check_position(photon):
                    self.photons.remove(photon)
                elif self.check_analysis_plan(photon):
                    pass
                    #self.save_photon(photon, ind, subind)
                    #self.photons.remove(photon)

mirror_focus = 0.3
yvertex = 0.5+0.3
thickness = 1.2

a = Simu(5, 5, 5, 0.03)
a.d2_source(0.0, -2.5, [0, 0, 1], 0.32, 21)

a.create_sphere_element([0.0, 0.0, 0.5], 1.0, 1.43) #from 2 to 3.
a.create_rectangle_element([-1.0, 1.0, -1.0, 1.0, 0.5, 2.0], 1.0, [0, 0, 1]) 

#a.create_parabolic_section_element([0.0, yvertex, 0.0], -1.0, 2*thickness, 3.0, 1.0) 
#a.create_rectangle_element([-2.45, 2.45, yvertex-mirror_focus, 2.0, -2.0, 0.0], 1.0, [0, 0, 1]) 

for val in numpy.linspace(0, 2.3, 3):
    a.create_analysis_plan(z=val)

a.show_elements(False, 'all')
a.run()
#a.show_elements(True, 'photons')
a.show_photons2D('xy')
