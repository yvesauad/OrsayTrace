import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider
import time
from tqdm import tqdm


class photon_list():
    def __init__(self, normal, value):
        self.normal = normal / numpy.linalg.norm(normal)
        self.value = value
        self.photons = numpy.asarray([])

    def add_photon(self, sphoton):
        new_photon = photon([0, 0, 0], [0, 0, 1])
        new_photon.set_attr( sphoton.get_attr() )
        self.photons = numpy.append(self.photons, new_photon)

    def distance_point_to_plane(self, pos):
        return numpy.abs(numpy.dot(self.normal, pos)-self.value)



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
        
        if sur_normal==[0, 0, 0]:
            print('***WARNING***: No surface normal.')
            return False

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
        #print(self.refraction_count, inc, sur_normal, self.n)
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
        self.ss = 2. #sub sampling for creating elements
        self.grid = (int(x/res), int(y/res), int(z/res))
        self.index = numpy.ones(self.grid)
        self.normal = numpy.asarray([numpy.zeros(self.grid), numpy.zeros(self.grid), numpy.zeros(self.grid)])
        self.photons = list()
        self.photon_lists=numpy.asarray([])
        self.x = numpy.linspace(-self.size[0]/2., self.size[0]/2., self.grid[0])
        self.y = numpy.linspace(-self.size[1]/2., self.size[1]/2., self.grid[1])
        self.x, self.y = numpy.meshgrid(self.x, self.y)

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

    def grid_to_pos(self, grid):
        pos = [ (2.0*grid[i]/self.grid[i]-1) * self.size[i]/2.0 for i in range(3)]
        return pos

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


    def show_photons2D(self, plan='xy'):
        if 'xz' not in plan and 'xz' not in plan and 'xy' not in plan:
            raise Exception('Please pick either xy, xz or yz as plan.')
        
        def unpack_2d_photons(photon_list, sindex):
            print(f'Unpacking 2D {len(photon_list.photons)} photons.')
            
            phx, phy, phz = list(), list(), list()
            nphx, nphy, nphz = list(), list(), list()
            intensity = list()

            for photon in photon_list.photons:
                ipos = (photon.pos)
                nor = photon.normal
                phx.append(ipos[0]); phy.append(ipos[1]); phz.append(ipos[2])
                nphx.append(nor[0]); nphy.append(nor[1]); nphz.append(nor[2])
                intensity.append(photon.intensity)
            
            if plan=='xy': axes[sindex].hist2d(phx, phy, 21, weights = intensity)
            if plan=='xz': axes[sindex].hist2d(phx, phz, 21, weights = intensity)
            if plan=='yz': axes[sindex].hist2d(phy, phz, 21, weights = intensity)

            #axes[sindex].set_title(str(photon_list.normal)+'_'+str(photon_list.value))
            axes[sindex].set_title(str(photon_list.value))

        fig, axes = plt.subplots(nrows=1, ncols=len(self.photon_lists), sharex=False, sharey=False)
            
        for index, photon_list in enumerate(self.photon_lists):
            unpack_2d_photons(photon_list, index)

        plt.show()
        
    def show_elements(self, saved_photons=False, mode='all'):
        
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        xrefl, yrefl, zrefl = self.grid_to_pos(numpy.where(self.index==-1.00))
        xrefr, yrefr, zrefr = self.grid_to_pos(numpy.where(self.index>1.00))

        def unpack_photons(photon_list):
            if hasattr(photon_list, 'photons'):
                title = 'PN: '+str(photon_list.normal)+'. d= '+format(photon_list.value, '.2f')
                photon_list = photon_list.photons
            else:
                title = 'Photons'
            
            print(f'Unpacking {len(photon_list)} photons.')
            
            phx, phy, phz = list(), list(), list()
            nphx, nphy, nphz = list(), list(), list()

            for photon in photon_list:
                ipos = photon.pos
                nor = photon.normal
                phx.append(ipos[0]); phy.append(ipos[1]); phz.append(ipos[2])
                nphx.append(nor[0]); nphy.append(nor[1]); nphz.append(nor[2])
            
            ax.scatter(phx, phz, phy, c='green', label=title)
            ax.quiver(phx, phz, phy, nphx, nphz, nphy, length=0.5)


        if saved_photons:
            for index, photon_list in enumerate(self.photon_lists):
                unpack_photons(photon_list)
        else:
            unpack_photons(self.photons)
        

        if mode=='all':
            if xrefl.any(): ax.scatter(xrefl, zrefl, yrefl, c='red', label='Refractive')
            if xrefr.any(): ax.scatter(xrefr, zrefr, yrefr, c='blue', label='Reflective')
            for index, photon_list in enumerate(self.photon_lists):
                ax.plot_surface(self.x, (photon_list.value-photon_list.normal[0]*self.x-photon_list.normal[1]*self.y)/photon_list.normal[2], self.y, color='yellow', alpha=0.2)

        ax.set_xlabel('X')
        ax.set_xlim(-self.size[0]/2.0, self.size[0]/2.0)
        
        ax.set_ylabel('Z')
        ax.set_ylim(-self.size[2]/2.0, self.size[2]/2.0)
       
        ax.set_zlabel('Y') 
        ax.set_zlim(-self.size[1]/2.0, self.size[1]/2.0)
        
        plt.legend()
        plt.show()

    def circular_source(self, r, c, normal=[0, 0, 1]):
        xc, yc, zc = c
        x = numpy.arange(xc-r, xc, self.res)
        y = numpy.arange(yc-r, yc, self.res)
        if not x.size>0:
            x=[0]
            y=[0]
        for xpos in tqdm(x, desc='Source'):
            for ypos in y:
                if xpos**2+ypos**2<=r**2:
                    self.photons.append(photon([xpos, ypos, zc], normal))
                    self.photons.append(photon([-xpos, ypos, zc], normal))
                    self.photons.append(photon([xpos, -ypos, zc], normal))
                    self.photons.append(photon([-xpos, -ypos, zc], normal))


    def d2_source(self, r, c=[0, 0, 0], normal=[0, 0, 1], na = 0.12, angles=1):
        xc, yc, zc = c
        x = numpy.arange(xc-r, xc+r, self.res)
        if not x.size>0:
            x=[0]
        naper = numpy.linspace(-na, na, angles)
        for xpos in tqdm(x, desc='Source'):
            for ypos in numpy.arange(-numpy.sqrt(r**2-(xpos-xc)**2), numpy.sqrt(r**2-(xpos-xc)**2)+self.res, self.res):
                for nax in naper:
                    for nay in naper:
                        normal2 = numpy.add(normal, numpy.multiply([1, 1, 0], [nax, nay, 0]))
                        self.photons.append(photon([xpos, ypos, zc], normal2))

    def create_parabolic_section_element(self, c, n, th, wid, pp):
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
        x = numpy.arange(xmin, xmax+self.res, self.res/self.ss)
        y = numpy.arange(ymin, ymax+self.res, self.res/self.ss)
        z = numpy.arange(zmin, zmax+self.res, self.res/self.ss)
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

        x = numpy.arange(c[0], c[0]+r+self.res, self.res/self.ss)
        y = numpy.arange(c[1], c[1]+r+self.res, self.res/self.ss)
        z = numpy.arange(c[2], c[2]+r+self.res, self.res/self.ss)
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


    def create_rectangle_element(self, val, n, normal, inclusive=True):
        if inclusive:
            fac = 1
        else:
            fac = 0
        xmin, xmax, ymin, ymax, zmin, zmax = val
        x = numpy.arange(xmin, xmax+fac*self.res, self.res/self.ss)
        y = numpy.arange(ymin, ymax+fac*self.res, self.res/self.ss)
        z = numpy.arange(zmin, zmax+fac*self.res, self.res/self.ss)
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

    def check_analysis_plan(self, photon):
        pos = photon.pos
        for index, planes in enumerate(self.photon_lists):
            if planes.distance_point_to_plane(pos)<=self.res/2.0:
                planes.add_photon(photon)

    def create_analysis_plan(self, normal, value):
        self.photon_lists = numpy.append(self.photon_lists, photon_list(normal, value))


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

#a = Simu(5, 5, 5, 0.03)
#a.d2_source(0.0, -2.5, [0, 0, 1], 0.32, 21)

#a.create_sphere_element([0.0, 0.0, 0.5], 1.0, 1.43) #from 2 to 3.
#a.create_rectangle_element([-1.0, 1.0, -1.0, 1.0, 0.5, 2.0], 1.0, [0, 0, 1]) 

#a.create_parabolic_section_element([0.0, yvertex, 0.0], -1.0, 2*thickness, 3.0, 1.0) 
#a.create_rectangle_element([-2.45, 2.45, yvertex-mirror_focus, 2.0, -2.0, 0.0], 1.0, [0, 0, 1]) 

#for val in numpy.linspace(0, 2.3, 3):
#    a.create_analysis_plan(z=val)

#a.show_elements(False, 'all')
#a.run()
#a.show_elements(True, 'photons')
#a.show_photons2D('xy')
