import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider
import time
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import ThreadPoolExecutor


class photon_list():
    def __init__(self, normal, value):
        self.normal = normal / numpy.linalg.norm(normal)
        self.value = value
        self.photons = numpy.asarray([])
        self.condition_dict = dict()

    def append_photon(self, photon):
        if self.condition_dict:
            for cond, val in self.condition_dict.items():
                for key, pval in photon.__dict__.items():
                    if cond==key:
                        if type(val) is tuple:
                            if pval>=val[0] and pval<=val[1]:
                                self.photons = numpy.append(self.photons, photon)
                        else:
                            if pval>=val:
                                self.photons = numpy.append(self.photons, photon)
        else:
            self.photons = numpy.append(self.photons, photon)


    def add_photon(self, sphoton):
        new_photon = photon([0, 0, 0], [0, 0, 1])
        new_photon.set_attr( sphoton.get_attr() )
        self.append_photon(new_photon)

    def distance_point_to_plane(self, pos):
        return numpy.abs(numpy.dot(self.normal, pos)-self.value)

    def add_symmetric_xphoton(self, sphoton):
        new_photon = photon([0, 0, 0], [0, 0, 1])
        new_photon.set_attr( sphoton.get_attr() )
        new_photon.pos = -new_photon.pos[0], new_photon.pos[1], new_photon.pos[2]
        new_photon.normal = -new_photon.normal[0], new_photon.normal[1], new_photon.normal[2]
        self.append_photon(new_photon)
    
    def add_symmetric_yphoton(self, sphoton):
        new_photon = photon([0, 0, 0], [0, 0, 1])
        new_photon.set_attr( sphoton.get_attr() )
        new_photon.pos = new_photon.pos[0], -new_photon.pos[1], new_photon.pos[2]
        new_photon.normal = new_photon.normal[0], -new_photon.normal[1], new_photon.normal[2]
        self.append_photon(new_photon)

    def avg_divergence(self, normal_ref):
        normal_ref = normal_ref / numpy.linalg.norm(normal_ref)
        vals = numpy.average([numpy.dot(photon.normal, normal_ref)**2 for photon in self.photons])
        return vals

    def get_positions(self):
        vals = numpy.asarray([photon.pos for photon in self.photons])
        return vals

    def get_relative_centroid_positions(self):
        #Relative position to centroid for each photon
        vals = self.get_positions() - self.avg_position()
        return vals

    def get_relative_centroid_distances(self):
        #Relative distance to centroid for each photon
        vals = numpy.sqrt(numpy.sum(numpy.power(self.get_relative_centroid_positions(), 2), axis=1))
        return vals

    def get_intensities(self):
        vals = numpy.asarray([photon.intensity for photon in self.photons])
        return vals

    def avg_position(self):
        vals = numpy.asarray([numpy.average([photon.pos[axis] for photon in self.photons]) for axis in range(3)])
        return vals
    
    def max_position(self):
        vals = [numpy.amax([photon.pos[axis] for photon in self.photons]) for axis in range(3)]
        return vals
    
    def min_position(self):
        vals = [numpy.amin([photon.pos[axis] for photon in self.photons]) for axis in range(3)]
        return vals
    
    def std_position(self):
        vals = [numpy.std([photon.pos[axis] for photon in self.photons]) for axis in range(3)]
        return vals

    def get_weighted_inverse(self):
        vals = numpy.divide(self.get_intensities(), self.get_relative_centroid_distances())
        return vals

    def get_average_weighted_inverse(self):
        vals = numpy.average(self.get_weighted_inverse())
        return vals

    def avg_distance_axis_z(self, c=[0, 0]):
        xc, yc = c[0], c[1]
        vals = numpy.average([numpy.sqrt((photon.pos[0]-xc)**2+(photon.pos[1]-yc)**2) for photon in self.photons])
        return vals
    
    def avg_distance_axis_y(self, c=[0, 0]):
        #Average relative distance to given point in a given axis for each photon
        xc, zc = c[0], c[1]
        vals = numpy.average([numpy.sqrt((photon.pos[0]-xc)**2+(photon.pos[2]-zc)**2) for photon in self.photons])
        return vals
    
    def get_weighted_inverse_axis_y(self, c=[0, 0]):
        assert (self.normal == [0, 1, 0]).all()
        pos = self.get_positions() - numpy.asarray([c[0], self.value, c[1]])
        dist = numpy.sqrt(numpy.sum(numpy.power(pos, 2), axis=1))
        wi = numpy.divide(self.get_intensities(), dist)
        avg_wi = numpy.average(wi)
        return avg_wi


class photon():
    def __init__(self, pos, normal, intensity=1.):
        self.pos = pos
        self.normal = normal / numpy.linalg.norm(normal)
        self.intensity = 1.
        self.n = 1.00
        self.last_surface = [0, 0, 0]
        self.refraction_count = 0
        self.reflection_count = 0

    def set_attr(self, values):
        self.pos, self.normal, self.intensity, self.n, self.last_surface, self.refraction_count, self.reflection_count = values
    
    def get_attr(self):
        return self.pos, self.normal, self.intensity, self.n, self.last_surface, self.refraction_count, self.reflection_count


    def reflection(self):
        '''
        Reflection photon object using last_surface attribute.

        Parameters
        ----------
        None
        

        Returns
        -------
        boolean
            True if reflection is sucessfull and false if you attempted a reflection in a grid which there was no surface normal vector.
        '''
        inc = self.normal
        inc = inc / numpy.linalg.norm(inc)
        sur_normal = self.last_surface
        
        if sur_normal==[0, 0, 0]:
            print('***WARNING***: No surface normal.')
            return False
        
        sur_normal = sur_normal / numpy.linalg.norm(sur_normal)
        cos_ang1 = -numpy.dot(inc, sur_normal)
        if cos_ang1<0:
            sur_normal = -sur_normal
            cos_ang1 = -numpy.dot(inc, sur_normal)
 
        refl = inc + 2*cos_ang1*sur_normal
        refl = refl / numpy.linalg.norm(refl)
        self.reflection_count+=1

        self.normal = refl
        #print(self.reflection_count, inc, sur_normal, refl, self.pos)
        return True
    
    def refraction(self, n2):
        '''
        Reflection photon object using last_surface attribute and a given index of refraction.

        Parameters
        ----------
        n2: float
            The index of refraction.
        

        Returns
        -------
        boolean
            True if refraction is sucessfull and false if you attempted a refraction in a grid which there was no surface normal vector
        '''
        
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
        '''
        Increment photon position

        Parameters
        ----------
        move: float
            Increment photon position of value
        

        Returns
        -------
        boolean
            True if is sucessfull.
        '''
        self.pos = self.pos + self.normal*value
        return True

    def update(self, value):
        sur_normal, index = value
        
        if sur_normal != [0, 0, 0]:
            self.last_surface = sur_normal
        
        if index!=self.n and index>0:
            self.refraction(index)
            return True
        elif index==-1:
            self.reflection()
            return True

        return False


class Simu:
    def __init__(self, x, y, z, res):
        self.size = [x, y, z]
        self.res = res
        self.ss = 2.01 #sub sampling for creating elements
        self.grid = (int(x/res), int(y/res), int(z/res))
        self.index = numpy.ones(self.grid)
        self.normal = numpy.asarray([numpy.zeros(self.grid), numpy.zeros(self.grid), numpy.zeros(self.grid)])
        self.photons = list()
        self.photon_lists=numpy.asarray([])
        self.x = numpy.linspace(-self.size[0]/2., self.size[0]/2., self.grid[0])
        self.y = numpy.linspace(-self.size[1]/2., self.size[1]/2., self.grid[1])
        self.x, self.y = numpy.meshgrid(self.x, self.y)

    def assign_normal(self, index, normal):
        if (normal!=numpy.zeros(3)).all(): 
            normal = normal / numpy.linalg.norm(normal)
        for i in range(3):
            self.normal[i][index[0], index[1], index[2]] = normal[i]
    
    def assign_block_normal(self, min_index, max_index, normal):
        if (normal!=numpy.zeros(3)).all():
            normal = normal / numpy.linalg.norm(normal)
        for i in range(3):
            self.normal[i][min_index[0]:max_index[0], min_index[1]:max_index[1], min_index[2]:max_index[2]] = normal[i]

    def assign_n(self, index, value):
        self.index[index[0], index[1], index[2]] = value
    
    def assign_block_n(self, min_index, max_index, value):
        self.index[min_index[0]:max_index[0], min_index[1]:max_index[1], min_index[2]:max_index[2]] = value

    def pos_to_grid_1D(self, value, index):
        return round((value/(self.size[index]/2.)+1)/2*self.grid[index])
        

    def pos_to_grid(self, pos):
        #this returns relative current photon position from [0:self.grid]
        indexes = [round((pos[i]/(self.size[i]/2.)+1)/2*self.grid[i]) for i in range(3)]
        return indexes

    def grid_to_pos(self, grid):
        pos = [ (2.0*grid[i]/self.grid[i]-1) * self.size[i]/2.0 for i in range(3)]
        return pos

    def index_from_pos(self, pos):
        #this returns current index of refraction of a given photon in a given position

        indexes = [round((pos[i]/(self.size[i]/2.)+1)/2*self.grid[i]) for i in range(3)]
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


    def show_photons2D(self, my_photon_lists, plan='xy'):
        if 'xz' not in plan and 'xz' not in plan and 'xy' not in plan:
            raise Exception('Please pick either xy, xz or yz as plan.')
        
        def unpack_2d_photons(photon_list, sindex):
            
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

            axes[sindex].set_title(str(photon_list.value))

        fig, axes = plt.subplots(nrows=1, ncols=len(my_photon_lists), sharex=False, sharey=False)
            
        for index, photon_list in enumerate(my_photon_lists):
            assert hasattr(photon_list, 'photons')
            unpack_2d_photons(photon_list, index)

        plt.show()

    def show_created_elements(self, mode):
        
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        xrefl, yrefl, zrefl = self.grid_to_pos(numpy.where(self.index==-1.00))
        xrefr, yrefr, zrefr = self.grid_to_pos(numpy.where(self.index>1.00))

        def unpack_photons(photon_list):
            print(f'Unpacking {len(photon_list)} photons.')
            
            phx, phy, phz = list(), list(), list()
            nphx, nphy, nphz = list(), list(), list()

            for photon in photon_list:
                ipos = photon.pos
                nor = photon.normal
                phx.append(ipos[0]); phy.append(ipos[1]); phz.append(ipos[2])
                nphx.append(nor[0]); nphy.append(nor[1]); nphz.append(nor[2])
            
            ax.scatter(phx, phz, phy, c='green', label='Photons')
            ax.quiver(phx, phz, phy, nphx, nphz, nphy, length=0.5)

        unpack_photons(self.photons)

        if 'all' in mode:
            if xrefl.any(): ax.scatter(xrefl, zrefl, yrefl, c='red', label='Reflective')
            if xrefr.any(): ax.scatter(xrefr, zrefr, yrefr, c='blue', label='Refractive')
            if '-noplan' not in mode:
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

    def show_elements(self, my_photon_lists, mode='all-noplan'):

        fig = plt.figure()
        ax = plt.axes(projection='3d')
        xrefl, yrefl, zrefl = self.grid_to_pos(numpy.where(self.index==-1.00))
        xrefr, yrefr, zrefr = self.grid_to_pos(numpy.where(self.index>1.00))

        def unpack_photons(photon_list):
            if hasattr(photon_list, 'photons'):
                title = 'PN: '+str(photon_list.normal)+'. d= '+format(photon_list.value, '.2f')
                photon_list = photon_list.photons
                title = ''
            
            #print(f'Unpacking {len(photon_list)} photons.')
            
            phx, phy, phz = list(), list(), list()
            nphx, nphy, nphz = list(), list(), list()

            for photon in photon_list:
                ipos = photon.pos
                nor = photon.normal
                phx.append(ipos[0]); phy.append(ipos[1]); phz.append(ipos[2])
                nphx.append(nor[0]); nphy.append(nor[1]); nphz.append(nor[2])
            
            ax.scatter(phx, phz, phy, c='green', label=title)
            ax.quiver(phx, phz, phy, nphx, nphz, nphy, length=0.5)


        for index, photon_list in enumerate(my_photon_lists):
            assert hasattr(photon_list, 'photons')
            unpack_photons(photon_list)
        

        if 'all' in mode:
            if xrefl.any(): ax.scatter(xrefl, zrefl, yrefl, c='red', label='Reflective', alpha=0.1)
            if xrefr.any(): ax.scatter(xrefr, zrefr, yrefr, c='blue', label='Refractive', alpha=0.1)
            if '-noplan' not in mode:
                for index, photon_list in enumerate(self.photon_lists):
                    if photon_list.normal[2]>0:
                        ax.plot_surface(self.x, (photon_list.value-photon_list.normal[0]*self.x-photon_list.normal[1]*self.y)/photon_list.normal[2], self.y, color='yellow', alpha=0.2)
                    elif photon_list.normal[1]>0:
                        ax.plot_surface(self.x, self.y, (photon_list.value-photon_list.normal[0]*self.x-photon_list.normal[2]*self.y)/photon_list.normal[1], color='yellow', alpha=0.2)
                    elif photon_list.normal[0]>0:
                        ax.plot_surface((photon_list.value-photon_list.normal[2]*self.x-photon_list.normal[1]*self.y)/photon_list.normal[0], self.x, self.y, color='yellow', alpha=0.2)

        ax.set_xlabel('X')
        ax.set_xlim(-self.size[0]/2.0, self.size[0]/2.0)
        
        ax.set_ylabel('Z')
        ax.set_ylim(-self.size[2]/2.0, self.size[2]/2.0)
       
        ax.set_zlabel('Y') 
        ax.set_zlim(-self.size[1]/2.0, self.size[1]/2.0)
        
        plt.legend()
        plt.show()

    def point_source(self, r, c, normal=[0, 0, 1]):
        xc, yc, zc = c
        x = numpy.arange(xc-r, xc+self.res, self.res)
        y = numpy.arange(yc-r, yc+self.res, self.res)
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
        x = numpy.arange(xc-r, xc+self.res, self.res)
        y = numpy.arange(yc-r, yc+self.res, self.res)
        if not x.size>0:
            x=[xc]
            y=[yc]
        naper = numpy.linspace(-na, na, angles)
        for xpos in tqdm(x, desc='Source'):
            for ypos in y:
                if (xpos-xc)**2+(ypos-yc)**2<=r**2:
                    for nax in naper:
                        for nay in naper:
                            normal2 = numpy.add(normal, numpy.multiply([1, 1, 0], [nax, nay, 0]))
                            self.photons.append(photon([xpos, ypos, zc], normal2))
                            self.photons.append(photon([-xpos+2*xc, ypos, zc], normal2))
                            self.photons.append(photon([xpos, -ypos+2*yc, zc], normal2))
                            self.photons.append(photon([-xpos+2*xc, -ypos+2*yc, zc], normal2))
    
    def rotate(self, ang, axis, origin, ROI = None):
        axis = axis / numpy.linalg.norm(axis)
        ux, uy, uz = axis

        xc, yc, zc = origin
        c = numpy.cos(ang)
        s = numpy.sin(ang)

        pos_origin = numpy.asarray([
            [xc],
            [yc],
            [zc],
            ])

        def rotate_pos(p):
            px, py, pz = p
            
            m = numpy.asarray([
                [c+ux**2*(1-c), ux*uy*(1-c)-uz*s, ux*uz*(1-c)+uy*s],
                [uy*ux*(1-c)+uz*s, c+uy**2*(1-c), uy*uz*(1-c)-ux*s],
                [uz*ux*(1-c)-uy*s, uz*uy*(1-c)+ux*s, c+uz**2*(1-c)],
                ])
            pos = numpy.asarray([
                    [px],
                    [py],
                    [pz],
                    ])
            return numpy.asarray((numpy.matmul(m, pos - pos_origin) + pos_origin).T)
        
        def rotate_normal(nor):
            nx, ny, nz = nor
            
            m = numpy.asarray([
                [c+ux**2*(1-c), ux*uy*(1-c)-uz*s, ux*uz*(1-c)+uy*s],
                [uy*ux*(1-c)+uz*s, c+uy**2*(1-c), uy*uz*(1-c)-ux*s],
                [uz*ux*(1-c)-uy*s, uz*uy*(1-c)+ux*s, c+uz**2*(1-c)],
                ])
            normal = numpy.asarray([
                    [nx],
                    [ny],
                    [nz],
                    ])
            return (numpy.matmul(m, normal).T)


        if ROI is not None:
            xmin, xmax, ymin, ymax, zmin, zmax = ROI
        else:
            xmin, xmax, ymin, ymax, zmin, zmax = -self.size[0]/2, self.size[0]/2.0-self.res, -self.size[0]/2.0, self.size[0]/2.0-self.res, -self.size[0]/2.0, self.size[0]/2.0-self.res
        

        min_index = self.pos_to_grid([xmin, ymin, zmin])
        max_index = self.pos_to_grid([xmax, ymax, zmax])


        points_to_rotate = None
        index_to_rotate = None
        normal_to_rotate = None

        x = numpy.arange(xmin, xmax+self.res, self.res/self.ss)
        y = numpy.arange(ymin, ymax+self.res, self.res/self.ss)
        z = numpy.arange(zmin, zmax+self.res, self.res/self.ss)
        
        for xpos in tqdm(x, desc='Rotate: '):
            for ypos in y:
                for zpos in z:
                    if self.index_from_pos([xpos, ypos, zpos]) != 1:
                        if points_to_rotate is None:
                            points_to_rotate = numpy.asarray(rotate_pos([xpos, ypos, zpos]))
                            index_to_rotate = numpy.asarray(self.index_from_pos([xpos, ypos, zpos]))
                            normal_to_rotate = numpy.asarray(rotate_normal(self.normal_from_pos([xpos, ypos, zpos])))
                        else:
                            points_to_rotate = numpy.append(points_to_rotate, rotate_pos([xpos, ypos, zpos]), axis=0)
                            index_to_rotate = numpy.append(index_to_rotate, self.index_from_pos([xpos, ypos, zpos]))
                            normal_to_rotate = numpy.append(normal_to_rotate, rotate_normal(self.normal_from_pos([xpos, ypos, zpos])), axis=0)
                                
                            index = self.pos_to_grid([xpos, ypos, zpos])
                            self.assign_n(index, 1.0)
                            self.assign_normal(index, [0, 0, 0])

        if points_to_rotate is None:
            raise Exception('There is nothing to rotate. Please check your elements in simulation cell or ROI.')

        for j, index_point in enumerate(points_to_rotate):
            xpos, ypos, zpos = index_point
            ind_refr = index_to_rotate[j]
            normal = normal_to_rotate[j]
            #if self.index_from_pos([xpos, ypos, zpos]) == ind_refr: 
            #    print('**WARNING***: Duplicate in rotation. Losing those points: ', self.pos_to_grid([xpos, ypos, zpos]), [xpos, ypos, zpos], ind_refr)
            grid_pos = self.pos_to_grid([xpos, ypos, zpos])
            self.assign_n(grid_pos, ind_refr)
            self.assign_normal(grid_pos, normal)

    def create_parabolic_section_element(self, c, n, th, wid, pp):
        
        def assign(pos):
            xpos, ypos, zpos = pos
            ind = self.pos_to_grid([xpos, ypos, zpos])
            self.assign_n(ind, n)
            self.assign_normal(ind, [(1/pp)*(xpos-x0), (1/pp)*(ypos-y0), 1])
        
        x0, y0, z0 = c
        ymax = y0 + 0.5*th; ymin = y0 - 0.5*th
        xmax = x0 + 0.5*wid; xmin = x0 - 0.5*wid
        zmax = z0; zmin = z0 - (1/2*pp)*(max(abs(ymax-y0), abs(ymin-y0))**2+max(abs(xmax-x0), abs(xmin-x0))**2)
        
        x = numpy.arange(x0, xmax, self.res/self.ss)
        y = numpy.arange(y0, ymax, self.res/self.ss)
        z = numpy.arange(zmin, z0, self.res/self.ss)
        for xpos in tqdm(x, desc='Parabolic'):
            for ypos in y:
                for zpos in z:
                    if abs(zpos-z0-(-1/(2*pp))*((ypos-y0)**2+(xpos-x0)**2))<=self.res:
                        assign([xpos, ypos, zpos])
                        assign([-xpos+2*x0, ypos, zpos])
                        assign([xpos, -ypos+2*y0, zpos])
                        assign([-xpos+2*x0, -ypos+2*y0, zpos])

        
    def create_sphere_element(self, c, r, n):
        
        def assign(pos):
            xpos, ypos, zpos = pos


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


    def create_rectangle_element(self, val, n, normal):
        
        xmin, xmax, ymin, ymax, zmin, zmax = val
        xc, yc, zc = (xmin+xmax)/2., (ymin+ymax)/2., (zmax+zmin)/2.
        xl, yl, zl = (xmax-xmin)/2., (ymax-ymin)/2., (zmax-zmin)/2.
        assert xl>=0 and yl>=0 and zl>=0

        x = [x for x in numpy.arange(xc-xl, xc+self.res, self.res/self.ss) if abs(x)<self.size[0]/2. and abs(-x+2*xc)<self.size[0]/2.]
        y = [y for y in numpy.arange(yc-yl, yc+self.res, self.res/self.ss) if abs(y)<self.size[1]/2. and abs(-y+2*yc)<self.size[1]/2.]
        z = [z for z in numpy.arange(zc-zl, zc+self.res, self.res/self.ss) if abs(z)<self.size[2]/2. and abs(-z+2*zc)<self.size[2]/2.]
        
        min_index = self.pos_to_grid([min(x), min(y), min(z)])
        max_index = self.pos_to_grid([-min(x)+2*xc+1, -min(y)+2*yc+1, -min(z)+2*zc+1])
        
        self.assign_block_n(min_index, max_index, n)
        self.assign_block_normal(min_index, max_index, normal)

    def create_cylinder_element(self, center, radius, length, n, normal):
        
        xc, yc, zc = center
        
        def assign(pos):
            xpos, ypos, zpos = pos
            ind = self.pos_to_grid([xpos, ypos, zpos])
            self.assign_n(ind, n)
            self.assign_normal(ind, normal)
        
        x = numpy.arange(xc-radius, xc+self.res, self.res/self.ss)
        y = numpy.arange(yc-radius, yc+self.res, self.res/self.ss)
        z = numpy.arange(zc - length/2., zc+self.res, self.res/self.ss)
        for xpos in tqdm(x, desc='Cylinder'):
            for ypos in y:
                if (xpos-xc)**2+(ypos-yc)**2<=radius**2:
                    for zpos in z:
                        assign([xpos, ypos, zpos])
                        assign([xpos, -ypos+2*yc, zpos])
                        assign([xpos, -ypos+2*yc, -zpos+2*zc])
                        assign([xpos, ypos, -zpos+2*zc])
                        
                        assign([-xpos+2*xc, ypos, zpos])
                        assign([-xpos+2*xc, -ypos+2*yc, zpos])
                        assign([-xpos+2*xc, -ypos+2*yc, -zpos+2*zc])
                        assign([-xpos+2*xc, ypos, -zpos+2*zc])




    def create_thin_lens(self, cplane, focus, aperture, index_refr, lens_type='plane-convex'):
        xc, yc, zp = cplane
        r = focus/2.0
        assert focus/2.>=aperture/2.
        if lens_type=='plane-convex': zc = zp - numpy.sqrt((focus/2.0)**2 - (aperture/2.)**2) #plane_convex
        if lens_type=='convex-plane': zc = zp + numpy.sqrt((focus/2.0)**2 - (aperture/2.)**2) #convex-plane

        self.create_sphere_element([xc, yc, zc], focus/2., index_refr)
        
        if lens_type=='plane-convex':
            self.create_rectangle_element([xc-r-self.res, xc+r, yc-r-self.res, yc+r, zc-r-self.res, zp], 1.0, [0, 0, 0])

        if lens_type=='convex-plane':
            self.create_rectangle_element([xc-r-self.res, xc+r, yc-r-self.res, yc+r, zp-self.res, zc+r], 1.0, [0, 0, 0])
        
        self.create_cylinder_element([xc, yc, zp], aperture/2., 0.0, index_refr, [0, 0, 1])

    def distance(self, vec1, vec2):
        return numpy.sum(numpy.power(numpy.subtract(vec1, vec2), 2))**0.5

    def is_photon_in_cell(self, photon):
        pos = photon.pos
        if -self.size[0]/2.<pos[0]<self.size[0]/2. and -self.size[1]/2.<pos[1]<self.size[1]/2. and -self.size[2]/2.<pos[2]<self.size[2]/2.:
           return True
        else:
            return False

    def check_analysis_plan(self, photon, rindex):
        pos = photon.pos
        for index, planes in enumerate(self.split_photon_lists[rindex]):
            if planes.distance_point_to_plane(pos)<=self.res/2.0:
                planes.add_photon(photon)

    def create_analysis_plan(self, normal, value, **kargs):
        plist = photon_list(normal, value)
        plist.condition_dict = kargs
        self.photon_lists = numpy.append(self.photon_lists, plist)
    
    def run_photon(self, photon_list, rindex):
        for photon in tqdm(photon_list[rindex], desc=f'Running'):
            while True:
                if photon.update(self.normal_and_index_from_pos(photon.pos)): photon.move(self.res)
                photon.move(self.res)
                if not self.is_photon_in_cell(photon):
                    self.photons[rindex] = [photons for photons in self.photons[rindex] if photons!=photon]
                    break
                self.check_analysis_plan(photon, rindex)


    def run(self, run_index=0, split=1, xsym=False, ysym=False):
        assert run_index<split
        n = split

        if xsym:
            self.photons = [photons for photons in self.photons if (photons.pos[0]>=0. and numpy.dot(photons.normal, [1, 0, 0])>=0.)]
        
        if ysym:
            self.photons = [photons for photons in self.photons if (photons.pos[1]>=0. and numpy.dot(photons.normal, [0, 1, 0])>=0.)]


        self.photons = numpy.asarray(self.photons)
        self.photons = numpy.array_split(self.photons, n)

        self.split_photon_lists = numpy.asarray([self.photon_lists for i in range(n)])

        self.run_photon(self.photons, run_index)

        for index, photon_list in enumerate(self.photon_lists):
            for photon in photon_list.photons:
                if xsym: self.photon_lists[index].add_symmetric_xphoton(photon)
        
        for index, photon_list in enumerate(self.photon_lists):
            for photon in photon_list.photons:
                if ysym: self.photon_lists[index].add_symmetric_yphoton(photon)

        return self.split_photon_lists[run_index]

    def merge_photon_lists(self, my_photon_lists):
        for index, photon_list in enumerate(my_photon_lists):
            for photon in photon_list.photons:
                self.photon_lists[index].append_photon(photon)
        return self.photon_lists

    def reset(self):
        pass

