import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm


class photon_list():
    '''
    photon_list is a class that bounds photons together by a plane and a condition. Plane is characterized by a normal vector and its value.

    Parameters
    ----------
    normal: list
        A 3 dimensional list stating the normal vector from plane
    value: float
        The plane independent variable in the form of x + y + z = value

    '''
    def __init__(self, normal: list, value: float):
        self.normal = normal / numpy.linalg.norm(normal)
        self.value = value
        self.photons = numpy.asarray([])
        self.condition_dict = dict()

    def append_photon(self, photon):
        '''
        Append a photon tp the photon_list given that the condition is True. 
    
        Parameters
        ----------
        photon: class.photon
            A photon coming from photon class.

        '''
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
        '''
        Instanteates a new photon before calling append_photon. Append photon checks if photon object satifies conditions.

        Parameters
        ----------
        sphoton: class.photon
            A photon coming from photon class.
        '''
        new_photon = photon([0, 0, 0], [0, 0, 1])
        new_photon.set_attr( sphoton.get_attr() )
        self.append_photon(new_photon)

    def distance_point_to_plane(self, pos: list):
        '''
        Calculates the distance from a given point to plane.

        Parameters
        ----------
        pos: A 3 dimensional list of a given position

        Returns
        -------
        float
            A float provenient from a numpy.abs based on point position and plane equation
        '''
        return numpy.abs(numpy.dot(self.normal, pos)-self.value)

    def add_symmetric_xphoton(self, sphoton):
        '''
        Instanteates a new X symmetric photon before callind append_photon. Append photon checks if photon object satifies conditions.

        Parameters
        ----------
        sphoton: class.photon
            A photon coming from photon class

        '''
        new_photon = photon([0, 0, 0], [0, 0, 1])
        new_photon.set_attr( sphoton.get_attr() )
        new_photon.pos = -new_photon.pos[0], new_photon.pos[1], new_photon.pos[2]
        new_photon.normal = -new_photon.normal[0], new_photon.normal[1], new_photon.normal[2]
        self.append_photon(new_photon)
    
    def add_symmetric_yphoton(self, sphoton):
        '''
        Instanteates a new Y symmetric photon before callind append_photon. Append photon checks if photon object satifies conditions.

        Parameters
        ----------
        sphoton: class.photon
            A photon coming from photon class
        '''

        new_photon = photon([0, 0, 0], [0, 0, 1])
        new_photon.set_attr( sphoton.get_attr() )
        new_photon.pos = new_photon.pos[0], -new_photon.pos[1], new_photon.pos[2]
        new_photon.normal = new_photon.normal[0], -new_photon.normal[1], new_photon.normal[2]
        self.append_photon(new_photon)

    def avg_divergence(self, vec_ref):
        '''
        Calculates the average photon divergence from given a vector reference.

        Parameters
        ----------
        vec: array_like
            A 3 dimensional array of given direction.

        Returns
        ------
        float
            A float provenient from a numpy.average based on vector direction and propagation direction for each photon

        '''
        vec_ref = vec_ref / numpy.linalg.norm(vec_ref)
        vals = numpy.average([numpy.dot(photon.normal, vec_ref)**2 for photon in self.photons])
        return vals

    def get_positions(self):
        '''
        Get positions of all photons in the photon_list.

        Returns
        -------
        numpy.asarray
            A numpy array
        '''
        vals = numpy.asarray([photon.pos for photon in self.photons])
        return vals

    def get_relative_centroid_positions(self):
        '''
        Get the relative position of all photons in the photon_list relative from centroid (calculted using avg_position)

        Returns
        -------
        numpy.asarray
        '''
        #Relative position to centroid for each photon
        vals = self.get_positions() - self.avg_position()
        return vals

    def get_relative_centroid_distances(self):
        '''
        Calculates the distances from a given point to plane.

        Returns
        -------
        numpy.asarray
        '''
        #Relative distance to centroid for each photon
        vals = numpy.sqrt(numpy.sum(numpy.power(self.get_relative_centroid_positions(), 2), axis=1))
        return vals

    def get_intensities(self):
        '''
        Gets the intensity of all photons in the photon_list

        Returns
        -------
        numpy.asarray
            A photons-dimensional numpy array
        '''
        vals = numpy.asarray([photon.intensity for photon in self.photons])
        return vals

    def avg_position(self):
        '''
        Calculates the average position of all photons in the photon_list.

        Returns
        ----------
        numpy.asarray
            A 3 dimensional numpy array

        '''
        vals = numpy.asarray([numpy.average([photon.pos[axis] for photon in self.photons]) for axis in range(3)])
        return vals
    
    def max_position(self):
        '''
        Calculates the max position in each axis for all photons in photon_list.

        Returns
        -------
        numpy.asarray
            A 3 dimensional numpy array
        '''
        vals = numpy.asarray([numpy.amax([photon.pos[axis] for photon in self.photons]) for axis in range(3)])
        return vals
    
    def min_position(self):
        '''
        Calculates the min position in each axis for all photons in photon_list.

        Parameters
        ----------
        pos: A 3 dimensional list of a given position

        Returns
        -------
        numpy.asarray
            A 3 dimensional numpy array
        '''
        vals = numpy.asarray([numpy.amin([photon.pos[axis] for photon in self.photons]) for axis in range(3)])
        return vals
    
    def std_position(self):
        '''
        Calculates the standard deviation of the position for all photons in photon_list.

        Returns
        -------
        numpy.asarray
            A 3 dimensional numpy array
        '''
        vals = numpy.asarray([numpy.std([photon.pos[axis] for photon in self.photons]) for axis in range(3)])
        return vals

    def get_weighted_inverse(self):
        '''
        Calculates the division of intensity and photon relative distance from centroid fro all photons in photon_list

        Returns
        -------
        numpy array
            A n dimensional numpy array
        '''
        vals = numpy.divide(self.get_intensities(), self.get_relative_centroid_distances())
        return vals

    def get_average_weighted_inverse(self):
        '''
        Calculates the average of get_weighted_inverse.
        
        Returns
        -------
        float
        '''
        vals = float(numpy.average(self.get_weighted_inverse()))
        return vals

    def avg_distance_axis_z(self, c=[0, 0]):
        '''
        Calculates the distance for all photons in photon_list relative to a (x, y) coordinate in  a Z normal plane.

        Parameters
        ----------
        c: The position in x, y plane

        Returns
        -------
        float
        '''
        xc, yc = c[0], c[1]
        vals = float(numpy.average([numpy.sqrt((photon.pos[0]-xc)**2+(photon.pos[1]-yc)**2) for photon in self.photons]))
        return vals
    
    def avg_distance_axis_y(self, c=[0, 0]):
        '''
        Calculates the distance for all photons in photon_list relative to a (x, z) coordinate in a Y nornal plane.

        Parameters
        ----------
        c: The position in x, z plane

        Returns
        -------
        float
        '''
        xc, zc = c[0], c[1]
        vals = float(numpy.average([numpy.sqrt((photon.pos[0]-xc)**2+(photon.pos[2]-zc)**2) for photon in self.photons]))
        return vals
    
    def get_average_weighted_inverse_axis_y(self, c=[0, 0]):
        '''
        Calculates the average of get_weighted_inverse relative to a (x, z) coordinate in Y normal plane.

        Parameters
        ----------
        c: The position in x, z plane

        Returns
        ------
        float
        '''
        assert (self.normal == [0, 1, 0]).all()
        pos = self.get_positions() - numpy.asarray([c[0], self.value, c[1]])
        dist = numpy.sqrt(numpy.sum(numpy.power(pos, 2), axis=1))
        wi = numpy.divide(self.get_intensities(), dist)
        avg_wi = float(numpy.average(wi))
        return avg_wi


class photon():

    '''
    A photon is the basic structure of this simulation. It contains the attributes pos, normal, intensity, n, last_surface, refraction_count and reflection_count

    Parameters
    ----------
    pos: list
        A 3 dimensional list stating each photon position
    normal: list
        A 3 dimensional list stating each photon normal (or wavevector)
    intensity: float
        A float for the initial photon intensity

    Notes
    -----

    Other attributes of photon available to user and used during simulation run:

    n: float
        Index of refraction
    last_surface: array_like
        As photons crosses the grid, it overwrites the last valid surface normal
    refraction_count: int
        Number os refractions photon has done
    reflection_count: int
        Number os reflections photon has done
    '''

    def __init__(self, pos, normal, intensity=1.):
        self.pos = pos
        self.normal = normal / numpy.linalg.norm(normal)
        self.intensity = 1.
        self.n = 1.00
        self.last_surface = [0, 0, 0]
        self.refraction_count = 0
        self.reflection_count = 0

    def set_attr(self, values):
        '''
        This sets all photon attributes. Can be used to clone photons.

        Parameters
        ----------
        values: array_like
            An array_like object as following: **[pos, normal, intensity, n, last_surface, refraction_count, reflection_count]**
        '''

        self.pos, self.normal, self.intensity, self.n, self.last_surface, self.refraction_count, self.reflection_count = values
    
    def get_attr(self):
        '''
        This gets all photon attributes. Can be used to clone photons.

        Returns
        -------
        array_like
            An array_like object as following: **[pos, normal, intensity, n, last_surface, refraction_count, reflection_count]**
        '''
        return self.pos, self.normal, self.intensity, self.n, self.last_surface, self.refraction_count, self.reflection_count


    def reflection(self):
        '''
        Reflection photon object using last_surface attribute.

        Returns
        -------
        boolean
            True if reflection is sucessfull and false if you attempted a reflection in a grid which there was no surface normal vector.
        '''
        inc = self.normal
        inc = inc / numpy.linalg.norm(inc)
        sur_normal = self.last_surface
        
        if all([sur_normal[i]==0 for i in range(3)]):
        #if all([sur_normal[i]!=0 for i in range(3)]):
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
        Refraction photon object using last_surface attribute and a given index of refraction.

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
        
        if all([sur_normal[i]==0 for i in range(3)]):
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

    def move(self, value: float):
        '''
        Increment photon position.

        Parameters
        ----------
        value: float
            Increment photon position by value and given by normal (photon direction)
        

        Returns
        -------
        boolean
            True if is sucessfull.
        '''
        self.pos = self.pos + self.normal*value
        return True

    def update(self, value):
        '''
        Update is a core function that updates photon last_surface and photon normal by means of refractions of reflections.

        Parameters
        ----------
        value: array_like
            A 2 dimensional array where first index is the surface normal and the second is the index of refraction.

        Returns
        -------
        boolean
            True if photon suffered a refraction or reflection. False if nothing happens with photon normal.

        
        Notes
        --------
        If return is True, photon moves a single step without refraction/reflection using move. This avoids non intentional double reflections/refractions

        '''

        sur_normal, index = value

        #if sur_normal != [0, 0, 0]:
        if not all([sur_normal[i]==0 for i in range(3)]):
            self.last_surface = sur_normal
        
        if index!=self.n and index>0:
            self.refraction(index)
            return True
        elif index==-1:
            self.reflection()
            return True

        return False


class Simu:
    '''
    A simulation class that will allow one to create sereral objects and perform spacial transformations, such as rotations. This class also contains a few convenient functions in order to display data easily to the user. 

    Parameters
    ----------
    x: float
        A float representing cell size in x
    y: float
        A float representing cell size in y
    z: float
        A float representing cell size in z
    res: float
        A float for simulation resolution.

    Notes
    -----

    Other attributes of simulation available to user and used during simulation run:

    ss: int
        A subsampling factor. Default is 2.01 and it says meshing will be done using sub pixels resolution.
    grid: array_like
        A 3 dimensional array for the number of grid points.
    index: numpy array
        A numpy array initialized as **numpy.ones(grid)**. This means initial cell is under vacuum for all wavelengths.
    normal: numpy array
        A numpy array initialized as **numpy.asarray([numpy.zeros(grid), numpy.zeros(grid), numpy.zeros(grid)])** to all surface normals. Initial cell has no valid surface normals because there is no surfaces.
    photons: list
        A list containing photon objects. Those photons will propagate and be saved along the simulation.
    photon_lists: numpy array
        A numpy array cointaining photon_lists objects. photons will be saved here. Array length is given by the number of analyzes plans.
    '''

    def __init__(self, x, y, z, res):
        self.size = [x, y, z]
        self.res = res
        self.ss = 2.01 #sub sampling for creating elements
        self.grid = (int(x/res), int(y/res), int(z/res))
        self.index = numpy.ones(self.grid)
        self.normal = numpy.asarray([numpy.zeros(self.grid), numpy.zeros(self.grid), numpy.zeros(self.grid)])
        self.photons = numpy.asarray([])
        self.photon_lists=numpy.asarray([])
        self.x = numpy.linspace(-self.size[0]/2., self.size[0]/2., self.grid[0])
        self.y = numpy.linspace(-self.size[1]/2., self.size[1]/2., self.grid[1])
        self.x, self.y = numpy.meshgrid(self.x, self.y)

    def assign_normal(self, index, normal):
        '''
        Assigns a surface normal in a given index

        Parameters
        ----------
        index: array_like
            A three dimensional array with grid position (in pixels).
        normal: array_like
            A three dimensional array with the correspondent surface normal.
        '''

        if (normal!=numpy.zeros(3)).all(): 
            normal = normal / numpy.linalg.norm(normal)
        for i in range(3):
            self.normal[i][index[0], index[1], index[2]] = normal[i]
    
    def assign_block_normal(self, min_index, max_index, normal):

        '''
        Assigns normal for multiple points at once.

        Parameters
        ----------
        min_index: array_like
            A three dimensional array with grid position (in pixels) for the botton left point.
        max_index: array_like
            A three dimensional array with grid position (in pixels) for the top right point.
        normal: array_like
            A three dimensiona array with the correspondent surface normal.
        '''
        if (normal!=numpy.zeros(3)).all():
            normal = normal / numpy.linalg.norm(normal)
        for i in range(3):
            self.normal[i][min_index[0]:max_index[0], min_index[1]:max_index[1], min_index[2]:max_index[2]] = normal[i]

    def assign_n(self, index, value):
        '''
        Assigns index of refraction for a given point

        Parameters
        ----------
        index: array_like
            A three dimensional array with grid position (in pixels).
        value: float
            The desired index of refraction
        '''
        self.index[index[0], index[1], index[2]] = value
    
    def assign_block_n(self, min_index, max_index, value):
        '''
        Assigns index of refraction for multiple points at once

        Parameters
        ----------
        min_index: array_like
            A three dimensional array with grid position (in pixels) for the botton left point.
        max_index: array_like
            A three dimensional array with grid position (in pixels) for the top right point.
        normal: array_like
            A three dimensiona array with the correspondent surface normal.
        '''
        self.index[min_index[0]:max_index[0], min_index[1]:max_index[1], min_index[2]:max_index[2]] = value

    def pos_to_grid_1D(self, value, index):
        '''
        Gets grid index from position for a single dimensional axis

        Parameters
        ----------
        value: float
            Desired position
        index:
            Desired axis (0 for 'x', 1 for 'y' and 2 for 'z')

        Returns
        -------
        int
            Grid position.

        Raises
        ------
        Exception
            If value is outside cell boundaries.
        '''

        assert abs(value)<=self.size[index]/2.0

        return round((value/(self.size[index]/2.)+1)/2*self.grid[index])
        

    def pos_to_grid(self, pos):
        '''
        Get grid index from position in the 3D space.

        Parameters
        ----------
        pos: float
            Desired position

        Returns
        -------
        array_like
            Grid position

        Raises
        -----
        Exception
            If value is outside cell boundaries
        '''
        assert all([abs(pos[i])<=self.size[i]/2.0 for i in range(3)])
        indexes = [round((pos[i]/(self.size[i]/2.)+1)/2*self.grid[i]) for i in range(3)]
        if any([indexes[i]>=self.grid[i] for i in range(3)]):
            indexes = [int((pos[i]/(self.size[i]/2.)+1)/2*self.grid[i]) for i in range(3)]
        return indexes

    def grid_to_pos(self, grid):
        '''
        Get position from grid index.

        Parameters
        ----------
        grid: array_like
            Grid Position.

        Returns
        -------
        array_like
            Position.
        
        Raises
        ------
        Exeception
            If grid is outside cell boundaries.
        '''
        pos = [ (2.0*grid[i]/self.grid[i]-1) * self.size[i]/2.0 for i in range(3)]
        return pos

    def index_from_pos(self, pos):
        '''
        Get index of refraction from position.

        Parameters
        ----------
        pos: array_like
            Desired Position.

        Returns
        -------
        float
            The correspondent index of refraction.

        Raises
        ------
        Exception
            If value is outside cell boundaries
        '''
        #this returns current index of refraction of a given photon in a given position
        
        assert all([abs(pos[i])<=self.size[i]/2.0 for i in range(3)])
        indexes = self.pos_to_grid(pos)
        return self.index[indexes[0], indexes[1], indexes[2]]

    def normal_from_pos(self, pos):
        '''
        Get normal direction from position.

        Parameters
        ----------
        pos: array_like
            Desired Position.

        Returns
        -------
        array_like:
            The correspondent normal direction.
        '''
        #this returns surface normal (if exists) of a given photon in a given position

        indexes = self.pos_to_grid(pos)
        return numpy.asarray([self.normal[0][indexes[0], indexes[1], indexes[2]],
            self.normal[1][indexes[0], indexes[1], indexes[2]], 
            self.normal[2][indexes[0], indexes[1], indexes[2]]])
    
    def normal_and_index_from_pos(self, pos):
        '''
        Get both normal surface and index of refraction from pos

        Parameters
        ----------
        pos: array_like
            Desired position.

        Returns
        -------
        tuple
            A tuple in which first element is the normal and the second is the refractive index.
        '''
        #this returns surface normal (if exists) and the refractive index of a given photon in a given position
        
        assert all([abs(pos[i])<=self.size[i]/2.0 for i in range(3)])

        indexes = self.pos_to_grid(pos)

        return ([self.normal[0][indexes[0], indexes[1], indexes[2]],
            self.normal[1][indexes[0], indexes[1], indexes[2]], 
            self.normal[2][indexes[0], indexes[1], indexes[2]]], self.index[indexes[0], indexes[1], indexes[2]])


    def show_photons2D(self, my_photon_lists, plan='xy'):
        '''
        See Also
        --------

        Deprecated function. Please access photon_list properties manually.
        '''
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
        '''
        Function shows in 3D elements created. Those are the initial conditions of the problem

        Parameters
        ----------
        mode: str
            Defines what is displayed. 'all' displays all, **'all-noplan'** or **'-noplan'** excludes analyzes plans and anything else display only photons.
        '''
        
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
        '''
        Shows in 3D elements at the end of simulation.

        Parameters
        ----------
        my_photon_lists: array_like
            Array_like containing **class.photon_list** objects.
        mode: str
            Defines what is displayed. 'all' displays all, **'all-noplan'** or **'-noplan'** excludes analyzes plans and anything else display only photons.
        '''

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


    def d2_source(self, r, c=[0, 0, 0], normal=[0, 0, 1], na = 0.0, angles=1):
        '''
        Creates a 2d source along axis z with a numerical aperture

        Parameters
        ----------
        r: float
            Source radius. If 0.0, creates a point source.
        c: array_like:
            Source center.
        normal: array_like
            Source normal. Photon propagation direction
        na: float
            Source numerical aperture along axis *[1, 1, 0]*.
        angles: int
            Numerical aperture discretization

        Examples
        -------
        >>> d2_source(0.1, [0, 0, 0], [0, 0, 1], na=0.0, angles=1)
        
        Creates a collimated source with radius 0.1 centered at [0, 0, 0].

        >>> d2_source(0.0, [0, 0, 0], [0, 0, 1], na=0.22, angles=11)

        Creates a point source centered at [0, 0, 0] with numerical aperture 0.22. Angles are discretized by 0.04 (-0.22 to 0.22) in 11 points.

        >>> d2_source(0.3, [0, 0, 0], [0, 0, 1], na=0.22, angles=11)

        Creates a diverging source centered at [0, 0, 0] with numerical aperture 0.22 and radius 0.3.

        Raises
        ------
        AssertionError
            If numerical aperture is higher than 0, propagation vector must be aligned with an axis. This
            function does not create an angle cone for an arbitrary propagation vector.
        '''


        if na>0:
            check_list = [normal[i] == 0 for i in range(3)]
            check_list.sort()
            assert check_list[1] #second element is true so we have at least two zeros

            rot_axis = numpy.subtract([1, 1, 1], normal)
            rot_vecs = []
            for i in range(3):
                if i!=numpy.where(rot_axis==0)[0]:
                    vec = [0, 0, 0]
                    vec[i] = 1
                    rot_vecs.append(vec)

            rot_vecs = numpy.asarray(rot_vecs)
            na_mesh = numpy.linspace(-na, na, angles)

            all_vecs = numpy.asarray([normal])
            for xna_vals in na_mesh:
                for yna_vals in na_mesh:
                    new_vec = numpy.add(
                        numpy.add(numpy.multiply(rot_vecs[0], xna_vals), numpy.multiply(rot_vecs[1], yna_vals)),
                    normal)
                #print(new_vec, vec)
                    if new_vec.tolist() not in all_vecs.tolist():
                        all_vecs = numpy.append(all_vecs, [new_vec], axis=0)

        else:
            all_vecs = numpy.asarray([normal])

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
                            #print(normal2)
                            self.photons = numpy.append(self.photons, photon([xpos, ypos, zc], normal2))
                            self.photons = numpy.append(self.photons, photon([-xpos+2*xc, ypos, zc], normal2))
                            self.photons = numpy.append(self.photons, photon([xpos, -ypos+2*yc, zc], normal2))
                            self.photons = numpy.append(self.photons, photon([-xpos+2*xc, -ypos+2*yc, zc], normal2))
    
    def rotate(self, ang, axis, origin, ROI = None):
        '''
        Rotate a selected ROI in a given direction and origin. This function only rotates grid points that have index of refraction different of 1.

        Parameters
        ----------
        ang: float
            Angle of rotation (in radians).
        axis: array_like
            A 3 dimensional axis of rotation.
        origin: array_like
            A 3 dimensional origin position.
        ROI: array_like
            A 6 dimensional ROI in the form of **[xmin, xmax, ymin, ymax, zmin, zmax]**.

        Raises
        ------
        Exception
            If there is no element to rotate.
        '''
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

    def create_parabolic_surface_element(self, c, n, th, wid, pp):
        '''
        Creates a parabolic element surface.

        Parameters
        ----------
        c: array_like
            A 3 dimensional position.
        n: float
            The index of refraction. -1 creates a reflective material.
        th: float
            Thickness of the element. Y direction total size is given by this value.
        wid: float
            Width of the element. X direction total size is given by this value.
        pp: float
            P parameter of the parabola.

        Notes
        -----
        Parabolic curve is defined as:

        .. math:: 2pp * (Z-  Z0) = (Y-Y0)^{2} - (X-X0)^{2}

        '''

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
        '''
        Creates a spherical element.

        Parameters
        ----------
        c: array_like
            A 3 dimensional sphere center.
        r: float
            Sphere radius.
        n: float
            Element index of refraction.

        '''
        
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
        '''
        Creates a rectangular element.

        Parameters
        ----------
        val: array_like
            A 6 dimensional array containing [xmin, xmax, ymin, ymax, zmin, zmax]
        n: float
            The index of refraction
        normal: array_like
            A 3 dimensional array containing the surface normal. 

        Notes
        -----
        You can use a rectangular element to destroy some of your active elements by setting index of refraction as 1.0 and normal vector as [0, 0, 0].
        '''
        
        xmin, xmax, ymin, ymax, zmin, zmax = val
        xc, yc, zc = (xmin+xmax)/2., (ymin+ymax)/2., (zmax+zmin)/2.
        xl, yl, zl = (xmax-xmin)/2., (ymax-ymin)/2., (zmax-zmin)/2.
        assert xl>=0 and yl>=0 and zl>=0

        x = [x for x in numpy.arange(xc-xl, xc+xl+self.res, self.res/self.ss) if abs(x)<self.size[0]/2.]
        y = [y for y in numpy.arange(yc-yl, yc+yl+self.res, self.res/self.ss) if abs(y)<self.size[1]/2.]
        z = [z for z in numpy.arange(zc-zl, zc+zl+self.res, self.res/self.ss) if abs(z)<self.size[2]/2.]
        
        min_index = self.pos_to_grid([min(x), min(y), min(z)])
        max_index = self.pos_to_grid([max(x), max(y), max(z)])
        
        self.assign_block_n(min_index, max_index, n)
        self.assign_block_normal(min_index, max_index, normal)

    def create_cylinder_element(self, center, radius, length, n, normal):
        '''
        Creates a cilindrical element.

        Parameters
        ----------
        center: array_like
            A 3 dimensional of center position.
        radius: float
            Cilinder radius.
        length: float
            Cilinder length.
        n: float
            Index of refraction.
        normal: array_like
            A 3 dimensional array of cilinder axis.
        '''
        
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
        '''
        Creates a thin lens based on other geometrical elements. Allow one to create either a plane-convex or a convex-plane lens.

        Parameters
        ----------
        cplane: array_like
            A 3 dimensional array representing the center of the lens. This point is always in the plane side of the lens.
        focus: float
            Lens focus.
        aperture: float
            Lens aperture. 
        index_refr: float
            Lens index of refraction.
        lens_type: str
            Must be either plane-convex or convex-plane

        Raises
        ------
        AssertionError
            If aperture is bigger than focus or if type is not 'plane-convex' or 'convex-plane', raises an exception error.
        '''
        xc, yc, zp = cplane
        r = focus/2.0
        assert focus/2.>=aperture/2.
        assert lens_type=='plane-convex' or lens_type=='convex-plane'
        if lens_type=='plane-convex': zc = zp - numpy.sqrt((focus/2.0)**2 - (aperture/2.)**2) #plane_convex
        if lens_type=='convex-plane': zc = zp + numpy.sqrt((focus/2.0)**2 - (aperture/2.)**2) #convex-plane

        u = 2 * self.res #margin of exclusion

        self.create_sphere_element([xc, yc, zc], focus/2., index_refr)
        if lens_type=='plane-convex':
            self.create_rectangle_element([xc-r-u, xc+r+u, yc-r-u, yc+r+u, zc-r-u, zp+u/2.], 1.0, [0, 0, 0])
        elif lens_type=='convex-plane':
            self.create_rectangle_element([xc-r-u, xc+r+u, yc-r-u, yc+r+u, zp-u/2., zc+r+u], 1.0, [0, 0, 0])
        
        self.create_cylinder_element([xc, yc, zp], aperture/2+self.res, self.res, index_refr, [0, 0, 1])

    def distance(self, vec1, vec2):
        '''
        Distance between two points.

        Parameters
        ----------
        vec1: array_like
            A 3 dimensional array representing first position.
        vec2: array_like
            A 3 dimensional array representing second position.


        '''
        return numpy.sum(numpy.power(numpy.subtract(vec1, vec2), 2))**0.5

    def is_photon_in_cell(self, photon):
        '''
        Check if photon is in the cell.

        Parameters
        ----------
        photon: class.photon

        Returns
        -------
        bool
        '''
        pos = photon.pos
        val = all([abs(pos[i])<self.size[i]/2. for i in range(3)])
        if val:
            return True
        else:
            return False

    def check_analysis_plan(self, photon):
        '''
        Check if photon is in a created analysis plan.

        Parameters
        ----------
        photon: class.photon
        rindex: int
            Run index. Used for multi process simulations.
        '''
        pos = photon.pos
        for index, planes in enumerate(self.photon_lists):
            if planes.distance_point_to_plane(pos)<=self.res/2.0:
                planes.add_photon(photon)

    def create_analysis_plan(self, normal, value, **kargs):
        '''
        Create an analysis plan.

        Parameters
        ----------
        normal: array_like
            A 3 dimensional array of plan normal vector.
        value: float
            Indenpendent value of a plan equation.
        **kargs: dict
            A dictionary that will be used as a condition dictionary. Photons will be added the plan only if satisfies this condition.

        Notes
        -----
        Condition must be an attribute of photon. It can take two forms: a single number or a tuple. Single number represents minimal values, while tuple represents a possible interval

        Examples
        --------
        >>> create_analysis_plan([0, 1, 0], a, reflection_count = 1)

        Example above creates a analyses plan with equation 
        
        .. math:: y=a

        and with a condition that reflection_count is higher than 1 (at least one reflection).
        
        >>> create_analysis_plan([1, 2, 3], 3, reflection_count = (1, 1))
        
        Example above creates a analyses plan with equation 
        
        .. math:: (x + y + z) / \sqrt{14} = 3

        and with a condition that reflection_count is exactly 1. Factor square root of 14 comes from vector normalization
        
        >>> create_analysis_plan([1, 0, 0], 0.4, refraction_count = (2, 30))
        
        Example above creates a analyses plan with equation 
        
        .. math:: x=0.4

        and with a condition that refraction_count is higher than 2 but smaller than 30.
        '''
        plist = photon_list(normal, value)
        plist.condition_dict = kargs
        self.photon_lists = numpy.append(self.photon_lists, plist)
    
    def run_photon(self, initial_photons, rindex):
        '''
        Simulatiom run for the initial set (or subset) of photons. It terminates when all photons leave the cell.

        Parameters
        ----------
        initial_photon_list: array_like
            The complete set (or subset) of the intial photons.
        rindex: int
            Run index. Used for multi process simulations. 

        See Also
        --------
        run: 
            Starts the simulation for each set or subset of photon. Can be used with multi processing.

        '''
        for photon in tqdm(initial_photons[rindex], desc=f'Running'):
            #while self.is_photon_in_cell(photon):
            while True:
                if photon.update(self.normal_and_index_from_pos(photon.pos)):
                    photon.move(self.res)
                    self.check_analysis_plan(photon)
                photon.move(self.res)
                if not self.is_photon_in_cell(photon):
                    self.photons[rindex] = [photons for photons in self.photons[rindex] if photons!=photon]
                    break
                self.check_analysis_plan(photon)

    def prepare_acquisition(self, split = 1):
        assert self.photons.size
        assert split>0

        self.photons = numpy.asarray(self.photons)
        self.photons = numpy.array_split(self.photons, split)

    def run(self, run_index=0, multiprocessing=False, return_dict=dict() , xsym=False, ysym=False):
        '''
        Simulation run main function.

        Parameters
        ----------
        run_index: int
            Run index. Used for multi process simulation
        split: int
            Split initial photons in equal parts of 'split'. Used for multi process simulation
        xsym: bool
            If simulation is symmetric with respect to x axis, you can reduce number of running photons by two.
        ysym: bool
            If simulation is symmetric with respect to y axis, you can reduce number of running photons by two.

        Returns
        -------
        array_like
            Returns an array of class.photon_list. Lenght is given by the name of planes created. Each element in the list has class.photon objects. Length for each element is the number of photons saved for the correspondent plane.

        Raises
        ------
        AssertionError
            Raises an assertionError if run_index is smaller than split. You cannot run subset 1 if you have not splitted your initial photon array. Maximum run_index is always split-1. Raises an AssertionError also if there is no photons in the initial set.
        '''

        if not multiprocessing:
            self.prepare_acquisition(1)

        assert len(self.photons)>run_index

        if xsym:
            self.photons = [photons for photons in self.photons if (photons.pos[0]>=0. and numpy.dot(photons.normal, [1, 0, 0])>=0.)]
        if ysym:
            self.photons = [photons for photons in self.photons if (photons.pos[1]>=0. and numpy.dot(photons.normal, [0, 1, 0])>=0.)]

        self.run_photon(self.photons, run_index)

        for index, photon_list in enumerate(self.photon_lists):
            for photon in photon_list.photons:
                if xsym: self.photon_lists[index].add_symmetric_xphoton(photon)
        
        for index, photon_list in enumerate(self.photon_lists):
            for photon in photon_list.photons:
                if ysym: self.photon_lists[index].add_symmetric_yphoton(photon)

        return_dict[run_index] = self.photon_lists
        return self.photon_lists

    def merge_photon_lists(self, my_photon_lists):
        '''
        Merge photon lists in a single array. This function is used for multi processing simulations.

        Parameters
        ----------
        my_photon_lists: array_like
            An array containing objects class.photon_list

        Returns
        -------
        array_like
            Return an array like containing photon_lists. 
        '''
        for index, photon_list in enumerate(my_photon_lists):
            for photon in photon_list.photons:
                if photon not in self.photon_lists[index].photons:
                    self.photon_lists[index].append_photon(photon)
                else:
                    #This does not work using the return dict. It is very strange that we have different objects
                    print('***WARNING***: Duplicating photon in the list. Photon will not be added.')
        return self.photon_lists

    def reset_photons(self):
        '''
        Reset initial photon lists in order to run simulation a second time.

        '''
        self.photons = numpy.asarray([])

        return True

    def reset_structures(self):
        '''
        Reset initial structures in order to run simulation a second time.
        '''
        self.index = numpy.ones(self.grid)
        self.normal = numpy.asarray([numpy.zeros(self.grid), numpy.zeros(self.grid), numpy.zeros(self.grid)])
        
        return True
    
