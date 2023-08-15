### Libraries
import numpy as np                   # NumpPy
from matplotlib import pyplot as plt # MatPlotLib
from numpy import pi                 # Pi

class Constants:
    STEFAN = 5.67e-8        # Stefan's Constant (Wm-2K-4)
    ALBEDO_LAND = 0.2       # average albedo of dry land                  
    ALBEDO_ICE = 0.85       # average albedo of ice / snow covered land
    ALBEDO_WATER = 0.4      # average albedo of open water (oceans / seas)
    K_AIR = 0.025           # thermal conductivity of air
    K_OCEAN = 0.592         # thermal conductivity of (sea)water

    def __init__(self) -> None:
        pass

class Planet(Constants):

    def __init__(self, luminosity: float, distance: float, radius: float, emissivity: float, pAIR: float, pOCEAN: float) -> None:
        super().__init__()
        self.luminosity = luminosity
        self.distance = distance
        self.radius = radius
        self.emissivity = emissivity
        self.P_AIR = pAIR
        self.P_OCEAN = pOCEAN

        self.global_intensity = (self.luminosity / (4 * pi * self.distance ** 2)) / 4

        self.latitudes = np.arange(87.5, -92.5, -5)

        self.composition = [[0.65, 0.19, 0.16] for i in range(int(len(self.latitudes)/3))]
        self.composition.extend([[0.81, 0.0, 0.19] for i in range(int(len(self.latitudes)/3))])
        self.composition.extend([[0.65, 0.19, 0.16] for i in range(int(len(self.latitudes)/3))])

        self.albedos = [(i[0] * self.ALBEDO_WATER + i[1] * self.ALBEDO_ICE + i[2] * self.ALBEDO_LAND) for i in self.composition]
        self.distances = [(self.radius * 5 * pi / 180 * np.cos(np.radians(theta * pi / 180)))/1000 for theta in self.latitudes]

    def local_intensity(self):
        latitudes = list(self.latitudes)
        local_intensity_array = [(self.global_intensity * np.cos(np.radians(theta)) * self.albedos[latitudes.index(theta)]) for theta in latitudes]

        return np.array(local_intensity_array)
    
    def surface_temperatures(self, intensity_array: list):
        surface_temperatures_array = [((j / self.STEFAN) ** 0.25) for j in intensity_array]
        
        return np.array(surface_temperatures_array)
    
    def delta_temperatures(self, temperatures_array: list):
        temperatures_array = list(temperatures_array)
        delta_temps = []
        one_to_eighteen = []
        nineteen_to_thirty_five = []

        for i in range(1, 19):
            delta_t = temperatures_array[i] - temperatures_array[i-1]
            one_to_eighteen.append(delta_t)

        for i in range(18, 36):
            delta_t = temperatures_array[i-1] - temperatures_array[i]
            nineteen_to_thirty_five.append(delta_t)
        
        delta_temps.append(one_to_eighteen)
        delta_temps.append(nineteen_to_thirty_five)

        return np.array(delta_temps)

    def meridional_heat_transport(self, delta_temperatures_array: list):
        delta_temperatures_array = list(delta_temperatures_array)
        meridional_heat_transport = []
        lst1 = []
        lst2 = []

        for t in delta_temperatures_array[0]:
            distance = self.distances[list(delta_temperatures_array[0]).index(t)]
            phi = (t * ((self.P_AIR * -1 * self.K_AIR) + (self.P_OCEAN * -1 * self.K_OCEAN)))/distance
            lst1.append(phi)

        for t in delta_temperatures_array[1]:
            distance = self.distances[list(delta_temperatures_array[1]).index(t)]
            phi = (t * ((self.P_AIR * -1 * self.K_AIR) + (self.P_OCEAN * -1 * self.K_OCEAN)))/distance
            lst2.append(phi)

        meridional_heat_transport.extend(lst1)
        meridional_heat_transport.extend(lst2)

        return np.array(meridional_heat_transport)
    
    def greenhouse_temperatures(self, temperatures_array: list):
        greenhouse_temperatures = [temp * ((1 / (1 - self.emissivity / 2)) ** 0.25) for temp in temperatures_array]
        
        return np.array(greenhouse_temperatures)
    
    def update_intensity(self, input_array: list, meridional_transfer_array: list):
        input_array = list(input_array)
        for i in range(len(input_array) - 1):
            if i == 1:
                input_array[i] -= meridional_transfer_array[i]
            elif i <= int(len(input_array) / 2):
                input_array[i] += meridional_transfer_array[i-1]
                input_array[i] -= meridional_transfer_array[i]
            elif i > int(len(input_array) / 2) and i < len(input_array):
                input_array[i] += meridional_transfer_array[i+1]
                input_array[i] -= meridional_transfer_array[i] 
            elif i == len(input_array):
                input_array[i] -= meridional_transfer_array[i] 

        return input_array

    def update_albedos(self, latitude_temperatures_dictionary: dict):
        ice_list = []
        for latitude in latitude_temperatures_dictionary:
            if latitude_temperatures_dictionary[latitude][-1] > 263.15:
                ice = False
            else:
                ice = True
            ice_list.append(ice)

        for i in self.composition:
            if not ice_list[self.composition.index(i)]:
                i[0] = 0.81
                i[1] = 0.0
            else:
                update = i[1] * 0.00005
                if latitude_temperatures_dictionary[latitude][-1] > latitude_temperatures_dictionary[latitude][-2]:
                    i[0] += update
                    i[1] -= update
                else:
                    i[0] -= update
                    i[1] += update

def plot_latitude_temp(time: list, latitude_temps: dict, latitude):
    plt.plot(time, latitude_temps[latitude])
    plt.xlabel("Time / years")
    plt.ylabel("Temperature / degrees K")
    plt.grid(True)
    plt.show() 

def plot_global_temp(time: list, average_temps: np.ndarray):
    plt.plot(time, average_temps)
    plt.xlabel("Time / years")
    plt.ylabel("Temperature / degrees K")
    plt.grid(True)
    plt.show()

def model():
    LUMINOSITY = 3.828e26       # luminosity of the star (W) - initialised to 1 L☉
    DISTANCE = 1495978707e2     # distance from planet to star (m) - initialised to 1 AU
    RADIUS = 63781e2            # radius of the planet (m)
    EMISSIVITY = 0.78           # radius of the planet (m)
    P_AIR = 0.7                 # atmospheric contribution to meridional heat transport
    P_OCEAN = 0.3               # oceanic contribution to meridional heat transport
    YEARS = 1000                # length of simulation in years

    earth = Planet(luminosity=LUMINOSITY,
                   distance=DISTANCE,
                   radius=RADIUS,
                   emissivity=EMISSIVITY,
                   pAIR=P_AIR,
                   pOCEAN=P_OCEAN)
    
    average_temps_array = []
    latitude_temps_dictionary = {i:[] for i in (np.arange(87.5, -92.5, -5))}

    local_insolation_array = earth.local_intensity()

    year = 0

    while year < YEARS:
        surface_temps_array = earth.surface_temperatures(intensity_array=local_insolation_array)
        greenhouse_temps_array = earth.surface_temperatures(intensity_array=local_insolation_array)
            
        average_temps_array.append(np.mean(greenhouse_temps_array))
            
        for i in earth.latitudes:
            latitude_temps_dictionary[i].append(greenhouse_temps_array[list(earth.latitudes).index(i)])

        for i in range(12):
            delta_temps_array = earth.delta_temperatures(temperatures_array=greenhouse_temps_array)

            meridional_heat_transport_array = earth.meridional_heat_transport(delta_temperatures_array=delta_temps_array)

            local_insolation_array = earth.update_intensity(input_array=local_insolation_array,
                                                            meridional_transfer_array=meridional_heat_transport_array)
                
        surface_temps_array = earth.surface_temperatures(intensity_array=local_insolation_array)
        greenhouse_temps_array = earth.greenhouse_temperatures(temperatures_array=surface_temps_array)

        for i in earth.latitudes:
            latitude_temps_dictionary[i].append(greenhouse_temps_array[list(earth.latitudes).index(i)])
                
        earth.update_albedos(latitude_temperatures_dictionary=latitude_temps_dictionary)

        change_intensity = 0.25 * earth.global_intensity

        if year % 50 == 0 and (year / 50) % 2 == 1:
            earth.emissivity *= 0.02
            earth.global_intensity += change_intensity
        elif year % 50 == 0 and (year / 50) % 2 == 0:
            earth.emissivity *= 0.02
            earth.global_intensity += change_intensity

        year += 1
    
    return np.array(average_temps_array), latitude_temps_dictionary

plot_global_temp(time=range(1000),
                 average_temps=model()[0])

plot_latitude_temp(time=range(2000),
                   latitude_temps=model()[1],
                   latitude=42.5)