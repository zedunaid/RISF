import pandas as pd
import math

class RISF:

    def __init__(self,d_risk):
        """
           The constructor initializes all the required constants that can be used for simulation
        """
        self.intial_depth=40
        self.d_risk=d_risk
        self.cols = {'date': 'Date',
                     'avg_air_tem_f': 'Average Air Temperature (F)',
                     'avg_air_tem_c': 'Average Air Temperature (C)',
                     'max_air_tem_F': 'Maximum Air Temperature (F)',
                     'max_air_tem_c': 'Maximum Air Temperature (C)',
                     'min_air_tem_f': 'Minimum Air Temperature (F)',
                     'min_air_tem_c': 'Minimum Air Temperature (C)',
                     'avg_rel_humidity_per': 'Average Relative Humidity (%)',
                     'max_rel_humidity_per': 'Maximum Relative Humidity (%)',
                     'min_rel_humidity_per': 'Minimum Relative Humidity (%)',
                     'total_per': 'Total Precipitation (in)',
                     'avg_solar_rad': 'Average Solar Radiation (W/m2)',
                     'avg_wind_speed_mph': 'Average Wind Speed (mph)',
                     'avg_wind_speed_mps': 'Average Wind Speed (mps)',
                     'max_wind_speed_mph': 'Maximum Wind Speed (mph)',
                     'max_wind_speed_mps': 'Maximum Wind Speed (mps)',
                     'min_wind_speed_mph': 'Minimum Wind Speed (mph)',
                     'min_wind_speed_mps': 'Minimum Wind Speed (mps)'
                     }
        self.constants_windVelocity = [4.87, 67.8, 5.42] #make constants
        self.constants_z2 = 10.0  #elevation constants
        self.constants_deltas = [4098, 0.6108, 17.27, 237.3]  #make constants
        self.constants_radiation_a = 0.65
        self.constants_radiation_b = -0.85
        self.constants_depth_calculate=[1.4235e+02, -1.5777e-05, 2.6079e-13]


        self.constants_lagoon_surface_area=[151254, -387.11]
        self.constants_lagoon_volume=[ 10994931.66,-94087.98,119.94]

        self.animal_count = 6120
        self.animalType='Feeder-finish'
        self.manure_generation_rate = {
                     'Farrow-wean': 4.39,
                     'Farrow-feeder': 5.30,
                     'Farrow-finish': 14.38,
                     'Wean-feeder': 0.30,
                     'Wean-finish': 1.17,
                     'Feeder-finish': 1.37
                     }




    def getDelta(self, average_air_tem_c):
        """
        Calculates delta for each average air temperature (C)
        :param self:
        :param average_air_tem_c:
        :return list of deltas calculated for average air temperature (C)
        """
        delta_average_air_tem_c = [(1000 * (
                self.constants_deltas[0] * self.constants_deltas[1] * math.exp((self.constants_deltas[2] * tem) / (tem + self.constants_deltas[3]))) / (
                                        pow(tem + self.constants_deltas[3], 2))) for tem in average_air_tem_c]
        return delta_average_air_tem_c

    # Todo: Get sensible function name


    def ea(self, min_air_tem_c, max_air_tem_c, max_rel_humidity_per, min_rel_humidity_per):
        """
        :param self:
        :param min_air_tem_c:
        :param max_air_tem_c:
        :param max_rel_humidity_per:
        :param min_rel_humidity_per:
        :return:
        """
        return 1000 * ((self.constants_deltas[1] * math.exp(
            (self.constants_deltas[2] * max_air_tem_c) / (max_air_tem_c + self.constants_deltas[3])) * (min_rel_humidity_per / 100))
                       + ((self.constants_deltas[1] * math.exp(
                    (self.constants_deltas[2] * min_air_tem_c) / (min_air_tem_c + self.constants_deltas[3]))) * (
                                  max_rel_humidity_per / 100))) / 2

    def es(self, min_air_tem_c, max_air_tem_c):
       """
       :param self:
       :param min_air_tem_c:
       :param max_air_tem_c:
       :return:
       """
       return 1000 * ((self.constants_deltas[1] * math.exp((self.constants_deltas[2] * max_air_tem_c) / (max_air_tem_c + self.constants_deltas[3])))
                       + (self.constants_deltas[1] * math.exp(
                   (self.constants_deltas[2] * min_air_tem_c) / (min_air_tem_c + self.constants_deltas[3])))
                       ) / 2

    def getNetRadiation(self, avg_solar_rad):
        """
        Calculated net radiation from average solar radiation
        :param self:
        :param avg_solar_rad:
        :return list of net radiation :
        """
        return [radiation * self.constants_radiation_a + self.constants_radiation_b for radiation in avg_solar_rad]



    def getWindSpeed(self, avg_wind_speed):
        """
        Calculates wind velocity form average wind speed (mps)
        :param self:
        :param avg_wind_speed:
        :return list of wind velocity:
        """
        return [velocity * (self.constants_windVelocity[0] / math.log(self.constants_windVelocity[1] * self.constants_z2 - self.constants_windVelocity[2])) for velocity
                in avg_wind_speed]


    def getAirDensity(self, average_air_tem_c):
        """
        Calculates air density from average air temperature (C)
        :param self:
        :param average_air_tem_c:
        :return list of air density:
        """
        p = 101325
        r = 287.5
        return [p / (r * (tem + 273)) for tem in average_air_tem_c]



    def calculateEvaporationRate(self, delta_air_tem_c, e_as, e_a, air_density, net_radiation,
                                 avg_wind_speed_at_two_meters):
        """
        Calculates final evaporation from given parameters
        :param self:
        :param delta_air_tem_c:
        :param e_as:
        :param e_a:
        :param air_density:
        :param net_radiation:
        :param avg_wind_speed_at_two_meters:
        :return list of evaporation rate for each day
        """
        evaporation = [0.0] * len(delta_air_tem_c)

# the below will remain constant thorughout
        lv = 2450000.0
        rho_w = 997.0
        p = 101325.0
        z0 = 0.00002
        k = 0.4
        gam = 67.4
        e_const = 0.622
        e_mult = 8.64e+7
        #  r = 287.5

        for i in range(len(delta_air_tem_c)):
            evaporation[i] = e_mult * (
                        (1 / (delta_air_tem_c[i] + gam)) * (((net_radiation[i] * delta_air_tem_c[i]) / (lv * rho_w)) +
                                                            ((e_const * (pow(k, 2)) * air_density[i] *
                                                              avg_wind_speed_at_two_meters[i] * gam) / (p * rho_w *
                                                                                                        (pow(math.log(self.constants_z2 / z0), 2)))) * (e_as[i] - e_a[i])))

        return evaporation



    def calculateLagoonSurfaceArea(self,depth):
        """
        Calculate surface area for lagoon
        :param self:
        :param depth:
        :return Surface area for lagoon from depth
        """
        return self.constants_lagoon_surface_area[0] + self.constants_lagoon_surface_area[1] * depth


    def calculateLagoonVolume(self,depth):
        """
        Calculate volume for lagoon
        :param self:
        :param depth:
        :return: Volume for lagoon from depth
        """

        return  self.constants_lagoon_volume[0] + self.constants_lagoon_volume[1]*depth + self.constants_lagoon_volume[2]*depth * depth


    def getDepthFromVol(self, volume):
        """
        Calculate depth from given volume
        :param self:
        :param volume:
        :return list of depths for each volumes
        """
        return self.constants_depth_calculate[0] + self.constants_depth_calculate[1] * volume + self.constants_depth_calculate[2] * volume*volume


    def calculateNewDepths(self, evaporation_rate, rainfall_rate):
        """
        Calculates new depth from evaportion_rate,rainfall_rate and animal_waste
        :param evaporation_rate:
        :param rainfall_rate:
        :return:
        """
        new_depth=[]
        animal_waste = self.animal_count*self.manure_generation_rate[self.animalType]  # might change later

        depth=self.intial_depth
        overflow_flag=[]


        for i in range(len(evaporation_rate)):

            lagoon_surface_area= self.calculateLagoonSurfaceArea(depth)
            lagoon_volume=self.calculateLagoonVolume(depth)

            evaporation_vol = evaporation_rate[i]*lagoon_surface_area
            rainfall_vol = rainfall_rate[i]*lagoon_surface_area

            lagoon_volume = lagoon_volume + rainfall_vol+animal_waste - evaporation_vol

            #Get new depth from the updated lagoon volume
            depth = self.getDepthFromVol(lagoon_volume)

            if depth <=1:
                overflow_flag.append("Lagoon overflow even")

            elif depth<= self.d_risk:
                overflow_flag.append("“overflow risk")
            else:
                overflow_flag.append("N/A")

            new_depth.append(depth)
        return new_depth,overflow_flag


#Main function
    def readInputFile(self):
        """
        This method gets the input file, invokes required methods to get the final depth for each day and decides about the appropriate flags
        :param self:
        :return:
        """
        print("Starting ...")
        workbook = pd.read_excel('wd_CLIN.xlsx', skiprows=12)
        try:
            workbook.fillna(0)
            average_air_tem_c = []

            # ToDo: find good variable name for this list

            e_a = []
            e_as = []

            net_radiation = self.getNetRadiation(workbook[self.cols['avg_solar_rad']])
            avg_wind_speed_at_two_meters = self.getWindSpeed(workbook[self.cols['avg_wind_speed_mps']])
            rainfall_rate = workbook[self.cols['total_per']]
            for index, row in workbook.iterrows():

                if row[self.cols['min_air_tem_c']] == '#VALUE!':
                    continue

                e_a.append(self.ea(float(row[self.cols['min_air_tem_c']]), float(row[self.cols['max_air_tem_c']]),
                                   float(row[self.cols['max_rel_humidity_per']]),
                                   float(row[self.cols['min_rel_humidity_per']])))
                e_as.append(self.es(row[self.cols['min_air_tem_c']], row[self.cols['max_air_tem_c']]))
                average_air_tem_c.append((row[self.cols['max_air_tem_c']] + row[self.cols['min_air_tem_c']]) / 2)
                avg_wind_speed_at_two_meters.append(row[self.cols['avg_wind_speed_mps']])

            air_density = self.getAirDensity(average_air_tem_c)
            delta_air_tem_c = self.getDelta(average_air_tem_c)

            # calculate evaporation rate
            evaporation_rate = self.calculateEvaporationRate(delta_air_tem_c, e_as, e_a, air_density, net_radiation,
                                                        avg_wind_speed_at_two_meters)

            print("\nPrinting average air temperature in Celcius")
            print(average_air_tem_c)
            print("\nPrint e_a")
            print(e_a)
            print("\nPrinting e_as")
            print(e_as)
            print("\nPrinting  net solar radiation")
            print(net_radiation)
            print("\nPrinting average wind velocity")
            print(avg_wind_speed_at_two_meters)
            print("\nPrinting air density")
            print(air_density)
            print("\nPrinting delta air temperature in celcius")
            print(delta_air_tem_c)
            print("\nPrinting evaporation")
            print(evaporation_rate)


            new_depth,overflow_flag= self.calculateNewDepths(evaporation_rate, rainfall_rate)

            print("\nPrinting new Depth")
            print(new_depth)

            print("\nPrinting overflow flags")
            print(overflow_flag)

        except:
            print("Error occurred while calculating evaporation")

