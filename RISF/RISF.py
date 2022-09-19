import pandas as pd
import math

class RISF:

    def __init__(self):

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
        self.windVelocity = [4.87, 67.8, 5.42]
        self.z2 = 10.0
        self.deltas = [4098, 0.6108, 17.27, 237.3]
        self.radiation_a = 0.65
        self.radiation_b = -0.85


        self.lagoon_surface_area=[151254,-387.11]

    def getDelta(self, average_air_tem_c):

        delta_average_air_tem_c = [(1000 * (
                self.deltas[0] * self.deltas[1] * math.exp((self.deltas[2] * tem)) / (tem + self.deltas[3])) / (
                                        pow(tem + self.deltas[3], 2))) for tem in average_air_tem_c]
        return delta_average_air_tem_c

    # Todo: Get sensible function name

    def ea(self, min_air_tem_c, max_air_tem_c, max_rel_humidity_per, min_rel_humidity_per):

        return 1000 * ((self.deltas[1] * math.exp(
            (self.deltas[2] * max_air_tem_c) / (max_air_tem_c + self.deltas[3])) * (min_rel_humidity_per / 100))
                       + ((self.deltas[1] * math.exp(
                    (self.deltas[2] * min_air_tem_c) / (min_air_tem_c + self.deltas[3]))) * (
                                  max_rel_humidity_per / 100))) / 2

    def es(self, min_air_tem_c, max_air_tem_c):

        return 1000 * ((self.deltas[1] * math.exp((self.deltas[2] * max_air_tem_c) / (max_air_tem_c + self.deltas[3])))
                       + (self.deltas[1] * math.exp(
                    (self.deltas[2] * min_air_tem_c) / (min_air_tem_c + self.deltas[3])))
                       ) / 2

    def getNetRadiation(self, avg_solar_rad):
        return [radiation * self.radiation_a + self.radiation_b for radiation in avg_solar_rad]

    def getWindSpeed(self, avg_wind_speed):

        return [velocity * (self.windVelocity[0] / math.log(self.windVelocity[1] * self.z2 - self.windVelocity[2])) for velocity
                in avg_wind_speed]

    def getAirDensity(self, average_air_tem_c):
        p = 101325
        r = 287.5
        return [p / (r * (tem + 273)) for tem in average_air_tem_c]

    def calculateEvaporation(self, delta_air_tem_c, e_as, e_a, air_density, net_radiation,
                             avg_wind_speed_at_two_meters):

        evaporation = [0.0] * len(delta_air_tem_c)

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
                                                            (pow(math.log(self.z2 / z0), 2)))) * (e_as[i] - e_a[i])))

        return evaporation

    def calculateLagoonSurfaceArea(self,depth):

        return self.lagoon_surface_area[0]+self.lagoon_surface_area[1]*depth




#Main function
    def readInputFile(self):
        print("Starting ...")
        workbook = pd.read_excel('wd_CLIN.xlsx', skiprows=12)

        try:

            workbook.replace(regex='#VALUE!', value=0.0)
            workbook.fillna(0)

            average_air_tem_c = []

            # ToDo: find good variable name for this list

            e_a = []
            e_as = []

            net_radiation = self.getNetRadiation(workbook[self.cols['avg_solar_rad']])
            avg_wind_speed_at_two_meters = self.getWindSpeed(workbook[self.cols['avg_wind_speed_mps']])

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

            # calculate evaporation

            evaporation = self.calculateEvaporation(delta_air_tem_c, e_as, e_a, air_density, net_radiation,
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
            print(evaporation)


            #########Lagoon surface Area #########
            depth=0
            lagoon_surface_area= self.calculateLagoonSurfaceArea(depth)
            evaporation_volume=[ evaporation_rate*lagoon_surface_area for evaporation_rate in evaporation]

            print("\nPrinting evaporation volume")
            print(evaporation_volume)


        except:
            print("Error occurred while calculating evaporation")
