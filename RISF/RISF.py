import pandas as pd
import math
import random
from datetime import datetime
from collections import defaultdict

class RISF:

    def __init__(self):
        """
           The constructor initializes all the required constants that can be used for simulation
        """
        self.d_irrigate=30

        self.mmToInchConversion = 0.0393701
        self.inchToFeetConversion=0.08333
        self.cuFtToGal=7.48052
        self.minVolPerField = 10000  #(gallon/acre)
        self.maxVolPerField = 27000  #(gallon/acre)

        self.cols = {'date': 'Date',
                     'max_air_tem_c': 'Maximum Air Temperature (C)',
                     'min_air_tem_c': 'Minimum Air Temperature (C)',
                     'max_rel_humidity_per': 'Maximum Relative Humidity (%)',
                     'min_rel_humidity_per': 'Minimum Relative Humidity (%)',
                     'total_per': 'Total Precipitation (in)',
                     'avg_solar_rad': 'Average Solar Radiation (W/m2)',
                     'avg_wind_speed_mps': 'Average Wind Speed (mps)',
                     }
        self.constants_windVelocity = [4.87, 67.8, 5.42]  # make constants
        self.constants_z2 = 10.0  # elevation constants
        self.constants_deltas = [4098, 0.6108, 17.27, 237.3]  # make constants

        self.constants_radiation_a = 0.65
        self.constants_radiation_b = -0.85

        # self.Lagoon_d_Coeffs = [1.4235e+02, -1.5777e-05, 2.6079e-13]
        #
        # self.Lagoon_A_Coeffs = [151254, -387.11]
        # self.Lagoon_V_Coeffs = [10994931.66, -94087.98, 119.94]

        # self.AnimalCount = 6120
        # self.AnimalType = 'Feeder-finish'
        self.wastageWater=0.25+1

        self.manure_generation_rate = {
            'Farrow-wean': 4.39,
            'Farrow-feeder': 5.30,
            'Farrow-finish': 14.38,
            'Wean-feeder': 0.30,
            'Wean-finish': 1.17,
            'Feeder-finish': 1.37
        }

        self.field_parameters = {
            "03-01": [1, "09-30", 3.0, "bermuda", 6.0, 46.0],
            "02-15": [2, "06-30", 4.0, "corn", 174.0, 0.78],
            "03-15": [3, "09-15", 8.0, "soybeans", 40.0, 3.91],
            "09-01": [4, "03-31", 5.0, "wheat", 100.0, 1.14],
        }



    def getFarmDetails(self):
        try:

            workbook = pd.read_excel('Input_Template_Farm_new.xlsx', skiprows=1,nrows=11,usecols=range(1,2))
            data=(workbook['Value'].values.tolist())

            #Assigning value to the variables from excel
            self.AnimalType=(data[0])
            self.AnimalCount=float(data[1])
            self.WaterLoss_Coeff= (data[2])
            self.Lagoon_V_Coeffs= [float(v.strip()) for v in data[3].split(',')]
            self.Lagoon_A_Coeffs=  [float(a.strip()) for a in data[4].split(',')]
            self.Lagoon_d_Coeffs=  [float(d.strip()) for d in data[5].split(',')]
            self.d_initial= float(data[6])
            self.d_start=float(data[7])    #ToDo When to use this
            self.d_stop = float(data[8])
            self.d_freeboard=float(data[9])
            self.Avg_N_lbkgal=float(data[10])
            print(data)
        except:
            print("Error while fetching farm data,please check the input file")

    def generateRandomVolume(self,acre):
        """

        :return:  random volume in the range 10,000- 27,000
        """
        return random.randint(self.minVolPerField*acre+1,self.maxVolPerField*acre-1)   #randomly generates volume from the range

    def getDelta(self, average_air_tem_c):
        """
        Calculates delta for each average air temperature (C)
        :param self:
        :param average_air_tem_c:
        :return list of deltas calculated for average air temperature (C)
        """
        delta_average_air_tem_c = [(1000 * (
                self.constants_deltas[0] * self.constants_deltas[1] * math.exp(
            (self.constants_deltas[2] * tem) / (tem + self.constants_deltas[3]))) / (
                                        pow(tem + self.constants_deltas[3], 2))) for tem in average_air_tem_c]
        return delta_average_air_tem_c

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
            (self.constants_deltas[2] * max_air_tem_c) / (max_air_tem_c + self.constants_deltas[3])) * (
                                    min_rel_humidity_per / 100))
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
        return 1000 * ((self.constants_deltas[1] * math.exp(
            (self.constants_deltas[2] * max_air_tem_c) / (max_air_tem_c + self.constants_deltas[3])))
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
        return [velocity * (self.constants_windVelocity[0] / math.log(
            self.constants_windVelocity[1] * self.constants_z2 - self.constants_windVelocity[2])) for velocity
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

        for i in range(len(delta_air_tem_c)):
            evaporation[i] = e_mult * (
                        (1 / (delta_air_tem_c[i] + gam)) * (((net_radiation[i] * delta_air_tem_c[i]) / (lv * rho_w)) +
                                                            ((e_const * (pow(k, 2)) * air_density[i] *
                                                              avg_wind_speed_at_two_meters[i] * gam) / (p * rho_w *
                                                                                                        (pow(math.log(self.constants_z2 / z0), 2)))) * (e_as[i] - e_a[i])))

        return evaporation

    def calculateLagoonSurfaceArea(self, depth):
        """
        Calculate surface area for lagoon
        :param self:
        :param depth:
        :return Surface area for lagoon from depth
        """
        return self.Lagoon_A_Coeffs[0] + self.Lagoon_A_Coeffs[1] * depth

    def calculateLagoonVolume(self, depth):
        """
        Calculate volume for lagoon
        :param self:
        :param depth:
        :return: Volume for lagoon from depth
        """

        return self.Lagoon_V_Coeffs[0] + self.Lagoon_V_Coeffs[1] * depth + self.Lagoon_V_Coeffs[2] *depth*depth + self.Lagoon_V_Coeffs[3] * depth*depth*depth + self.Lagoon_V_Coeffs[4] * depth*depth*depth*depth


    def getDepthFromVol(self, volume):
        """
        Calculate depth from given volume
        :param self:
        :param volume:
        :return list of depths for each volumes
        """
        return self.Lagoon_d_Coeffs[0] + self.Lagoon_d_Coeffs[1] * volume + self.Lagoon_d_Coeffs[2] *volume*volume + self.Lagoon_d_Coeffs[3] *volume*volume*volume + self.Lagoon_d_Coeffs[4] *volume*volume*volume*volume


    def isIrrigationReq(self,irrigate_fields,lagoon_volume):
        fields_volumes=[]

        for key,val in irrigate_fields.items():
            for i in range(len(val)):
                fields_volumes.append([val[i][0]/val[i][1], key, i])

        # sort on ratio of current/total

        fields_volumes.sort(key=lambda x:x[0])
        # print(fields_volumes)
        lbsTogalConversion = 1000/self.Avg_N_lbkgal
        irrigate_vol=0
        for values in fields_volumes:
                if irrigate_fields[values[1]][values[2]][0]*lbsTogalConversion < self.minVolPerField*irrigate_fields[values[1]][values[2]][2]:
                     continue
                if irrigate_fields[values[1]][values[2]][0]*lbsTogalConversion>self.maxVolPerField*irrigate_fields[values[1]][values[2]][2]:
                      volume_alloted =  self.generateRandomVolume(irrigate_fields[values[1]][values[2]][2])
                      if volume_alloted> lagoon_volume:
                          continue

                else:
                    volume_alloted = min(lagoon_volume,irrigate_fields[values[1]][values[2]][0]*lbsTogalConversion)

                # print("lg",lagoon_volume,irrigate_fields[values[1]][values[2]][0],volume_alloted/lbsTogalConversion)
                lagoon_volume-=volume_alloted
                irrigate_vol+=volume_alloted
                irrigate_fields[values[1]][values[2]][0]-=(volume_alloted/lbsTogalConversion)

        return  irrigate_vol

    def calculateNewDepths(self, evaporation_rate, rainfall_rate,dates):
        """
        Calculates new depth from evaportion_rate,rainfall_rate and animal_waste
        :param evaporation_rate:
        :param rainfall_rate:
        :return:
        """
        new_depth = []
        animal_waste = self.AnimalCount * self.manure_generation_rate[self.AnimalType]*self.wastageWater # might change later

        depth = self.d_initial
        overflow_flag = []
        irrigate_fields= defaultdict(list)
        invent_irri_vol=[]
        invent_lagoon_vol=[]
        delta_change=[]
        daily_evap=[]
        daily_rainfall=[]
        day=0
        for i in range(len(evaporation_rate)):

            lagoon_surface_area = self.calculateLagoonSurfaceArea(depth)
            if day==0:
                lagoon_volume = self.calculateLagoonVolume(depth)
            day+=1
            evaporation_vol = evaporation_rate[i] * lagoon_surface_area * self.mmToInchConversion*self.cuFtToGal*self.inchToFeetConversion
            rainfall_vol = rainfall_rate[i] * lagoon_surface_area*self.inchToFeetConversion*self.cuFtToGal
            daily_rainfall.append(rainfall_vol)
            daily_evap.append(evaporation_vol)
            data = dates[i].split("-")
            cur_date = data[1]+'-'+data[2]


            if cur_date in self.field_parameters:
                end_date = self.field_parameters[cur_date][1]
                irrigate_fields[end_date].append([float(self.field_parameters[cur_date][2])*float(self.field_parameters[cur_date][4])*float(self.field_parameters[cur_date][5]),
                                                  float(self.field_parameters[cur_date][2])*float(self.field_parameters[cur_date][4])*float(self.field_parameters[cur_date][5]),self.field_parameters[cur_date][2]])
                # and ((depth < self.d_irrigate) or (rainfall_vol==0 and depth<self.d_stop)):   #condition to check if it rained or not


            irrigation_volume=0
            if i%7==0 and depth<self.d_stop and rainfall_vol==0:  #irrigation decision per week
                irrigation_volume=self.isIrrigationReq(irrigate_fields,lagoon_volume)

            invent_irri_vol.append(irrigation_volume)
            incrementDelta= rainfall_vol + animal_waste - evaporation_vol - irrigation_volume
            #toDo wastewater
            delta_change.append(incrementDelta)
            lagoon_volume = lagoon_volume + incrementDelta


            # Get new depth from the updated lagoon volume
            depth = self.getDepthFromVol(lagoon_volume)


            if depth <= 1:
                overflow_flag.append("Lagoon overflow event")

            elif depth <= self.d_freeboard:
                overflow_flag.append("overflow risk")
            else:
                overflow_flag.append("N/A")

            new_depth.append(depth)
            invent_lagoon_vol.append(lagoon_volume)

            print(cur_date,irrigate_fields)
            irrigate_fields.pop(cur_date,None)

        print(len(dates),len(invent_irri_vol),len(invent_lagoon_vol),len(new_depth))

        cols=[dates,invent_irri_vol,new_depth,invent_lagoon_vol,overflow_flag,delta_change,daily_rainfall,daily_evap]
        df1= pd.DataFrame(cols).transpose()
        df1.columns=["Dates","Vol used for irrigation","New depths","Lagoon Volumes","overFlow flag","Delta change","rainfall","evaporation"]
        df1.to_excel("output.xlsx",    sheet_name='Sheet_name_1')

        return new_depth, overflow_flag









    # Main function
    def readInputFile(self):
        # """
        # This method gets the input file, invokes required methods to get the final depth for each day and decides about the appropriate flags
        # :param self:
        # :return:
        # """
            print("Starting ...")
            workbook = pd.read_excel('wd_CLIN.xlsx', skiprows=12)
        # try:
            workbook.fillna(0)
            average_air_tem_c = []

            # ToDo: find good variable name for this list

            e_a = []
            e_as = []

            net_radiation = self.getNetRadiation(workbook[self.cols['avg_solar_rad']])
            avg_wind_speed_at_two_meters = self.getWindSpeed(workbook[self.cols['avg_wind_speed_mps']])
            rainfall_rate = workbook[self.cols['total_per']]
            dates=[]
            for index, row in workbook.iterrows():

                if row[self.cols['min_air_tem_c']] == '#VALUE!':
                    continue

                e_a.append(self.ea(float(row[self.cols['min_air_tem_c']]), float(row[self.cols['max_air_tem_c']]),
                                   float(row[self.cols['max_rel_humidity_per']]),
                                   float(row[self.cols['min_rel_humidity_per']])))
                e_as.append(self.es(row[self.cols['min_air_tem_c']], row[self.cols['max_air_tem_c']]))
                average_air_tem_c.append((row[self.cols['max_air_tem_c']] + row[self.cols['min_air_tem_c']]) / 2)
                avg_wind_speed_at_two_meters.append(row[self.cols['avg_wind_speed_mps']])
                dates.append(row[self.cols['date']])
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

            new_depth,overflow_flag= self.calculateNewDepths(evaporation_rate, rainfall_rate,dates)

            print("\nPrinting new Depth")
            print(new_depth)
            print(dates)
            print("\nPrinting overflow flags")
            print(overflow_flag)
        # except:
        #     print("Error occurred while calculating evaporation")

