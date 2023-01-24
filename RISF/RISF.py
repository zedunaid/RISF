import pandas as pd
import math
import os
import glob
import random
from datetime import datetime
from collections import defaultdict

"""
“Impact of Future Climate Events on NC Animal Agriculture Systems
"""

class RISF:

    def __init__(self):
        """
           The constructor initializes all the required constants that can be used for simulation
        """
        self.d_irrigate= 30 # Modeling conservative operator

        self.mmToInchConversion = 0.0393701
        self.inchToFeetConversion=0.08333
        self.cuFtToGal=7.48052
        self.minVolPerField = 10000  #(gallon/acre)
        self.maxVolPerField = 27000  #(gallon/acre)

        self.cols = {'date': 'Date',
                     'max_air_tem_f': 'Maximum Air Temperature (F)',
                     'min_air_tem_f': 'Minimum Air Temperature (F)',
                     'max_rel_humidity_per': 'Maximum Relative Humidity (%)',
                     'min_rel_humidity_per': 'Minimum Relative Humidity (%)',
                     'total_per': 'Total Precipitation (in)',
                     'avg_solar_rad': 'Average Solar Radiation (W/m2)',
                     'avg_wind_speed_ms': 'Average Wind Speed (ms)',
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
        # self.wastageWater=0.25+1

        self.manure_generation_rate = {
            'Farrow-wean': 4.39,
            'Farrow-feeder': 5.30,
            'Farrow-finish': 14.38,
            'Wean-feeder': 0.30,
            'Wean-finish': 1.17,
            'Feeder-finish': 1.37
        }

        # self.field_parameters = {
        #     "03-01": [1, "09-30", 3.0, "bermuda", 6.0, 46.0],
        #     "02-15": [2, "06-30", 4.0, "corn", 174.0, 0.78],
        #     "03-15": [3, "09-15", 8.0, "soybeans", 40.0, 3.91],
        #     "09-01": [4, "03-31", 5.0, "wheat", 100.0, 1.14],
        # }



    def getFarmDetails(self,farmFile):
        try:

            workbook = pd.read_excel(farmFile, skiprows=1,nrows=11,usecols=range(1,2))
            data=(workbook['Value'].values.tolist())

            #Assigning value to the variables from excel
            self.AnimalType=(data[0])
            self.AnimalCount=float(data[1])
            self.wastageWater= float(data[2])+1
            self.Lagoon_V_Coeffs= [float(v.strip()) for v in data[3].split(',')]
            self.Lagoon_A_Coeffs=  [float(a.strip()) for a in data[4].split(',')]
            self.Lagoon_d_Coeffs=  [float(d.strip()) for d in data[5].split(',')]
            self.d_initial= float(data[6])
            self.d_start=float(data[7])    
            self.d_stop = float(data[8])
            self.d_freeboard=float(data[9])
            self.Avg_N_lbkgal=float(data[10])
            self.d_irrigate = self.d_stop    # Conservative operator
           # print(data)
        except:
            print("Error while fetching farm data,please check the input file")


    def getFieldDetails(self,fieldFile):
        workbook = pd.read_excel(fieldFile, skiprows=1,nrows=5,usecols=range(6,12))

        # print(workbook)
        self.field_parameter={}
        self.crop_mapper={}
        for index,row in workbook.iterrows():
            window_start_date= str(row['Start Appl. Window Date'])[5:11].strip() # Extracting only year and month from the date string
            window_end_date=str(row['End Appl. Window Date'])[5:11].strip()
            crop_code = row['Crop Code']
            crop_name = row['Crop name']
            n_removal_per_yield = float(row['N removal per unit yield (lb/yield)'])
            self.field_parameter[window_start_date]= [crop_code, window_end_date, crop_name, n_removal_per_yield] # Append for same window
            self.crop_mapper[crop_code]= window_start_date
       
        workbook = pd.read_excel(fieldFile)
        self.number_of_fields= workbook.iloc[0][1]
        self.field_input={}

        for i in range(3,3+self.number_of_fields):
            field_id = workbook.iloc[i][0]
            field_area = float(workbook.iloc[i][1])
            field_crop_code = workbook.iloc[i][2]
            crop_yield_per_acre = float(workbook.iloc[i][3])

            window_start_date = self.crop_mapper[field_crop_code]

            crop_code = self.field_parameter[window_start_date][0]
            window_end_date = self.field_parameter[window_start_date][1]
            crop_name = self.field_parameter[window_start_date][2]
            n_removal_per_yield = self.field_parameter[window_start_date][3]
            
           #Taking cumulative of nitrogen required incase there are more than one occurence of same field
            if window_start_date not in self.field_input:
                self.field_input[window_start_date] = [[crop_code]+ [window_end_date]+ [field_area]+ [crop_name]+ [crop_yield_per_acre]+ [n_removal_per_yield]+ [field_id]]
            else:
                self.field_input[window_start_date].append([crop_code]+ [window_end_date]+ [field_area]+ [crop_name]+ [crop_yield_per_acre]+ [n_removal_per_yield]+ [field_id])

        #print("printing final field")
        #print(self.field_input) # [CC,CN,FieldArea,NRemovel,CropYield,NRemoval,FieldId]

        """
           self.field_parameters = {
            "03-01": [2, "09-30", 3.0, "corn", 6.0, 46.0],
           
            "03-15": [3, "09-15", 8.0, "soybeans", 40.0, 3.91],
            "09-01": [4, "03-31", 5.0, "wheat", 100.0, 1.14],
        }"""


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


    def isIrrigationReq(self,irrigate_fields,lagoon_volume,vol_per_field):
        fields_volumes=[]

        for key,val in irrigate_fields.items():
            for i in range(len(val)):
                current_N = val[i][0]
                total_N = val[i][1]
                window_end_date = key
                fields_volumes.append([current_N/total_N, window_end_date, i])

        # sort on ratio of current/total
        fields_volumes.sort(key=lambda x:x[0])

        # print(fields_volumes)
        lbsTogalConversion = 1000/self.Avg_N_lbkgal

        irrigate_vol=0
        for values in fields_volumes:
                window_end_date = values[1]
                index = values[2]
                irr_list = irrigate_fields[window_end_date][index]
                current_N_credits = irr_list[0]
                field_area = irr_list[2]
                field_id = irr_list[3]
                N_concentration = self.Avg_N_lbkgal

                # Calculating Nitrogen Volume in gallons
                field_N_Vol = (current_N_credits / N_concentration) * 1000
                min_field_N_Vol = self.minVolPerField * field_area
                max_field_N_Vol = self.maxVolPerField * field_area
                
                # Calculating volume allocated to the field
                if field_N_Vol < min_field_N_Vol:
                     continue
                if field_N_Vol > max_field_N_Vol:
                      volume_alloted =  self.generateRandomVolume(field_area)
                      if volume_alloted> lagoon_volume:
                          continue

                else:
                    volume_alloted = min(lagoon_volume,field_N_Vol)

                vol_per_field[field_id]+=volume_alloted
                lagoon_volume = lagoon_volume - volume_alloted
                irrigate_vol = irrigate_vol + volume_alloted

                # Conversion of Volume units (Gallons) to N units (lbs)
                amount_N_removed = (volume_alloted * N_concentration / 1000)

                # Reducing N Credits from the field
                irrigate_fields[window_end_date][index][0]-= amount_N_removed

        return  irrigate_vol

    def calculateNewDepths(self, evaporation_rate, rainfall_rate,dates):
        """
        Calculates new depth from evaportion_rate,rainfall_rate and animal_waste
        :param evaporation_rate:
        :param rainfall_rate:
        :return:
        """
        new_depth = []

        # Animal Waste Calculation (AnimalCount * Manure Generation rate of that Animal * % Wastage Water)
        animal_waste = self.AnimalCount * self.manure_generation_rate[self.AnimalType]*self.wastageWater # might change later

        depth = self.d_initial
        overflow_flag = []
        irrigate_fields= defaultdict(list)
        invent_irri_vol=[]
        invent_lagoon_vol=[]
        delta_change=[]
        daily_evap=[]
        daily_rainfall=[]
        exceedance_lagoon_volume=[]
        day=0
        volume_allocation_per_field={}
        
        # Initializing volume allocation of each field as empty 
        for i in range(1,self.number_of_fields+1):
            volume_allocation_per_field[i]=[]

        for i in range(len(evaporation_rate)):

            lagoon_surface_area = self.calculateLagoonSurfaceArea(depth)
            if day==0:
                # Calculating initial lagoon volume based on init depth
                lagoon_volume = self.calculateLagoonVolume(depth)
            day+=1

            # Evaporation Volume in gallons = Evap rate * lagoon surface area
            evaporation_vol = evaporation_rate[i] * lagoon_surface_area * self.mmToInchConversion * self.cuFtToGal * self.inchToFeetConversion

            # Rainfall Volume in gallons = Rainfall rate * lagoon surface area
            rainfall_vol = rainfall_rate[i] * lagoon_surface_area * self.inchToFeetConversion * self.cuFtToGal
            daily_rainfall.append(rainfall_vol)
            daily_evap.append(evaporation_vol)
            data = dates[i].split("-")
            cur_date = data[1]+'-'+data[2]

            if cur_date in self.field_input:
                for window in self.field_input[cur_date]:
                    window_end_date = window[1]
                    field_area = float(window[2])
                    crop_yield_per_acre = float(window[4])
                    n_removal_per_yield = float(window[5])
                    field_id = window[6]
                    current_N = field_area * crop_yield_per_acre * n_removal_per_yield
                    total_N = current_N
                    if window_end_date not in irrigate_fields:
                        irrigate_fields[window_end_date] = [[current_N, total_N, field_area, field_id ]]
                    else:
                        irrigate_fields[window_end_date].append([current_N, total_N,field_area, field_id])



            irrigation_volume=0
            vol_per_field=[0]*(self.number_of_fields+1)
            if i%7==0 and depth<=self.d_irrigate and rainfall_vol==0:  #irrigation decision per week
                irrigation_volume=self.isIrrigationReq(irrigate_fields,lagoon_volume,vol_per_field)



            #Get volume allocated for each field
            for i in range(1,self.number_of_fields+1):
                volume_allocation_per_field[i].append(vol_per_field[i])

            # toDo find a new variable name for invent_irri_vol 
            invent_irri_vol.append(irrigation_volume)
            incrementDelta= rainfall_vol + animal_waste - evaporation_vol - irrigation_volume

            #toDo wastewater
            delta_change.append(incrementDelta)
            lagoon_volume = lagoon_volume + incrementDelta

            #toDo allocation of volume to each field during the weekly cycle

            prev_day_depth = depth
            # Get new depth from the updated lagoon volume
            depth = self.getDepthFromVol(lagoon_volume)


            if depth <= 1:
                overflow_flag.append("Lagoon overflow event")

            elif depth <= self.d_freeboard:
                overflow_flag.append("overflow risk")
            else:
                overflow_flag.append("N/A")

            if depth<0:
                exceedance_lagoon_volume.append(incrementDelta)
                #The depth remains as previous day as we move the extra volume to exceedance
                depth=prev_day_depth
                #Removing the added increment volume from lagoon
                lagoon_volume-=incrementDelta
            else:
                exceedance_lagoon_volume.append(0)

            new_depth.append(depth)
            invent_lagoon_vol.append(lagoon_volume)

         #   print(cur_date,irrigate_fields)
            irrigate_fields.pop(cur_date,None)

        #print(len(dates),len(invent_irri_vol),len(invent_lagoon_vol),len(new_depth))

        self.file_name =  str(datetime.now())+".xlsx"
        directory_output = os.getcwd()+'/Output_Files/'

        if not os.path.exists(directory_output):
            os.makedirs(directory_output)
        else:
            files = glob.glob(directory_output+"*")
            for f in files:
                os.remove(f)

        cols=[dates,invent_irri_vol,new_depth,invent_lagoon_vol,overflow_flag,delta_change,daily_rainfall,daily_evap,exceedance_lagoon_volume]
        for i in range(1, self.number_of_fields+1):
            tem=volume_allocation_per_field[i]
            cols.append(tem)

        df1= pd.DataFrame(cols).transpose()

        col_labels=["Dates","Vol used for irrigation","New depths","Lagoon Volumes","overFlow flag","Delta change","rainfall","evaporation","exceedance Lagoon Volume"]
        for i in range(1, self.number_of_fields+1):
            col_labels.append("Field "+str(i))

        df1.columns=col_labels
       # df1.to_excel(directory_output+"Report-"+self.file_name, sheet_name='Report')
        df2 = df1.copy()

        df2['YearMonth'] = pd.to_datetime(df1['Dates']).apply(lambda x: '{year}-{month}'.format(year=x.year, month=x.month))

        df2 = df2.groupby('YearMonth',sort=False).sum()
        col_labels=["Delta change","rainfall","evaporation","exceedance Lagoon Volume"]
        for i in range(1, self.number_of_fields+1):
            col_labels.append("Field "+str(i))

        df2 = df2[col_labels]
#        print(df2)
        #df2.to_excel(directory_output+"Aggregated-Report-"+self.file_name, sheet_name='aggregation')
        with pd.ExcelWriter(directory_output+"Report-"+self.file_name) as writer:
   
    # use to_excel function and specify the sheet_name and index
    # to store the dataframe in specified sheet
            df1.to_excel(writer, sheet_name="Report")
            df2.to_excel(writer, sheet_name="Aggregation")

        self.file_name = directory_output+"Report-"+self.file_name
        return new_depth, overflow_flag



    # Main function
    def readInputFile(self,climateFile):
            """
            This method gets the input file, invokes required methods to get the final depth for each day and decides about the appropriate flags
            :param self:
            :return:
            """
            print("Starting ...")
            workbook = pd.read_excel(climateFile, skiprows=12)
        # try:
            workbook.fillna(0)
            workbook.replace(to_replace='QCF',value=0,inplace=True)
            average_air_tem_c = []

            # ToDo: find good variable name for this list

            e_a = []
            e_as = []

            net_radiation = self.getNetRadiation(workbook[self.cols['avg_solar_rad']])
            avg_wind_speed_at_two_meters = self.getWindSpeed(workbook[self.cols['avg_wind_speed_ms']])
            rainfall_rate = workbook[self.cols['total_per']]
            dates=[]
            for index, row in workbook.iterrows():
                min_air_tem_f = row[self.cols['min_air_tem_f']]
                min_air_tem_c = 0.5556*(float(min_air_tem_f)-32)
                max_air_tem_f = row[self.cols['max_air_tem_f']]
                max_air_tem_c = 0.5556*(float(max_air_tem_f)-32)
                max_rel_humidity_per = float(row[self.cols['max_rel_humidity_per']])
                min_rel_humidity_per = float(row[self.cols['min_rel_humidity_per']])

                if row[self.cols['min_air_tem_f']] == '#VALUE!':
                    continue

                e_a.append(self.ea( min_air_tem_c, max_air_tem_c, max_rel_humidity_per, min_rel_humidity_per))
                e_as.append(self.es(min_air_tem_c, max_air_tem_c))
                average_air_tem_c.append((max_air_tem_c + min_air_tem_c)/2)
                avg_wind_speed_at_two_meters.append(row[self.cols['avg_wind_speed_ms']])
                dates.append(row[self.cols['date']])
            air_density = self.getAirDensity(average_air_tem_c)
            delta_air_tem_c = self.getDelta(average_air_tem_c)

            # calculate evaporation rate
            evaporation_rate = self.calculateEvaporationRate(delta_air_tem_c, e_as, e_a, air_density, net_radiation,
                                                        avg_wind_speed_at_two_meters)

            new_depth,overflow_flag= self.calculateNewDepths(evaporation_rate, rainfall_rate,dates)

