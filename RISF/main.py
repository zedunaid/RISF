from RISF import  *
if __name__ == "__main__":

    obj = RISF()
    #Get Farm Details from Excel
    obj.getFarmDetails()

    # #Get Field Details from Excel
   # obj.getFieldDetails()

    #Get New Depths and Irrigate fields
    obj.readInputFile()



    