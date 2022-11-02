from crypt import methods
# from urllib import request
from flask import Flask , render_template, request , send_file
from RISF.RISF import *

app = Flask(__name__)

def getSimulationReport(farmFile,fieldFile):
    obj = RISF()
    #Get Farm Details from Excel
    obj.getFarmDetails(farmFile)

    # #Get Field Details from Excel
    obj.getFieldDetails(fieldFile)

    #Get New Depths and Irrigate fields
    obj.readInputFile()
    return obj.file_name

@app.route('/submit',methods=['POST'])
def submitFile():
   print("inside submit")
   farmFile = request.files['farmFile']
   fieldFile= request.files['fieldFile']
   print(farmFile,fieldFile)
   file_name = getSimulationReport(farmFile,fieldFile)

   return send_file('./RISF/Output_Files/Report-'+file_name,as_attachment=True, download_name="farm_File.xlsx")


@app.route('/')
def hello_world():
   #  getSimulationReport()

    return render_template('index.html')

if __name__ == '__main__':
    app.run(host='0.0.0.0',port='8080',debug=True)