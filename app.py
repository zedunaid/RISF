from crypt import methods
from fileinput import filename
# from urllib import request
from flask import Flask , render_template, request , send_file
from RISF.RISF import *

app = Flask(__name__)
var={}
def getSimulationReport(farmFile,fieldFile):
    obj = RISF()
    #Get Farm Details from Excel
    obj.getFarmDetails(farmFile)

    # #Get Field Details from Excel
    obj.getFieldDetails(fieldFile)

    #Get New Depths and Irrigate fields
    obj.readInputFile()
    return obj.file_name

@app.route('/download',methods=['GET'])
def download():
    
    return send_file('RISF/Output_Files/Report-'+var['FileName'],as_attachment=True, download_name="Simulation-Reports.xlsx")



@app.route('/submit',methods=['POST'])
def submitFile():
   print("inside submit")
   farmFile = request.files['farmFile']
   fieldFile= request.files['fieldFile']
   print(farmFile,fieldFile)
   var['FileName'] = getSimulationReport(farmFile,fieldFile)

   return render_template('download.html')


@app.route('/')
def hello_world():
   #  getSimulationReport()

    return render_template('index.html')

if __name__ == '__main__':
    app.run(host='0.0.0.0',port='8080',debug=True)