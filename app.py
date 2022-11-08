from crypt import methods
from fileinput import filename
# from urllib import request
from flask import Flask , render_template, request , send_file
from RISF.RISF import *

app = Flask(__name__)
var={}
def getSimulationReport(farmFile,fieldFile,climateFile):
    obj = RISF()
    #Get Farm Details from Excel
    obj.getFarmDetails(farmFile)

    # #Get Field Details from Excel
    obj.getFieldDetails(fieldFile)

    #Get New Depths and Irrigate fields
    obj.readInputFile(climateFile)
    return obj.file_name

@app.route('/download',methods=['GET'])
def download():
    print(var['FileName'],"pring")    
    return send_file(var['FileName'],as_attachment=True, download_name="Simulation-Reports.xlsx")



@app.route('/submit',methods=['POST'])
def submitFile():
   print("inside submit")
   farmFile = request.files['farmFile']
   fieldFile= request.files['fieldFile']
   climateFile= request.files['climateFile']
   print(farmFile,fieldFile)
   var['FileName'] = getSimulationReport(farmFile,fieldFile,climateFile)
   sleep(5)
   print(var['FileName']
   return render_template('download.html')


@app.route('/')
def hello_world():
   #  getSimulationReport()

    return render_template('index.html')

if __name__ == '__main__':
    app.run(debug=True)
