from flask import Flask, render_template, request, send_file
from RISF.RISF import *

app = Flask(__name__)
var = {}


def getSimulationReport(farmFile, fieldFile, climateFile):
    obj = RISF()
    # Get Farm Details from Excel
    obj.getFarmDetails(farmFile)

    # #Get Field Details from Excel
    obj.getFieldDetails(fieldFile)

    # Get New Depths and Irrigate fields
    obj.readInputFile(climateFile)
    return obj.file_name


@app.route('/submit', methods=['POST'])
def submitFile():
    farmFile = request.files['farmFile']
    fieldFile = request.files['fieldFile']
    climateFile = request.files['climateFile']
    result_file = climateFile.filename.split('.')[0] + "_" + farmFile.filename.split('.')[0] + "_" + \
                  fieldFile.filename.split('.')[0] + ".xlsx"
    FileName = getSimulationReport(farmFile, fieldFile, climateFile)
    return send_file(FileName, as_attachment=True, download_name=result_file)


@app.route('/')
def hello_world():
    #  getSimulationReport()

    return render_template('index.html')


if __name__ == '__main__':
    app.run(debug=True)
