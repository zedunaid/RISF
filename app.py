from flask import Flask, render_template, request, send_file
from RISF.RISF import *
import time
import calendar
import os
import glob
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

def getManualSimulationReport(inputdata, isDeveloper = False):
    file_path = {'GOLD': os.path.join('Input_files','GOLD_2015-2020.xlsx'),
                 'KINS': os.path.join('Input_files','KINS_2015_2020.xlsx'),
                 'CLIN': os.path.join('Input_files','CLIN_2015_2020.xlsx'),
                 'WILD': os.path.join('Input_files','WILD_2015-2020.xlsx')
                }
    station = inputdata['station']
    obj = RISF()
    obj.saveInputData(inputdata)
    if isDeveloper:
        directory_output = os.getcwd()+'/Output_Files/'
        dev_file_path = directory_output + 'Dev-Report.xlsx';
        if not os.path.exists(directory_output):
            os.makedirs(directory_output)
        else:
            print("Deleting Dev-Report.xlsx")
            files = glob.glob(directory_output+"Dev-Report.xlsx")
            for f in files:
                os.remove(f)

        no_of_simulations = int(inputdata['noofsimulations'])
        start_time = time.time()
        for i in range(0,no_of_simulations):
            # Get Farm Details from Manual
            obj.getManualFarmDetails(inputdata,isDeveloper)
            # Get Field Details from Manual
            obj.getManualFieldDetails(inputdata, isDeveloper)
            # Get New Depths and Irrigate fields
            obj.readInputFile(file_path[station],isDeveloper, i+1, no_of_simulations)
        print("---Total execution time is %s seconds ---" % (time.time() - start_time))
    else:
        # Get Farm Details from Manual
        obj.getManualFarmDetails(inputdata,isDeveloper)
        # Get Field Details from Manual
        obj.getManualFieldDetails(inputdata, isDeveloper)
        # Get New Depths and Irrigate fields
        obj.readInputFile(file_path[station],isDeveloper)
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
    return render_template('home.html')

@app.route('/upload', methods=['POST'])
def upload_data():
    return render_template('index.html')

@app.route('/manualhome', methods=['POST'])
def manual_home():
    return render_template('manual.html')

@app.route('/developerhome', methods=['POST'])
def developer_home():
    return render_template('developer.html')

@app.route('/manual', methods=['POST'])
def manual_data():
    if request.method == "POST":
        inputdata = request.form
        result_file = inputdata['station']+ "_"+ str(calendar.timegm(time.gmtime()))+ ".xlsx"
        FileName = getManualSimulationReport(inputdata)
    return render_template('report.html')
    #return send_file(FileName, as_attachment=True, download_name=result_file)

@app.route('/developer', methods=['POST'])
def developer_data():
    if request.method == "POST":
        inputdata = request.form
        result_file = "DevReport_"+ str(calendar.timegm(time.gmtime()))+ ".xlsx"
        FileName = getManualSimulationReport(inputdata, True)
    #return render_template('report.html')
    return send_file(FileName, as_attachment=True, download_name=result_file)

@app.route('/download', methods=['POST'])
def download_report():
    FileName = glob.glob(os.path.join('Output_Files','*.xlsx'))[0]
    out_file_name = FileName.split('/')[1]
    return send_file(FileName, as_attachment=True, download_name=out_file_name)

@app.route('/visualize', methods=['POST'])
def visualize_data():
    FileName = glob.glob(os.path.join('Output_Files','*.xlsx'))[0]
    InputFileName = glob.glob(os.path.join('Input_Files','Input-*.xlsx'))[0]
    obj = RISF()
    string, vol_field_data, nitro_field_data = obj.visualize(FileName)
    animaltype, avgnumofanimals, length, width, depth, sideslope,inp_table = obj.fetchInputData(InputFileName)
    vol_field_data = vol_field_data.applymap(lambda x: "{:,}".format(x) if isinstance(x, int) else x)
    nitro_field_data = nitro_field_data.applymap(lambda x: "{:.2f}".format(x) if isinstance(x, float) else x)
    return render_template('visualize.html',plot=string,  table=vol_field_data.reset_index().to_html(classes="data",index=False, header=True), table2=nitro_field_data.reset_index().to_html(classes="data",index=False, header=True),animaltype=animaltype, avgnumofanimals=avgnumofanimals, length=length, width=width, depth=depth, sideslope=sideslope,inp_table=inp_table.to_html(classes="data",index=False, header=True))

if __name__ == '__main__':
    app.run(host='0.0.0.0',port=5000,debug=True)
