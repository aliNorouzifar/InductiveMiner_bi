import PySimpleGUI as sg

def input():
    sg.theme("DarkGrey7")
    layout = [  [sg.Text('Settings', font='Any 20')],
                    [sg.Text('support ', justification='center', font='Any 12'), sg.Slider(range=(0,1), resolution=0.1, orientation='h', border_width =2, s=(100,20), key='-sup-')],
                    [sg.Text('ratio      ', font='Any 12') , sg.Slider(range=(0,1), resolution=0.1, orientation='h', border_width =2,s=(100,20), key='-r-')],
                    [sg.Text('Desirable Log(.xes)   ', font='Any 12'), sg.FileBrowse(key="-Desirable Log-")],
                    [sg.Text('Undesirable Log(.xes)', font='Any 12'), sg.FileBrowse(key="-Undesirable Log-")],
                    [sg.Button('OK')]]

    window = sg.Window('Inputs', layout, size=(600, 300))

    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED or event == "OK":
            break

    window.close()
    return float(values["-sup-"]), float(values["-r-"]), values["-Desirable Log-"], values["-Undesirable Log-"]


def output(acc,F1,acc_s,F1_s,fitp,prc,time):
    layout = [[sg.Text("Results", font='Any 20')],
              [sg.Text("acc: " + acc, font='Any 12')],
              [sg.Text("F1: " + F1, font='Any 12')],
              [sg.Text("acc_s: " + acc_s, font='Any 12')],
              [sg.Text("F1_s: " + F1_s, font='Any 12')],
              [sg.Text("fit+: " + fitp, font='Any 12')],
              [sg.Text("prc: " + prc, font='Any 12')],
              [sg.Text("time: " + time, font='Any 12')],
              [sg.Button("Exit")]]

    window = sg.Window("Outputs", layout, size=(600, 300))
    while True:
        event, values = window.read()
        if event == "Exit" or event == sg.WIN_CLOSED:
            break
    window.close()

