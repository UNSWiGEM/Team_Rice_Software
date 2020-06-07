import PySimpleGUI as sg
import create_thermometers as tm

layout = [[sg.Text("---Please enter your sequences---")],
          [sg.Text("5' context: "), sg.Input()],
          [sg.Text("3' context: "), sg.Input()],
          [sg.Text("Variable region: "), sg.Input()],
          [sg.Text("---Please enter your evoluation parameters---")],
          [sg.Text("Number of generations: "), sg.Input()],
          [sg.Text("Crossover probability: "), sg.Input()],
          [sg.Text("Shuffle mutation probability: "), sg.Input()],
          [sg.Text("Point mutation probability: "), sg.Input()],
          [sg.Text("---Please enter your temperature bounds---")],
          [sg.Text("Upper bound (°C): "), sg.Input()],
          [sg.Text("Lower bound (°C): "), sg.Input()],
          [sg.Text("---Please enter your ribosome binding site location---")],
          [sg.Text("Start location: "), sg.Input()],
          [sg.Text("End location: "), sg.Input()],
          [sg.Text("---Please enter a random seed---")],
          [sg.Text("Random seed: "), sg.Input()],
          [sg.Button("Run"), sg.Exit()]]

def main():
	window = sg.Window("Thermometer Evolution", layout)

	while True:
		event, values = window.Read()

		if event is None or event == 'Exit':
			break

		context_before = str(values[0])
		context_after = str(values[1])
		variable_region = str(values[2])

		ngen = int(values[3])
		cxpb = float(values[4])
		mutpb = float(values[5])
		otmut = float(values[6])

		upper_temp = int(values[7])
		lower_temp = int(values[8])

		rbs_start = int(values[9])
		rbs_end = int(values[10])

		seed = int(values[11])

		if event == 'Run':
			sg.Popup("started thermometer running.")
			tm.run(context_before, context_after, variable_region,
				   ngen, cxpb, mutpb, otmut,
				   upper_temp, lower_temp,
				   rbs_start, rbs_end,
				   seed)

			sg.Popup("Thermometer running finished!")

			break

	window.Close()

if __name__ == "__main__":
	main()
