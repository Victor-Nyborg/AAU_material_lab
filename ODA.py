"""
Author: Victor Ørsøe Nyborg von@build.aau.dk

Description:
    Takes the raw data file from the Oxygen Diffusion Apparatus (ODA) LabVIEW program and calculate the oxygen diffusivity.

    Features:
        - If a leak of the ODA is known, the script can also correct the oxygen diffusion.
        - Can calculate the oxygen diffusion based on the Taylor method and the Currie method.
"""
from typing import Any
from dataclasses import dataclass
import pandas as pd
import numpy as np
from scipy import  stats
from scipy.optimize import brentq
import tkinter as tk
from tkinter import filedialog, ttk

@dataclass
class PrintWithMuteFunction:
    """
    This print text, with the possibility to mute the print function into the console.
    """
    muted: bool = False # Mute function, change to True for mute any text print.
    base_length: int = 50 # Base character length of headers

    def __call__(self, *text: Any, **kwargs) -> None:
        if not self.muted:
            print(*text, **kwargs)
        else:
            pass

    def double(self, *text: Any) -> None:
        """
        print the text twice in console.
        """
        self(*text)
        self(*text)

    def title(self, *text: Any, width: int = base_length, symbol: str = '=', end: str = '') -> None:
        """
        print a title in the console surrounded by the defined symbol.
        :param text: Text which should be printed.
        :param width: The width span of the header (n characters).
        :param symbol: The character to use for the header, around the text.
        :param end: print any after the title (could be a new line).
        :return: None
        """
        text = ' '.join(text)
        text = ' ' + text + ' '
        text = f'{text:{symbol}^{width}}'
        self(symbol*width)
        self(text)
        self(symbol * width, end=end)

    def header(self, *text: Any, width: int = base_length, symbol: str = '=') -> None:
        """
        print a header in console, where the section name is centered in the middle of the header.
        :param text: Text which should be printed.
        :param width: The width span of the header (n characters).
        :param symbol: The character to use for the header, around the text.
        :return: None
        """
        text = ' '.join(text)
        text = ' ' + text + ' '
        text = f'{text:{symbol}^{width}}'
        self(text)

    def line(self, symbol: str, length: int = base_length):
        text = symbol*length
        self(text)


@dataclass
class OxygenDiffusion:
    """
    This is the object that calculates the oxygen diffusion, based on a file path to the Oxygen Diffusion Apparatus
    LabVIEW program result file.

    INPUT:
        filepath: str | File path to the data file.
    """
    filepath: str # File path to raw data file

    muted: bool = False # Mute any print in console
    step: int = 10 # Number of data point to increase for each iteration, in the best fit algorithm

    def __post_init__(self):
        self.Print: PrintWithMuteFunction = PrintWithMuteFunction(self.muted)
        self.Print.title('Running ODA data treatment', end='\n\n')
        ## Sample info
        self.Print.header('Load', self.filepath.split('\\')[-1])
        self.data = self._read_data_file()
        self.Print('File loaded')

        self.max_voltage: float = self.data['Sensor signal'].max() # Ambient oxygen voltage signal, sample has reached equilibrium.
        self.init_voltage: float = self.data['Sensor signal'].min() # Initial condition of oxygen voltage signal inside oxygen gas chamber.
        self.Print('Measurement conditions determined',
              ('Ambient', float(self.max_voltage)),
              ('Initial', float(self.init_voltage)),
              sep='\n\t')

        denominator = self.init_voltage - self.max_voltage
        arg = (self.data['Sensor signal'] - self.max_voltage) / denominator
        temporary_array = np.full(len(self.data), np.nan) # Temporary array which is going to be filled with the calculated natural log of the relative oxygen concentration inside the oxygen gas chamber.
        mask = (self.data['Sensor signal'] >= self.max_voltage).to_numpy()
        temporary_array[mask] = 0
        temporary_array[~mask] = np.log(arg.loc[~mask])
        self.data['Ln(C)'] = temporary_array
        del denominator, arg, temporary_array, mask
        self.Print('Natural log calculated')

        best_fit_r2: float = 0
        start_index: int =  0
        slope: float = 0
        min_fraction = 0.4 # TODO consider other forms for criteria (n amount of time, or n number of data points, voltage relative span)
        # At least 40% of the data should be included in the fit, the start index increase by 10 between each iteration

        while start_index < len(self.data) * (1 - min_fraction):
            end_index: int = start_index + int(len(self.data) * min_fraction)
            while end_index < len(self.data):
                temporary_data = self.data.loc[start_index: end_index, ['Elapsed time', 'Ln(C)']]
                s, intercept, r, p, std_err = stats.linregress(temporary_data['Elapsed time'], temporary_data['Ln(C)'])
                r2 = r**2

                if r2 > best_fit_r2:
                    best_fit_r2 = r2
                    slope = s

                end_index += self.step
            start_index += self.step
        self.Print('Best slope calculated',
              ('r2:', float(round(best_fit_r2, 2))),
              ('Slope:', float(round(slope, 6))),
              sep='\n\t')
        self.slope: float = slope
        self.Print.line('=')


    def _read_data_file(self) -> pd.DataFrame:
        raw_data = pd.read_csv(self.filepath, sep='\t') \
            .drop(columns='Voltage oxygen sensor ambient [V]')
        raw_data['Time stamp'] = pd.to_datetime(raw_data['Time stamp'], format='%Y/%m/%d - %H:%M:%S.%f')
        col = raw_data.columns.tolist()
        assert len(col) == 2, 'Number of columns do not fit the requirement.'
        raw_data.rename(columns={col[1]: 'Sensor signal'}, inplace=True)
        raw_data['Elapsed time'] = (raw_data['Time stamp'] - raw_data.at[0, 'Time stamp']).dt.total_seconds()
        return raw_data

    def plot_raw_data(self): # TODO plot the raw data, might be a helpful feature in the GUI
        pass

    def plot_best_fit(self): # TODO plot best fit, against the raw data, might be a helpful feature in the GUI
        pass

    def taylor(self, l: float) -> float:
        """
        Uses the Taylor method to calculate the oxygen diffusivity.

        slope = -Dp * A / (l * V) -> Dp = -slope * l * V / A

        where:
            slope: best fitted slope for the Ln(Cr)
            Dp: oxygen diffusivity
            A: exposed area
            l: length of sample
            V: Volume of the oxygen chamber

        :input:
        l: length of sample in meters

        :return:
        The oxygen diffusivity [m2/s]
        """
        self.Print.header('Taylor method', symbol=' ')
        ODA_chamber_height: float = 0.103  # Height in meters of oxygen gas chamber.
        ODA_chamber_diameter: float = 0.0532  # Diameter in meters of oxygen gas chamber.
        exposed_surface: float = 0.0502 ** 2 / 4 * np.pi  # Exposed surface of the sample in square meters.
        ODA_chamber_volume = ODA_chamber_height * ODA_chamber_diameter ** 2 / 4 * np.pi  # Volume of the oxygen gas chamber.
        dp = -self.slope * ODA_chamber_volume * l / exposed_surface
        self.Print(f'Dp={dp:.8f} m2/s')
        self.Print.line('-')
        return dp

    def currie(self, e: float, l: float) -> float:
        """
        Uses the Currie method to calculate the oxygen diffusivity.

        slope = -Dp * a1^2 / e -> Dp = -slope * e / a1^2

        where:
            slope: best fitted slope for the Ln(Cr)
            Dp: oxygen diffusivity
            a1: positive roots for the general equation
            e: the air filled porosity of the sample

        :input:
        e: air filled porosity
        l: length of sample in meters

        :return:
        The oxygen diffusivity [m2/s]
        """
        self.Print.header('Currie method', symbol=' ')
        a = 0.103
        h = e/a

        f = lambda x: x * np.tan(x) - (h * l) # Function
        a1l, r = brentq(f, a=0, b=np.pi/2, maxiter=200, full_output=True) # Find the first positive root
        a12 = (a1l / l)**2
        dp = -self.slope * e / a12
        self.Print(f'Dp={dp:.8f} m2/s')
        self.Print.line('-')
        return dp


### GUI
class Inputs(tk.Frame):
    filepath: str = None
    oda_obj: None | OxygenDiffusion = None
    def __init__(self, root: tk.Frame, master):
        h = root.winfo_height()
        w = root.winfo_width()
        super().__init__(root, bg='lightblue', bd=3, relief='raised', height=int(h/2)-10, width=w-10)
        self.root = root
        self.master = master

        self.L_inputs = tk.Label(self, text='Input manager', bg='lightblue', font=14)
        self.L_inputs.grid(row=0, column=0)

        self.L_status = tk.Label(self, text='Select file', bg='white', bd=3, height=3, width=35)
        self.L_status.grid(row=1, column=0, pady=10)

        self.B_filepath = tk.Button(self, text='Select data file', command=self.select_file)
        self.B_filepath.grid(row=2, column=0, pady=5)

        self.L_method = tk.Label(self, text='Select calculation method', bg='lightblue')
        self.L_method.grid(row=3, column=0, pady=20)
        methods = ['Taylor', 'Currie']
        self.CB_methods = ttk.Combobox(self, values=methods)
        self.CB_methods.grid(row=4, column=0)
        self.CB_methods.bind('<<ComboboxSelected>>', self.show_inputs)

    def select_file(self):
        self.filepath: str = filedialog.askopenfilename(title='Select ODA data file')
        self.L_status.config(text='File selected :)\nPlease wait')
        self.update()
        self.oda_obj = OxygenDiffusion(self.filepath)
        self.root.oda_obj = self.oda_obj
        self.L_status.config(text='File loaded.\nSelect calculation method')

    def show_inputs(self, event):
        method = self.CB_methods.get()
        print(method)
        self.L_status.config(text=f'The {method} method is chosen.\nPlease insert inputs')
        if method == 'Taylor':
            self.master.Properties.taylor()
        elif method == 'Currie':
            self.master.Properties.currie()


# noinspection PyUnresolvedReferences
class Properties(tk.Frame):
    taylor_shown: bool = False
    currie_shown: bool = False
    def __init__(self, root: tk.Frame, master: tk.Tk):
        h = root.winfo_height()
        w = root.winfo_width()
        super().__init__(root, bg='lightblue', bd=3, relief='raised', height=int(h/2)-10, width=w-10)
        self.root = root
        self.master = master

        self.L_text = tk.Label(self, text='Calculation input', bg='lightblue', font=14)
        self.L_text.grid(row=0, column=0)

        self.B_run = tk.Button(self, text='Run calculations!', command=self.run_calc)
        self.B_run.grid(row=5, column=0)

    def taylor(self):
        if self.currie_shown:
            self.L_currie1.grid_forget()
            self.E_currie1.grid_forget()
            self.L_currie2.grid_forget()
            self.E_currie2.grid_forget()
            self.update()
            self.currie_shown = False
        self.L_taylor = tk.Label(self, text='Length of sample in meters')
        self.L_taylor.grid(row=1, column=0)
        self.E_taylor = tk.Entry(self)
        self.E_taylor.grid(row=2, column=0)
        self.taylor_shown = True

    def currie(self):
        if self.taylor_shown:
            self.L_taylor.grid_forget()
            self.E_taylor.grid_forget()
            self.update()
            self.taylor_shown = False
        self.L_currie1 = tk.Label(self, text='Length of sample in meters')
        self.L_currie1.grid(row=1, column=0)
        self.E_currie1 = tk.Entry(self)
        self.E_currie1.grid(row=2, column=0)
        self.L_currie2 = tk.Label(self, text='Air filled porosity of sample')
        self.L_currie2.grid(row=3, column=0, pady=10)
        self.E_currie2 = tk.Entry(self)
        self.E_currie2.grid(row=4, column=0)
        self.currie_shown = True

    def run_calc(self):
        dp: float | None = None
        if self.taylor_shown:
            l = float(self.E_taylor.get())
            dp = self.root.oda_obj.taylor(l)
        elif self.currie_shown:
            l = float(self.E_currie1.get())
            e = float(self.E_currie2.get())
            dp = self.root.oda_obj.currie(e, l)
        self.root.dp = dp

class Outputs(tk.Frame):
    def __init__(self, root: tk.Tk):
        h = root.winfo_height()
        w = root.winfo_width()
        super().__init__(root, bg='lightblue', bd=3, relief='raised', height=h-10, width=int(w/2)-10)
        self.root = root

class App(tk.Tk):
    title_name: str = 'ODA data treatment' # Name of the app
    width: int = 800
    height: int = 700

    def __init__(self):
        super().__init__()

        self.title(self.title_name)
        size = f'{self.width}x{self.height}'
        self.geometry(size)

        self.Frame1 = tk.Frame(self, height=self.height-10, width=int(self.width/2)-10)
        self.Frame1.grid(row=0, column=0)
        self.update()
        self.Inputs: Inputs = Inputs(self.Frame1, self)
        self.Inputs.grid(row=0, column=0)
        self.update()
        self.Properties = Properties(self.Frame1, self)
        self.Properties.grid(row=1, column=0)
        self.update()

        self.Outputs: Outputs = Outputs(self)
        self.Outputs.grid(row=0, column=1, padx=10)
        self.update()

        # TODO complete input frame and work on output frame, could be a interactive 2x1 plot, with the Dp in a labelbox below.
    def show_prop(self, method: str) -> None:
        if method == 'Taylor':
            self.Properties.taylor()
        elif method == 'Currie':
            self.Properties.currie()



if __name__ == '__main__':
    # filepath: str = r"C:\Users\AC03LH\Aalborg Universitet\Biobased building materials workshop - Documents\General\Student projects\2024 IE10 - Victor, Asbjørn and Lotte\Data\Experiment\Wood_wool\Wood_wool_test_cell_1_sample_3"

    # obj = OxygenDiffusion(filepath=filepath)
    # obj.taylor(l=0.051)
    # obj.currie(e=0.97, l=0.051)

    app = App()
    app.mainloop()