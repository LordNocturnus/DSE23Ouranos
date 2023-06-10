import numpy as np
import pandas as pd
import os
import platform
import subprocess
import re

# Check what OS we are running in to select correct version
if platform.system() == "Windows":
    gram_path = "GRAM_Suite_1.5/Windows/UranusGRAM.exe"
elif platform.system() == "Linux":
    gram_path = "GRAM_Suite_1.5/Linux/UranusGRAM.exe"
else:
    raise OSError("Please ues a proper OS! Gram only supports Windows and Linux")
if not os.path.exists(__file__[:-7]+"GRAM_Suite_1_5"):
    raise ImportError("Gram not found please contact Felix for the files")

_path = os.path.dirname(__file__)


class GRAM(object):

    def __init__(self):
        self.gram_path = gram_path
        self.spice_path = _path+"/GRAM_Suite_1_5/SPICE"
        self.list_path = _path+"/atmos_LIST"
        self.col_path = _path+"/atmos_OUTPUT"
        self.month = 3
        self.day = 25
        self.year = 2020
        self.hour = 12
        self.minute = 30
        self.seconds = 0.0
        self.seed = 1001
        self.densityperturbation = 1.0
        self.stepsize = 0.0
        self.trajectorypath = _path+"/traj_data.txt"
        self.runs = 1
        self.altitudes = np.arange(7000, -291, -1, dtype=float)
        self.lat = np.zeros_like(self.altitudes)
        self.long = np.zeros_like(self.altitudes)
        self.time = np.arange(0, len(self.altitudes), 1, dtype=float)
        self.data = None

    def compile_config(self):
        """ compiles the necessary data for GRAM into a configuration file that GRAM can read"""
        txt = " $INPUT\n"
        txt += f"  SpicePath             = '{self.spice_path}'\n"
        txt += f"  ListFileName          = '{self.list_path}'\n"
        txt += f"  ColumnFileName        = '{self.col_path}'\n"

        txt += "  EastLongitudePositive = 1\n\n"

        txt += f"  TimeFrame = 1\n"
        txt += f"  TimeScale = 1\n"
        if 1 <= self.month <= 12 and type(self.month) == int:
            txt += f"  Month = {self.month}\n"
        else:
            raise ValueError("month needs to be integer between 1 and 12")
        if type(self.day) == int:
            txt += f"  Day = {self.day}\n"
        else:
            raise ValueError("day needs to be and integer")
        if 1970 <= self.year <= 2069 and type(self.year) == int:
            txt += f"  Year = {self.year}\n"
        else:
            raise ValueError("year needs to be integer between 1970 and 2069")
        if 0 <= self.hour <= 23 and type(self.hour) == int:
            txt += f"  Hour = {self.hour}\n"
        else:
            raise ValueError("hour needs to be integer between 0 and 23")
        if 0 <= self.minute <= 59 and type(self.minute) == int:
            txt += f"  Minute = {self.minute}\n"
        else:
            raise ValueError("minute needs to be integer between 0 and 59")
        if 0.0 <= self.seconds < 60.0:
            txt += f"  Seconds = {self.seconds}\n\n"
        else:
            raise ValueError("seconds needs to be between 0.0 and 60.0")

        if 0 <= self.seed <= 29999 and type(self.seed) == int:
            txt += f"  InitialRandomSeed = {self.seed}\n"
        else:
            raise ValueError("initial random seed needs to be integer between 0 and 29999")

        if 0.0 <= self.densityperturbation <= 2.0:
            txt += f"  DensityPerturbationScale = {self.densityperturbation}\n"
        else:
            raise ValueError("Density perturbation needs to be between 0.0 and 2.0 (3 sigma ~= 1.0")
        txt += f"  MinimumRelativeStepSize  = 0.0\n"

        txt += f"  TrajectoryFileName = '{_path}/traj_data.txt'\n"

        if 1 <= self.runs and type(self.runs) == int:
            txt += f"  NumberOfMonteCarloRuns = {self.runs}\n"

        txt += f"  FastModeOn        = 0\n"
        txt += f"  ExtraPrecision    = 0\n"

        txt += f" $END"
        with open(_path+"/gram_config.txt", "w") as file:
            file.write(txt)

    def compile_trajectory(self):
        """compiles the trajectory data into the required format"""
        if not len(self.altitudes) == len(self.lat) or not len(self.altitudes) == len(self.long) or not \
             len(self.altitudes) == len(self.time):
            raise IndexError("length of the altitude latitude longitude and timestamps for the GRAM trajectory should be the same")
        if any(self.time < 0):
            raise ValueError("time should be positive")
        if any(self.altitudes > 7000) or any(self.altitudes < -290):
            raise ValueError("GRAM altitude should be between 7000 and -290")
        if any(np.abs(self.lat) > 90.0):
            raise ValueError("GRAM latitude should be +- 90° with negative being south")
        if any(self.long < 0.0) or any(self.long > 360.0):
            raise ValueError("GRAM longitude should be between 0° and 360°")

        with open(_path+"/traj_data.txt", "w") as file:
            for k, _ in enumerate(self.altitudes):
                file.write(f"{self.time[k]} {self.altitudes[k]} {self.lat[k]} {self.long[k]}\n")

    def run(self):
        self.compile_config()
        self.compile_trajectory()
        out = subprocess.check_output([_path+"/GRAM_Suite_1_5/Windows/UranusGRAM.exe", "-file",
                                       _path+"/gram_config.txt"])

        if not re.search("Files output: ", str(out)):
            print(str(out))
            raise RuntimeError("GRAM failed to run please check input files")
        self.read_data()

    def read_data(self):
        """
        loads the GRAM generated data into a pandas dataframe and filters out useless data look at GRAM documentation in
        /GRAM_Suit_1_5/Documentation/Uranus-Gram User Guide.pdf for full ist of parameters"""

        self.data = pd.read_csv(_path+"/atmos_OUTPUT.csv")
        #
        #self.data.drop(columns=["AverageMolecularWeight", "LocalSolarTime_hr",
         #                      "NSWindPerturbation_ms", "PerturbedEWWind_ms", "PerturbedNSWind_ms",
          #                     "EWStandardDeviation_ms", "NSStandardDeviation_ms", "EWWind_ms", "NSWind_ms",
           #                    "EWWindPerturbation_ms", "LongitudeOfTheSun", "SubsolarLatitude_deg",
            #                   "SubsolarLongitudeE_deg", "SolarZenithAngle_deg", "OneWayLightTime_min",
             #                  "OrbitalRadius_AU", "SecondsPerSol", "PressureAtSurface_Pa", "TotalNumberDensity_m3",
              #                 "H2nd_m3", "H2mass_pct", "H2mole_pct", "H2amw",
               #                "Hend_m3", "Hemass_pct", "Hemole_pct", "Heamw",
                #               "CH4nd_m3", "CH4mass_pct", "CH4mole_pct", "CH4amw"], inplace=True)
        self.data.drop(columns=self.data.columns[-1], inplace=True)
        #print(self.data.head())
        for c in self.data.columns:
            #print(c)
            if self.data[c].isnull().any():
                raise ValueError(f"GRAM returned NaN values in {c} please investigate")


if __name__ == "__main__":
    test = GRAM()
    test.run()
