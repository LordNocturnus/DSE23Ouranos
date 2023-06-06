import os
import platform

# Check what OS we are running in to select correct version
if platform.system() == "Windows":
    gram_path = "GRAM_Suite_1.5/Windows/UranusGRAM.exe"
elif platform.system() == "Linux":
    gram_path = "GRAM_Suite_1.5/Linux/UranusGRAM.exe"
else:
    raise OSError("Please ues a proper OS! Gram only supports Windows and Linux")


if __name__ == "__main__":
    print("Hello World")