import os
import subprocess

dirpath = os.path.dirname(os.path.abspath(__file__))
os.chdir(dirpath)

for i in range(0,50):
    print(i , "/50")
    subprocess.check_call(['generateknot.exe'])