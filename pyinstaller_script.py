# Script used to create executable binary file, using the Pyinstaller package (6.17.0)
# Tip: hook-tkinterdnd2.py and hook-scipy.py must be in a folder named 'hook'.

import PyInstaller.__main__

PyInstaller.__main__.run(
    [
        "main.py",
        "--onefile",
        "--windowed",
        "--additional-hooks-dir=hook",
        "--icon=SpectroIBIS_icon.ico",
    ]
)
