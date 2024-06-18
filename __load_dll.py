import os

STATUS = "Develop"
# STATUS = "Release"

# Load the shared library
current_dir = os.path.dirname(os.path.abspath(__file__))
if STATUS == "Develop":
    if os.name == "nt":
        # dll_path = os.path.join(current_dir, "build/bin/Debug/rsf.dll")
        dll_path = os.path.join(current_dir, "build/bin/Release/rsf.dll")
    else:
        # dll_path = os.path.join(current_dir, "build/lib/Debug/librsf.so")
        dll_path = os.path.join(current_dir, "build/lib/Release/librsf.so")
elif STATUS == "Release":
    if os.name == "nt":
        dll_path = os.path.join(current_dir, "bin/rsf.dll")
    else:
        dll_path = os.path.join(current_dir, "lib/librsf.so")
else:
    raise AssertionError('''The variable named STATUS in "__load_dll.py" is not correctly set.''')