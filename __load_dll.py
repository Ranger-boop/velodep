import os
import ctypes

STATUS = "Develop"
# STATUS = "Release"

# Load the shared library
current_dir = os.path.dirname(os.path.abspath(__file__))
if STATUS == "Develop":
    if os.name == "nt":
        # dll_path = os.path.join(current_dir, "build/bin/Debug/velodep.dll")
        dll_path = os.path.join(current_dir, "build/bin/Release/velodep.dll")
    else:
        # dll_path = os.path.join(current_dir, "build/lib/Debug/libvelodep.so")
        dll_path = os.path.join(current_dir, "build/lib/Release/libvelodep.so")
elif STATUS == "Release":
    if os.name == "nt":
        dll_path = os.path.join(current_dir, "bin/velodep.dll")
    else:
        dll_path = os.path.join(current_dir, "lib/libvelodep.so")
else:
    raise AssertionError('''The STATUS variable in "__load_dll.py" is not correctly set.''')
dll = ctypes.cdll.LoadLibrary(dll_path)