import os
import ctypes

STATUS = "Develop"
# STATUS = "Release"

# Load the shared library
current_dir = os.path.dirname(os.path.abspath(__file__))
if STATUS == "Develop":
    if os.name == "nt":
        # dll_path = os.path.join(current_dir, "build/Debug/bin/velodep.dll")
        dll_path = os.path.join(current_dir, "build/Release/bin/velodep.dll")
    else:
        # dll_path = os.path.join(current_dir, "build/Debug/lib/libvelodep.so")
        dll_path = os.path.join(current_dir, "build/Release/lib/libvelodep.so")
elif STATUS == "Release":
    if os.name == "nt":
        dll_path = os.path.join(current_dir, "bin/velodep.dll")
    else:
        dll_path = os.path.join(current_dir, "lib/libvelodep.so")
else:
    raise AssertionError('''The STATUS variable in "__load_dll.py" is not correctly set.''')
dll = ctypes.cdll.LoadLibrary(dll_path)