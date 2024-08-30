# velodep
A python, C, and C++ package for simulating the slip behavior and aperture/permeability evolution of a fracture using the rate-and-state friction (RSF) law and the displacement- and velocity-dependent aperture (DVA) model.

## Usage

In the ***Example*** folder, there is an example source code file and an example data.

The file named "example.ipynb" is the example source code file, which includes the source code of the inversion of consitutive parameters of both the rate-and-state friction law and the displacement- and velocity-dependent aperture model, as well as the source code of the calculation of slip behavior and permeability evolution based on these two models. The excel file named "example_monotonic_pressurization.xlsx" is the example data needed by the source code.

To run the example source code, please do the following steps:

1. After download this package, move the ***Example*** folder to its parent folder. Now you are in the parent folder of ***velodep***.
  
2. Consecutively run the following commands to compile the package based on cmake:
     
    - On windows platform using VisutalStudio compiler:
     
        ```powershell
        cd velodep
        mkdir build
        cd buid
        cmake -A x64 ..
        cmake --build . --config Release
        ```
     
    - On Linux platform using gcc compiler:
     
        ```bash
        cd velodep
        mkdir build
        cd buid
        cmake -A x64 -DCMAKE_BUILD_TYPE=Release ..
        cmake --build .
        ```

3. After the package has been compiled, open the "example.ipynb" file and run the code in the file.

For a published work using this package, please refer to https://doi.org/10.6084/m9.figshare.26871748.v1.