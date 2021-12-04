# SimpleMPM3
2D and 3D material points method for:
1. Large deformation
2. Water - soil interaction
3. Soil - rigid structure interaction

# External dependencies
Eigen, Hdf5, Qt5 (5.15.0), Intel Threading Building Block (TBB), OpenGL

Eigen libraries are all the header files, and put into Vendors folder

Hdf5, Qt5, TBB comes in vcpkg:
• tbb:x64-windows
• hdf5:x64-windows
• freetype:x64-windows

Packages needed in Python:
• numpy
• matplotlib
• h5py

Graphic cards and associating drivers are needed for OpenGL

System variables are needed:
• VCPKG_ROOT for vcpkg install folder .../vcpkg
• QT5_ROOT as qt installation folder .../Qt
• VENDORS_ROOT as .../Vendors

Add ...\Qt\5.15.0\msvc2019_64\bin into path or run "windeployqt.exe" in the program folder for Qt libraries  
