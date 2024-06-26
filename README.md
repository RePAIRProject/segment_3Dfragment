﻿# segment_3Dfragment
This repository is designed to segment a fresco fragment into meaningful surfaces, including the painted (intact) surface, the opposite surface, and the "sidewalls" surfaces.

# 1) Description
This project employs a region-growing algorithm based on the curvedness values of the vertices of the mesh, leveraging differential geometry principles. The algorithm processes a 3D mesh or point cloud of a fresco fragment to identify and segment these distinct surfaces.
![segment 3d image](https://github.com/RePAIRProject/segment_3Dfragment/assets/38216201/2fcb68ab-6570-42ab-8162-8fc619753285)


# 2) Installation
This project was developed on a Windows machine using CMake and Visual Studio compiler. The project heavily dependent on [libigl library](https://github.com/libigl/libigl).

1. Please clone the project 
```bash
git clone https://github.com/RePAIRProject/segment_3Dfragment.git
cd segment_3Dfragment
```
2. Compile and build using cmake
```bash
mkdir build
cd build
cmake ..
cmake --build .  # or 'make' if you're using Makefiles as the generator
sudo cmake --install .  # Use sudo if you need superuser permissions to install
```


## 3) Usage
To run the segmentation algorithm to extract the intact surface, use the following command:
```bash
./segment_3Dfragment --input-File <PATH>
```
Where *`--input-File`* specifies the path to the OBJ file of the fragment to be segemented. 
The output file would be save in the same directory of the inputted file, using its name and the postfix *intact*.

Other optional parameters are the following:
- *`--script`*: and then specify `opposite-surface` or `sidewalls-surface`
- *`--disable-save`*: if used, then the segmented mesh won't be outputted, and in particular, a new obj file would not be saved.
- *`--enable-debug`*: if used, a visualization of the fragemnts and the segmentation process would popped out to the screen. Press 1 for detecting the results of the region-growing process, press 2 for the outputted segment (the only one with color, the rest would be white), and press 3 for viewing the sufficiently large enough segments.
- *`--intact-Normal-Std-Thershold`*: Choose this value in the range (0,0.3) (you can exceed 0.3 in extreme cases). It indicates the maximum std of the normals of the desired surface. Higher value means that the surface is accepped as less smooth.
- *`--intact-Similarity-Fracture`*: Choose this value roughly in the range (0.4,0.85). Lower values would lead to faster spread of a segment in the region growing process. 


# 4) Known Issues
Bug descriptions.

# 5) Relevant publications
Some publications.

