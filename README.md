# FasterCap

## Copyright

Copyright (c) 2019
FastFieldSolvers S.R.L. - http://www.fastfieldsolvers.com

## Description

FasterCap is a powerful three- and two-dimensional capactiance extraction program. 

For pre-compiled binaries, support, consultancy and additional information please visit http://www.fastfieldsolvers.com/
Access to the download pages is free, and you may access anonymously if you want.

The source code is released under LGPL 2.1 or any later LGPL version, see the file LICENSE.txt for the details.

## Compiling

FasterCap can be compiled for MS Windows and for *nix. For multiple platform compilation support, FasterCap uses CMake, plus Code::Blocks for the GUI version.
The following software tools and complier chains are required. Used versions are indicated, higher version may work, but were not tested.

###  MS Windows 64 bits

- CMake 3.6.0
- Code::Blocks, version 13.12
- TDM-GCC 64 bits, version 4.8.1
- wxWidgets, version 3.0

#### Notes

- Use a fresh TDM-GCC package, not the one shipped with Code::Blocks, as the latter is missing the TDM OpenMP package.
      
- You need to pre-compile wxWidgets with TDM-GCC. Do not use pre-compiled versions, as the used compiler switch configurations may be very different when generating the binaries. Do NOT compile as monolithic, so .exe are smaller. Do two compiles, debug + release, using MSDOS makefiles and TDM-GCC (no need to use MSYS).

- Run CMake until you configure all the required parameters (you need some knowledge of how CMake works). You may want to set the switch FASTFIELDSOLVERS_HEADLESS as ON for a DOS-only version ("headless") of the software, however under Windows this is useless; just run FasterCap in shell-only mode with the appropriate -b switch (see FasterCap documentation)

###  Linux 64 bits

- CMake 2.8.12
- Code::Blocks, version 13.12
- GCC, version 4.8.1
- wxWidgets (wxGTK), version 3.0

#### Notes

- Run CMake until you configure all the required parameters (you need some knowledge of how CMake works), then generate a CodeBlock project. You may want to set the switch FASTFIELDSOLVERS_HEADLESS as ON for a shell-only version ("headless") of the software, or you can keep this OFF and just run FasterCap in shell-only mode with the appropriate -b switch (see FasterCap documentation)

Example (using a GUI-less CMake):

Move to the build directory you choose (different from the source files directory) and type the following command, where "../FasterCap" is the path to the base directory containing the FasterCap source code:
    
`cmake -G"CodeBlocks - Unix Makefiles" -DCMAKE_BUILD_TYPE=Release ../FasterCap`

Then you can open Code::Blocks, open the FasterCap.cbp project file created by CMake under the build diretory, and compile it.

###  Linux 64 bits headless

- CMake 3.5.1
- GCC, version 4.8.1 or higher (higher version tested: on 5.4.0)
- wxWidgets, version 3.0.2 or higher (higher version tested: 3.0.2)

#### Notes

- Run CMake until you configure all the required parameters (you need some knowledge of how CMake works). You **must** set the switch FASTFIELDSOLVERS_HEADLESS as ON, as in a headless Linux distro you have no drivers at all for video.

Example: 

Move to the build directory you choose (different from the source files directory) and type the following command, where "../FasterCap" is the path to the base directory containing the FasterCap source code:
    
`cmake -G"Unix Makefiles"  -DCMAKE_BUILD_TYPE=Release -DFASTFIELDSOLVERS_HEADLESS=ON ../FasterCap`

Then you can launch the build process with:
    
`make`
  
Remark: at run time you will see the an "Assert failure" message. This is a wxWidgets issue, see ["debug message when running without session manager"](http://trac.wxwidgets.org/ticket/16024).
The assert is harmless, and it is fixed starting from wxWidgets 3.1.1

## Additional packages

FasterCap also requires two additional source code packages, that are available through the same official repositories you can access from http://www.fastfieldsolvers.com/ or directly from GitHub.

The packages are:

- LinAlgebra
- Geometry

Both package directories, with the above names, must be at the same hierarchy level in the folder structure of the FasterCap source code directory, and are handled by CMake.

## Additional information

For any additional information please visit [FastFieldSolvers](https://www.fastfieldsolvers.com/), write on the [FastFieldSolvers Forum](https://www.fastfieldsolvers.com/forum) or [contact us](https://www.fastfieldsolvers.com/contact.htm).

