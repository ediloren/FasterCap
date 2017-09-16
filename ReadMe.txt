==================================================================================

  FasterCap

  Copyright (c) 2017  FastFieldSolvers S.R.L. -  http://www.fastfieldsolvers.com

==================================================================================

FasterCap is a powerful three- and two-dimensional capactiance extraction program. 

For pre-compiled binaries, support, consultancy and additional information please visit http://www.fastfieldsolvers.com/
Access to the download pages is free, and you may access anonymously if you want.

The source code is released under LGPL 2.1 or any later LGPL version, see the file LICENSE.txt for the details.

---------------------------
How-to-compile instructions
---------------------------

FasterCap can be compiled for MS Windows and for *nix. To support multiple platform compilation,
FasterCap uses a multi-target Code::Blocks project.
The following software tools and compilator chains are required. Used versions are indicated, higher
version may work, but were not tested. 

** MS Windows 64 bits:

  - Code::Blocks, version 13.12
  - TDM-GCC 64 bits, version 4.8.1
  - wxWidgets, version 3.0

  Note that in Code::Blocks: 

    - you need to proper configure the compiler chain to point to TDM-GCC.
      Follow the Code::Blocks help to understand how to do it.
      Use a fresh TDM-GCC package, not the one shipped with Code::Blocks,
      as the latter is missing the TDM OpenMP package.

    - you need to define the global variable 'wx', pointing to the wxWidgets installation folder.
      This is done under 'Settings'->'Global variable editor'.
      Important remark: no spaces in the wxWidgets path, otherwise gcc called from Code::Block
      will fail (starts the path at the first space) 

  You also need to compile wxWidgets with TDM-GCC. Do not use pre-compiled versions, as the 
  used compiler switch configurations may be very different when generating the binaries.
  Do NOT compile as monolithic, so .exe are smaller. Do two compiles, debug + release,
  using MSDOS makefiles and TDM-GCC (no need to use MSYS).
 

** MS Windows 32 bits:

  Not maintained any more


** Linux 64 bits:

  - Code::Blocks, version 13.12
  - GCC, version 4.8.1
  - wxWidgets (wxGTK), version 3.0

  Note that in Code::Blocks:
  
    you should set up include and lib paths for wxWidgets. However looking in how
    wxSmith does it, we don't need to define an environment variable #wx and include library directories
    etc, because there is a shortcut: under Project->Properties, button 'Project's build options..', we configure for
    'release' and 'debug' the proper compiler AND linker setting using for both the 'Other options'
    settings (TWO different pages). These allow to call an external program, in this case 'wx-config'. 
    This upon 'make install' of wxGTK is installed under /usr/local/bin (so wx-config is in the standard path).
    Run wx-config from a shell, with no parameters, and you'll see it returns the commands and options needed
    to use wxWidgets.
    Warning: the supported configurations are related to the builds you did. So if you built with no static support,
    if you select --static=yes you get an error (and the compiler will receive as output of wx-config an error message,
    instead of the gcc options to compile the sources, that include library and include paths)
    Therefore, if you need to compile linking different versions of wxWidgets, e.g. release and debug,
    you must make and install the two different versions. The wx-conf command is updated with the configs
    of the builds. Then if you have different versions, you should pass the relative parameters to wx-config
    to choose the desired one. For instance, for unicode you pass wx-conf --unicode, for ANSI you pass
    wx-config --unicode=no, for static vs. dynamic linking --static etc. Select wx-config --help for a list
    of the options. Note that wx-config has a default selection, that cannot be changed with the options.
    This is because wx-config is only a symbolic link to the corresponding shell script. See where wx-config is
    ('which wx-config') go there and see ('ls -l wx-config'). Then you can 'rm wx-config' (from root) and
    create another link ('ln -s <<mylink>> wx-config').
    When installing wxWidgets3.0, wx-config options change, for instance you must not specify --unicode=no,
    since this option is not supported / not meaningful and returns an error.
  
  
** Linux 32 bits:

  Not maintained any more


FasterCap also requires two additional source code packages, that are available through the same official repositories
you can access from http://www.fastfieldsolvers.com/

The packages are:

- LinAlgebra
- Geometry

Both package directories, with the above names, must be at the same hierarchy level in the folder structure 
of the FasterCap source code directory.

