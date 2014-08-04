# Getting Started with Crux

In this short tutorial series you will learn the basics of Crux. The tutorial is presented as a demo project, in which we analyze mass spectra data from a [SAMPLE DESCRIPTION]. Each section in the tutorial steps through the process of utilizing a specific Crux command. A normal project workflow encompasses chaining together multiple Crux commands. We will show you which commands can be chained together and in what order. This tutorial is by no means comprehensive, but will give you a firmer grasp on the Crux workflow for tandem mass spectrometry analysis.

### Prerequisites

* Command Line Basics

### Table of Contents

1. [Installing Crux](#1)
* [Setting up the demo project](#2)
* [Understanding common conventions in Crux](#3)
* [Assigning high-resolution precursor masses using Bullseye](#4)
* [Assigning peptides to spectra using Comet](#5)
* [Assigning peptides to spectra using Tide](#6)
* [Post-processing with Percolator](#7)
* [Post-processing with Barista](#8)
* [Quantification using spectral-counts](#9)


## <a href="#1" name="1">#1</a> Installing Crux

Before we really get started with Crux, we are going to need to install it first. You have the choice of installing precompiled binaries or building Crux from the source code. The compiled binaries offer the easiest installation option.

#### Linux
* [Installing the Pre-compiled Binary](#linux-binary)
* [Installing from source](#linux-source)

#### Mac OS X
* [Installing the Pre-compiled Binary](#mac-os-x-binary)
* [Installing from source](#mac-os-x-source)

#### Windows
* [Installing the Pre-compiled Binary](#windows-binary)
* [Installing from source](#windows-source)

> Note:  
> All Crux downloads, including previous versions of Crux, can be found at [DOWNLOAD URL].

### <a href="#linux-binary" name="linux-binary">#</a> Linux: Installing the Pre-compiled Binary

####Download the files
Crux is freely available for download under an Apache 2.0 license. On the download page click on the link "Linux [architecture] Build (zipped)" to download them to your computer. If your computer is running a 32-bit version of Linux select the i686 version. If your computer is running a 64-bit version of Linux select the x86_64 version. If a window pops up, select the "Save" button. To check that the file was transfered correctly, run this command

	md5sum [downloaded file]
and compare the output to the "MD5 checksum" link next to the file you downloaded.

####Unpack and install the files
The file you downloaded should be named something like crux-x.xx-XXXX-GNU-Linux.zip where x.xx is replaced with the most current release version. and XXXX is either i686 or x86_64. Move this file to a convenient location and unpack the file with this command.

	unzip crux-x.xx-XXXX-GNU-Linux.zip
You should now have a new directory, crux-x.xx-XXXX-GNU-Linux, with two subdirectories bin, and doc. The crux programs are in the bin directory. You may want to move these files somewhere more convenient for you, but otherwise, that completes the installation! 

####Set the PATH
On Linux and OS X systems you can always run any of the programs by specifying its full path. For example,

	<install directory>/bin/crux search-for-matches
In this final step of the installation, we'll set things up so that you can also run any of the programs by just typing its name. This is done with an environment variable called `$PATH`. `$PATH` is a list of places that the computer looks for executable programs. Find out what is in the current list with this command:

	echo $PATH
The value returned might look something like this: `/usr/local/bin:/usr/bin:/bin:/usr/bin/X11:/usr/games`

The first directory in this colon-separated list, `/usr/local/bin`, is the first place the computer looks for a program. If it doesn't find the program there it tries /usr/bin and so on down the list. You can either move the crux executable files into one of these directories already in your PATH, or you can add <install directory> to PATH with this command. (Remember to replace the word <install directory> with the actual location of your installation)

	export PATH=$PATH:<install directory>/bin
Try `echo $PATH` again to see that the change was made. `<install directory>/bin` should be at the end of the list. Now you can run the crux programs by just typing their names. Try this.

	crux search-for-matches
The change to `$PATH` is only temporary. As soon as you close the window, the change will disappear. In order to make the change every time you log on, add the line `export PATH=$PATH:<install directory>/bin` to the file `.bashrc` located in your home directory. You will have to log on again for the change to take effect.

You're done, continue on to step 2, [Setting up the demo project](#2).

### <a href="#linux-source" name="linux-source">#</a> Linux: Installing from Source

####Prerequisites
To build Crux from source on Linux you must have installed the following software:

* gcc
* g++
* cmake
* make
* svn
* wget (must support SSL)
Building crux also requires the Proteowizard and Percolator libraries. These will automatically be copied from the Proteowizard and Percolator repositories and statically linked into crux.

####Download the files
The Crux source is freely available for download under an Apache 2.0 license>. On the download page click on the link "Source (tarred, gzipped)", and choose save to disk if prompted. To check that the file was transfered correctly, run this command

	md5sum [downloaded file]
and compare the output to the "MD5 checksum" link next to the file you selected on the Crux download page.

The file you downloaded should be named something like `crux-x.xx.source.tar.gz`, where x.xx is replaced with the most current release version.

####Unpack the files
Move the downloaded file to a convenient location and unpack it with this command.

	tar -zxvf crux-x.xx.Source.tar.gz
A new directory named `crux-x.xx.Source` will be created and a long list of files will stream by as they are being put into the new directory.
Configure the build
Move into the `crux-x.xx.Source` directory from step 1. Run the command

	cmake -DCMAKE_INSTALL_PREFIX:PATH=<install directory> .
where <install directory> is the location where you wish to install the crux programs. If you want crux to go in `$HOME/bin`, then `<install directory>` would be `$HOME`. The installation process will automatically put the programs in a directory called bin.
####Build and install
To complete the process, run these two commands.

	make
	make install
The crux programs are in the bin subdirectory of the `<install directory>` That completes the installation!

####Set the PATH
On Linux and OS X systems you can always run any of the programs by specifying its full path. For example,

	<install directory>/bin/crux search-for-matches
In this final step of the installation, we'll set things up so that you can also run any of the programs by just typing its name. This is done with an environment variable called `$PATH`. `$PATH` is a list of places that the computer looks for executable programs. Find out what is in the current list with this command:

	echo $PATH
The value returned might look something like this: `/usr/local/bin:/usr/bin:/bin:/usr/bin/X11:/usr/games`

The first directory in this colon-separated list, `/usr/local/bin`, is the first place the computer looks for a program. If it doesn't find the program there it tries /usr/bin and so on down the list. You can either move the crux executable files into one of these directories already in your PATH, or you can add <install directory> to PATH with this command. (Remember to replace the word <install directory> with the actual location of your installation)

	export PATH=$PATH:<install directory>/bin
Try `echo $PATH` again to see that the change was made. `<install directory>/bin` should be at the end of the list. Now you can run the crux programs by just typing their names. Try this.

	crux search-for-matches
The change to `$PATH` is only temporary. As soon as you close the window, the change will disappear. In order to make the change every time you log on, add the line `export PATH=$PATH:<install directory>/bin` to the file `.bashrc` located in your home directory. You will have to log on again for the change to take effect.

You're done, continue on to step 2, [Setting up the demo project](#2).

### <a href="#mac-os-x-binary" name="mac-os-x-binary">#</a> Mac OS X: Installing the Pre-compiled Binary
####Download the files
Crux is freely available for download under an Apache 2.0 license. On the download page click on the link "OS X x86_64 Build (zipped)" to download them to your computer. If a window pops up, select the "Save" button. To check that the file was transfered correctly, run this command

	md5 [downloaded file]
and compare the output to the "MD5 checksum" link next to the file you downloaded.

####Unpack and install the files
The file you downloaded should be named something like `crux-x.xx.Darwin.x86_64.zip` where x.xx is replaced with the most current release version. Move this file to a convenient location and unpack the file with this command.

	unzip crux-x.xx.Darwin.x86_64.zip
You should now have a new directory, `crux-x.xx.Darwin.x86_64`, with two subdirectories bin, and doc. The crux programs are in the bin directory. You may want to move these files somewhere more convenient for you, but otherwise, that completes the installation! 

####Set the PATH
On Linux and OS X systems you can always run any of the programs by specifying its full path. For example,

	<install directory>/bin/crux search-for-matches
In this final step of the installation, we'll set things up so that you can also run any of the programs by just typing its name. This is done with an environment variable called `$PATH`. `$PATH` is a list of places that the computer looks for executable programs. Find out what is in the current list with this command:

	echo $PATH
The value returned might look something like this: `/usr/local/bin:/usr/bin:/bin:/usr/bin/X11:/usr/games`

The first directory in this colon-separated list, `/usr/local/bin`, is the first place the computer looks for a program. If it doesn't find the program there it tries /usr/bin and so on down the list. You can either move the crux executable files into one of these directories already in your PATH, or you can add <install directory> to PATH with this command. (Remember to replace the word <install directory> with the actual location of your installation)

	export PATH=$PATH:<install directory>/bin
Try `echo $PATH` again to see that the change was made. `<install directory>/bin` should be at the end of the list. Now you can run the crux programs by just typing their names. Try this.

	crux search-for-matches
The change to `$PATH` is only temporary. As soon as you close the window, the change will disappear. In order to make the change every time you log on, add the line `export PATH=$PATH:<install directory>/bin` to the file `.bashrc` located in your home directory. You will have to log on again for the change to take effect.

You're done, continue on to step 2, [Setting up the demo project](#2).

### <a href="#mac-os-x-source" name="mac-source">#</a> Mac OS X: Installing from Source

####Prerequisites
To build Crux from source on OS X you must have installed the following software:

* gcc
* g++
* cmake
* make
* svn
* wget (must support SSL)

Building crux also requires the Proteowizard and Percolator libraries. These will automatically be copied from the Proteowizard and Percolator repositories and statically linked into crux.

The C++ compiler that comes with Apples's XCode is partially incompatible with the GNU C++ compilers and can't compile some parts of crux and the supporting libraries. If you are building on OS X you will need to install a copy of the GNU C and C++ compilers from a repository such as MacPorts. We have used version 4.7 of the GNU g++ compiler downloaded from MacPorts. Once you've installed MacPorts the commands to install the GNU C and C++ compilers is

	sudo port install gcc47
MacPorts allows you to install multiple versions of the GNU compilers. You'll want to set your prefered version as the default for MacPorts. Fo example:

	sudo port select --set gcc mp-gcc47
You will also have to configure your environment to use the MacPorts compilers. This is done using the environment variables PATH, which sets the list of directories that are searched for executables, CC, which identifies the C compiler, and CXX, which identifies the C++ compiler. Before staring your build set these environment variables to point to the MacPorts compilers. For example:

	export PATH=/opt/local/bin:$PATH 
	export CC=/opt/local/bin/gcc 
	export CXX=/opt/local/bin/g++

You may want to include these commands in your .profile or .bashrc.

cmake, svn, and wget, can also be installed from MacPorts. Note that you must install the ssl variant of wget:

	sudo port install wget +ssl
	
####Download the files
The Crux source is freely available for download under an Apache 2.0 license. On the download page click on the link "Source (tarred, gzipped)", and choose save to disk if prompted. To check that the file was transfered correctly, run this command

	md5sum [downloaded file]
	
and compare the output to the "MD5 checksum" link next to the file you selected on the Crux download page.

The file you downloaded should be named something like `crux-x.xx.source.tar.gz`, where x.xx is replaced with the most current release version.

####Unpack the files
Move the downloaded file to a convenient location and unpack it with this command.

	tar -zxvf crux-x.xx.Source.tar.gz
	
A new directory named `crux-x.xx.Source` will be created and a long list of files will stream by as they are being put into the new directory.
Configure the build
Move into the `crux-x.xx.Source` directory from step 1. Run the command

	cmake -DCMAKE_INSTALL_PREFIX:PATH=<install directory> .

where `<install directory>` is the location where you wish to install the crux programs. If you want crux to go in `$HOME/bin`, then `<install directory>` would be `$HOME`. The installation process will automatically put the programs in a directory called bin.
Build and install
To complete the process, run these two commands.

	make
	make install
The crux programs are in the bin subdirectory of the `<install directory>` That completes the installation! 

####Set the PATH
On Linux and OS X systems you can always run any of the programs by specifying its full path. For example,

	<install directory>/bin/crux search-for-matches
In this final step of the installation, we'll set things up so that you can also run any of the programs by just typing its name. This is done with an environment variable called `$PATH`. `$PATH` is a list of places that the computer looks for executable programs. Find out what is in the current list with this command:

	echo $PATH
The value returned might look something like this: `/usr/local/bin:/usr/bin:/bin:/usr/bin/X11:/usr/games`

The first directory in this colon-separated list, `/usr/local/bin`, is the first place the computer looks for a program. If it doesn't find the program there it tries /usr/bin and so on down the list. You can either move the crux executable files into one of these directories already in your PATH, or you can add <install directory> to PATH with this command. (Remember to replace the word <install directory> with the actual location of your installation)

	export PATH=$PATH:<install directory>/bin
Try `echo $PATH` again to see that the change was made. `<install directory>/bin` should be at the end of the list. Now you can run the crux programs by just typing their names. Try this.

	crux search-for-matches
The change to `$PATH` is only temporary. As soon as you close the window, the change will disappear. In order to make the change every time you log on, add the line `export PATH=$PATH:<install directory>/bin` to the file `.bashrc` located in your home directory. You will have to log on again for the change to take effect.

You're done, continue on to step 2, [Setting up the demo project](#2).

### <a href="#windows-binary" name="windows-binary">#</a> Windows: Installing the Pre-compiled Binary

####Prerequisites
To run Crux on Windows you must have installed the following software:

* Thermo 32-bit MSFileReader

Be sure to install the 32-bit version of the library.

####Download the files
Crux is freely available for download under an Apache 2.0 license. On the download page click on the link "Windows x86 Build (zipped)" to download them to your computer. If a window pops up, select the "Save" button.

####Unpack and install the files
Browse to the downloaded file in Windows Explorer. Right click on the file and choose "Expand All ...". You should now have a new directory, `crux-x.xx`, with two subdirectories bin, and doc. The crux programs are in the bin directory. You may want to move these files somewhere more convenient for you, but otherwise, that completes the installation! Skip ahead to Setting the PATH on Windows for notes on how to run the executables.

####Set the PATH
On Windows systems you can always run any of the programs by specifying its full path. For example,

	C:\ <install directory>\bin\crux search-for-matches
	
In this final step of the installation, we'll set things up so that you can also run any of the programs by just typing its name. This is done with an environment variable called PATH. PATH is a list of directoires that the computer looks for executable programs. Find out what is in the current list with this command:

	C:\ echo %PATH%
	
The value returned might look something like this:

	C:\Windows;C:\Windows\system32;c:\Program Files\Microsoft Visual Studio 10.0\Common7\IDE\;c:\Program Files\Microsoft Visual Studio 10.0\VC\BIN;C:\Program Files\Common Files\Microsoft Shared\ Windows Live;C:\Windows\System32\Wbem;C:\Windows\S ystem32\WindowsPowerShell\v1.0\;c:\Program files\Notepad++;C:\Program Files\Tort oiseSVN\bin;C:\Program Files\Windows Live\Shared;C:\Program Files\CMake 2.8\bin; C:\Program Files\PuTTY;C:\Program Files\Subversion\bin;C:\GnuWin32\bin;
	
The first directory in this semicolon-separated list, `C:\Windows`, is the first place the computer looks for a program. If it doesn't find the program there it tries `C:\WIndows\system32` and so on down the list. You can either move the crux executable files into one of these directories already in your PATH, or you can add `<install directory>` to PATH with this command. (Remember to replace the word `<install directory>` with the actual location of your installation)

	set PATH=%PATH%;<install directory>\bin;
The change to PATH is only temporary. To make a permanent change to the path follow the following steps:

1. Open the "Control Panel".
* Click on "System and Security".
* Click on "Advanced System Settings" in the list on the left (you may be asked to give an administrative password).
* You'll presented with a dialog box titled "System Properties". Click on the button labeled "Environment Variables".
* You'll be presented with a new dialog box containing two list boxes: one for user environment variables, and one for system environment variables. Select the PATH variable in the system variables list, and click on the "Edit ..." button.
* Add the full path to the directory containing the Crux executables to the list of directories. Remember that each PATH entry should end with a semi-colon.

You're done, continue on to step 2, [Setting up the demo project](#2).

### <a href="#windows-source" name="windows-source">#</a> Windows: Installing from Source

####Prerequisites
To build Crux from source on Windows you must have installed the following software:

* Microsoft Visual C++ 2010
* Microsoft Visual C++ 2010 Service Pack 1
* Thermo 32-bit MSFileReader
* CMake
* SVN
* wget (must support SSL)

Other verions of Visual Studio will not be able to build Crux! You must have Visual Studio 2010, and you must install Service Pack 1. The Express version of Visual Studio 2010 does work.

wget is available as part of the GnuWin32 package. The easiest way to install this is via the GetGnuWin32 package manager. wget must support SSL, so install the optional packages provided by GetGnuWin32 .

wget and svn must be included in the execution PATH for your system. See the section Setting the PATH on Windows for more information.

####Download the file
The Crux source is freely available for download under an Apache 2.0 license>. On the download page click on the link "Source (zipped)", and choose save to disk if prompted.

The file you downloaded should be named something like `crux-x.xx.Source.zip`, where x.xx is replaced with the most current release version.

####Unpack the files
Unpack the downloaded file using Windows Explorer. Right click on the file icon and click on "Extract All ...". A new directory named `crux-x.xx.Source` will be created containing the source code and documentation for Crux.

Note Windows has a fixed maximum path length of 270 characters. The Crux build creates several levels of directories so it is very easy to exceed this limit. This will cause the build to fail. To avoid this problem we recommend keeping the Crux source and build trees very near the disk root (typically `C:\`)

Open the command line development environment.

Under the "Start" button open the "Visual Studio 2010" folder, then open the "Visual Studio Tools" folder, and double-click on "Visual Studio x64 Win64 Command Prompt". This will launch a cmd shell configured to use the Microsoft development tools.

####Configure
In the cmd shell you just opened, move into the `crux-x.xx.Source` directory from step 1. Run the command

	cmake -DCMAKE_INSTALL_PREFIX:PATH=<install directory> -G "Visual Studio 10" .
where `<install directory>` is the location where you wish to install the crux programs.

####Build and install
To complete the process run the command

	cmake --build . --config Release --target ALL_BUILD

Once you've run the configuration step you can also build and debug Crux from within the Visual Studio IDE. CMake will generate a file named `crux.sln`. Open this file in Visual Studio. You can select the build configuration: Release or Debug in the combo box under the menu bar.

Note The Crux build will not automatically generate debugging symbols even if you've selected the Debug configuration. Right click on the Crux folder in the "Solution Explorer" and select "Properties" from the pop-up menu. In the "Properties" page select the "Debugging" item under the "Linker" properties. Set "Generate Debug Info" to "Yes".

Changes made to the Visual Studio solution and project files will be lost the next time CMake regenerates the project.

####Set the PATH
On Windows systems you can always run any of the programs by specifying its full path. For example,

	C:\ <install directory>\bin\crux search-for-matches
	
In this final step of the installation, we'll set things up so that you can also run any of the programs by just typing its name. This is done with an environment variable called PATH. PATH is a list of directoires that the computer looks for executable programs. Find out what is in the current list with this command:

	C:\ echo %PATH%
	
The value returned might look something like this:

	C:\Windows;C:\Windows\system32;c:\Program Files\Microsoft Visual Studio 10.0\Common7\IDE\;c:\Program Files\Microsoft Visual Studio 10.0\VC\BIN;C:\Program Files\Common Files\Microsoft Shared\ Windows Live;C:\Windows\System32\Wbem;C:\Windows\S ystem32\WindowsPowerShell\v1.0\;c:\Program files\Notepad++;C:\Program Files\Tort oiseSVN\bin;C:\Program Files\Windows Live\Shared;C:\Program Files\CMake 2.8\bin; C:\Program Files\PuTTY;C:\Program Files\Subversion\bin;C:\GnuWin32\bin;
	
The first directory in this semicolon-separated list, `C:\Windows`, is the first place the computer looks for a program. If it doesn't find the program there it tries `C:\WIndows\system32` and so on down the list. You can either move the crux executable files into one of these directories already in your PATH, or you can add `<install directory>` to PATH with this command. (Remember to replace the word `<install directory>` with the actual location of your installation)

	set PATH=%PATH%;<install directory>\bin;
The change to PATH is only temporary. To make a permanent change to the path follow the following steps:

1. Open the "Control Panel".
* Click on "System and Security".
* Click on "Advanced System Settings" in the list on the left (you may be asked to give an administrative password).
* You'll presented with a dialog box titled "System Properties". Click on the button labeled "Environment Variables".
* You'll be presented with a new dialog box containing two list boxes: one for user environment variables, and one for system environment variables. Select the PATH variable in the system variables list, and click on the "Edit ..." button.
* Add the full path to the directory containing the Crux executables to the list of directories. Remember that each PATH entry should end with a semi-colon.

You're done, continue on to step 2, [Setting up the demo project](#2).


## <a href="#2" name="2">#2</a> Setting up the demo project

### Project Structure
In an effort to make the purpose of the various files more transparent, our demo project will take on the following structure.

	crux-demo
	|-- input
	`-- output

The root of our project will be called `crux-demo`, and once you've created it the rest of the tutorial will assume that is your current working directory when issuing commands.  The `input` directory will contain our datasets. The `output` directory is where we will store all of the output from Crux.
 
Let's give our demo project some structure.

	cd <DIR>         # change directory to desired project location
	mkdir crux-demo  # create the root directory
	cd crux-demo
	mkdir input      # create the directory for our input data
	mkdir output     # create the directory for our output data


> Note:  
> We've created the crux-demo root directory.  
> That means all blocks of commands are intended to be run from the root directory!  

### Downloading the Sample Input Data

Crux requires two types of input data, the mass spectra and a protein database. For this project, our mass spectra will be in the [ms1]() and [ms2]() file formats, and our protein database will be in the [fasta]() file format.

> Note:  
> If you wish to use other file formats later on, see the page on [file formats]() for information on acceptable file formats as well as information on converting to and from various file formats.

Let's gather a sample dataset to work with.
	
	cd input
	
	# Make a directory for the human dataset
	mkdir human
	cd human
	
	# Download the human dataset
	curl -O [DOWNLOAD URL]/human.ms1.gz     # ms1 mass spectra file
	curl -O [DOWNLOAD URL]/human.ms2.gz     # ms2 mass spectra file
	curl -O [DOWNLOAD URL]/human.fasta      # protein database file 
		
	# Extract gunzipped files
	gunzip human.ms1.gz
	gunzip human.ms2.gz
	
	cd ../..
	
	# Make a directory for Crux's output corresponding to the human dataset
	mkdir -p output/human

Your project should now have the following directory structure. Only newly added files are shown.

	crux-demo
	├── input
	│   └── human
	│       ├── human.fasta
	│       ├── human.ms1
	│       └── human.ms2
	└── output
	    └── human

## <a href="#3" name="3">#3</a> Understanding common conventions in Crux

Each Crux command has a set of parameters for customizing its usage. Typical parameters include options for specifying the formats of input and output data and for fine-tuning aspects of the algorithms used during spectral analysis. Crux allows you to specify parameters as command line arguments and as a [Crux parameter file]() with the `--parameter-file` command line option. Command specific parameters can be found in the [Crux documentation](). By default, Crux commands save used parameters in a [Crux parameter file](). This allows you to easily reproduce a customized analysis even when specifying parameters via the command line.

Two common options for controlling the output from a Crux command are the `--output-dir` and `--overwrite` options. `--output-dir` is used to specify where you want to save the output from a command. `--overwrite` is used to specify whether you wish to to allow the command to overwrite files if they already exist. By default, `--output-dir` is set to `crux-ouput` and `--overwrite` is set to `F`, for false.

During command execution, Crux outputs alot of useful information to the console such as performance details, warnings, and errors. For convenience, Crux saves all of theses messages in [Crux log files]().

Analysis results are saved in a [Crux tab-delimited file](), which is a Crux-specific file format that is simply a human readable text file containing a table. Crux tab-delimited files contain columns for all possible results from all Crux commands. Individual Crux commands only fill in the columns relevant to their part of the analysis, so don't be surprised if some columns are blank. More details on the Crux tab-delimited file format can be found in the [documentation]()

As you may be aware, there are many different file formats for MS/MS data and protein databases. Within Crux, support for MS/MS file formats varies between specific commands, all of which are covered in the [Crux documentation](). Crux is much more restrictive when it comes to the format of protein databases, exclusively supporting the FASTA file format. For converting between file formats when working on your own projects after the tutorial, both [proteowizard]() and the [Trans Proteomic Pipline]() have the conversion tools needed.



## <a href="#4" name="4">#4</a> Assigning high-resolution precursor masses using Bullseye

Bullseye assigns high-resolution precursor masses to mass spectra. Under the hood, bullseye first uses [Hardklör]() to identify persistent peptide isotope distributions in MS1 spectra. Then the monoisotopic mass from each persistent peptide isotope distribution is assigned to the  MS2 spectra corresponding to a specific time and mass range.

### From the [Bullseye Documentation]()

#####Usage

	crux bullseye [options] <MS1 spectra> <MS2 spectra>
	
#####Input

File|Description|File Fomats
-|-|-
`<MS1 spectra>` | MS1 spectra | .ms1 .bms1 .cms1 .mzXML
`<MS2 spectra>` | MS2 spectra | .ms2 .bms2 .cms2 .mzXML

#####Output

File|Description|File Fomat
-|-|-
`bullseye.params.txt` | [Crux parameter file]() | .txt
`bullseye.pid.<format>` | MS2 spectra with high resolution masses successfully assigned | .ms2 
`bullseye.no-pid.<format>`| MS2 spectra that could not be assigned high resolution masses | .ms2 
`bullseye.log.txt`| [Crux log file]() | .txt


### Using bullseye in our demo project

Let's assign high resolution precursor masses to our mass spectra. We'll output the results to the new directory, `output/human/bullseye`.

	crux bullseye --output-dir output/human/bullseye input/human/human.ms1 input/human/human.ms2

Your project should now have the following directory structure. Only newly added files are shown.

	crux-demo
	├── input
	│   └── human
	└── output
	    └── human
	        └── bullseye
	            ├── ISOTOPE.DAT
	            ├── bullseye.log.txt
	            ├── bullseye.no-pid.ms2
	            ├── bullseye.params.txt
	            ├── bullseye.pid.ms2
	            └── hardklor.mono.txt
	           
Note that bullseye outputs two MS2 spectra files, `bullseye.pid.ms2` and `bullseye.no-pid.ms2`. Bullseye only assigns precursor masses to MS2 spectra that correspond to a persistent peptide isotope distribution from an MS1 spectra. Consequently, not all MS2 spectra can be assigned precursor masses. `bullseye.pid.ms2` contains the spectra for which precursor masses could be assigned and `bullseye.no-pid.ms2` contains the rest. Bullseye's output links in with the rest of Crux by providing a new smaller subset of MS2 spectra with high-resolution precursor masses (`bullseye.pid`) that can be used as input to any of the other Crux commands that require MS2 spectra as input.

> Todo: explain how to use get-ms-spectrum to verify results



## <a href="#5" name="5">#5</a> Assigning peptides to spectra using Comet

Comet is the most simple and easy to use of the Crux peptide identification search engines. Comet can also be used to create decoy search results that can be used for post processing with percolator later on. 

### From the [Comet Documentation]()

#####Usage

	crux comet [options] <spectra> <protein input>

#####Input

File | Description | File Formats
-|-|-
`<spectra>` | MS2 spectra | .ms2 .cms2 .bms2 mzML1.1 mzML1.0 mzXML MGF Agilent Bruker FID/YEP/BAF ThermoRAW WatersRAW MGF mzIdentML
`<protein input>` | Protein database | .fasta

#####Output

File | Description | File Formats
-|-|-
`comet.params.txt` | [Crux parameter file]() | .txt
`comet.target.txt` | [Crux tab-delimited file]() | .txt
`comet.log.txt` | [Crux log file]() | .txt

### Using comet in our demo project

Let's assign peptides to our mass spectra. We now have the choice of using the MS2 spectra output from bullseye or the original MS2 spectra from our dataset. Let's just use the original MS2 spectra from our dataset for simplicity. Later on, we are going to want to be able to post-process our search results, which requires that we also perform a decoy search with comet now.  We can tell comet to perform a decoy search with the `--decoy_search` option, with `0` for no search, `1` for a concatenated search (outputs decoy and target results in a single [Crux tab-delimited file]()), or `2` for separate searches (outputs decoy and target results in two separate [Crux tab-delimited files]()). Separate searches will be depreciated in future versions of Crux, so let's perform a concatenated search. 

	crux comet --output-dir output/human/comet --decoy_search 1 input/human/human.ms2 input/human/human.fasta

Your project should now have the following directory structure. Only newly added files are shown.

	crux-demo
	├── input
	│   └── human
	└── output
	    └── human
	        ├── bullseye
	        └── comet
	            ├── comet.log.txt
	            ├── comet.params.txt
	            ├── comet.target.pep.xml
	            └── comet.target.txt
	         
Congratulations! You now have your peptide spectrum matches stored in `comet.target.txt`. The only thing left todo is post-process the results with one of the Crux post-processing tools for further analysis.


## <a href="#6" name="6">#6</a> Assigning peptides to spectra using Tide

Tide is the fastest of the Crux search engines, but also requires indexing the input protein database before searching for matches. Indexing a database only needs to be done once, so tide is split into two separate commands, `crux tide-index` and `crux tide-search`. This make performing multiple searches on the same data much faster and more convenient. Naturally, the output from `tide-index` is the input of `tide-search`.

### From the [tide-index documentation]()

#####Usage

	crux tide-index [options] <protein input file> <index name>

#####Input

File | Description | File Formats
-|-|-
`<protein input file>` | Protein database | .fasta
`<index name>` | Name of the index to be created  | -

#####Output

File | Description | File Formats
-|-|-
`tide-index.params.txt` | [Crux parameter file]() | .txt
`tide-index.log.txt` | [Crux log file]() | .txt


### From the [tide-search documentaion]()

#####Usage

	crux tide-search [options] <spectra> <peptide index>

#####Input

File | Description | File Formats
-|-|-
`<spectra>` | MS2 spectra. | .ms2 .cms2 .bms2 mzML1.1 mzML1.0 mzXML MGF Agilent Bruker FID/YEP/BAF ThermoRAW WatersRAW MGF mzIdentML
`<peptide index>` | Index created with tide-index. | index directory name

#####Output

File | Description | File Formats
-|-|-
`tide-search.target.txt` | [Crux tab-delimited file]() containing the target PSMs | .txt
`tide-search.decoy.txt` | [Crux tab-delimited file]() containing the decoy PSMs | .txt
`tide-search.params.txt` | [Crux parameter file]() | .txt
`tide-search.log.txt` | [Crux log file]() | .txt

	
### Using tide in our demo project

Let's begin by indexing our protein database. Don't be alarmed by warnings about invalid sequences, they are expected and will be silenced in future versions of Crux.

	crux tide-index --output-dir output/human/tide input/human/human.fasta output/human/tide/tide-index

Now lets assign peptides to our mass spectra. Like searching with comet, we'll want to perform a concatenated search (decoy and target results in a single [Crux tab-delimited file]()), which is the recommended decoy search method. Unlike comet, tide-search performs a separate search by default, which is why `tide-search.target.txt` and `tide-search.decoy.txt` are listed in the `tide-search` output table above. We can override this default with the `--concat` option, specifying `T` for true. Instead of the two separate output files, `tide-search.target.txt` and `tide-search.decoy.txt`, only the single file, `tide-search.txt`, containing the concatenated decoy and target search results will be output.

	crux tide-search --output-dir output/human/tide --concat T input/human/human.ms2 output/human/tide/tide-index

Your project should now have the following directory structure. Only newly added files are shown.

	crux-demo
	├── input
	│   └── human
	└── output
	    └── human
	        ├── bullseye
	        ├── comet
	        └── tide
	            ├── tide-index
	            │   ├── auxlocs
	            │   ├── pepix
	            │   ├── pepix.nopeaks.tmp
	            │   └── protix
	            ├── tide-index.decoy.fasta
	            ├── tide-index.log.txt
	            ├── tide-index.params.txt
	            ├── tide-search.log.txt
	            ├── tide-search.params.txt
	            └── tide-search.txt

Go ahead and give yourself another pat on the back, you've sucessfully identified potential peptide spectrum matches with tide. Just like the output from comet, the only thing left todo is post-process these results with one of Crux's post-processing tools.

## <a href="#7" name="7">#7</a> Post-processing with Percolator
Percolator is a semi-supervised learning algorithm that dynamically learns to separate target from decoy peptide-spectrum matches. Percolator also exists as a standalone and can be found at [per-colator.com](http://per-colator.com). The version of Percolator that is included in Crux differs slightly with the original, for more information on the differences, consult the [Percolator documentation]()

### From the [Percolator documentation]()

#####Usage

	crux percolator [options] <search results>

#####Input

File | Description | File Formats
-|-|-
`<search results>` | A collection of target and decoy peptide spectrum matches. | pin.xml, SQT, PepXML, Crux tab-delimited file, a list of files, or a tab-delimited table of features.

#####Output

File | Description | File Formats
-|-|-
`percolator.target.pout.xml` | [percolator xml file]() | .xml
`percolator.proteins.txt` | [Crux tab-delimited file]() containing the protein matches | .txt
`percolator.peptides.txt` | [Crux tab-delimited file]() containing the peptide matches | .txt
`percolator.psms.txt` | [Crux tab-delimited file]() containing the spectrum matches | .txt
`percolator.pep.xml` | [pepXML file]() containing the spectrum matches | .xml
`percolator.mzid` |  [mxIdent file]() containing the protein, peptide, and spectrum matches | .mzid
`percolator.params.txt` | [Crux parameter file]() | .txt
`percolator.log.txt` | [Crux log file]() | .txt

### Using percolator in our demo project

Let's post-process the peptide spectrum matches from comet with percolator. Because we performed a concatenated search with comet, we only need to provide `comet.target.txt` as input.

	crux percolator --output-dir output/human/comet/percolator output/human/comet/comet.target.txt 
	
Similarly, when post-processing the search results from tide, we only need to provide the single file `tide-search.txt` as input, since it contains the target and decoy search results.
	
	crux percolator --output-dir output/human/tide/percolator output/human/tide/tide-search.txt 

Your project should now have the following directory structure. Only newly added files are shown.

	crux-demo
	├── input
	│   └── human
	└── output
	    └── human
	        ├── bullseye
	        ├── comet
	        │   └── percolator
	        │       ├── make-pin.pin
	        │       ├── percolator.decoy.peptides.txt
	        │       ├── percolator.decoy.psms.txt
	        │       ├── percolator.log.txt
	        │       ├── percolator.params.txt
	        │       ├── percolator.target.peptides.txt
	        │       └── percolator.target.psms.txt
	        └── tide
	            ├── percolator
	            │   ├── make-pin.pin
	            │   ├── percolator.decoy.peptides.txt
	            │   ├── percolator.decoy.psms.txt
	            │   ├── percolator.log.txt
	            │   ├── percolator.params.txt
	            │   ├── percolator.target.peptides.txt
	            │   └── percolator.target.psms.txt
	            └── tide-index

As it runs, percolator reports a summary of the analysis to the console (the entire output also gets saved to `percolator.log.txt` for convenience). During the first 10 iterations of training the learning algorithm, the output displays the number of PSMs over a q-value of 0.01 are estimated by cross-validation. After training, a table of the set of weights found will be shown. The table will contain 2 rows, the top row containing the normalized weights, and the bottom row containing the raw weights. The details of the interpretation of the individual weights can be found in the [Percolator documentation](). Another notable piece of the analysis summary is the line containing `Selecting pi_0=`, which indicates the percent estimate of number of PSMs that are incorrect matches, given in decimal notation. Next you will want to quantify this analysis with spectral-counts.

## <a href="#8" name="8">#8</a> Post-processing with Barista

Barista performs both peptide spectrum match verification and protein inference. 

###From the [barista documentation]()

#####Usage

	crux barista [options] <protein-database> <spectra> <search results>

#####Inputs

File | Description | File Formats
-|-|-
`<protein-database>` | Target and a decoy protein databases, concatenated or by providing the decoy database separately | .fasta
`<spectra>` | MS2 spectra | .ms2
`<search results>` | [Crux tab-delimited file]() | .txt

#####Outputs

File | Description | File Formats
-|-|-
`barista.xml`| [barista xml file]() containing the protein, subset proteins, peptides, and spectrum matches | .xml 
`barista.target.proteins.txt` | [Crux tab-delimited file]() containing a ranked list of groups of indistinguishable target proteins with associated Barista scores and q-values and with peptides that contributed to the identification of the protein group. | .txt
`barista.target.subset-proteins.txt` | [Crux tab-delimited file]() containing a subset of indistinguishable proteins from `barista.target.proteins.txt` | .txt
`barista.target.peptides.txt` | [Crux tab-delimited file]() containing a ranked list of target peptides with the associated Barista scores and q-values. | .txt
`barista.target.psm.txt` | [Crux tab-delimited file]() containing a ranked list of target peptide-spectrum matches with the associated Barista scores and q-values. | .txt
`barista.params.txt` | [Crux parameter file]() | .txt
`barista.log.txt` | [Crux log file]() | .txt


###Using barista in our demo project

As shown in the input table above, barista requires a concatenated target-decoy database, so we'll have to create that now.

	crux create-index --decoys protein-shuffle input/human/human.fasta output/human/temp-index
	mkdir output/human/protein-databases
	mv output/human/temp-index/human-random.fasta output/human/protein-databases/decoy-protein-shuffle.fasta
	rm -rf output/human/temp-index
	cat input/human/human.fasta output/human/protein-databases/decoy-protein-shuffle.fasta > output/human/protein-databases/target-decoy-protein-shuffle.fasta

Now, lets post-process the peptide spectrum matches from comet with barista.
	
	crux barista --output-dir output/human/comet/barista output/human/protein-databases/target-decoy-protein-shuffle.fasta input/human/human.ms2 output/human/comet/comet.target.txt

And again for the tide matches.

	crux barista --output-dir output/human/tide/barista output/human/protein-databases/target-decoy-protein-shuffle.fasta input/human/human.ms2 output/human/tide/tide-search.txt

Your project should now have the following directory structure. Only newly added files are shown.
	
	crux-demo
	├── input
	│   └── human
	└── output
	    └── human
	        ├── bullseye
	        ├── comet
	        │   ├── barista
	        │   │   ├── barista.log.txt
	        │   │   ├── barista.params.txt
	        │   │   ├── barista.target.peptides.txt
	        │   │   ├── barista.target.proteins.txt
	        │   │   ├── barista.target.psms.txt
	        │   │   ├── barista.target.subset-proteins.txt
	        │   │   └── barista.xml
	        │   └── percolator
	        ├── protein-databases
	        │   ├── decoy-protein-shuffle.fasta
	        │   └── target-decoy-protein-shuffle.fasta
	        └── tide
	            ├── barista
	            │   ├── barista.log.txt
	            │   ├── barista.params.txt
	            │   ├── barista.target.peptides.txt
	            │   ├── barista.target.proteins.txt
	            │   ├── barista.target.psms.txt
	            │   ├── barista.target.subset-proteins.txt
	            │   └── barista.xml
	            └── percolator

> Todo: Explain output

## <a href="#9" name="9">#9</a> Quantification using spectral-counts
 Spectral-counts ranks a list of peptides or proteins from a collection of scored PSMs by one of four types of quantification scores. 
 
###From the [spectral-counts documentation]()

#####Usage

	crux spectral-counts [options] <input PSMs>

#####Inputs

File | Description | File Formats
-|-|-
`<input PSMs>` | [Crux tab-delimited file]() containing PSMs | .txt

#####Output

File | Description | File Formats
-|-|-
`spectral-counts.target.txt` | [Crux tab-delimited file]() containing ranked peptides or proteins | .txt
`spectral-counts.log.txt` | [Crux log file]() | .txt

###Using spectral-counts in our demo project
	crux spectral-counts --output-dir output/human/comet/percolator/spectral-counts --protein-database input/human/human.fasta output/human/comet/percolator/percolator.target.psms.txt
	
	#crux spectral-counts --output-dir output/human/tide/percolator/spectral-counts --protein-database input/human/human.fasta output/human/tide/percolator/percolator.target.psms.txt
	#Outputs:
	#	ERROR: Flanking AA count did not match protein count!
	#	ERROR: Failed to parse peptide src.
	#	ERROR: Failed to parse peptide (tab delimited)
	#	ERROR: Failed to parse tab-delimited PSM match
	
	crux spectral-counts --output-dir output/human/comet/barista/spectral-counts --protein-database input/human/human.fasta output/human/comet/barista/barista.target.psms.txt
	
	#crux spectral-counts --output-dir output/human/tide/barista/spectral-counts --protein-database input/human/human.fasta output/human/tide/barista/barista.target.psms.txt
	#Outputs:
	#	ERROR: Flanking AA count did not match protein count!
	#	ERROR: Failed to parse peptide src.
	#	ERROR: Failed to parse peptide (tab delimited)
	#	ERROR: Failed to parse tab-delimited PSM match

Your project should now have the following directory structure. Only newly added files are shown.

	crux-demo
	├── input
	│   └── human
	└── output
	    └── human
	        ├── bullseye
	        ├── comet
	        │   ├── barista
	        │   │   └── spectral-counts
	        │   │       ├── spectral-counts.log.txt
	        │   │       ├── spectral-counts.params.txt
	        │   │       └── spectral-counts.target.txt
	        │   └── percolator
	        │       └── spectral-counts
	        │           ├── spectral-counts.log.txt
	        │           ├── spectral-counts.params.txt
	        │           └── spectral-counts.target.txt
	        ├── protein-databases
	        └── tide
	            ├── barista
	            │   └── spectral-counts
	            │       ├── spectral-counts.log.txt
	            │       ├── spectral-counts.params.txt
	            │       └── spectral-counts.target.txt
	            ├── percolator
	            │   └── spectral-counts
	            │       ├── spectral-counts.log.txt
	            │       ├── spectral-counts.params.txt
	            │       └── spectral-counts.target.txt
	            └── tide-index
	            
> Todo: Explain output

---

>This tutorial was written by Dylan Holmes on June 20th, 2014.