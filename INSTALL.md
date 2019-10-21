# Installation and running instructions

## Downloading our software

The latest version can be either downloaded as a release, or you can run:

```
$ git clone https://github.com/igemsoftware2019/Team_Rice_Software.git
```

## Prerequisites

### Installing NuPack

[NuPack](http://www.nupack.org/) command line software suite must be installed and on your PATH. NuPack must be manually built and installed on your machine. To do so, [register for downloads on NuPack](http://www.nupack.org/downloads/register). Download NuPack v3.2.2 source code. This should be a tarball. Untar the source code and change directory into the source code directory with the following commands:

```
$ tar -xzf nupack3.2.2.tar.gz
$ cd nupack3.2.2
```

NuPack requires CMake and a C/C++ compiler to be built. Please install such a compiler appropriate to your operating system if it is not already installed. You can use a package manager such as <code>apt</code> for Debian-based operating systems or [Homebrew](https://brew.sh) for macOS. CMake can either be installed using the same package manager or found [here](https://cmake.org/).  

If you're using a newer version of NuPack, please consult the [current documentation](http://www.nupack.org/downloads/documentation) for build instructions. Otherwise if using NuPack v3.2.2, build and install NuPack using the following commands: 

```
$ mkdir build && cd build
$ cmake ../
$ make
$ sudo make install
```

If you wish to perform a custom install of NuPack, please consult [the documentation](http://www.nupack.org/downloads/documentation).

**Note**: At this time, it does not appear that NuPack works properly on macOS 10.14.x or 10.15.x (on the Rice iGEM 2019 software team's personal machine). NuPack should work on any version of Ubuntu. NuPack on Windows Subsystem for Linux (WSL) is untested to the best of our knowledge, but should work in theory.

### Python packages

Our software requires Python 3 and <code>pip3</code> to be installed and on your PATH (Python 3 must be executable using the command <code>python3</code>. Again, both of these can be installed using a package manager (e.g. <code>apt</code>, Homebrew) or following instructions from [the offical Python website](https://www.python.org/).

After Python 3 and <code>pip3</code> are installed (and our GitHub source code is cloned/downloaded), change your current directory to our source code and install Python package dependencies with the following command:

```
$ pip3 install -r requirements.txt
```
This should install [DEAP](https://deap.readthedocs.io/en/master/), [SCOOP](https://scoop.readthedocs.io/), [PySimpleGUI](https://pysimplegui.readthedocs.io/en/latest/), and their respective dependencies into your global Python packages. Please use a [virtualenv](https://docs.python.org/3/library/venv.html) if you don't feel comfortable doing so.

## Running our software

You can run our software by double clicking on <code>Run.sh</code> or typing the following command at your terminal:

```
$ ./Run.sh
```

If there is a permissions error, make the file executable:

```
$ sudo chmod a+x Run.sh
$ ./Run.sh
```
