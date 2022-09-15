# Rcrust version 2022-09
Rcrust is a program for phase stability calculations with path dependence.
View the Rcrust manual for more information.
To launch Rcrust double click the Rcrust workspace "Rcrust.RData" file in the code folder.

Its main references are:
Mayne, M. J., Moyen, J. F., Stevens, G., & Kaislaniemi, L. (2016). Rcrust: a tool for calculating path‐dependent open system processes and application to melt loss. Journal of Metamorphic Geology.

## Setup
*Git Tutorial - Public Repo.md* provides basic setup instructions for the
repository if you are not yet comfortable with Git.

R must be installed in order to use Rcrust. For installation of R,
see https://cloud.r-project.org/ - There are instructions for both Linux
and Windows. (There is also a Mac version; however Rcrust has not been
adapted to Mac yet). For Linux, there are different versions depending on your
exact OS (Debian, Ubuntu, etc). Probably the easiest way to install it is with
the following two commands:
```bash
sudo apt-get update
sudo apt-get install r-base
sudo apt install libgeos-dev
```

## Running Via GUI

For setup and execution of Rcrust, see *Getting Started (Rcrust Manual).pdf* -
This PDF was written primarily for Windows: On Linux, it is not possible to
double-click the RCrust.RData and immediately execute it (since R GUI is
not available on Linux).

There are two ways to execute the GUI on Linux:

The first method is similar to above: If you have **Rstudio** or similar
installed, double-click the *Rcrust.RData* file from a file explorer; this
will open it in the editor. Then from the terminal inside the editor, execute
`.First()`. This should launch the GUI.

The second method (which I suggest) is to use the terminal to load it via a script. This can be
done by executing the `RUN_GUI_LINUX.sh` script. Open a terminal (shell)
instance, and ensure you are in the `code/` directory. Then
execute the script via the following command:
```bash
bash ./RUN_GUI_LINUX.sh
```

For the initial first run, the program may fail and complain that it is unable to
 install the packages. If this is the case, see the
**Initial Setup / Package Install Errors** heading below.
All future executions should run correctly.

If there is an error saying that it is unable to execute the `meemum` file
(often due to permission errors), this can be a bug that occurs when
a zip/tar.bz2 file removes execution permissions on Linux. To fix this,
go to the `data/` directory, find the relevant meemum file you are using
in your input for your project, right click on it and tick the box saying
"allow execution". The alternative is to use the terminal and the command
```bash
chmod +x <meemum_file_name>
```
For example, if you are using the `meemum_686` executible, the command will
be
```Bash
chmod +x meemum_686
```

## Running From Command Line
This applies to both the Linux version as well as the Windows version.
On Windows, all testing was done using *Git Bash*, since this would need
to have been installed to clone the repository originally. It *should*
work with the command console of Windows, however I have not done much
testing on it.

#### Initial Setup / Package Install Errors
Rscript sometimes struggles to install the required packages when the
script is run (usually has to do with permission errors). The suggested
way around this is to either manually install the packages via the
interpreter.

This can be done by the following: either just `R` or `R --vanilla`.
Once the interpreter is running, execute the following commands:
```bash
install.packages("shiny")
install.packages("raster")
install.packages("rgeos")
install.packages("RColorBrewer")
```
If it asks to use a personal library, say 'yes'. This should bypass some other permission errors later.

<!--OLD: Ignore this line: install.packages("grDevices") -->

If there are still issues with the libraries not being found, it is
possible it installed them in the wrong location. You can manually
move them across by doing the following:

+ Open the interpreter (`R --vanilla`)
+ in it, run `.libPaths()`
+ If this returns two paths, then go to the two locations, and copy
the directories of the libraries from one to the other
(for example, this may be from documents/r/library to
program_files/r/library)
+ If it only returns one path, try run R in admin mode, and repeat the step above.

#### Execution
Move into the `code/` directory, then execute the script with an input file.
```bash
cd code
```
```bash
Rscript main.r <input_file>
```
For example, `Rscript main.r Example1` to run project Example1 and
input file *Example1.txt*

This will create a `<input_file>.RData` file which can be loaded
to view the results and create the desired outputs.

## Linux Variation Of Files
Only the new (Linux) `meemum` file in `data/` is needed.
Currently, the following versions of the Linux Meemum are already included
in Rcrust:
+ Version 686
+ Version 689

For future releases: The entire archive can be downloaded from
https://petrol.natur.cuni.cz/~ondro/Perple_X_Z.Z.Z_Linux_64_gfortran.tar.gz
Where *Z_Z_Z* is the version corresponding to the version of Perple_X used
in Windows. All the other files are common between the two archives (Windows
and Linux), and thus only the `meemum` file is needed as an alternative
to the executable. It can be renamed to `meemum_ZZZ`(with ZZZ being the
version number) to match the format of the Windows Meemum naming.

Another change is that in `data/hp11ver_xxx.dat` data file, the
comment/reference is edited so that the *ó* is replaced by *o*. While this
is not strictly needed, `grep` will throw multiple warnings every run if this is
not done.

For further Perple_X documentation, which Meemum is a part of, see:
http://www.perplex.ethz.ch/perplex/ibm_and_mac_archives/

### Modifying A Project to Run on Linux
To modify a Project to work on Linux, only 1 modification must be made to the
input file:
+ `meemum_path` must be changed from `meemum_xxx.exe` to `meemum_xxx`, with xxx
being the version number.

### Magic Number Issue
Sometimes when a Project is loaded from an input file, then run
(producing a .RData file) on Windows, and then re-loaded on Linux,
Linux will have difficulty loading the file due to a
*bad restore file magic number*. This seems to be an issue with
how the different R versions save the file. Do not worry, however -
if you delete the .RData file, you can load the input file as normal,
and can just re-run the calculation.

There doesn't seem to be an issue with loading the file in Windows
after being calculated on Linux, however.

## Additional Comments
+ Ensure that in your Project *input.txt* file the projects_directory
and working_file variables are correct.
+ Inputs and Outputs directories for projects are not auto-created
if they are missing.
