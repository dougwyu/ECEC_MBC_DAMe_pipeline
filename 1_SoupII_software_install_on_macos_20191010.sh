#!/bin/bash
set -e
set -u
set -o pipefail

####################################################################################################
# Software installation on macOS
####################################################################################################
# The script uses quite a bit of open-source software and has been run on macOS.  Much of this software can be installed using Homebrew, which is a package manager for macOS. The whole installation process can require several hours, depending on internet speed and how much you have previously installed

## Command line tools for macOS
# installs compilers and other Unix tools for programming
# if you have Xcode installed, you do not need to install CLT. (CLT is much smaller than Xcode)
xcode-select --install

## Homebrew for macOS
# go to http://brew.sh and follow the instructions for installing Homebrew on macOS

## after Homebrew  is installed, run these brew installations
brew update; brew upgrade; brew cleanup # get the latest version homebrew
brew tap brewsci/bio # a "tap" is a repository of "installation formulas" of specialist software, here, bioinformatics
brew install git
brew install coreutils
brew install wget
brew install gnu-sed # (gsed == GNU version of sed == Linux version of sed)
brew install grep # gnu-grep
brew install gawk # gnu-awk
brew install brewsci/bio/seqkit
brew install brewsci/bio/vsearch # https://github.com/torognes/vsearch
brew install parallel # GNU parallel
brew install spades # http://cab.spbu.ru/software/spades/
brew install seqtk # https://github.com/lh3/seqtk
brew install sickle # https://github.com/najoshi/sickle
brew install python@2
brew install python@3
brew install samtools
brew install bwa
brew cleanup # remove unneeded files


## for software installed from github, I have a dedicated folder at ~/src/
mkdir ~/src/

## AdapterRemoval
# https://github.com/MikkelSchubert/adapterremoval
cd ~/src
git clone https://github.com/MikkelSchubert/adapterremoval.git
cd adapterremoval
make
sudo make install # requires your password
AdapterRemoval -h

## pandaseq
# https://github.com/neufeld/pandaseq/releases/download/v2.11/PANDAseq-2.11.pkg
cd ~/src/
wget https://github.com/neufeld/pandaseq/releases/download/v2.11/PANDAseq-2.11.pkg
# go to ~/src in the macOS Finder, double click on the PANDAseq-2.11.pkg, and follow the installation instructions
pandaseq -h

# python numpy and matplotlib
# https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html
cd ~/src/
wget https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh
bash Miniconda2-latest-MacOSX-x86_64.sh  # this installs miniconda, which is a package manager for python
conda install numpy
conda install matplotlib

## Begum
# http://github.com/shyamsg/Begum
cd ~/src
git clone https://github.com/shyamsg/Begum.git
# no building needed because these are python2 scripts

## Sumatra
cd ~/src/
wget https://git.metabarcoding.org/obitools/sumatra/uploads/251020bbbd6c6595cb9fce6077e29952/sumatra_v1.0.20.tar.gz
tar -zxvf sumatra_v1.0.20.tar.gz
cd sumatra_v1.0.20/
make CC=clang # in macOS, disables OpenMP, which isn't on macOS
mv sumatra /usr/local/bin # will need to use sudo mv on Ubuntu

## Sumaclust
cd ~/src/
wget https://git.metabarcoding.org/obitools/sumaclust/uploads/69f757c42f2cd45212c587e87c75a00f/sumaclust_v1.0.20.tar.gz
tar -zxvf sumaclust_v1.0.20.tar.gz
cd sumaclust_v1.0.20/
make CC=clang # in macOS, disables OpenMP, which isn't on macOS
mv sumaclust /usr/local/bin/


## GUI programs

## Atom:  Text editor
# download binary from https://atom.io # (currently 1.38.2)

# In Atom, go to Atom/Preferences (or type cmd+, (command comma)). This opens the Settings tab.  Click on the "Install" button on the left and type 'platformio-ide-terminal' in the 'Search packages' window, and when the package appears, click on the blue install button.
# Platformio-ide-terminal installs a terminal inside Atom. You can open terminal windows by clicking on the + sign in the lower left of the Atom window.

# This allows you to send commands from an Atom text window to a terminal window (just like in RStudio), using 'ctrl-enter'.  In some cases, 'ctrl-enter' does not work. This is because the 'key map' is missing.
# following the instructions on:  https://github.com/platformio/platformio-atom-ide-terminal/issues/67
# open Keymap... under the Atom menu
# add this text to the end of the keymap.cson file:
'atom-workspace atom-text-editor:not([mini])': 'ctrl-enter': 'platformio-ide-terminal:insert-selected-text'



## R:  Statistical analysis
# download binary from https://cran.r-project.org

## RStudio:  GUI to R
# download binary from https://rstudio.org

## R packages:  additional functionality in R
# Launch RStudio and run these commands in R. This step can take hours the first time

install.packages(c("tidyverse", "data.table", "vegan", "car", "RColorBrewer", "devtools", "BiocManager", "metacoder"), dependencies = TRUE)
library(BiocManager)
BiocManager::install(c("GenomicRanges", "Biobase", "IRanges", "AnnotationDbi", "dada2", "phyloseq")) # install Bioconductor packages
BiocManager::valid() # identify packages that are out of date or unexpected versions. Creates a command for updating out of date bioconductor packages
