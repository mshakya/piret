#!/usr/bin/env bash

set -e # Exit as soon as any line in the bash script fails

ROOTDIR=$( cd $(dirname $0) ; pwd -P ) # path to main PiReT directory


echo
exec &> >(tee -a  install.log)
exec 2>&1 # copies stderr onto stdout

# create a directory where all dependencies will be installed
cd $ROOTDIR
mkdir -p thirdParty
cd thirdParty


# create a directory to add short cuts to dependencies
mkdir -p $ROOTDIR/bin
# export "PATH=$PATH:$ROOTDIR/bin/"
export "PATH=$ROOTDIR/thirdParty/miniconda/bin/:$ROOTDIR/bin/:$PATH"
if [[ "$OSTYPE" == "darwin"* ]]
then
{
export PERL5LIB="$ROOTDIR/ext/lib/perl5:$ROOTDIR/ext/lib/perl5/darwin-thread-multi-2level:$PERL5LIB"
}
else
{  
export PERL5LIB="$ROOTDIR/ext/lib/perl5:$ROOTDIR/lib/perl5/auto/:$PERL5LIB"
}
fi

# Add pythonpath
# export PYTHONPATH="$ROOTDIR/thirdParty/miniconda/lib/python2.7/site-packages/:$PYTHONPATH"

#Add R path
export R_LIBS="$ROOTDIR/ext/lib/R:$R_LIBS:$R_LIBS_USER"

# Minimum Required versions of dependencies
miniconda_VER=4.4.11
samtools_VER=1.3.1
bamtools_VER=2.4.0
R_VER=3.4.2

hisat2_VER=2.0.5
htseq_VER=0.6.1
#jbrowse_VER=1.12.1
stringtie_VER=1.3.3
subread_VER=1.5.0

# minimum required version of Scripting languages
# perl5_VER=5.16.3
python3_VER=3.5.2

#minimum required version of R modules
R_edgeR_VER=3.14.0
R_DESeq2_VER=1.12.4

#minimum required version of python modules
python_pandas_VER=0.19.2
python_Bio_VER=1.68
python_luigi_VER=2.7.3
python_plumbum_VER=1.6.3

################################################################################
#                           Installation recipes
################################################################################

install_python()
{
echo "--------------------------------------------------------------------------
                           installing python v$python3_VER
--------------------------------------------------------------------------------
"
conda install --yes -n piret python=$python3_VER
ln -sf $ROOTDIR/thirdParty/miniconda/bin/python $ROOTDIR/bin/python

echo "
------------------------------------------------------------------------------
                           python v$python3_VER installed
------------------------------------------------------------------------------
"
}

install_python_pandas()
{
echo "--------------------------------------------------------------------------
                installing Python module pandas $python_pandas_VER
--------------------------------------------------------------------------------
"
conda install --yes -n piret -c anaconda pandas=$python_pandas_VER
echo "
--------------------------------------------------------------------------------
                           pandas installed
--------------------------------------------------------------------------------
"
}

install_python_Bio()
{
echo "--------------------------------------------------------------------------
                installing Python module Bio $python_Bio_VER
--------------------------------------------------------------------------------
"
conda install --yes -n piret -c anaconda biopython=$python_Bio_VER
echo "
--------------------------------------------------------------------------------
                           Biopython installed
--------------------------------------------------------------------------------
"
}

install_python_luigi()
{
echo "--------------------------------------------------------------------------
                installing Python module luigi
--------------------------------------------------------------------------------
"
conda install --yes -n piret -c anaconda luigi=$python_luigi_VER
# pip install luigi
echo "
--------------------------------------------------------------------------------
                           luigi installed
--------------------------------------------------------------------------------
"
}


install_python_plumbum()
{
echo "--------------------------------------------------------------------------
                installing Python module plumbum
--------------------------------------------------------------------------------
"
conda install --yes -n piret -c conda-forge plumbum=1.6.3
echo "
--------------------------------------------------------------------------------
                           plumbum installed
--------------------------------------------------------------------------------
"
}


install_hisat2()
{
echo "--------------------------------------------------------------------------
                           installing hisat2 v$hisat2_VER
--------------------------------------------------------------------------------
"
conda install --yes -n piret -c bioconda hisat2=$hisat2_VER
ln -sf $ROOTDIR/thirdParty/miniconda/bin/hisat2 $ROOTDIR/bin/hisat2
ln -sf $ROOTDIR/thirdParty/miniconda/bin/hisat2-build $ROOTDIR/bin/hisat2-build
echo "
------------------------------------------------------------------------------
                           hisat2 v$hisat2_VER installed
------------------------------------------------------------------------------
"
}

install_stringtie()
{
echo "--------------------------------------------------------------------------
                           installing stringtie v$stringtie_VER
--------------------------------------------------------------------------------
"
conda install --yes -n piret -c bioconda stringtie=$stringtie_VER
ln -sf $ROOTDIR/thirdParty/miniconda/bin/stringtie $ROOTDIR/bin/stringtie
echo "
------------------------------------------------------------------------------
                           stringtie v$stringtie_VER installed
------------------------------------------------------------------------------
"
}


install_featureCounts()
{
echo "--------------------------------------------------------------------------
                           installing subread v$subread_VER
--------------------------------------------------------------------------------
"
conda install --yes -n piret -c bioconda subread=$subread_VER
ln -sf $ROOTDIR/thirdParty/miniconda/bin/featureCounts $ROOTDIR/bin/featureCounts
echo "
------------------------------------------------------------------------------
                           featureCounts v$subread_VER installed
------------------------------------------------------------------------------
"
}

install_samtools()
{
echo "--------------------------------------------------------------------------
                           Downloading samtools v $samtools_VER
--------------------------------------------------------------------------------
"
conda install --yes -n piret -c bioconda samtools=$samtools_VER
ln -sf $ROOTDIR/thirdParty/miniconda/bin/samtools $ROOTDIR/bin/samtools
echo "
--------------------------------------------------------------------------------
                           samtools v $samtools_VER installed
--------------------------------------------------------------------------------
"
}
install_bamtools()
{
echo "--------------------------------------------------------------------------
                           Downloading bamtools v $bamtools_VER
--------------------------------------------------------------------------------
"
conda install --yes -n piret -c bioconda bamtools=$bamtools_VER
ln -sf $ROOTDIR/thirdParty/miniconda/bin/bamtools $ROOTDIR/bin/bamtools
echo "
--------------------------------------------------------------------------------
                           bamtools v $bamtools_VER installed
--------------------------------------------------------------------------------
"
}

install_R()
{
echo "--------------------------------------------------------------------------
                           Installing R v $R_VER
--------------------------------------------------------------------------------
"
conda install --yes -n piret -c r r-base=$R_VER
ln -sf $ROOTDIR/thirdParty/miniconda/bin/R $ROOTDIR/bin/R
ln -sf $ROOTDIR/thirdParty/miniconda/bin/Rscript $ROOTDIR/bin/Rscript
echo "
--------------------------------------------------------------------------------
                           R v $R_VER installed
--------------------------------------------------------------------------------
"
}

install_miniconda()
{
echo "--------------------------------------------------------------------------
                           downloading miniconda
--------------------------------------------------------------------------------
"

if [[ "$OSTYPE" == "darwin"* ]]
then
{

  curl -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
  chmod +x miniconda.sh
  ./miniconda.sh -b -p $ROOTDIR/thirdParty/miniconda -f
  conda update --yes -n base conda
  # export PATH=$ROOTDIR/thirdParty/miniconda/bin:$PATH
  ln -sf $ROOTDIR/thirdParty/miniconda/bin/conda $ROOTDIR/bin/conda
  ln -sf $ROOTDIR/thirdParty/miniconda/bin/pip $ROOTDIR/bin/pip


}
else
{  

  wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
  chmod +x miniconda.sh
  ./miniconda.sh -b -p $ROOTDIR/thirdParty/miniconda -f
  conda update --yes -n base conda
  # export PATH=$ROOTDIR/thirdParty/miniconda/bin:$PATH
  ln -sf $ROOTDIR/thirdParty/miniconda/bin/conda $ROOTDIR/bin/conda
  ln -sf $ROOTDIR/thirdParty/miniconda/bin/pip $ROOTDIR/bin/pip

}
fi
echo "--------------------------------------------------------------------------
                miniconda installed and and bin added to PATH
--------------------------------------------------------------------------------
"
}


install_FaQCs()
{
echo "------------------------------------------------------------------------------
                      Installing FaQC
------------------------------------------------------------------------------
"
conda install --yes -n piret -c anaconda zlib=1.2.11
export CPLUS_INCLUDE_PATH=$ROOTDIR/thirdParty/miniconda/include/:$CPLUS_INCLUDE_PATH
cd $ROOTDIR/thirdParty
tar xvzf FaQCs-2.06.tar.gz
cd FaQCs
make
ln -sf $ROOTDIR/thirdParty/FaQCs/FaQCs $ROOTDIR/bin/FaQCs
cd $ROOTDIR
echo "
------------------------------------------------------------------------------
                        FaQC-2.01 Installed
------------------------------------------------------------------------------
"
}

install_gffread()
{
echo "--------------------------------------------------------------------------
                      installing gffread
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda gffread -n piret
echo "
--------------------------------------------------------------------------------
                           gffread installed
--------------------------------------------------------------------------------
"
}

# install_jbrowse()
# {
# echo "--------------------------------------------------------------------------
#                       installing Jbrowse
# --------------------------------------------------------------------------------
# "
# cd $ROOTDIR/thirdParty
# tar xvzf JBrowse-1.11.6.tar.gz
# cd JBrowse-1.11.6
# ./setup.sh
# mkdir -p -m 775 data
# cd $ROOTDIR/thirdParty
# ln -sf $ROOTDIR/thirdParty/JBrowse-1.11.6 $ROOTDIR/bin/JBrowse
# echo "
# --------------------------------------------------------------------------------
#                       Jbrowse installed
# --------------------------------------------------------------------------------
# "
# }

install_R_edgeR()
{
echo "--------------------------------------------------------------------------
                installing latest R package edgeR $R_edgeR_VER
--------------------------------------------------------------------------------
"
mkdir -p $ROOTDIR/ext/lib/R
echo "source('https://bioconductor.org/biocLite.R') 
      biocLite('edgeR', lib='$ROOTDIR/ext/lib/R')" | Rscript -
echo "
--------------------------------------------------------------------------------
                           latest version of edgeR installed
--------------------------------------------------------------------------------
"
}

install_R_DESeq2()
{
echo "--------------------------------------------------------------------------
                installing latest R package DESeq2 $R_DESeq2_VER
--------------------------------------------------------------------------------
"
mkdir -p $ROOTDIR/ext/lib/R
echo "source('https://bioconductor.org/biocLite.R') 
      biocLite('DESeq2', lib='$ROOTDIR/ext/lib/R')" | Rscript -
echo "
--------------------------------------------------------------------------------
                           DESeq2 installed
--------------------------------------------------------------------------------
"
}

checkSystemInstallation()
{
    IFS=:
    for d in $PATH; do
      if test -x "$d/$1"; then return 0; fi
    done
    return 1
}

checkLocalInstallation()
{
    IFS=:
    for d in $ROOTDIR/bin; do
      if test -x "$d/$1"; then return 0; fi
    done
    return 1
}


checkRpackages()
{

  echo "if(\"$1\" %in% rownames(installed.packages()) == FALSE) {0} else {1}"| Rscript - 

}


checkPythonModule()
{
  python -c "import $1" 2>&1
  return $?
}

containsElement () {
  local e
  for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
  return 1
}

print_usage()
{
cat << EOF
usage: $0 options
    If no options, it will check existing installation and run tools installation for those uninstalled.
    options:
    <help | help| -h>          show this help
    list                       show available tools for updates
    tools_name                 install/update individual tool
    force                      force to install all list tools locally
    
    ex: To update bowtie2 only
        $0 bowtie2
    ex: To update bowtie2 and bwa
        $0 bowtie2 bwa
        
EOF
}


### Main ####

mkdir -p  $ROOTDIR/thirdParty/miniconda/bin

if [ "$#" -ge 1 ]
then
  for f in $@
  do
    case $f in
      -h|help|-help)
        print_usage
        exit 0;;
      list)
        print_tools_list
        exit 0;;
      Alignment)
        for tool in "${alignments_tools[@]}"
        do
            install_$tool
        done
        echo -e "Alignment tools installed.\n"
        exit 0;;
      Utility)
        for tool in "${utility_tools[@]}"
        do
            install_$tool
        done
        echo -e "Utility tools installed.\n"
        exit 0;;
      force)
        for tool in "${all_tools[@]}"
        do
            install_$tool
        done
        ;;
      *)
        if ( containsElement "$f" "${alignments_tools[@]}" || containsElement "$f" "${utility_tools[@]}" )
        then
            install_$f
        else
            echo "$f: no this tool in the list"
            print_tools_list
        fi
        exit 0;;
    esac
  done
fi

###############################################################################
if ( checkSystemInstallation conda )
then
  conda_installed_VER=`conda --version 2>&1 | perl -nle 'print $& if m{conda \d+\.\d+\.\d+}'`;
  if ( echo $conda_installed_VER $miniconda_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then 
    echo " - found conda $conda_installed_VER"
    if [ -d "$ROOTDIR/thirdParty/miniconda" ]; then
      echo "conda is pointed to right environment"
      echo "creating piret environment"
      conda create --name piret --yes
      source activate piret
    else
      echo "Creating a separate conda enviroment ..."
      conda create --name piret --yes
      source activate piret
    fi
  else
    echo "Required version of conda ($miniconda_VER) was not found"
    install_miniconda
    conda create --name piret --yes
    source activate piret
  fi
else
  echo "conda was not found"
  install_miniconda
  conda create --name piret --yes
  source activate piret
fi

###############################################################################
if ( checkSystemInstallation python )
then
  python_installed_VER=`python -V 2>&1 | perl -nle 'print $& if m{Python \d+\.\d+\.\d+}'`;
  if ( echo $python_installed_VER $python3_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then 
    echo " - found python $python_installed_VER"
  else
    echo "Required version of python ($python3_VER) was not found"
    install_python
  fi
else
  echo "Python was not found"
  install_python
fi

###############################################################################
if ( checkSystemInstallation R )
    then
      R_installed_VER=`R --version 2>&1 | grep "version"| perl -nle 'print $& if m{version \d+\.\d+\.\d+}'`;
      if ( echo $R_installed_VER $R_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
      then 
        echo " - found R $R_installed_VER"
      else
        echo "Required version of R version $R_VER was not found"
        install_R
      fi
    else
      echo "R was not found"
      install_R
fi

################################################################################
if ( checkSystemInstallation hisat2 )
then
  hisat2_installed_VER=`hisat2 --version 2>&1 | grep 'hisat2-align-s version' | perl -nle 'print $& if m{version \d+\.\d+\.\d+}'`;
  if ( echo $hisat2_installed_VER $hisat2_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then 
    echo " - found hisat2 $hisat2_installed_VER"
  else
    echo "Required version of hisat2 was not found"
    install_hisat2
  fi
else
  echo "hisat2 was not found"
  install_hisat2
fi

################################################################################
if ( checkSystemInstallation stringtie )
then
  stringtie_installed_VER=`stringtie --version 2>&1 | perl -nle 'print $& if m{\d\.\d\.\d}'`;
  if ( echo $stringtie_installed_VER $stringtie_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then 
    echo " - found stringtie $stringtie_installed_VER"
  else
    echo "Required version of stringtie was not found"
    install_stringtie
  fi
else
  echo "stringtie was not found"
  install_stringtie
fi
################################################################################
# check if required bioconductor R packages are installed
################################################################################


# this is required for edgeR dependency installations (Rcpp)
conda install --yes gxx_linux-64 -n piret
conda install --yes gfortran_linux-64 -n piret
# install optparse
Rscript --no-init-file -e "if('optparse' %in% rownames(installed.packages()) == TRUE){packageVersion('optparse');}"  | awk '{print " - found optparse "$2}'
Rscript --no-init-file -e "if('optparse' %in% rownames(installed.packages()) == FALSE){install.packages('optparse',repos='https://cran.r-project.org')}";
# install dplyr
Rscript --no-init-file -e "if('dplyr' %in% rownames(installed.packages()) == TRUE){packageVersion('dplyr');}"  | awk '{print " - found dplyr "$2}'
Rscript --no-init-file -e "if('dplyr' %in% rownames(installed.packages()) == FALSE){install.packages('dplyr',repos='https://cran.r-project.org')}";
# install R ggplot2 packages
Rscript --no-init-file -e "if('ggplot2' %in% rownames(installed.packages()) == TRUE){packageVersion('ggplot2');}"  | awk '{print " - found ggplot2 "$2}'
Rscript --no-init-file -e "if('ggplot2' %in% rownames(installed.packages()) == FALSE){install.packages('ggplot2',repos='https://cran.r-project.org')}";
# install R reshape2 packages
Rscript --no-init-file -e "if('reshape2' %in% rownames(installed.packages()) == TRUE){packageVersion('reshape2');}"  | awk '{print " - found reshape2 "$2}'
Rscript --no-init-file -e "if('reshape2' %in% rownames(installed.packages()) == FALSE){install.packages('reshape2',repos='https://cran.r-project.org')}";
# install R pheatmap packages
Rscript --no-init-file -e "if('pheatmap' %in% rownames(installed.packages()) == TRUE){packageVersion('pheatmap');}"  | awk '{print " - found pheatmap "$2}'
Rscript --no-init-file -e "if('pheatmap' %in% rownames(installed.packages()) == FALSE){install.packages('pheatmap',repos='https://cran.r-project.org')}";
# install R edgeR packages
Rscript --no-init-file -e "if('edgeR' %in% rownames(installed.packages()) == TRUE){packageVersion('edgeR');}"  | awk '{print " - found edgeR "$2}'
Rscript --no-init-file -e "if('edgeR' %in% rownames(installed.packages()) == FALSE){source('https://bioconductor.org/biocLite.R');biocLite('edgeR')}";
# install R deseq2 packages
Rscript --no-init-file -e "if('DESeq2' %in% rownames(installed.packages()) == TRUE){packageVersion('DESeq2');}"  | awk '{print " - found DESeq2 "$2}'
Rscript --no-init-file -e "if('DESeq2' %in% rownames(installed.packages()) == FALSE){source('https://bioconductor.org/biocLite.R');biocLite('DESeq2')}";
# install R pathview package
Rscript --no-init-file -e "if('pathview' %in% rownames(installed.packages()) == TRUE){packageVersion('pathview');}"  | awk '{print " - found pathview "$2}'
Rscript --no-init-file -e "if('pathview' %in% rownames(installed.packages()) == FALSE){source('https://bioconductor.org/biocLite.R');biocLite('pathview')}";
# install R gage package
Rscript --no-init-file -e "if('gage' %in% rownames(installed.packages()) == TRUE){packageVersion('gage');}"  | awk '{print " - found gage "$2}'
Rscript --no-init-file -e "if('gage' %in% rownames(installed.packages()) == FALSE){source('https://bioconductor.org/biocLite.R');biocLite('gage')}";

################################################################################
if ( checkSystemInstallation featureCounts )
then
  echo " - found featureCounts"
else
  echo "featureCounts was not found"
  install_featureCounts
fi

################################################################################
if ( checkSystemInstallation samtools )
then
  samtools_installed_VER=`samtools 2>&1| grep 'Version'|perl -nle 'print $& if m{\d+\.\d+.\d+}'`;
    if ( echo $samtools_installed_VER $samtools_VER| awk '{if($2>=$3) exit 0; else exit 1}' )
    then
    echo " - found samtools $samtools_installed_VER"
    else
    echo "Required version of samtools $samtools_VER was not found"
    install_samtools
    fi
else
  echo "samtools was not found"
  install_samtools
fi
################################################################################
if ( checkSystemInstallation bamtools )
then
  bamtools_installed_VER=`bamtools -v 2>&1| grep 'bamtools'|perl -nle 'print $& if m{\d+\.\d+.\d+}'`;
    if ( echo $bamtools_installed_VER $bamtools_VER| awk '{if($2>=$3) exit 0; else exit 1}' )
    then
    echo " - found bamtools $bamtools_installed_VER"
    else
    echo "Required version of bamtools $bamtools_VER was not found"
    install_bamtools
    fi
else
  echo "bamtools was not found"
  install_bamtools
fi


################################################################################
if ( checkSystemInstallation gffread )
then
 echo "gffread is found"
else
 echo "gffread is not found"
 install_gffread
fi

################################################################################
if ( checkSystemInstallation FaQCs )
then
 echo "FaQCs is found"
else
 echo "FaQCs is not found"
 install_FaQCs
fi

################################################################################
#                        Python Modules
################################################################################
if ( checkPythonModule pandas)
  then
  python_pandas_installed_VER=`python -c "import pandas; print pandas.__version__" | perl -nle 'print $& if m{\d+\.\d+\.\d+}'`
  if (echo $python_pandas_installed_VER $python_pandas_VER | awk '{if($1>=$2) exit 0; else exit 1}' )
  then
    echo " - found Python module pandas $python_pandas_installed_VER"
  else
    echo "Required version of pandas $python_pandas_VER was not found" 
    install_python_pandas
  fi
else
    echo "pandas was not found"
    install_python_pandas
fi

################################################################################
if ( checkPythonModule luigi)
  then
    echo " - found Python module luigi"
  else
    install_python_luigi
fi

################################################################################

if ( checkPythonModule plumbum)
  then
    echo " - found Python module plumbum"
else
    echo "plumbum was not found"
    install_python_plumbum
fi

################################################################################
if ( checkPythonModule Bio)
  then
  python_Bio_installed_VER=`python -c "import Bio; print Bio.__version__" | perl -nle 'print $& if m{\d\.\d+}'`
  if (echo $python_Bio_installed_VER $python_Bio_VER | awk '{if($1>=$2) exit 0; else exit 1}' )
  then
    echo " - found Python module Bio $python_Bio_installed_VER"
  else
    echo "Required version of Bio $python_Bio_VER was not found" 
    install_python_Bio
  fi
else
    echo "Bio was not found"
    install_python_Bio
fi


#installing some lightweight developer tools
echo "Installing some developer tools codecov, coverage, pytest-cov"
conda install --yes -c conda-forge codecov -n piret
conda install --yes -c conda-forge coverage -n piret
conda install --yes -c conda-forge pytest-cov -n piret

echo "
All done! Please Restart the Terminal Session.
Run
runPiReT -h
for usage.
Read the README for more information!
To run a test data set
sh tests/test_pipeline_linux.sh
Thanks!
	"
