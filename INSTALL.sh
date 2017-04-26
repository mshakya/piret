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
export "PATH=$PATH:$ROOTDIR/bin/"
printenv
if [[ "$OSTYPE" == "darwin"* ]]
then
{
export PERL5LIB="$ROOTDIR/ext/lib/perl5:$ROOTDIR/ext/lib/perl5/darwin-thread-multi-2level:$PERL5LIB"
}
else
#TODO: update this specific to unix
{  
export PERL5LIB="$ROOTDIR/ext/lib/perl5:$ROOTDIR/lib/perl5/darwin-thread-multi-2level/:$PERL5LIB"
}
fi
printenv

# Add pythonpath
# export PYTHONPATH="$ROOTDIR/thirdParty/miniconda/lib/python2.7/site-packages/:$PYTHONPATH"

#Add R path
export R_LIBS="$ROOTDIR/ext/lib/R:$R_LIBS:$R_LIBS_USER"

# Minimum Required versions of dependencies
bowtie2_VER=2.2.8
bwa_VER=0.7.15
cpanm_VER=1.7039
miniconda_VER=4.2.12
samtools_VER=1.3.1
jellyfish_VER=2.2.6
bedtools_VER=2.26.0
R_VER=3.3.1
hisat2_VER=2.0.5
htseq_VER=0.6.1
jbrowse_VER=1.12.1

# minimum required version of Scripting languages
perl5_VER=5.8.0
python2_VER=2.7.12

#minimum required version of Perl modules
perl_String_Approx_VER=3.27
perl_Parllel_ForkManager_VER=1.17

#minimum required version of Python modules
python_numpy_VER=1.1.12
python_matplotlib_VER=1.5.3

#minimum required version of Python modules
R_edgeR_VER=3.14.0
R_DESeq2_VER=1.12.4

################################################################################
#                           Installation recipes
################################################################################

install_python()
{
echo "--------------------------------------------------------------------------
                           installing python v$python2_VER
--------------------------------------------------------------------------------
"
conda install --yes python=$python2_VER -p $ROOTDIR/thirdParty/miniconda
ln -sf $ROOTDIR/thirdParty/miniconda/bin/python $ROOTDIR/bin/python

echo "
------------------------------------------------------------------------------
                           python v$python2_VER installed
------------------------------------------------------------------------------
"
}


install_hisat2()
{
echo "--------------------------------------------------------------------------
                           installing hisat2 v$hisat2_VER
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda hisat2=$hisat2_VER -p $ROOTDIR/thirdParty/miniconda
ln -sf $ROOTDIR/thirdParty/miniconda/bin/hisat2 $ROOTDIR/bin/hisat2
ln -sf $ROOTDIR/thirdParty/miniconda/bin/hisat2-build $ROOTDIR/bin/hisat2-build
echo "
------------------------------------------------------------------------------
                           hisat2 v$hisat2_VER installed
------------------------------------------------------------------------------
"
}

install_jellyfish()
{
echo "--------------------------------------------------------------------------
                           installing jellyfish v$jellyfish_VER
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda jellyfish=$jellyfish_VER -p $ROOTDIR/thirdParty/miniconda
ln -sf $ROOTDIR/thirdParty/miniconda/bin/jellyfish $ROOTDIR/bin/jellyfish
echo "
------------------------------------------------------------------------------
                           jellyfish v$jellyfish_VER installed
------------------------------------------------------------------------------
"
}

install_perl_Parallel_ForkManager()
{
echo "--------------------------------------------------------------------------
  Installing Perl Module Parallel-ForkManager v$perl_Parllel_ForkManager_VER
--------------------------------------------------------------------------------
"

cpanm Parallel::ForkManager@$perl_Parllel_ForkManager_VER -l $ROOTDIR/ext
# conda install --yes -c bioconda perl-parallel-forkmanager=$perl_Parllel_ForkManager_VER
echo "
--------------------------------------------------------------------------------
      Parallel-ForkManager-$perl_Parllel_ForkManager_VER Installed
--------------------------------------------------------------------------------
"
}

install_bwa()
{
echo "--------------------------------------------------------------------------
                           Downloading bwa v$bwa_VER
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda bwa=$bwa_VER -p $ROOTDIR/thirdParty/miniconda
ln -sf $ROOTDIR/thirdParty/miniconda/bin/bwa $ROOTDIR/bin/bwa
echo "
--------------------------------------------------------------------------------
                           bwa v$bwa_VER installed
--------------------------------------------------------------------------------
"
}

install_htseq()
{
echo "------------------------------------------------------------------------------
                           downloading htseq v $htseq_VER
------------------------------------------------------------------------------
"
conda install --yes -c bioconda htseq=$htseq_VER -p $ROOTDIR/thirdParty/miniconda
ln -sf $ROOTDIR/thirdParty/miniconda/bin/htseq-count $ROOTDIR/bin/htseq-count
echo "
------------------------------------------------------------------------------
                           htseq v $htseq_VER installed
------------------------------------------------------------------------------
"
}

install_samtools()
{
echo "--------------------------------------------------------------------------
                           Downloading samtools v $samtools_VER
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda samtools=$samtools_VER -p $ROOTDIR/thirdParty/miniconda
ln -sf $ROOTDIR/thirdParty/miniconda/bin/samtools $ROOTDIR/bin/samtools
echo "
--------------------------------------------------------------------------------
                           samtools v $samtools_VER installed
--------------------------------------------------------------------------------
"
}

install_cpanm()
{
echo "--------------------------------------------------------------------------
                           Installing cpanm
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda perl-app-cpanminus=$cpanm_VER -p $ROOTDIR/thirdParty/miniconda
ln -sf $ROOTDIR/thirdParty/miniconda/bin/cpanm $ROOTDIR/bin/cpanm
echo "
--------------------------------------------------------------------------------
                           cpanm installed
--------------------------------------------------------------------------------
"
}


install_bedtools()
{
echo "--------------------------------------------------------------------------
                           Compiling bedtools v $bedtools_VER
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda bedtools=$bedtools_VER -p $ROOTDIR/thirdParty/miniconda
ln -sf $ROOTDIR/thirdParty/miniconda/bin/genomeCoverageBed $ROOTDIR/bin/genomeCoverageBed
ln -sf $ROOTDIR/thirdParty/miniconda/bin/bedtools $ROOTDIR/bin/bedtools
echo "
--------------------------------------------------------------------------------
                           bedtools v $bedtools_VER compiled
--------------------------------------------------------------------------------
"
}

install_R()
{
echo "--------------------------------------------------------------------------
                           Installing R v $R_VER
--------------------------------------------------------------------------------
"
conda install --yes -c r r-base=$R_VER -p $ROOTDIR/thirdParty/miniconda
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

  curl -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda2-4.2.12-MacOSX-x86_64.sh
  chmod +x miniconda.sh
  ./miniconda.sh -b -p $ROOTDIR/thirdParty/miniconda -f
  # export PATH=$ROOTDIR/thirdParty/miniconda/bin:$PATH
  ln -sf $ROOTDIR/thirdParty/miniconda/bin/conda $ROOTDIR/bin/conda
  ln -sf $ROOTDIR/thirdParty/miniconda/bin/pip $ROOTDIR/bin/pip

}
else
{  

  wget https://repo.continuum.io/miniconda/Miniconda2-4.2.12-Linux-x86_64.sh -O miniconda.sh
  chmod +x miniconda.sh
  ./miniconda.sh -b -p $ROOTDIR/thirdParty/miniconda -f
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



install_perl_string_approx()
{
echo "------------------------------------------------------------------------------
                 Installing Perl Module String-Approx-3.27
------------------------------------------------------------------------------
"
#TODO: Figure out how to download this using curl
# curl -k -l http://search.cpan.org/CPAN/authors/id/J/JH/JHI/String-Approx-3.27.tar.gz -o String-Approx-3.27.tar.gz

cd $ROOTDIR/thirdParty
tar xvzf String-Approx-3.27.tar.gz
cd String-Approx-3.27
perl Makefile.PL 
make
mkdir -p $ROOTDIR/ext/lib/perl5/auto
mkdir -p $ROOTDIR/ext/lib/perl5/auto/String
mkdir -p $ROOTDIR/ext/lib/perl5/auto/String/Approx
cp -fR blib/lib/* $ROOTDIR/ext/lib/perl5
cp -fR blib/arch/auto/String/Approx/Approx.* $ROOTDIR/ext/lib/perl5/auto/String/Approx/
cd $ROOTDIR/thirdParty
echo "
------------------------------------------------------------------------------
                        String-Approx-3.27 Installed
------------------------------------------------------------------------------
"
}

install_gffread()
{
echo "--------------------------------------------------------------------------
                      installing gffread
--------------------------------------------------------------------------------
"
cd $ROOTDIR/thirdParty/miniconda/bin
git clone https://github.com/gpertea/gclib
git clone https://github.com/gpertea/gffread gffread_git
cd gffread_git
make
cp gffread ../
ln -sf $ROOTDIR/thirdParty/miniconda/bin/gffread $ROOTDIR/bin/gffread
cd $ROOTDIR/thirdParty
echo "
--------------------------------------------------------------------------------
                           gffread installed
--------------------------------------------------------------------------------
"
}

install_jbrowse()
{
echo "--------------------------------------------------------------------------
                      installing Jbrowse
--------------------------------------------------------------------------------
"
cd $ROOTDIR/thirdParty
tar xvzf JBrowse-1.11.6.tar.gz
cd JBrowse-1.11.6
./setup.sh
mkdir -p -m 775 data
cd $ROOTDIR/thirdParty
ln -sf $ROOTDIR/thirdParty/JBrowse-1.11.6 $ROOTDIR/bin/JBrowse
echo "
--------------------------------------------------------------------------------
                      Jbrowse installed
--------------------------------------------------------------------------------
"
}

install_python_numpy()
{
echo "--------------------------------------------------------------------------
                installing Python module numpy $python_numpy_VER
--------------------------------------------------------------------------------
"
conda install --yes -c anaconda numpy=$python_numpy_VER
echo "
--------------------------------------------------------------------------------
                           numpy installed
--------------------------------------------------------------------------------
"
}

install_python_matplotlib()
{
echo "--------------------------------------------------------------------------
                installing Python module matplotlib $python_matplotlib_VER
--------------------------------------------------------------------------------
"
conda install --yes -c anaconda matplotlib=$python_matplotlib_VER
echo "
--------------------------------------------------------------------------------
                           matplotlib installed
--------------------------------------------------------------------------------
"
}

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

checkPerlModule()
{
   # perl -e "use lib \"$ROOTDIR/lib/lib/perl5\"; use $1;"
   perl -e "use $1";
   return $?
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

print_tools_list()
{


   echo "Available tools for updates/re-install"
   echo -e "\nAlignment"
   for i in "${alignments_tools[@]}"
   do
       echo "* $i"
   done
   echo -e "\nUtility"
   for i in "${utility_tools[@]}"
   do
       echo "* $i"
   done
   echo -e "\nCount"
   for i in "${count_tools[@]}"
   do
       echo "* $i"
   done
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
    else
      echo "Creating a separate conda enviroment ..."
      conda create --yes -p $ROOTDIR/thirdParty/miniconda
    fi
  else
    echo "Required version of conda ($miniconda_VER) was not found"
    install_miniconda
  fi
else
  echo "conda was not found"
  install_miniconda
fi

###############################################################################
if ( checkSystemInstallation python )
then
  python_installed_VER=`python -V 2>&1 | perl -nle 'print $& if m{Python \d+\.\d+\.\d+}'`;
  if ( echo $python_installed_VER $python_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then 
    echo " - found python $python_installed_VER"
  else
    echo "Required version of python ($python2_VER) was not found"
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
# check if required bioconductor R packages are installed
################################################################################

if ( checkRpackages edgeR )
  then
  R_edgeR_installed_VER=`echo "cat(unname(installed.packages()[,3][\"edgeR\"]))" | Rscript - | sed "s/\"//g"`;
  echo $R_edgeR_installed_VER 
  if (echo $R_edgeR_installed_VER $R_edgeR_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then
  echo " - found edgeR $R_edgeR_installed_VER"
  else
    echo "Required version of edgeR $R_edgeR_VER was not found"
    install_R_edgeR
  fi
else
  echo "edgeR is not found"
    install_R_edgeR
fi

#------------------------------------------------------------------------------#

if ( checkRpackages DESeq2 )
  then
  R_DESeq2_installed_VER=`echo "cat(unname(installed.packages()[,3][\"DESeq2\"]))" | Rscript - | sed 's/\"//g'`;
  echo $R_DESeq2_installed_VER 
  if (echo $R_DESeq2_installed_VER $R_DESeq2_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then
  echo " - found DESeq2 $R_DESeq2_installed_VER"
  else
    echo "Required version of DESeq2 $R_DESeq2_VER was not found"
    install_R_DESeq2
  fi
else
  echo "DESeq2 is not found"
    install_R_DESeq2
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
if ( checkSystemInstallation htseq-count )
then
  htseq_installed_VER=`htseq-count -h | grep 'version'| perl -nle 'print $& if m{version \d+\.\d+\.\d+}'`
  if ( echo $htseq_installed_VER $htseq_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then
    echo " - found htseq $htseq_installed_VER"
  else
    echo "Required version of htseq was not found"
    install_htseq
  fi
else
  echo "htseq was not found"
  install_htseq
fi
################################################################################
if ( checkSystemInstallation jellyfish )
then
  jellyfish_installed_VER=`jellyfish --version | perl -nle 'print $& if m{jellyfish \d+\.\d+\.\d+}'`
  if ( echo $jellyfish_installed_VER $jellyfish_VER | awk '{if($2>=$3) exit 0; else exit 1}')
  then
    echo " - found $jellyfish_installed_VER"
  else
    echo "Required version of jellyfish was not found"
  fi
else
  echo "jellyfish was not found"
  install_jellyfish
fi
################################################################################
if ( checkSystemInstallation bwa )
then
bwa_installed_VER=`bwa 2>&1| grep 'Version'  | perl -nle 'print $& if m{Version: \d+\.\d+\.\d+}'`;
  if  ( echo $bwa_installed_VER $bwa_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then
    echo " - found BWA $bwa_installed_VER"
  else
    install_bwa
  fi
else
  echo "bwa is not found"
  install_bwa
fi
################################################################################
if ( checkSystemInstallation samtools )
then
  samtools_installed_VER=`samtools 2>&1| grep 'Version'|perl -nle 'print $& if m{Version: \d+\.\d+.\d+}'`;
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
if ( checkSystemInstallation bedtools )
then
  bedtools_installed_VER=`bedtools --version | perl -nle 'print $& if m{\d+\.\d+\.\d+}'`;
  if ( echo $bedtools_installed_VER $bedtools_VER | awk '{if($1>=$2) exit 0; else exit 1}' )
  then
    echo " - found bedtools $bedtools_installed_VER"
  else
    echo "Required version of bedtools $bedtools_VER was not found"
    install_bedtools
    fi
else
  echo "bedtools is not found"
  install_bedtools
fi


################################################################################
#TODO: add a way to check version here as well
if [ -x $ROOTDIR/bin/JBrowse/bin/prepare-refseqs.pl ]
then
  echo "JBrowse is found"
else
  echo "JBrowse is not found"
  install_jbrowse
fi

################################################################################
if ( checkSystemInstallation cpanm )
then
  cpanm_installed_VER=`cpanm -V 2>&1| head -n 1 | grep 'version' | perl -nle 'print $& if m{version \d+\.\d+}'`;
  if  ( echo $cpanm_installed_VER $cpanm_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then 
    echo " - found cpanm $cpanm_installed_VER"
  else
  echo "Required version of cpanm was not found"
  install_cpanm
  fi
else 
  echo "cpanm was not found"
  install_cpanm
fi


# commenting this out until i figure out what is this actually used for
#if ( checkSystemInstallation gffread )
#then
#  echo "gffread is found"
#else
#  echo "gffread is not found"
#  install_gffread
#fi

################################################################################
#                        Perl Modules
################################################################################

if ( checkPerlModule Parallel::ForkManager )
then
  perl_Parallel_ForkManager_installed_VER=`perl -MParallel::ForkManager -e 'print $Parallel::ForkManager::VERSION ."\n";'`
  if ( echo $perl_Parallel_ForkManager_installed_VER $perl_Parallel_ForkManager_VER | awk '{if($1>=$2) exit 0; else exit 1}')
  then
    echo " - found Perl module Parallel::ForkManager $perl_Parallel_ForkManager_installed_VER"
  else
    echo "Required version of Parallel::ForkManager $perl_Parallel_ForkManager_VER was not found"
    install_perl_Parallel_ForkManager
  fi
else
  echo "Perl Parallel::ForkManager is not found"
  install_perl_Parallel_ForkManager
fi

################################################################################

if ( checkPerlModule String::Approx )
then
  perl_String_Approx_installed_VER=`cpan -D String::Approx 2>&1 | grep 'Installed' | perl -nle 'print $& if m{Installed: \d+\.\d+}'`
  if ( echo $perl_String_Approx_installed_VER $perl_String_Approx_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then
    echo " - found Perl module String::Approx $perl_String_Approx_installed_VER"
  else
    echo "Required version of String::Approx $perl_String_Approx_VER was not found"
    install_perl_string_approx
  fi
else
  echo "Perl String::Approx was not found"
  install_perl_string_approx
fi
################################################################################
#                        Python Modules
################################################################################
if ( checkPythonModule numpy)
  then
  python_numpy_installed_VER=`python -c "import numpy; print numpy.__version__" | perl -nle 'print $& if m{\d+\.\d+\.\d+}'`
  if (echo $python_numpy_installed_VER $python_numpy_VER | awk '{if($1>=$2) exit 0; else exit 1}' )
  then
    echo " - found Python module numpy $python_numpy_installed_VER"
  else
    echo "Required version of numpy $python_numpy_VER was not found" 
    install_python_numpy
  fi
else
    echo "numpy was not found"
    install_python_numpy
fi

################################################################################

if ( checkPythonModule matplotlib)
  then
  python_matplotlib_installed_VER=`python -c "import matplotlib; print matplotlib.__version__" | perl -nle 'print $& if m{\d+\.\d+\.\d+}'`
  if (echo $python_matplotlib_installed_VER $python_matplotlib_VER | awk '{if($1>=$2) exit 0; else exit 1}' )
  then
    echo " - found Python module matplotlib $python_matplotlib_installed_VER"
  else
    echo "Required version of matplotlib $python_matplotlib_VER was not found" 
    install_python_matplotlib
  fi
else
    echo "matplotlib was not found"
    install_python_matplotlib
fi


echo "
All done! Please Restart the Terminal Session.
Run
./runPiReT.pl
for usage.
Read the README for more information!
To run a test data set
cd test_data
sh test_pipeline_linux.sh
OR
sh test_pipeline_MacOSX.sh
Thanks!
	"
