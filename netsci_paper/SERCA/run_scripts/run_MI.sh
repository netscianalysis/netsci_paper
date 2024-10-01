#!/bin/bash
#
# Bashscript that calls the run_simple_MI.py script
# Created by Marcus Hock 
# Updated 05/23/2024
# Not necessary to run MI calcs, but can be useful to modify variables
# specifically if running batch jobs on a super computer. 
#
#
# ENVIRONMENT SETUP Variables 
# These will need to be changed for different running locations (Mojave, TSCC, etc)
CONDA_PATH=/home/marcus/anaconda3
NETSCI_ENV=netsci
PATH_TO_MI=~/Documents/SERCA/serca-MI/

# JOB SETUP Variables 
START=0 # Starting Frame Value 
STOP=999 # Stopping value (do not exceed number of frames. For 1000 frame simulation enter 999)
STRIDE=5 # Stride for analysis. Likely do not need every frame 

# Path to PDB 'topology' file 
TOPOLOGY_FILE=INTSERT_TOPOLOGY_FILEPATH

# Path to DCD 'trajectory' file 
TRAJECTORY_FILE=INTSERT_TRAJECTORY_FILEPATH

# Printing working directory
echo $(pwd)


# Setup conda
# This path will need to be modified on different systems
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('$CONDA_PATH/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "$CONDA_PATH/etc/profile.d/conda.sh" ]; then
        . "$CONDA_PATH/etc/profile.d/conda.sh"
    else
        export PATH="$CONDA_PATH/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

# May have a different conda name for the environment
conda activate $NETSCI_ENV


SCRIPTS_DIR=${PATH_TO_MI}'analysis_scripts'


echo Running with python located at...
echo $(which python)
echo Running the command:
echo python ${SCRIPTS_DIR}/run_simple_MI.py -pdb $TOPOLOGY_FILE -dcd $TRAJECTORY_FILE -start $START -stop $STOP -stride $STRIDE 
python ${SCRIPTS_DIR}/run_simple_MI.py -pdb $TOPOLOGY_FILE -dcd $TRAJECTORY_FILE -start $START -stop $STOP -stride $STRIDE 



