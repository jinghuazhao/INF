# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# Uncomment the following line if you don't like systemctl's auto-paging feature:
# export SYSTEMD_PAGER=

# User specific aliases and functions

alias h='history'
alias ls='ls -hF'
alias df='df -h'
alias du='du -h'
export PATH=$HOME/bin:/data/jinhua/ImageMagick-7.0.8-22/bin:$PATH

function R.3.5.3()
{
  export R_LIBS=/data/$USER/R:$HOME/R:/services/tools/R/3.5.3-ICC-MKL/lib64/R/library
  module load intel/redist/2019 intel/perflibs/64/2019 gcc/5.4.0 R/3.5.3-ICC-MKL
  source /data/jinhua/parallel-20190222/bin/env_parallel.bash
  alias R='/services/tools/R/3.5.3-ICC-MKL/bin/R -q $@'
}

function R.3.3.1
{
  export R_LIBS=/home/$USER/R:/services/tools/R/3.3.1-ICC-MKL/lib64/R/library
  module load intel/redist/2019 intel/perflibs/64/2019 gcc/5.4.0 R/3.3.1-ICC-MKL
}
