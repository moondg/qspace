#!bin/bash

# update ubuntu and install requirements
sudo apt update
sudo apt upgrade -y
sudo apt install git libgmp-dev libmpfr-dev -y

# specify folders
QSPACE_DIR=$HOME/Documents/MATLAB/qspace-v4-pub
echo "Default installation path is $QSPACE_DIR"

# prepare before compile
cd $(dirname $QSPACE_DIR)
git clone https://bitbucket.org/qspace4u/qspace-v4-pub.git
mkdir $QSPACE_DIR/build

user_config() {
  cat <<'EOF' > $QSPACE_DIR/system/matlab_setup_user.sh
export MATLAB_ROOT=/opt/matlab/R2025a
export MSLOTS=$(grep -c processor /proc/cpuinfo)
EOF
}

bash_config() {
  cat <<EOF > $HOME/.qspace_config
# =============================================================================
# =                             QSpace v4.0 setup                             =
# =============================================================================
# Assumed: Matlab R2025a container image based on Ubuntu
# 1) From qspace user-guide.pdf p.96, define alias
MYMATLAB=$QSPACE_DIR
EOF
  cat <<'EOF' >> $HOME/.qspace_config
alias mlsetup="source $MYMATLAB/system/matlab_setup.sh"
# 2) From qspace user-guide.pdf p.96-97, add path
export PATH=$PATH:$MYMATLAB/Source
export PATH=$PATH:$MYMATLAB/system
# `matlab` is in PATH "/usr/local/bin", however, `mex` is not.
# Actually, both `matlab` and `mex` are in "/opt/matlab/R2025a/bin".
# `matlab` in the "/usr/local/bin" is just symbolic link to actual path
# Therefore, add actual path to PATH
export PATH=/opt/matlab/R2025a/bin:$PATH
# 3) Default path of LMA, RC_STORE, RC_SYNC is not proper in this case
export LMA=$MYMATLAB/build
export RC_STORE=$LMA/RCStore
export RC_SYNC=$LMA/RCSync
EOF
}

# update config
user_config
bash_config

. ~/.bashrc # necessary
if [ -z "$QSPACE_INSTALLER_HISTORY" ]; then
  echo "export QSPACE_INSTALLER_HISTORY=1" >> ~/.bashrc
  echo "source ~/.qspace_config" >> ~/.bashrc
fi