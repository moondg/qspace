# QSpace v4.1 install guide

See [QSpace](https://bitbucket.org/qspace4u/qspace-v4-pub) for more information.

> [!NOTE]
> This installer assumes using [Mathworks Matlab image](https://hub.docker.com/r/mathworks/matlab)

Automated installation can be done by

1. Copy `installer.sh` to computer.
2. Modify `QSPACE_DIR` in `installer.sh`.
3. Run command below
    ```bash
    bash installer.sh
    . ~/.bashrc
    mlsetup
    cd $MYMATLAB/Source
    make -B all
    ```

# Description
+ Default installation path is `$HOME/Documents/MATLAB/qspace`, which can be changed by modifying `QSPACE_DIR` in `installer.sh`.
+ `LMA` and `RC_STORE`, `RC_SYNC` are set to be in `QSPACE_DIR/build`.
