# QSpace v4.0 install guide

See [QSpace](https://bitbucket.org/qspace4u/qspace-v4-pub) for more information.

> [!NOTE]
> This installer assumes using [Mathworks Matlab image](https://hub.docker.com/r/mathworks/matlab)

Automated installation can be done by

```
bash installer.sh
. ~/.bashrc
mlsetup
cd $MYMATLAB/Source
make -B all
```

# Description
+ Default installation path is `$HOME/Documents/MATLAB/qspace-v4-pub`, which can be changed by modifying `QSPACE_DIR` in `installer.sh`.
+ `LMA` and `RC_STORE`, `RC_SYNC` are set to be in `QSPACE_DIR/build`.