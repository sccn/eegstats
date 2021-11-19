This plugin compute frequency band power, alpha peak frequency, and alpha asymetry. It uses the [restingIAF](https://github.com/corcorana/restingIAF) MATLAB code for some of the computation.

# Do not download Github zip file

Downloading the zip file, you will be missing the restingIAF code dependency. Instead download a released version or check out the repository with dependencies.

# Checkout repository

Make sure to copy submodule when you clone the repository

```
git clone --recurse-submodules https://github.com/arnodelorme/eegstats.git
```

# Graphic interface

The plugin may be used from the command line or from its GUI. See the [pop_eegstats.m](https://github.com/arnodelorme/eegstats/blob/master/pop_eegstats.m) header for more information.

![](eegstats_gui.png)

# Version history

v1.0 - initial version
