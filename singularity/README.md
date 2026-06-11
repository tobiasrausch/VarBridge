You can build a [VarBridge](https://github.com/tobiasrausch/VarBridge) singularity container (SIF file) using

`sudo singularity build varbridge.sif varbridge.def`

Once you have built the container you can run analysis using

`singularity exec varbridge.sif varbridge --help`
