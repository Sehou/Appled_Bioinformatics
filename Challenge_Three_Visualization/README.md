# Admixture_Q_estimate_with_altitude

Admixture_Q_estimate_with_altitude is a Python script for plotting and visualizing admixture Q estimates along with the altitude. The output plot of this script gives an idea of the admixture among populations, which should be accounted in the functional diversity study within species.

The script first load both Q file and meta data file. Afterwards, the two data frames are merged and samples (row) containing NAN values are dropped. In the next step, it performs a K-mean clustering with k equal to the number of column in the Q file to separate the sample into k clusters.  Using this cluster assignment, the sample are sort per cluster. The sorted data is finally then plot.

## Installation

This script does not require any installation. However, its repository must be cloned before usage. This can be done as follow:

```bash
git clone https://github.com/Sehou/Applied_Bioinformatics.git
```

## Usage

```python
python Admixture_Q_estimate_with_altitude.py <options>
```

The options -m for the meta data and -q for Q file are mandatory. All the available options can be see by running: 

```bash
python Admixture_Q_estimate_with_altitude.py --help
```

## Output

Once the program finish executing, the plot pop-ups. However, the plot is also saved into the specified output direction if any or in the current working directory under the name Admixture_altitude.png.