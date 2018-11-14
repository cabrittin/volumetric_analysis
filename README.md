# volumetric_analysis
Companion code for generating the figures from 

Brittin, Cook, Hall, Emmons and Cohen. 'Volumetric reconstruction of 
*Caenorhabditis elegans* nerve ring supports combinatorial CAM expression 
model for synaptic specificity'. (2018) Under review. 

If you use any part of this code please site this paper. 

## Installation
Clone or download respository. Make sure to maintain the relative paths.

You will need to install and set up the appropriate mysql databases.

### Prerequisites
Python3 will need to be installed along with the following python packages 
```
numpy v1.14.3
scipy v1.0.1
matplotlib v2.2.2
matplotlib-venn 0.11.5
mysqlclient 1.3.12
python-igraph 0.7.1.post6
networkx 2.1

```

## Usage
To generate figures use the following 
```
python paper_figures figNum -o /path/to/outputfigure
```
where figNum is the figure number. For example, for Figure 3a use f3a. For Figure S3a use fs3a. 

## Author

* **Christopher Brittin** - *Initial work* - [cabrittin](https://github.com/cabrittin)


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
