# Installation

Please follow the guide below to install CGMega and its dependent software.

## System Requirements

Operating system: Linux or are highly recommended. CGMega was developed and tested on Linux. Windows Subsystem for Linux (WSL) is not recommended.

Memory: Memory usage depends on the size of your Hi-C and PPI dataset and it may cost 

CPU: Intel(R) Xeon(R) Silver 4210R CPU @ 2.40GHz. This is NOT the least requirement.

## Software Requirements

CGMega was developed using python 3.8. The interpretation part adopted GNNExplainer. 

- python=3.8.12
- pytorch=1.9.1
- networkx=2.6.3
- neoloop=0.4.3
- transformers=4.19.4
- wandb=0.12.18(optional)
- pyg=2.0.3

## Installation

We recommend installing CGMega using Git from our Github repository through the following command:

```
git clone https://github.com/NBStarry/CGMega.git
```

```note
The software package of CGMega is still under testing. We are going to release a stable version in the near future.
```

source: `{{ page.path }}`
