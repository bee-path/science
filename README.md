Beepath Science Module
========================================================================

Copyright 2013 Oleguer Sagarra and Mario Gutiérrez Roig. All rights reserved. Code under License GPLv3.
______________________________________________________________________________________


## Contents of the package

This python package contains the class definition and functions used in the analysis of the data from the Bee-path experiment.
It may be used to analyze point-like GPS general mobility data and to aggregate such points into displacements and paused states.
For more details check out the webpage of the project  [Bee-path](http://bee-path.net/?lang=en) and the related publications.

The module is fully documented in a *quasi-standard*, *pythonic* way (but could be improved).

## References 

[1] Gutiérrez-Roig, M et altri
	[Active and reactive behaviour in human mobility: the influence of attraction points on pedestrians](http://arxiv.org/abs/1511.03604)

[2] Sagarra, O et altri
    [Citizen Science Practices for Computational Social Science Research: The Conceptualization of Pop-Up Experiments](http://journal.frontiersin.org/article/10.3389/fphy.2015.00093/full)
## Installation 

Being a usual python package, to install just type:

```
	$ python setup.py install
```	

And to import,

```
    $ import beepath_science 
```

and use as needed.

Note that some dependencies are needed. Mainly:

- [Numpy](http://www.numpy.org/)
- [Pyproj](https://github.com/jswhit/pyproj)
- [Shapely](http://toblerity.org/shapely/)


## License

Copyright 2013 Oleguer Sagarra (osagarra@ub.edu) and Mario Gutiérrez Roig (mariogutierrezroig@gmail.com)
Code under License GPLv3.
All rights reserved. 


## Disclaimer

The software is provided as it is, with absolutely NO WARRANTY. If you detect any bugs, please let us know.





