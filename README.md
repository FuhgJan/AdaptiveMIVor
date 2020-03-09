# AdaptiveMIVor
The Monte Carlo-intersite Voronoi (MiVor) adaptive scheme is an adaptive sampling technqiue for ordinary Kriging. 

The following library provides a Matlab implementation of the MiVor algorithm. For tests of the method three different test cases with different degrees of complexity are provided within this framework.


## Getting Started

We encourage those who are interested in using this code to run the main file and pick a test case.

### Prerequisites

Matlab version R2017a or higher.

## Examples 

The three following working examples similarly to what is found in the paper are included in this code for which the classification boundary limit is set to 0.0 for each case.

---

<p align="center">
  <img align="middle" src="./docs/TestCase1_Image.png" alt="Example 1" width="250" height="250" />
  <img align="middle" src="./docs/TestCase2_Image.png" alt="Example 1" width="250" height="250" />
  <img align="middle" src="./docs/TestCase3_Image.png" alt="Example 1" width="250" height="250" />
</p>

---

### Example 1
The following Gifs show an exemplary sampling process for example 1. 
<p align="center">
 <kbd><img align="middle" src="./docs/TestCase1_MetaVor.gif" alt="ODE Demo" width="300" height="300" border="50"  /></kbd>
  <kbd><img align="middle" src="./docs/TestCase1_Vor.gif" alt="ODE Demo" width="400" height="300" border="50" /></kbd>

</p>

### Adaptive sampling process Examples 2 and 3
The following Gifs show an exemplary sampling process for example 2 and 3 respectively. 
<p align="center">
 <kbd><img align="middle" src="./docs/TestCase2_Meta.gif" alt="ODE Demo" width="350" height="350" border="50" /></kbd>
  <kbd><img align="middle" src="./docs/TestCase3_Meta.gif" alt="ODE Demo" width="350" height="350" border="50" /></kbd>
</p>


---

# References

Please cite this code with:

Fuhg, Jan N., and Amelie Fau. "An innovative adaptive kriging approach for efficient binary classification of mechanical problems." arXiv preprint arXiv:1907.01490 (2019).



