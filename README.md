# Microbiota-to-Recon (m2r)
Microbiota-to-Recon (M2R) program allows to modify your model (a human genome-scale metabolic model like RECON3D) by incorporating information from selected microbiota.

## Dependencies
* `python2` or `python3`
* `cobrapy 0.18.1` [Link](https://cobrapy.readthedocs.io/)
* `os`
* `pickle`

## Requirements
The program assumes the following:
* All necessary packages are installed (see dependencies).
* Microbiota collection is downloaded (e.g., see [Link](https://www.vmh.life/#downloadview)).
* A model (e.g., RECON3D, see [Link](https://www.vmh.life/#downloadview)) is prepared.

## Running the program
In the console or in a Python IDE run the following:

`python m2r.py`

or

`python3 m2r.py`

Then, please provide all necessary information (i.e., paths to directories, values of hyperparameters). The output of the program is a new model with information from microbiota.

## Citing the code
If you find our program interesting or you use it in your work, please cite the following paper:
```
@article{m2r2020,
  title={M2R: Modifying human genome-scale metabolic models with information from microbiota,
  author={Weglarz-Tomczak, Ewelina and Tomczak, Jakub M. and Brul, Stanley},
  journal={(submitted)},
  year={2020}
}
```
