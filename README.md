## dce_pseudo_test_retest

This package contains source code for the following paper:

*Yannick Bliesener, R. Marc Lebel, Jay Acharya, Richard Frayne, Krishna S. Nayak, Pseudo Test-Retest Evaluation of millimeter-resolution whole-brain Dynamic Contrast-enhanced MRI in Patients with High-Grade Glioma, Radiology, under review.*

Code is provided in Matlab programming language.

## Example Usage

* Please modify `GRtools/GRoptset.m` and `TKESTIMATION_optset.m` to your data. 
* Run `demo_recon.m`.

## Dependencies

* This software packages uses the Orchestra SDK for Matlab to load GE raw data (p-files). The Orchestra SDK is copyright GE Healthcare (https://collaborate.mr.gehealthcare.com/).
* (optional) Fitting of tracer kinetic models can be accelerated with GPUs. This option requires the GPUfit toolbox: https://github.com/ironictoo/Gpufit

## Citation

If you use it, please refer to our work.

```
@article{author = {Bliesener, Yannick and Lebel, R. Marc and Acharya, Jay and Frayne, Richard, and Nayak, Krishna S.},
title = {Pseudo Test-Retest Evaluation of millimeter-resolution whole-brain Dynamic Contrast-enhanced MRI in Patients with High-Grade Glioma},
}
```

## Contact

Please contact bliesene@usc.edu if you have any questions.

Yannick Bliesener, University of Southern California, MREL (Magnetic Resonance Engineering Laboratory, PI: Krishna S. Nayak, https://mrel.usc.edu/) 2020.
