## Reaction-Diffusion-Strain (RDE) model for cranial vault development


### BC

Matlab files to calculate boundary condition and sample data.

```markdown
1. disp_field.m //estimate displacement field between surfaces in .suf file format
2. interpolation.m //estimate displacement vecotrs on each queried points that coincide with computational mesh elements
```


### RDE_DM

Source code of RDE model with Dynamic Mesh. 

Execute commands below after loading foam-extend-3.1

```markdown
wclean
wmake
```

Solver named RDE_DM will be registered in $(FOAM_USER_APPBIN)


### cranial_vault_example

Example of simulation of cranial vault development using RDE model.

```markdown
./clean
./run
```


### Contact

Chanyoung Lee (lcys0914@gmail.com)
