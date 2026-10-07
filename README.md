# Generalized Learning in MNNs

The following repository includes the data and code used in a study on generalized learning and elliptical interpolation in mechanical neural networks (MNNs).

## MNNSimulation

The code used to simulate MNNs using the finite truss method assigning local stiffness matrices for each beam element. It is an iteration of the simulation tool also available on another repository[https://github.com/PietroSainaghi/MNNsimulation].

The file `Main_MNN` is the main code to be run for simulating an MNN of desired width and depth, where target behaviors, beam properties, and optimization settings can be controlled.

The behavior sets used in the optimization study are found in `SavedBehaviors` under the `generalizedbehavior` naming scheme.

# GeneralizationPostProcessing

The code to compare each simulation result, and the data used for publication are found here.

The file `Processgeneralized` comparutes overall and instantaneous mean squared error for training and testing datasets.

The files in the `Results` folder do not include all structures that are generally exported by `Main_MNN` due to file size constraints on GitHub. This includes the history of stiffness combinations every optimizer iteration.
