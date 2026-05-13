# ibm-enkf
Test code for hybrid assimilation scheme for Individual Based Models

Running the simulation:
- Model and assimilation settings are set in runIBM.jl main() function
- To load code into Julia REPL, evaluate runIBM.jl
- To run the code, call the following functions:
    - recordtwin() runs model with minimal ensemble and records the run to use as twin.
    - dryrun() runs ensemble and computes corrections using RS method, but corrections are not applied
    - resamplerun() runs ensemble with corrections computed using the RS method
    - normrun() runs ensemble with corrections computed using the OT method
