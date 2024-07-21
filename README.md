# mexeig - LAPACK's Eigensolvers Interface for MATLAB

## Installation

1. Download all files.
2. Add the path of "mexeig-main" folder to MATLAB's search path.
3. Execute the following command to compile the mex files:
```
mexeig('compile');
```

## Usage

```
[Outputs,...] = mexeig(FunctionName, Matrix, Options...);

%
% FunctionName: 'syev', 'syevd', 'syevr', 'syevx', 'geev', 'geevx',
%                'sygv', 'sygvd', 'sygvx', 'ggev', 'ggev3', 'ggevx'
% Matrix      : single or double matrix
% Option1     : single or double positive definite matrix (for generalized EVD),
%               'matrix', 'vector' (for specifying the output form),
%               'N', 'P', 'S', 'B' (for BALANC for 'geevx' and 'ggevx')
%               or abstol
% Option2     : 'matrix', 'vector' (for specifying the output form),
%               'N', 'P', 'S', 'B' (for BALANC for 'geevx' and 'ggevx')
%               or abstol
% Option3     : 'matrix', 'vector' (for specifying the output form), or abstol
%
```
 
