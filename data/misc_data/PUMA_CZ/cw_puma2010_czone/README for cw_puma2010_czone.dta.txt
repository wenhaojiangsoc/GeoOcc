README for crosswalk file cw_puma2010_czone
David Dorn


Notes: 
Create a variable "puma2010" whose first two digits are the state FIPS code and whose last five digits are the PUMA code within a state. Use the Stata command "joinby" to combine the Census microdata with the crosswalk file using the variable "puma2010". The variable "afactor" in the crosswalk indicates which fraction of a PUMA's population maps to a given CZ. To analyze weighted Census data, the person weight from the Census needs to be multiplied with "afactor".

Please cite as source for this file:
David Autor, David Dorn and Gordon Hanson. "When Work Disappears: Manufacturing Decline and the Falling Marriage-Market Value of Young Men." American Economic Review: Insights, 1(2): 161-178, 2019.
