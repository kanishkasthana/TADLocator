TADLocator
==========

The script and associated data for my rotation in the Ren lab in the Fall Quarter of 2014 at UC San Diego.

The report describing the method used can be seen <a href='https://www.dropbox.com/s/8d9i9d9v2ikkcw0/An%20Independent%20Method%20to%20Identify%20Topological%20Domain%20Boundaries.pdf?dl=0'>here </a>.
The presentation can also be seen <a href='https://www.dropbox.com/s/1hc38hf9z8c3ogx/An%20Independent%20Method%20to%20Identify%20Topological%20Domains.pptx?dl=0'> here </a>.
You don't need to download dropbox to view this file, just close the signup window to view the file.

<b> Update: getBins function addded </b>
This function allows the user to get the Bin numbers for Domain boundaries in a MATRIX file for a specific stringency level.
Here stringency is defined as the number of standard deviation from the mean of the distribution of rates of the windowed Bias Score.
So for example a stringency of 1 will give you all the bins that have a rate lower than -1*(standard deviation) of the histogram
of rates described in the report.

