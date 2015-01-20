TADLocator
==========

The script and associated data for my rotation in the Ren lab in the Fall Quarter of 2014 at UC San Diego.

The report describing the method used can be seen <a href='https://drive.google.com/file/d/0B6aRfRgg95OGemNGcVU0SDRDaGs/view?usp=sharing'>here </a>.
The presentation can also be seen <a href='https://drive.google.com/file/d/0B6aRfRgg95OGUjBXa1lMbkpnOEU/view?usp=sharing'> here </a>.

<b> Update: getBins function addded </b>
This function allows the user to get the Bin numbers for Domain boundaries in a MATRIX file for a specific stringency level.
Here stringency is defined as the number of standard deviation from the mean of the distribution of rates of the windowed Bias Score.
So for example a stringency of 1 will give you all the bins that have a rate lower than -1*(standard deviation) of the histogram
of rates described in the report.

