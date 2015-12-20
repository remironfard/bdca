# Introduction #

BDCA for EEG aims to identify discriminating components in the EEG. Each component is represented by a spatial and a temporal activity profile. Ambiguities are resolved using independence across trials.

  * Read the [Tutorial](http://code.google.com/p/bdca/wiki/Tutorial) to see how you should tune the parameters.

# How to obtain the code #

**The best way to get the code is this:**

  1. Make sure SVN is installed on your machine ([why?](http://en.wikipedia.org/wiki/Subversion_(software))):
    * (Mac) To install SVN on a Mac follow steps 1,2,3 given [here](http://www.wikihow.com/Install-ubversion-on-Mac-OS-X)
    * (Linux) To install SVN on Linux is very easy using apt-get or find it in the synaptic package manager.
    * (Windows / Redhat / Solaris) go [here](http://www.collab.net/downloads/subversion/)
  1. Click the "Source" tab or click [here](http://code.google.com/p/bdca/source/checkout) to see exactly how you use SVN to receive the up-to-date source code. It will tell you to run the following command from a terminal:
```
svn checkout http://bdca.googlecode.com/svn/trunk/ bdca-read-only
```
the source will then be downloadet to a directory called _bdca-read-only_ under your current directory.

**Alternatively, (NOT RECOMMENDED)**

  * Click the "Downloads" tab, and download a release which will NOT be up to date with bug fixes and new features.

# Dependencies #

The BDCA toolbox relies (currently) on other toolboxes. You need to download and install the following toolboxes in order to use the BDCA toolbox:

  * [IMMOPTIBOX](http://www2.imm.dtu.dk/~hbn/immoptibox/) - Get the file: [download](http://www2.imm.dtu.dk/~hbn/immoptibox/immoptibox.zip)

# Staying up to date with bugfixes and new features #

If the source code on this page changes you can get the updates in a simple step:
  1. In a terminal, go to the _bdca-read-only_ directory where the toolbox is located on your machine.
  1. run the following command from the terminal:
```
svn update
```
that's it!
# References #

  * Dyrholm, M., Christoforou, C., Parra, L. C., ''Bilinear Discriminant Component Analysis'', _Journal of Machine Learning Research_, 8(May):1097-1111, 2007

  * Dyrholm, M., Parra, L. C., ''Smooth bilinear classification of EEG'', In _proceedings of the 28th Annual International Conference of the IEEE Engineering in Medicine and Biology Society_, New York, August 2006.