# sim-inspect

Evaluation code of inspect journal article.

## Installation

Clone repository into "scratch" directory of ns-3 distribution.
This code is designed to work with ns-3 version 3.26.
You can clone that version with the following commands:

    git clone --branch ns-3.26 https://github.com/nsnam/ns-3-dev-git.git ns-3.26

Next, enter the directory and checkout this repository at "scratch":

    cd ns-3.26/scratch
    git clone https://github.com/inspect-paper/simulation.git
    cd ..

Finally, build the whole ns-3 distribution (including the simulation) with:

    ./waf configure
    ./waf build
  
You can run the simulation then with...

    ./waf --run 'sim-inspect --help'

...replacing --help with more useful arguments.
