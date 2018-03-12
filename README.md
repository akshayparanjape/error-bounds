# error_bounds_Hiwi_i6

This repository contains code used for error_bounds 

To compile the program -follow the steps below

1. Change the mods (in /beck2) file according to what is available on the system
2. unloads modules by using the command    . beck2/mods
3. making bin folder, if it doesn't create on its own. Sometimes you may also require to create sub folder (like components, distribution)
4. create a soft link while running the experiment   
	ln -s beck2/bin bin
5. Run the experiment file, python experiment.py [config|run|plot]
