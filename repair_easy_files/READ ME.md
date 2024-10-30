#Repairing easy files

Most of the time when calling load() or pop_easy() in Matlab it will call some version of load() that is expecting data to be
uniform in nature (same number of columns and rows). However, if the enobio recording is stopped abruptly then this will
interupt the saving process and the number of values saved may not match the number of columns required. Similarly, if the 
recording is stopped and then started again it may be stopped with an incomplete sample. This is definitely a bug within the 
enobio recording software.

##Fixing files manually
To fix this right click and open the .easy file in a text editor (such as notepad), this will open the whole file and you can
look at the data. Often the issue is on the last line so you can skip to the end and see if the last line matches up with
everything else. If not delete the last line and try to load the file again.

If this doesn't work then you will need to look for errors within the .easy file. You can do this by importing the data into
excel. Import the data as if it were a csv file with newline parameters. Then go to data > filter. Click on the down filter
arrow for the last column. You're expecting to see "(Blanks). Filter for just the blanks and note down which rows they are in.
Next unfilter the data and take a look at the rows that have been noted down. You need to replace any bad data within the 
datastream with 0s or Nans (probably 0s though). If this is a line or two it may be worth doing this manually. 

If however, there is a lot of time that needs to be replaced with 0s you can use the scripts in this folder. To do this you will
need to separate the good data from the data that needs to be replaced withs 0s. Copy the good data above the break (bad data)
into a separate excel spreadsheet. Save this as a .easy file with something appropriate like the filename plus "_1", then save
the data after the break (bad data) in a separate file plus "_2". If there are multiple breaks to repair then repeat as necessary.

Next run the scripts.

##Stitching scripts automatically
The scripts simply load the data in the order that they are found in the folder. For each pair of files it takes the last 
timestamp from the first file and the first timestampe from the second one. It calculates the gap of missing data (assuming a
sampling rate of 500Hz) and then concatenates the files as file_1 + blanks + file_2. It does this for all files in the folder
until they are all concatenated.

The only thing to keep in mind is to make sure that only files you would like to concatenate together are in the folder that these
scripts are run over.
