import csv

## Here the chimeras.txt is the list that is created from a grep command in unix that prints outs all the positive chimeras
lol = list(csv.reader(open('chimeras.txt', 'rb'), delimiter = '\t'))
f = open('chimlist.txt','w')
for i in range(0,len(lol)):
    f.write(lol[i][1])
    f.write("\n")


