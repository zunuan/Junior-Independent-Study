#!/usr/bin/env python
# coding: utf-8

# ## Deinterleave FASTQ

# In[ ]:


fileName = 'SRR8619000.fastq'
parseName = fileName.split('.')[0]
inputFile = open(fileName, 'r')
outputFile_1 = open(parseName+'_1.fastq', 'w')
outputFile_2 = open(parseName+'_2.fastq', 'w')

flag = 1
cnt = 1
for inputLine in inputFile.readlines():
    if cnt == 5:
        cnt = 1
        flag *= -1
    if cnt == 1 or cnt == 3:
        parseIndex = inputLine.find('.', inputLine.find('.') + 1)  # find second delimiter
        inputLine = inputLine[:parseIndex] + inputLine[parseIndex+2:]  # remove different pair read names
    if flag == 1:
        outputFile_1.write(inputLine)
    else:
        outputFile_2.write(inputLine)
    cnt += 1

inputFile.close()
outputFile_1.close()
outputFile_2.close()

