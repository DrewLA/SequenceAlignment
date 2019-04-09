Author: Andrew Lewis
Date: 11-26-2018
Program: Alignment.java

Compilation: javac Alignment.java
Execution: Java Alignment 
Dependencies: Java HashMap, Java File.

This implementation of a sequence alignment class utlizes a hashmap to store the endoded characters against indices in a similarity matrix.


The similarity matrix is read from test.txt and contains the pair costs of the corresponding DNA Nucleotide bases A,C,G,T.
The input order must always correspond in the specidifed order inorder to be correctly processed by the hashmap. 

getMin() takes as input the last two lines of test.txt which are the DNA sequences to be compared.
getMin prints out a paired sequence(top-bottom).
The pair sequence, and their corresponding costs are the minimum possible total cost of an alignment between the two sequences.

The alogrithm behind getMin is based on the recurrence relation:
            OPT(i,j) = min[cost(si,ji)+OPt(i-1,j-1), gapCost+OPT(i,j-1), gapCost+OPT(i-1,j)]
    
    The optiminal alignment at each point of the sequence can be thought of as the best case between the next three possible alignments. 

the dashweight is read in as the first integer in the test.txt file. 

Project 2, CPSC 320 - Trinity College
<a href="http://www.cs.trincoll.edu/~miyazaki/cpsc320/project2.html">