How to run:

1) SW algorithm:

Copy paste all the source code in SW_SP.py into an IDE like PyCharm and press option+shift+E to run them in the console and load them into memory so they are global. Alternatively, you can copy paste into Jupiter notebook. You can then run fasta files like test_input.fa (without the comments and uneven spacing) by calling read_do_sw(file)), it's pairwise, so all the sequences should be paired with headers like '>seq1/2'. Can run on individual sequences by calling smith_waterman(seq1, seq2). 

2) KMeans algorithm:
Copy paste all the source code in KM_SP.py into an IDE like PyCharm and press option+shift+E to run them in the console and load them into memory so it is global. The file line at the top of the file should point to the location on your computer that the Biase_2014.csv or gene-expression data, the lines below should be run afterward and are used to standardize the data for all future  functions and operations. They are in main so they are globally available. Following the commented line "#################################CHECK WITH KMEANS##############################" is the code to run, kmeans(file_scaled, num_clus, rand), where file_scaled is the standardized file, num_clus is the number of clusters to create and rand is the random number initializer to seed the random centroid function. Following lines are to perform the intercluster comparisons and p-value calculations. 
