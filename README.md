# Fall-2018-Project-3

## Group members
Ian Alvarez ian_alvarez23@csu.fullerton.edu

Francis Nguyen lqtruongnguyen@csu.fullerton.edu

Antonio Lopez antonio_lopez@csu.fullerton.edu

## The Problem

Longest common subsequence

Input: Two strings: S1 and S2 and a scoring matrix P.

Output: The score of the substring of S2 which has the maximum score of all possible substrings of S1 and the strings representing the alignment of the two strings. 

## The Algorithm

    for i = 1 to n

      for j = 1 to m
    
         Get values from up, left and diag and add appropriate penalties
         
         Find maximum of three values (or 0, whichever is greater – important!)
         
         Set the maximum value in the array
         
         Set the value for the back pointer according to which direction was chosen
         
    Find the best score from the bottom row

    Starting at location of best score, follow backtrack pointers until you get to “?” translating the decisions back into the two alignment strings (including markers for deletion and insertion).

    Reverse the strings (because we do it backwards).

## Results

![alt text](https://github.com/CSUF-CPSC-335-FA18-Bernstein/project-three-every_villians_is_lemons/blob/master/project3_results.PNG)
