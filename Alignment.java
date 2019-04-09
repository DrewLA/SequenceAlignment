/**
 * This program prints out the minimum cost alignment of two given DNA sequences. 
 *  
 * It reads data into a similarity matrix and maps the corresponding DNA codes to the indices of the matrix.
 * 
 * It implements a getMin() function which prints out the paired aligned sequence with the minimum cost. 
 * getMin() runs in O(mn) time.
 * 
 * It is designed specifically for Project 2, CPSC 320 - Trinity College
 * <a href="http://www.cs.trincoll.edu/~miyazaki/cpsc320/project2.html">
 * 
 * Compilation: javac Alignment.java
 * Execution: Java Alignment 
 * Dependencies:  Java HashMap, Java File.
 * 
 * Package with a test file document "test.txt" which can be modified under the specified encoding.
 * 
 * @author Andrew Lewis
 * @date 26-11-2018
 */

import java.util.*;
import java.util.HashMap;
import java.util.Scanner;
import java.io.File;
import java.io.FileNotFoundException;
import java.lang.Integer;

/**
 * This class 'Alignment' takes input from test.txt and utilizes a global matrix to store the pair-costs
 * representation of respective sequence elements. 
 * 
 * The helper function getMin is specifically taiored to use a mapping of the DNA codes 'A''C''G''T' to an
 * inputed pair-costs (similarity) matrix. 
 * 
 * getMin prints outs the pair sequence alignment with the minimum cost based on the prescirbed gap weight 
 * and the data from the costs matrix.  
 *    */
public class Alignment {

    /**
     * similarity matrix to store corresponding penalities
     */
    static int matrix[][];

    /**
     * Variable for the gap penalty
     */
    static int dashWeight;

    /**
     * This method prints out the least cost alignment of two DNA endoded string sequences. 
     * 
     * @param s DNA sequence A
     * @param t DNA sequence B
     */
    public static void getMin(String s, String t) {
        HashMap<Character, Integer> codes = new HashMap<>();
        //Mapping the DNA codes to indices of the similarity matrix
        codes.put('A', 0);codes.put('C', 1);codes.put('G', 2);codes.put('T', 3);
        int i, j; // intialising variables 
        int cost;
        int m = s.length(); // length of dna1 
        int n = t.length(); // length of dna2 
        
        // graph table for ieterated solutions to the  optimal subproblems  
        int g[][] = new int[n + m + 1][n + m + 1]; 
        //Initialize g
        for (int[] x1 : g) 
            Arrays.fill(x1, 0); 
    
        //Case for empty sequence
        for (i = 0; i <= (n + m); i++){ 
            g[i][0] = i * dashWeight; 
            g[0][i] = i * dashWeight; 
        }  

        // calcuting the  minimum penalty 
        for (i = 1; i <= m; i++) {       
             for (j = 1; j <= n; j++){ 
                cost = matrix[codes.get(s.charAt(i-1))][codes.get(t.charAt(j-1))];
                    //System.out.println(x.charAt(i-1)+" "+y.charAt(j-1)+" "+cost);
                if (s.charAt(i - 1) == t.charAt(j - 1)) { 
                    g[i][j] = g[i - 1][j - 1]; 
                } 
                else{
                    //System.out.println("gotten");
                    g[i][j] = Math.min(Math.min(g[i - 1][j - 1] + cost ,  
                                                g[i - 1][j] + dashWeight),  
                                                g[i][j - 1] + dashWeight); 
                } 
            } 
        } 
        
        /**
         * <-----------------------Traversing G for the valid pairs------------------------->
         */
        
        int l = n + m; //max g index
        
        i = m; j = n; //Traversal variables
        
        int spos = l; //index of sequence s
        int tpos = l; //index of sequence t

        // Final sequences for the DNA alignment 
        int dnaA[] = new int[l + 1];  
        int dnaB[] = new int[l + 1]; 
        //traverse G and select the nodes that satisify the minimum costs path
        while ( !(i == 0 || j == 0)){ 
            cost = matrix[codes.get(s.charAt(i-1))][codes.get(t.charAt(j-1))];
            if (s.charAt(i - 1) == t.charAt(j - 1)) { 
                dnaA[spos--] = (int)s.charAt(i - 1); 
                dnaB[tpos--] = (int)t.charAt(j - 1); 
                i--; j--; 
            }  
            else if (g[i - 1][j - 1] + cost == g[i][j]) { 
                dnaA[spos--] = (int)s.charAt(i - 1); 
                dnaB[tpos--] = (int)t.charAt(j - 1); 
                i--; j--; 
            } 
            else if (g[i - 1][j] + dashWeight == g[i][j]) { 
                dnaA[spos--] = (int)s.charAt(i - 1); 
                dnaB[tpos--] = (int)'-'; 
                i--; 
            } 
            else if (g[i][j - 1] + dashWeight == g[i][j]) { 
                dnaA[spos--] = (int)'-'; 
                dnaB[tpos--] = (int)t.charAt(j - 1); 
                j--; 
            } 
        } 

        while (spos > 0) { 
            if (i > 0) dnaA[spos--] = (int)s.charAt(--i); 
            else dnaA[spos--] = (int)'-'; 
        } 
        while (tpos > 0) { 
            if (j > 0) dnaB[tpos--] = (int)t.charAt(--j); 
            else dnaB[tpos--] = (int)'-'; 
        } 
        //Counting gaps 
        int cnt = 1; 
        for (i = l; i >= 1; i--) { 
            if ((char)dnaB[i] == '-' && (char)dnaA[i] == '-') { 
                cnt = i + 1; 
                break; 
            } 
        } 

        // Printing the final sequence 
        char a,b; //varaibles for storing pair values
        for (i = cnt; i <= l; i++) {   
            System.out.print((char)dnaA[i]); 
        } 
        System.out.print("\n"); 
        for (i = cnt; i <= l; i++) {   
            System.out.print((char)dnaB[i]);  
        } 
        System.out.println("");
        //Printing the respective costs in sequence pairs
        for (i = cnt; i <= l; i++) {
            a = ((char)dnaA[i]); b = ((char)dnaB[i]);
            if (a == '-' || b == '-') {
                System.out.print(dashWeight);
            }
            else{
               System.out.print(matrix[codes.get(a)][codes.get(b)]);
            }
        }
        System.out.println("\n");
        System.out.print("with the minimum edit distance of "); 
        System.out.print(g[m][n] + ".\n");
        return; 
    }

    /**
     * Program driver
     * @param arg
     */
	public static void main(String[] args) {
        
        String dna1;    //Variable for sequence 1
        String dna2;    //Variable for sequence 2
        //Initialize similarity matrix
        matrix = new int[4][4]; 

        try {
            File test = new File("test.txt");   //"test.txt" is the test file found in the relative root of main class
            Scanner sc = new Scanner(test);
            dashWeight = sc.nextInt();
            sc.nextLine();
            //extract data from each line and add to similarity matrix.
            for (int i = 0; i < 4; i++){
                String[] lines = sc.nextLine().split(" ");
                for (int j = 0; j < 4; j++) {                    
                    matrix[i][j]= Integer.parseInt(lines[j]);
                }
           
            }
            dna1 = sc.next();   //Read and store first sequence
            dna2 = sc.next();   //Read and store second sequence
            System.out.println("The best alignment is: "+ "\n"); 
            getMin(dna1, dna2);
           
            sc.close();
        } catch (Exception e) {
            //TODO: handle exception
            System.out.println(e);
        }         
		
	}

}