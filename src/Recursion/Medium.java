package Recursion;

import java.util.List;

public class Medium {
    public static void main(String[] args) {

    }
    /* SubSequences and SubSets

    A contiguous/non contiguous sequence which follows the order
    An empty {}  array is also a subsequence
   total no. of subsequences --> 2^n

Subsets --> order does not matter but only considers single pair example
{2,1} is the same subset as {1,2} since order doesn't matter in subsets.

     */

    //    Print all subsequences of an array using recursion

    //    Using Recursion
    public static void print_ss(int[] arr, List<Integer> ds, int idx) {
/*
Approach --> Take it/not take it
//T.C --> 2^N * N
//S.C -->  N  (Stack space )
 */
//Base case
        if (idx == arr.length) {
            System.out.println(ds);
            return;
        }
        ds.add(arr[idx]);
        print_ss(arr, ds, idx + 1); //take it
        ds.remove(arr[idx]);
        print_ss(arr, ds, idx + 1); //not take it

    }
}