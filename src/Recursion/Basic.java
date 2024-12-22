package Recursion;

import java.util.Arrays;

public class Basic {
    public static void main(String[] args) {
//print_n(5);
        System.out.println(sum_n(5666));
        System.out.println(sumTillN(5666));
        int[] arr = {1, 2, 3, 4};
        System.out.println(Arrays.toString(arr));
        rev_arr(arr, 0, arr.length - 1);
        System.out.println(Arrays.toString(arr));
        System.out.println(fib(16));
    }

    //    print 1 to n
    public static void print_n(int n) {
        //Base Case
        if (n <= 0) {
            return;
        }
        print_n(n - 1);
        System.out.print(n + " ");
    }

    //    print 1 to n
    public static int sum_n(int n) {
        //Base Case
        if (n == 0) {
            return n;
        }
        return sum_n(n - 1) + n;
    }

    //    get sum at constant time O(1)
    public static int sumTillN(int n) {
        return n * (n + 1) / 2;
    }

    //Reverse an array using recursion
    public static void rev_arr(int[] arr, int i, int j) {
        //base case
        if (i >= j) {
            return;
        }
        //swap
        int tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
        rev_arr(arr, i + 1, j - 1);
    }

    //    CHECK PALINDROME
    public static boolean is_pal(String str, int i, int j) {
//        Base Case l
        if (i >= j) {
            return true;
        }
        if (str.charAt(i) != str.charAt(j)) return false;
        return is_pal(str, i + 1, j - 1);
    }

    //    MULTIPLE RECURSION CALLS
//fibonacci using loop T.C -> O(2pow(N))
    public static int fib(int n) {
//0 1 1 2 3 5 8 13 . . .  . . . . . . . .. . .
        if (n < 2) return n;
        return fib(n - 1) + fib(n - 2);
    }

    // Recursion Based Sorting Algorithms

    //    MERGE SORT
    public static int[] merge(int[] arr1, int[] arr2) {
//Merges two sorted arrays
        int i = 0, j = 0, k = 0, len1 = arr1.length, len2 = arr2.length;
        int[] res = new int[len1 + len2];
        while (i < len1 && j < len2) {
            if (arr1[i] <= arr2[j]) {
                res[k] = arr1[i];
                i++;
            } else {
                res[k] = arr2[j];
                j++;
            }
            k++;
        }

        if (i < len1) {
            while (i < len1) {
                res[k] = arr1[i];
                i++;
                k++;
            }
        }
        if (j < len2) {
            while (j < len2) {
                res[k] = arr2[j];
                j++;
                k++;
            }
        }
        return res;
    }

    public static int[] merge_sort(int[] arr, int i, int j) {
//        Base Case
        if (i == j) {//array contains single element
            return new int[]{arr[i]};
        }

        int mid = i + (j - i) / 2;
        int[] left = merge_sort(arr, i, mid);
        int[] right = merge_sort(arr, mid + 1, j);

        return merge(left, right);
//    T.C --> O(N LOG2(N))  logarithmic bcoz problem is dividing by 2 everytime
//    S.C -> O(n) due to the additional space required for merging subarrays
    }

    public static void swap(int[] arr, int i, int j) {
        int tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
    }

    //    Recursive Bubble sort
    public static void rec_Bub_sort(int[] arr, int a, int b, int n) {


    }
//    Recursive Selection  sort


//    Recursive Insertion sort


    //    Recursively check is array sorted
    public static boolean is_sorted(int[] arr, int i) {
//Base case
        if (i == arr.length - 1) {
            return true;//single element is always sorted
        }
        if (arr[i] > arr[i + 1]) return false;

        return is_sorted(arr, i + 1);
    }

    /*    QUICK SORT  --> DIVIDE & CONQUER
    1. Pick an element and place it at its correct place in an sorted array
        pivot can be 1st , last , median , random doesn't matter

    2. smaller on the left & larger on the right
       after 1st iteration that pivot will at its correct place

    3.place the pivot at its correct place get the partition index

    4. make call for the left and for the right

    */

    public static int getPartidx(int[] arr, int lo, int hi) {

        int pivot_ele = arr[lo];
        int i = lo, j = hi;
        while (i < j) {
            while (arr[i] <= pivot_ele && i < hi) i++;
            while (arr[j] > pivot_ele && j > lo) j--;
            if (i < j)
                swap(arr, i, j);
        }
        swap(arr, pivot_ele, arr[j]);
        return j;
    }

    public static void quick_sort(int[] arr, int lo, int hi) {
        if (lo < hi) { //single element itself sorted
            int partition_idx = getPartidx(arr, lo, hi);
            quick_sort(arr, lo, partition_idx - 1);
            quick_sort(arr, partition_idx + 1, hi);
        }

//        T.C--> O(N LOG N )

    }
}